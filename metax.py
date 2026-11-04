#!/usr/bin/env python3
"""MangroveMetaX Ultra workflow controller (Python stdlib; core uses nf-core/mag)."""
import argparse,csv,json,os,shutil,subprocess,sys,platform,shlex
from pathlib import Path
from scripts.host_paths import conda_root
BASE=Path(__file__).resolve().parent

def load(project):
    path=Path(project).resolve(); conf=json.loads(path.read_text());conf['_base']=str(BASE);conf['_config_path']=str(path)
    def absolutize(p):
        q=Path(p).expanduser()
        return q.resolve() if q.is_absolute() else (BASE/q).resolve()
    conf['sample_path']=absolutize(conf['samples_tsv'])
    conf['root']=absolutize(conf['results_dir']);conf['raw']=conf['root']/'nfcore_mag'
    conf['asset_path']=conf['root']/'assets.json';conf['ext']=conf['root']/'extensions'
    conf['mag_review_path']=absolutize(conf.get('mag_review_tsv','config/approved_mags.tsv'))
    return conf

def rows_for(conf):
    with open(conf['sample_path']) as f:return list(csv.DictReader(f,delimiter='\t'))
def spawn(cmd,dry_run=False,env=None):
    print('+',subprocess.list2cmdline([str(x) for x in cmd]),flush=True)
    if not dry_run:subprocess.run([str(x) for x in cmd],check=True,env=env)
def require_file(p):
    if not Path(p).is_file():raise FileNotFoundError(f'Expected file: {p}')
def validate_run(conf,full=False,dry=False):
    from scripts.validate_inputs import read_sheet,validate_fastq_pair
    records=read_sheet(conf['sample_path'],check_files=not dry)
    for r in records:
        if not dry:
            summary=validate_fastq_pair(r['r1'],r['r2'],0 if full else 10000)
            print(r['sample'],summary)
    return records

def present_db(conf, key):
    raw=(conf.get('reference_databases') or {}).get(key) or ''
    if not str(raw).strip():return ''
    path=Path(str(raw)).expanduser()
    path=path if path.is_absolute() else (BASE/path)
    return str(path) if path.exists() else ''

def gtdb_memory_problem(conf):
    """Return a blocker string when pplacer cannot spill to disk, else ''."""
    need=float(conf.get('resources',{}).get('gtdbtk_min_ram_plus_swap_gb',96))
    try:
        info={}
        for line in Path('/proc/meminfo').read_text().splitlines():
            name,value=line.split(':',1)
            info[name]=int(value.strip().split()[0])/(1024*1024)
        total=info.get('MemTotal',0)+info.get('SwapTotal',0)
    except (OSError,ValueError):
        return ''
    if total+0.5>=need:return ''
    return (f'GTDB-Tk pplacer needs about {need:.0f} GB of RAM plus swap. This host has '
            f'{info.get("MemTotal",0):.0f} GB RAM and {info.get("SwapTotal",0):.0f} GB swap. '
            'The pplacer index is loaded whole, outside the Nextflow memory cap. '
            'Add a swap file until RAM plus swap reaches that total, then enable it with mkswap and swapon.')

def build_core(conf):
    if conf['core_engine']!='nfcore_mag':raise ValueError('Use nfcore_mag production engine; SqueezeMeta is benchmark-only, never simultaneously in production.')
    mode=conf.get('resource_mode','research')
    if mode not in ('mac64','research','hpc','linux64'):
        raise ValueError('resource_mode must be mac64, research, hpc or linux64')
    if mode=='mac64' and conf.get('core_enable_gtdbtk'):
        raise ValueError('GTDB-Tk pplacer is prohibited in mac64 core; use taxonomy-ani or explicit taxonomy-hpc on a high-RAM machine')
    if mode=='mac64' and conf.get('core_benchmark_spades'):
        raise ValueError('metaSPAdes benchmark is not enabled in mac64 mode; use research/HPC mode to benchmark it')
    from scripts.validate_inputs import read_sheet
    rows=read_sheet(conf['sample_path'],False)
    if conf['assembly_mode'] not in ('per_sample','coassembly_group'):raise ValueError('Unsupported assembly_mode')
    if conf['assembly_mode']=='coassembly_group':raise ValueError('coassembly_group requires group-aware asset bridging, not yet supported. Use per_sample or curate coassembly manually.')
    work=conf['root'];work.mkdir(parents=True,exist_ok=True)
    sheet=work/'nfcore_samples.csv'
    with open(sheet,'w',newline='') as f:
        writer=csv.writer(f);writer.writerow(['sample','group','short_reads_1','short_reads_2','long_reads','short_reads_platform'])
        for r in rows:writer.writerow([r['sample'],r['group'],r['r1'],r['r2'],'','ILLUMINA'])
    command=['nextflow','run','nf-core/mag','-r',conf['nfcore_mag_version'],'-profile',conf['profile'],
             '--input',str(sheet),'--outdir',str(conf['raw']),'--save_clipped_reads','--save_assembly_mapped_reads',
             '--refine_bins_dastool','--postbinning_input','refined_bins_only','--bin_domain_classification','--run_checkm2','--skip_spades']
    if conf.get('core_benchmark_spades'):
        command=[x for x in command if x!='--skip_spades']
    # MEGAHIT stays the production assembler. metaSPAdes is an explicit benchmark flag, not a second result table.
    if mode=='mac64':
        r=conf.get('resources',{})
        cpu=int(r.get('max_cpus',8));ram=int(r.get('max_process_memory_gb',40))
        if cpu<1 or ram<4 or cpu>8 or ram>40 or int(conf['threads'])>8:
            raise ValueError('mac64 hard safety limits: threads <=8, max_cpus <=8, max process RAM <=40 GB')
        command.extend(['-c',str(BASE/'profiles'/'mac64.config'),'--max_cpus',str(cpu),'--max_memory',f'{ram}.GB'])
        if not r.get('enable_semibin_on_mac',False):command.append('--skip_semibin')
        if not conf.get('core_enable_gtdbtk',False):command.append('--skip_gtdbtk')
    elif mode=='linux64':
        r=conf.get('resources',{})
        cpu=int(r.get('max_cpus',24));ram=int(r.get('max_process_memory_gb',48))
        if cpu<1 or ram<8 or cpu>32 or ram>56:
            raise ValueError('linux64 limits: 1-32 CPUs and at most 56 GB per container on a 64 GB host')
        command.extend(['-c',str(BASE/'profiles'/'linux64.config'),'--max_cpus',str(cpu),'--max_memory',f'{ram}.GB'])
        # pplacer stays outside this 48 GB cgroup. See taxonomy_gtdb.
        command.append('--skip_gtdbtk')
        if conf.get('core_enable_busco',False):
            command+=['--run_busco','--busco_db_lineage','auto']
        splits=str(r.get('genomad_splits',8))
        # Local paths are used when they already exist. Otherwise nf-core/mag downloads them.
        if conf.get('core_enable_gunc'):
            command.append('--run_gunc')
            if present_db(conf,'gunc'):command.extend(['--gunc_db',present_db(conf,'gunc')])
        if conf.get('core_enable_metaeuk'):
            if present_db(conf,'metaeuk'):command.extend(['--metaeuk_db',present_db(conf,'metaeuk')])
            else:command.extend(['--metaeuk_mmseqs_db',str(conf.get('metaeuk_mmseqs_db') or 'UniRef50')])
        if conf.get('core_enable_genomad'):
            command+=['--run_virus_identification','--genomad_splits',splits]
            if present_db(conf,'genomad'):command.extend(['--genomad_db',present_db(conf,'genomad')])
    if conf.get('binning_mode')=='resource_aware':
        command+=['--skip_concoct','--skip_comebin','--skip_metabinner']
        if '--skip_semibin' not in command:command.append('--skip_semibin')
    elif '--skip_semibin' not in command and conf.get('semibin_environment'):
        command+=['--semibin_environment',str(conf['semibin_environment'])]
    db=conf['reference_databases']
    omitted=[]
    def attach(enabled,key,flag):
        path=present_db(conf,key) if mode=='linux64' else (db.get(key) or '')
        if enabled and path:
            command.extend([flag,path] if flag else [])
            if key=='gunc' and '--run_gunc' not in command:command.append('--run_gunc')
            return True
        if enabled and mode!='linux64' and key=='metaeuk':
            raise ValueError('MetaEuk database required')
        if enabled and mode!='linux64' and key=='gtdbtk':
            raise ValueError('GTDB-Tk database required when enabled')
        if enabled and mode!='linux64' and key=='gunc':
            raise ValueError('GUNC database required when enabled')
        if enabled:omitted.append(key)
        return False
    # GTDB-Tk pplacer is never attached to the linux64 Nextflow command.
    if mode!='linux64':
        attach(conf.get('core_enable_gtdbtk'),'gtdbtk','--gtdb_db')
        attach(conf.get('core_enable_gunc'),'gunc','--gunc_db')
        attach(conf.get('core_enable_metaeuk'),'metaeuk','--metaeuk_db')
        if conf.get('core_enable_genomad'):
            path=db.get('genomad') or ''
            if path:
                splits=str(conf.get('resources',{}).get('genomad_splits',8))
                command+=['--run_virus_identification','--genomad_db',path,'--genomad_splits',splits]
            else:
                raise ValueError('geNomad DB required')
    if present_db(conf,'cat_db'):
        command+=['--cat_db',present_db(conf,'cat_db')]
    if '--skip_gtdbtk' not in command and not conf.get('core_enable_gtdbtk'):command.append('--skip_gtdbtk')
    conf['_core_omitted']=omitted
    command+=['-resume']
    return command

def shell_status(conf):
    root=conf['root'];core=root/'nfcore_mag'; assets=conf['asset_path']
    d={'core_result_exists':core.is_dir(),'core_multiqc':(core/'multiqc'/'multiqc_report.html').is_file(),
       'asset_bridge':assets.is_file(),'extensions':sorted(x.name for x in conf['ext'].glob('*')) if conf['ext'].exists() else []}
    print(json.dumps(d,indent=2))

def main():
    p=argparse.ArgumentParser(description='MangroveMetaX Ultra v3.1: auditable, resumable run-all research workflow')
    p.add_argument('--config',default=str(BASE/'config/project.json'))
    sp=p.add_subparsers(dest='action',required=True)
    sp.add_parser('plan');sp.add_parser('status');sp.add_parser('doctor');sp.add_parser('audit');sp.add_parser('bridge');sp.add_parser('mag-review')
    ar=sp.add_parser('run-all',help='Run all computational branches and write execution report')
    ar.add_argument('--dry-run',action='store_true'); ar.add_argument('--allow-partial',action='store_true');ar.add_argument('--force',action='store_true')
    ai=sp.add_parser('taxonomy-ani');ai.add_argument('--dry-run',action='store_true')
    ix=sp.add_parser('skani-index');ix.add_argument('--reference-list',required=True);ix.add_argument('--out',required=True);ix.add_argument('--dry-run',action='store_true')
    g=sp.add_parser('taxonomy-hpc');g.add_argument('--run',action='store_true');g.add_argument('--acknowledge-high-memory',action='store_true')
    sp.add_parser('taxonomy-gtdb')
    v=sp.add_parser('validate');v.add_argument('--full',action='store_true')
    c=sp.add_parser('core');c.add_argument('--dry-run',action='store_true');c.add_argument('--skip-input-validation',action='store_true')
    from scripts.integrations import RUNNERS
    m=sp.add_parser('module');m.add_argument('name',choices=['taxonomy','ssu','euk_contig','eggnog','dbcan','dram2','checkv','fegenie','amrfinder','antismash','metabolic_c',*sorted(RUNNERS)]);m.add_argument('--dry-run',action='store_true')
    a=p.parse_args();conf=load(a.config)
    if a.action=='run-all':
        from scripts.autorun import run_all
        run_all(conf,dry_run=a.dry_run,allow_partial=a.allow_partial,force=a.force)
    elif a.action=='plan':
        with open(BASE/'docs'/'STEP_REGISTRY.tsv') as f:
            for r in csv.DictReader(f,delimiter='\t'):print(f"{r['ID']:>3}  {r['STATUS']:<18} {r['OWNER']:<22} {r['STEP']}")
    elif a.action=='status':shell_status(conf)
    elif a.action=='doctor':
        print('Host:',platform.system(),platform.machine(),'| profile:',conf.get('resource_mode'))
        for x in ['python3','nextflow','docker','skani','kraken2','bracken','phyloFlash.pl','tiara','whokaryote.py','emapper.py','run_dbcan','exec_annotation','mmseqs','coverm','dRep','eukcc','iphop','minced','hmmsearch','humann','nonpareil','FeGenie.py','amrfinder','antismash','checkv','gtdbtk','SqueezeMeta.pl','atlas']:
            print(f'{x:20} {shutil.which(x) or "NOT FOUND (stage will be BLOCKED until installed)"}')
        if conf.get('resource_mode')=='mac64':
            print('Mac64: Docker memory target 48GB; per-task max 40GB; GTDB-Tk pplacer DISABLED. Request caps are NOT physical usage guarantees.')
            if platform.machine() in ('arm64','aarch64'):
                print('Apple Silicon: Linux amd64 containers may require Rosetta/QEMU and can be slow or incompatible.')
        if conf.get('resource_mode')=='linux64':
            r=conf.get('resources',{})
            print(f"Linux64: {r.get('max_cpus',24)} CPUs, {r.get('max_process_memory_gb',48)} GB container cap, one heavy task at a time, GPU containers via nvidia runtime.")
            print('GTDB-Tk pplacer is a separate stage (taxonomy_gtdb), not inside the Nextflow 48 GB cap.')
            print('Core downloads geNomad, GUNC (progenomes), and MetaEuk UniRef50 when those paths are empty.')
            problem=gtdb_memory_problem(conf)
            if problem:print('GTDB-Tk:',problem)
            missing=[k for k,v in conf.get('reference_databases',{}).items() if k!='bracken_read_length' and not present_db(conf,k)]
            if missing:print('Database paths not installed yet:',', '.join(missing))
        print('Reference DBs must be installed separately; doctor does not claim DB/tool compatibility.')
    elif a.action=='validate':validate_run(conf,a.full)
    elif a.action=='core':
        cmd=build_core(conf)
        if not a.skip_input_validation and not a.dry_run:validate_run(conf)
        env=os.environ.copy()
        if conf.get('resource_mode')=='mac64':
            env.setdefault('NXF_OPTS','-Xms1g -Xmx4g')
        elif conf.get('resource_mode')=='linux64':
            heap=int(conf.get('resources',{}).get('nxf_java_heap_gb',6))
            env.setdefault('NXF_OPTS',f'-Xms1g -Xmx{heap}g')
            omitted=conf.get('_core_omitted') or []
            if omitted:print('Core optional databases not installed, omitted from this Nextflow command:',', '.join(omitted))
        spawn(cmd,a.dry_run,env=env)
    elif a.action=='bridge':
        from scripts.bridge_nfcore import build
        data=build(conf['raw'],conf['asset_path'],rows_for(conf),conf['assembler'].upper(),conf['mag_source_glob'])
        print(f"Bridged {len(data['samples'])} samples and {len(data['mags'])} MAGs; warnings={data['warnings']}")
    elif a.action=='mag-review':
        from scripts.mag_ani import make_review
        require_file(conf['asset_path'])
        data=json.loads(conf['asset_path'].read_text())
        path=make_review(data,conf['mag_review_path'])
        print('Review CheckM2 results, fill QC and approved=true only after review:',path)
    elif a.action=='skani-index':
        cmd=['skani','sketch','-l',str(Path(a.reference_list).resolve()),'-o',str(Path(a.out).resolve()),'-t',str(conf['threads'])]
        print('skani DB construction requires GTDB reference FASTAs and substantial disk space. Run separately from core.');spawn(cmd,a.dry_run)
    elif a.action=='taxonomy-ani':
        from scripts.mag_ani import run_ani
        run_ani(conf,a.dry_run)
    elif a.action=='taxonomy-gtdb':
        from scripts.mag_ani import hpc_gtdb_command
        if not conf.get('gtdbtk_accept_swap'):
            raise ValueError('Set gtdbtk_accept_swap=true before running pplacer on a 64 GB host.')
        problem=gtdb_memory_problem(conf)
        if problem:raise ValueError(problem)
        db=present_db(conf,'gtdbtk')
        if not db:raise FileNotFoundError('Set reference_databases.gtdbtk to an existing GTDB-Tk data directory')
        r=conf.get('resources',{})
        cmd=hpc_gtdb_command(conf,out_name='taxonomy_gtdb',cpus=int(r.get('gtdbtk_classify_cpus',4)),pplacer_cpus=int(r.get('gtdbtk_pplacer_cpus',1)))
        (conf['root']/'taxonomy_gtdb'/'scratch').mkdir(parents=True,exist_ok=True)
        env=os.environ.copy();env['GTDBTK_DATA_PATH']=db
        env['PATH']=str(conda_root()/'envs'/'metax-gtdbtk'/'bin')+os.pathsep+env.get('PATH','')
        wrapped=[]
        if shutil.which('nice'):wrapped+=['nice','-n','19']
        if shutil.which('ionice'):wrapped+=['ionice','-c','3']
        spawn(wrapped+cmd,env=env)
    elif a.action=='taxonomy-hpc':
        from scripts.mag_ani import hpc_gtdb_command
        cmd=hpc_gtdb_command(conf,dry=not a.run)
        if a.run:
            if conf.get('resource_mode')!='hpc' or not a.acknowledge_high_memory:
                raise ValueError('Only runs with resource_mode=hpc AND --acknowledge-high-memory on an adequately provisioned Linux host.')
            (conf['root']/'taxonomy_gtdb_hpc'/'scratch').mkdir(parents=True,exist_ok=True)
            if not (conf.get('reference_databases',{}).get('gtdbtk')):
                raise ValueError('Set reference_databases.gtdbtk to the installed R232 path')
            env=os.environ.copy();env['GTDBTK_DATA_PATH']=conf['reference_databases']['gtdbtk']
            spawn(cmd,env=env)
        else:
            print('Optional high-memory command (NOT executed):',shlex.join(cmd))
    elif a.action=='audit':
        require_file(conf['asset_path']);spawn([sys.executable,str(BASE/'scripts/audit.py'),'--assets',str(conf['asset_path']),'--out',str(conf['root']/'audit')])
    elif a.action=='module':
        from scripts.advanced_modules import run_module
        run_module(conf,a.name,a.dry_run)

if __name__=='__main__':
    try:main()
    except (FileNotFoundError,FileExistsError,ValueError,RuntimeError,subprocess.CalledProcessError,OSError) as err:sys.exit('ERROR: '+str(err))
