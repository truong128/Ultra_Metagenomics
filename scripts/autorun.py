"""Fail-transparent, resumable, dependency-aware end-to-end computational scheduler.

Never marks a stage successful from a zero return code alone if expected output is absent.
Does not claim to perform wet-lab analyses or infer ecological activity from DNA.
"""
import csv,datetime,hashlib,json,os,shlex,shutil,subprocess,sys,traceback
from pathlib import Path
from dataclasses import dataclass

BASE=Path(__file__).resolve().parents[1]

def _nonpareil_ready(prefix):
    """Nonpareil 3.5.5 writes the curve to -o itself. Older builds append .npo or .npa."""
    prefix=Path(prefix)
    named=[Path(str(prefix)+'.npo'), Path(str(prefix)+'.npa')]
    if any(p.is_file() and p.stat().st_size for p in named):
        return True
    if prefix.is_file() and prefix.stat().st_size:
        with prefix.open() as handle:
            return handle.readline().startswith('# @impl: Nonpareil')
    return False

def _nonpareil_done(folder, conf):
    from scripts.validate_inputs import read_sheet
    rows=read_sheet(conf['sample_path'], False)
    return bool(rows) and all(_nonpareil_ready(folder/r['sample']) for r in rows)

@dataclass(frozen=True)
class Stage:
    id:str
    deps:tuple=()
    executable:str=''
    db:str=''
    mag_required:bool=False
    description:str=''

STAGES=[
 Stage('validate',description='FASTQ integrity'),
 Stage('core',('validate',),'nextflow',description='nf-core/mag QC, assembly, mapping, Prodigal, binning, CheckM2, optional geNomad/MetaEuk'),
 Stage('bridge',('core',),description='Link core FASTAs, MAGs, sample identities'),
 Stage('audit',('bridge',),description='Compute actual assembly and protein metrics'),
 Stage('mag_auto',('bridge',),description='Select MAGs by CheckM2 QC thresholds, not manual quality approval'),
 Stage('taxonomy_ani',('mag_auto',),'skani','skani_representatives',True,'GTDB representative ANI candidate screen'),
 Stage('taxonomy',('bridge',),'kraken2','kraken2',False,'Read-level taxonomy, independent of assembly'),
 Stage('bracken',('taxonomy',),'bracken','kraken2',False,'Read abundance correction'),
 Stage('ssu',('bridge',),'phyloFlash.pl','phyloflash_silva',False,'SSU/18S reconstruction'),
 Stage('euk_contig',('bridge',),'tiara',description='Eukaryotic contig evidence'),
 Stage('eggnog',('bridge',),'emapper.py','eggnog_data',False,'eggNOG annotation'),
 Stage('dbcan',('bridge',),'run_dbcan','dbcan',False,'CAZyme functions'),
 Stage('dram2',('mag_auto',),'nextflow','',True,'Genome-resolved metabolic reconstruction'),
 Stage('fegenie',('mag_auto',),'FeGenie.py','',True,'Iron metabolism'),
 Stage('functional_summary',('eggnog',),description='Per-sample KEGG Ortholog gene annotation counts (not abundance)'),
 Stage('figures',('audit',),description='Descriptive FASTA-based SVG figures'),
 Stage('amrfinder',('bridge',),'amrfinder',description='AMR candidate screening'),
 Stage('antismash',('mag_auto',),'antismash','',True,'BGC discovery'),
 Stage('metabolic_c',('mag_auto',),'perl','',True,'Community metabolism, legacy METABOLIC-C Perl dependencies'),
 Stage('checkv',('bridge',),'checkv','checkv',False,'Viral genome completeness'),
 Stage('taxonomy_gtdb',('mag_auto',),'gtdbtk','gtdbtk',True,'GTDB-Tk classify_wf with one pplacer CPU and disk scratch'),
 Stage('whokaryote',('bridge',),'whokaryote.py',description='Gene-structure eukaryotic versus prokaryotic contigs'),
 Stage('nonpareil',('validate',),'nonpareil',description='Read-set coverage and redundancy'),
 Stage('kofam',('bridge',),'exec_annotation','kofam',False,'KOfamScan diagnostic KOs'),
 Stage('catalogue',('bridge',),'mmseqs',description='Nonredundant protein catalogue'),
 Stage('coverm',('mag_auto',),'coverm','',True,'MAG abundance from existing BAMs'),
 Stage('drep',('mag_auto',),'dRep','',True,'MAG dereplication across samples'),
 Stage('eukcc',('mag_auto',),'eukcc','eukcc',True,'Eukaryotic MAG completeness'),
 Stage('iphop',('checkv',),'iphop','iphop',False,'Viral host prediction from geNomad FASTA'),
 Stage('crispr',('mag_auto',),'minced','',True,'CRISPR arrays on QC-threshold MAGs'),
 Stage('taxprofiler',('validate',),'nextflow','taxprofiler_sheet',False,'Kaiju, MetaPhlAn and DIAMOND read taxonomy'),
 Stage('humann',('validate',),'humann','humann_nucleotide',False,'Read-level pathway abundance; reference bias applies'),
 Stage('markers_hmm',('bridge',),'hmmsearch','cnps_hmms',False,'Biogeochemical HMM marker screen'),
 Stage('benchmark_squeezemeta',('bridge',),'SqueezeMeta.pl',description='Separate SqueezeMeta benchmark, not merged into production contigs'),
 Stage('benchmark_atlas',('bridge',),'atlas',description='Separate metagenome-atlas benchmark, not merged into production contigs'),
 Stage('report',description='Comprehensive execution audit and HTML summary')
]
_report_at=next(i for i,s in enumerate(STAGES) if s.id=='report')
STAGES.insert(_report_at, Stage('taxonomy_full', tuple(s.id for s in STAGES if s.id!='report'), 'kraken2', '', False,
    'Full PlusPFP Kraken2, memory-mapped, only after every other stage'))

MODULE={s:s for s in ['taxonomy','ssu','euk_contig','eggnog','dbcan','dram2','fegenie','amrfinder','antismash','metabolic_c','checkv',
 'whokaryote','nonpareil','kofam','catalogue','coverm','drep','eukcc','iphop','crispr','taxprofiler','humann','markers_hmm',
 'benchmark_squeezemeta','benchmark_atlas']}


def utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()

def path_signature(conf):
    params=Path(conf['sample_path']);inputs=[]
    if params.is_file():
        with params.open() as handle:
            data=list(csv.DictReader(handle,delimiter='\t'))
        for r in data:
            for k in ('r1','r2'):
                p=Path(r.get(k,''));inputs.append([str(p),p.stat().st_size,p.stat().st_mtime_ns] if p.is_file() else [str(p),'MISSING'])
    config={k:v for k,v in conf.items() if k not in ('_base','root','raw','asset_path','ext','sample_path','mag_review_path')}
    body=json.dumps([config,inputs],default=str,sort_keys=True)
    return hashlib.sha256(body.encode()).hexdigest()


def path_status(conf,stage):
    root=conf['root']
    if stage=='core':return (conf['raw']/'pipeline_info').is_dir() and (conf['raw']/'Assembly').is_dir()
    if stage=='bridge':return conf['asset_path'].is_file()
    if stage=='audit':return (root/'audit'/'audit.json').is_file()
    if stage=='mag_auto':return conf['mag_review_path'].is_file() or (conf['asset_path'].is_file() and not json.loads(conf['asset_path'].read_text()).get('mags'))
    if stage=='taxonomy_ani':return (root/'taxonomy_ani'/'ani_candidates.tsv').is_file()
    if stage=='report':return (root/'execution'/'report.json').is_file()
    if stage=='functional_summary':return (root/'summaries'/'eggnog_ko_gene_counts.tsv').is_file()
    if stage=='figures':return len(list((root/'figures').glob('*.svg')))==3
    if stage=='taxonomy':
        return bool(list((conf['ext']/stage).glob('*/report.txt')))
    if stage=='taxonomy_full':
        from scripts.validate_inputs import read_sheet
        rows=read_sheet(conf['sample_path'], False)
        return all((conf['ext']/stage/r['sample']/'report.txt').is_file() for r in rows)
    if stage=='bracken':return bool(list((conf['ext']/stage).glob('*.tsv')))
    if stage=='ssu':
        from scripts.validate_inputs import read_sheet
        rows=read_sheet(conf['sample_path'], False)
        folder=conf['ext']/stage
        return bool(rows) and all(any((folder/r['sample']).glob(r['sample']+'.phyloFlash*')) for r in rows)
    if stage=='euk_contig':
        tables=[p for p in (conf['ext']/stage).glob('*_tiara.tsv') if not p.name.startswith('log_')]
        return len(tables) == len(json.loads(conf['asset_path'].read_text()).get('samples',[]))
    if stage=='eggnog':return bool(list((conf['ext']/stage).glob('**/*.emapper.annotations')))
    if stage=='dbcan':return bool(list((conf['ext']/stage).glob('**/*overview*')))
    if stage=='amrfinder':return bool(list((conf['ext']/stage).glob('*.tsv')))
    if stage=='checkv':return bool(list((conf['ext']/stage).glob('**/quality_summary.tsv')))
    if stage=='dram2':return (conf['ext']/stage/'DRAM2').is_dir()
    if stage=='fegenie':return (conf['ext']/stage/'FeGenie').is_dir()
    if stage=='antismash':return bool(list((conf['ext']/stage).glob('**/index.html')))
    if stage=='metabolic_c':return (conf['ext']/stage/'METABOLIC_C').is_dir()
    if stage=='taxonomy_gtdb':return any((root/'taxonomy_gtdb').glob('**/*.summary.tsv'))
    if stage=='nonpareil':return _nonpareil_done(conf['ext']/stage, conf)
    if stage=='kofam':return bool(list((conf['ext']/stage).glob('*.tsv')))
    if stage=='catalogue':return bool(list((conf['ext']/stage).glob('*_cluster.tsv')))
    if stage=='coverm':return (conf['ext']/stage/'mag_abundance.tsv').is_file()
    if stage=='drep':return (conf['ext']/stage/'data_tables').is_dir() or (conf['ext']/stage/'NOT_ENOUGH_MAGS.txt').is_file()
    if stage=='eukcc':return bool(list((conf['ext']/stage).glob('**/eukcc.csv'))) or (conf['ext']/stage/'no_eukaryotic_bins.txt').is_file()
    if stage=='iphop':return (conf['ext']/stage).is_dir() and any((conf['ext']/stage).iterdir()) if (conf['ext']/stage).exists() else False
    if stage=='crispr':return bool(list((conf['ext']/stage).glob('*.gff')))
    if stage=='whokaryote':return bool(list((conf['ext']/stage).glob('*/whokaryote_predictions_*.tsv')))
    if stage=='taxprofiler':return (conf['ext']/stage/'results').is_dir()
    if stage=='humann':return bool(list((conf['ext']/stage).glob('*/*pathabundance*.tsv')))
    if stage=='markers_hmm':
        from scripts.validate_inputs import read_sheet
        rows=read_sheet(conf['sample_path'], False)
        return bool(rows) and all((conf['ext']/stage/(r['sample']+'.tbl')).is_file() for r in rows)
    if stage=='benchmark_squeezemeta':return (conf['ext']/stage/'project').exists()
    if stage=='benchmark_atlas':return (conf['ext']/stage).is_dir()
    if stage in MODULE:return False
    return True


def refresh_config(conf):
    """Pick up databases and tool paths that finish downloading while an earlier stage is still running."""
    path=Path(conf.get('_config_path',''))
    if not path.is_file():return
    fresh=json.loads(path.read_text())
    if isinstance(fresh.get('reference_databases'),dict):
        conf['reference_databases']=fresh['reference_databases']
    for key in ('metabolic_c_script','metaeuk_mmseqs_db','core_enable_genomad','core_enable_gunc','core_enable_metaeuk','resources','gtdbtk_accept_swap'):
        if key in fresh:conf[key]=fresh[key]


def missing_reqs(stage,conf):
    problems=[]
    if stage.executable and not shutil.which(stage.executable):problems.append('Executable missing: '+stage.executable)
    if stage.db:
        p=conf.get('reference_databases',{}).get(stage.db,'')
        if not p or not Path(p).exists():problems.append('Reference DB missing: reference_databases.'+stage.db)
    if stage.id=='core' and 'docker' in str(conf.get('profile','')) and not shutil.which('docker'):
        problems.append('Docker CLI unavailable for core profile=docker')
    if stage.id=='checkv' and not conf.get('viral_fasta') and not conf.get('core_enable_genomad'):
        problems.append('geNomad disabled in core and no viral_fasta given')
    if stage.id=='taxonomy_full':
        raw=(conf.get('reference_databases') or {}).get('kraken2_full') or ''
        index=Path(raw)/'hash.k2d' if raw else Path()
        if not index.is_file():
            problems.append('Full PlusPFP index is queued; classification waits until the other stages finish and the 270 GB index is extracted')
    if stage.id=='taxonomy_gtdb':
        from metax import gtdb_memory_problem
        problem=gtdb_memory_problem(conf)
        if problem:problems.append(problem)
    if stage.id=='metabolic_c':
        script=conf.get('metabolic_c_script','')
        if not script or not Path(script).is_file():problems.append('Set metabolic_c_script path to upstream METABOLIC-C.pl')
    return problems


def mag_state(conf):
    if not conf['asset_path'].exists():return None
    data=json.loads(conf['asset_path'].read_text())
    if not data.get('mags'):return 'NO_MAGS'
    if not conf['mag_review_path'].exists():return None
    from scripts.mag_ani import reviewed_mags
    try:return 'READY' if reviewed_mags(conf['mag_review_path']) else 'NO_QUALIFIED_MAGS'
    except ValueError as e:
        if 'No approved MAGs' in str(e):return 'NO_QUALIFIED_MAGS'
        raise


def execute_stage(s,conf,log):
    cmd=[sys.executable,str(BASE/'metax.py'),'--config',str(conf['_config_path'])]
    if s.id=='validate':cmd+=['validate','--full']
    elif s.id in ('core','bridge','audit'):cmd += [s.id]
    elif s.id=='mag_auto':
        from scripts.auto_mag_qc import auto_review
        outcome=auto_review(conf['asset_path'],conf['mag_review_path'])
        log.write_text(json.dumps(outcome,indent=2)+'\n')
        return outcome['status']
    elif s.id=='taxonomy_ani':cmd+=['taxonomy-ani']
    elif s.id=='taxonomy_gtdb':cmd+=['taxonomy-gtdb']
    elif s.id=='functional_summary':
        from scripts.functional_summary import summarize_eggnog
        summarize_eggnog(conf);log.write_text('Generated TSV from eggNOG annotations; gene counts are NOT functional abundance\n');return 'SUCCEEDED'
    elif s.id=='figures':
        from scripts.scientific_figures import write_svg
        write_svg(conf);log.write_text('Generated 3 descriptive assembly metric SVGs\n');return 'SUCCEEDED'
    elif s.id=='bracken':
        from scripts.downstream import run_bracken
        run_bracken(conf,log)
        return 'SUCCEEDED'
    elif s.id=='taxonomy_full':
        from scripts.queue_kraken_full import run_full_kraken
        run_full_kraken(conf, log)
        return 'SUCCEEDED'
    elif s.id in MODULE:cmd+=['module',s.id]
    else:raise ValueError('Invalid stage: '+s.id)
    env=os.environ.copy()
    if conf.get('resource_mode')=='mac64':env.setdefault('NXF_OPTS','-Xms1g -Xmx4g')
    elif conf.get('resource_mode')=='linux64':
        heap=int(conf.get('resources',{}).get('nxf_java_heap_gb',6))
        env.setdefault('NXF_OPTS',f'-Xms1g -Xmx{heap}g')
    with log.open('w') as handle:
        handle.write('$ '+shlex.join(cmd)+'\n');handle.flush()
        subprocess.run(cmd,stdout=handle,stderr=subprocess.STDOUT,check=True,env=env,cwd=BASE)
    return 'SUCCEEDED'


def write_json(path,obj):
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(obj,indent=2,default=str)+'\n')
    tmp.replace(path)


def summary(conf,records,signature):
    from scripts.downstream import write_report
    write_report(conf,records,signature)


def run_all(conf,*,dry_run=False,allow_partial=False,force=False):
    root=conf['root']/'execution'
    root.mkdir(parents=True,exist_ok=True)
    status_path=root/'status.json'
    signature=path_signature(conf)
    before=json.loads(status_path.read_text()) if status_path.exists() else {}
    rec=before.get('stages',{}) if before.get('signature')==signature and not force else {}
    # Do not retain any previously failed/skipped results as successful when re-running.
    rec={k:v for k,v in rec.items() if v.get('status') in ('SUCCEEDED','CACHED')}
    if dry_run:
        plan=[]
        for s in STAGES:
            problems=missing_reqs(s,conf) if s.id not in ('report',) else []
            plan.append({'stage':s.id,'dependencies':list(s.deps),'missing':problems,'may_require_mags':s.mag_required})
        print(json.dumps({'mode':conf.get('resource_mode'),'stages':plan,'note':'Dry run: no FASTQ, tool or biological stage executed'},indent=2))
        return plan
    for s in STAGES:
        if s.id=='report':continue
        log=root/(s.id+'.log')
        status='SUCCEEDED';reason='';started=utc()
        engine=s.id.split('_',1)[1] if s.id.startswith('benchmark_') else ''
        if engine and engine not in set(conf.get('benchmark_engines') or []):
            status='HELD';reason='Duplicate full engine stays out of the production tables. Add "'+engine+'" to benchmark_engines to run it in results/extensions/'+s.id
        dep_bad=[d for d in s.deps if rec.get(d,{}).get('status') not in ('SUCCEEDED','CACHED','HELD')]
        if status=='SUCCEEDED' and dep_bad:status='BLOCKED';reason='Upstream stages not complete: '+', '.join(dep_bad)
        if status=='SUCCEEDED' and s.mag_required:
            state=mag_state(conf)
            if state in ('NO_MAGS','NO_QUALIFIED_MAGS'):status='NOT_APPLICABLE';reason=state
            elif state is None:status='BLOCKED';reason='MAG review evidence unavailable'
        if status=='SUCCEEDED':
            refresh_config(conf)
            missing=missing_reqs(s,conf)
            if missing:status='BLOCKED';reason='; '.join(missing)
        old=rec.get(s.id,{})
        if status=='SUCCEEDED' and old.get('status') in ('SUCCEEDED','CACHED') and path_status(conf,s.id):
            status='CACHED';reason='Resume: same input configuration and expected output'
        elif status=='SUCCEEDED':
            try:
                result=execute_stage(s,conf,log)
                if result=='NO_MAGS':status='SUCCEEDED'
                if not path_status(conf,s.id):
                    raise RuntimeError('Stage command exited but expected outputs not found; inspect '+str(log))
            except Exception as e:
                status='FAILED';reason=type(e).__name__+': '+str(e)
                with log.open('a') as handle:handle.write('\nERROR '+reason+'\n'+traceback.format_exc())
        rec[s.id]={'status':status,'reason':reason,'started_utc':started,'finished_utc':utc(),'log':str(log)}
        print(f'[{status:<15}] {s.id}: {reason}',flush=True)
        write_json(status_path,{'signature':signature,'updated_utc':utc(),'stages':rec})
        if status=='FAILED':
            # Keep collecting transparent BLOCKED dependencies and generate final report, but don't run unsafe downstream.
            continue
    try:summary(conf,rec,signature)
    except Exception as e:
        print('Report error:',e,file=sys.stderr)
        rec['report']={'status':'FAILED','reason':str(e),'finished_utc':utc()}
    else:rec['report']={'status':'SUCCEEDED','finished_utc':utc()}
    write_json(status_path,{'signature':signature,'updated_utc':utc(),'stages':rec})
    unfinished={k:v for k,v in rec.items() if v['status'] not in ('SUCCEEDED','CACHED','HELD')}
    if unfinished and not allow_partial:raise RuntimeError(f'End-to-end INCOMPLETE: {len(unfinished)} stages not executed. Details: {root/"report.html"}. Use --allow-partial to explicitly accept incomplete coverage.')
    return rec
