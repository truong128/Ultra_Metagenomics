"""Post-core adapters: each module only adds new evidence, never repeats core assembly/QC/Prodi​gal/bins."""
import gzip,json,os,shutil,subprocess,tempfile
from pathlib import Path
from scripts.host_paths import conda_root, database_root

def tool_env(extra_path=None, **values):
    """Copy the current environment and put a tool's own bin directory first."""
    env=os.environ.copy()
    if extra_path:
        env['PATH']=str(extra_path)+os.pathsep+env.get('PATH','')
    env.update(values)
    return env
def execute(argv,dry,folder=None,env=None):
    import shlex
    print('+',shlex.join([str(x) for x in argv]),flush=True)
    if not dry:
        exe=shutil.which(str(argv[0]))
        if not exe:raise FileNotFoundError('Required executable: '+str(argv[0]))
        subprocess.run([str(x) for x in argv],check=True,cwd=folder,env=env)
def must(db,key,dry):
    p=db.get(key,'')
    if not p:raise ValueError('Set reference_databases.'+key+' in config/project.json')
    if not dry and not Path(p).exists():raise FileNotFoundError(p)
    return p
def products(conf,dry):
    p=conf['asset_path']
    if not p.exists():
        if dry:raise FileNotFoundError('Run bridge after core first. Dry-run modules need assets.json too.')
        raise FileNotFoundError(str(p))
    return json.loads(p.read_text())
def unpack(src,tmp):
    src=Path(src);dst=tmp/src.name.replace('.gz','')
    if str(src).endswith('.gz'):
        with gzip.open(src,'rb') as inp,open(dst,'wb') as out:shutil.copyfileobj(inp,out)
        return dst
    return src

def _dram2_env():
    """Nextflow 26 needs Java 17 through 26. Prefer the Ubuntu OpenJDK 21 binary when it is installed."""
    env=tool_env(NXF_SYNTAX_PARSER='v1')
    for candidate in (os.environ.get('JAVA_CMD', ''),
                      '/usr/lib/jvm/java-21-openjdk-amd64/bin/java',
                      '/usr/lib/jvm/java-17-openjdk-amd64/bin/java'):
        if candidate and Path(candidate).is_file():
            env['JAVA_CMD']=candidate
            break
    return env


def write_dram2_config(dest):
    """Nextflow config for the DRAM2 databases that are present. KEGG is omitted."""
    root = database_root()
    dest = Path(dest)
    dest.mkdir(parents=True, exist_ok=True)
    path = dest / 'dram2_databases.config'
    path.write_text(
        "// Local DRAM2 references. KEGG is omitted; it requires a license.\n"
        "params {\n"
        "    max_cpus   = 24\n"
        "    max_memory = '48.GB'\n"
        "    max_time   = '240.h'\n"
        f"    kofam_db   = '{root}/kofam/profiles'\n"
        f"    kofam_list = '{root}/kofam/ko_list'\n"
        f"    dbcan_db           = '{root}/dbcan'\n"
        "    dbcan_version      = 'db_v5-2-9_5-5-2026'\n"
        f"    dbcan_version_file = '{root}/dbcan/version.txt'\n"
        f"    merops_db = '{root}/dram2/merops'\n"
        f"    vog_db   = '{root}/dram2/vogdb/hmm'\n"
        f"    vog_list = '{root}/dram2/vogdb/vog_annotations_latest.tsv.gz'\n"
        f"    camper_hmm_db      = '{root}/dram2/camper/hmm'\n"
        f"    camper_hmm_list    = '{root}/dram2/camper/hmm/camper_hmm_scores.tsv'\n"
        f"    camper_mmseqs_db   = '{root}/dram2/camper/mmseqs'\n"
        f"    camper_mmseqs_list = '{root}/dram2/camper/mmseqs/camper_scores.tsv'\n"
        f"    fegenie_db   = '{root}/dram2/fegenie'\n"
        f"    fegenie_list = '{root}/dram2/fegenie/fegenie_iron_cut_offs.txt'\n"
        f"    uniref_db = '{root}/dram2/uniref'\n"
        "}\n"
        "docker { fixOwnership = true }\n"
        "process {\n"
        "    maxForks = 1\n"
        "    resourceLimits = [cpus: 24, memory: 48.GB, time: 240.h]\n"
        "}\n"
    )
    return path


def dram2_anno_dbs():
    """Databases installed for DRAM2. KEGG is omitted because it requires a license."""
    root=database_root()
    dram=root/'dram2'
    chosen=[]
    if (root/'kofam'/'profiles').is_dir() and (root/'kofam'/'ko_list').is_file():
        chosen.append('kofam')
    if (root/'dbcan'/'dbCAN.hmm').is_file() and (root/'dbcan'/'version.txt').is_file():
        chosen.append('dbcan')
    if (dram/'merops'/'merops.mmsdb').is_file():
        chosen.append('merops')
    if any((dram/'vogdb'/'hmm').glob('*.hmm')):
        chosen.append('vog')
    if (dram/'fegenie'/'fegenie.hmm').is_file():
        chosen.append('fegenie')
    if (dram/'camper'/'hmm'/'camper.hmm').is_file() and (dram/'camper'/'mmseqs'/'camper.mmsdb').is_file():
        chosen.append('camper')
    if (dram/'uniref'/'uniref.mmsdb').is_file():
        chosen.append('uniref')
    if not chosen:
        raise FileNotFoundError('DRAM2 reference databases are not installed')
    return chosen

def mag_folder(conf,out,dry):
    from scripts.mag_ani import reviewed_mags
    mags=[m['path'] for m in reviewed_mags(conf['mag_review_path'], require_files=not dry)]
    if not mags:raise ValueError('No approved MAGs in MAG review TSV.')
    p=out/'mags_input'
    if not dry:
        p.mkdir(parents=True,exist_ok=True)
        for m in mags:
            src=Path(m)
            if not src.exists():raise FileNotFoundError(src)
            if src.name.endswith('.gz'):
                link=p/src.name[:-3]
                if not link.exists():
                    with gzip.open(src,'rb') as reader,open(link,'wb') as writer:shutil.copyfileobj(reader,writer)
            else:
                link=p/src.name
                if not link.exists():link.symlink_to(src)
    return p

def run_module(conf,name,dry=False):
    data=products(conf,dry);db=conf['reference_databases'];threads=str(conf['threads'])
    dest=conf['ext']/name
    if not dry:dest.mkdir(parents=True,exist_ok=True)
    if name=='taxonomy':
        kraken=must(db,'kraken2',dry)
        # Input is raw reads for independent read-level profiling. Dedicated quality-controlled
        # taxonomy should use the saved fastp outputs with documented read identifiers/denominators.
        from scripts.validate_inputs import read_sheet
        for r in read_sheet(conf['sample_path'],not dry):
            output=dest/r['sample'];
            if not dry:output.mkdir(parents=True,exist_ok=True)
            execute(['kraken2','--db',kraken,'--paired']+(['--gzip-compressed'] if r['r1'].endswith('.gz') else [])+['--threads',threads,'--report',output/'report.txt',
                     '--output',output/'kraken.output',r['r1'],r['r2']],dry)
    elif name=='ssu':
        silva=must(db,'phyloflash_silva',dry)
        from scripts.validate_inputs import read_sheet
        for r in read_sheet(conf['sample_path'],not dry):
            output=dest/r['sample']
            if not dry:output.mkdir(parents=True,exist_ok=True)
            root=Path(__file__).resolve().parents[1]
            reads=[str(Path(r[k]) if Path(r[k]).is_absolute() else (root/r[k]).resolve()) for k in ('r1','r2')]
            # phyloFlash calls BBMap scripts such as readlength.sh by name. They live in
            # the helper environment, not on the default PATH.
            phylo_path=str(conda_root()/'envs'/'metax-phyloflash'/'bin')+os.pathsep+str(conda_root()/'envs'/'metax-phylo-tools'/'bin')
            execute(['phyloFlash.pl','-dbhome',silva,'-lib',r['sample'],'-CPUs',threads,'-read1',reads[0],'-read2',reads[1],'-almosteverything'],dry,folder=output,
                    env=tool_env(phylo_path))
    elif name=='euk_contig':
        for r in data['samples']:
            src=r['assembly'];output=dest/(r['sample']+'_tiara.tsv')
            if dry:print('+ tiara -i',src,'-o',output,'-t',threads,'--probabilities (requires decompression)')
            else:
                with tempfile.TemporaryDirectory(dir=dest) as td:
                    fa=unpack(src,Path(td))
                    execute(['tiara','-i',fa,'-o',output,'-t',threads,'--probabilities'],False)
    elif name=='eggnog':
        egg=must(db,'eggnog_data',dry)
        for r in data['samples']:
            out=dest/r['sample']
            if not dry:out.mkdir(parents=True,exist_ok=True)
            if dry:print('+ emapper.py -i',r['proteins'],'--itype proteins --data_dir',egg,'-o',r['sample'])
            else:
                with tempfile.TemporaryDirectory(dir=dest) as td:
                    fa=unpack(r['proteins'],Path(td))
                    execute(['emapper.py','-i',fa,'--itype','proteins','--data_dir',egg,
                             '-o',r['sample'],'--output_dir',out,'--cpu',threads],False)
    elif name=='dbcan':
        dbc=must(db,'dbcan',dry)
        for r in data['samples']:
            out=dest/r['sample']
            if not dry:out.mkdir(parents=True,exist_ok=True)
            if dry:print('+ run_dbcan CAZyme_annotation --input_raw_data',r['proteins'],'--db_dir',dbc,'--output_dir',out)
            else:
                with tempfile.TemporaryDirectory(dir=dest) as td:
                    fa=unpack(r['proteins'],Path(td))
                    execute(['run_dbcan','CAZyme_annotation','--input_raw_data',fa,'--mode','protein',
                             '--db_dir',dbc,'--output_dir',out,'--threads',threads],False)
    elif name=='dram2':
        mags=mag_folder(conf,dest,dry)
        anno=','.join(dram2_anno_dbs())
        execute(['nextflow','run','BortonWrightonLabs/DRAM','-r','dev',
                 '-c',write_dram2_config(dest),
                 '--input_fasta',mags,'--outdir',dest/'DRAM2',
                 '--anno_dbs',anno,
                 '--rename','--annotate','--qc','--summarize','--sum_ecos','aquatic,bgc',
                 '--visualize','--traits','-profile',conf['profile'],'-resume'],dry,
                env=_dram2_env())
    elif name=='checkv':
        chk=must(db,'checkv',dry)
        explicit=conf.get('viral_fasta','')
        if explicit:
            candidates=[Path(explicit)]
        else:
            root_v=Path(data['nfcore_root'])
            candidates=[]
            for pattern in ('VirusIdentification/geNomad/**/*_viruses.fna','VirusIdentification/geNomad/**/*_viruses.fna.gz','VirusIdentification/geNomad/**/*_virus.fna','VirusIdentification/geNomad/**/*_virus.fna.gz'):
                candidates+=sorted(root_v.glob(pattern))
        if not candidates:raise FileNotFoundError('No viral FASTAs. Activate core geNomad and provide its DB, or explicitly set viral_fasta.')
        for i,viral in enumerate(dict.fromkeys(candidates)):
            if not dry and not viral.is_file():raise FileNotFoundError(viral)
            output=dest/(viral.stem+'_'+str(i))
            execute(['checkv','end_to_end',viral,output,'-d',chk,'-t',threads],dry)
    elif name=='fegenie':
        mags=mag_folder(conf,dest,dry)
        fegenie_share=conda_root()/'envs'/'metax-fegenie'/'share'/'fegenie-1.2'
        # FeGenie overwrites HMMlib.txt from $iron_hmms and rscripts.txt from $rscripts.
        execute(['FeGenie.py','-bin_dir',mags,'-bin_ext','fa','-t',threads,'-out',dest/'FeGenie'],dry,folder=dest,
                env=tool_env(conda_root()/'envs'/'metax-fegenie'/'bin',
                             rscripts=str(fegenie_share/'rscripts'),
                             iron_hmms=str(fegenie_share/'hmms'/'iron')))
    elif name=='amrfinder':
        for r in data['samples']:
            output=dest/(r['sample']+'.tsv')
            if dry:print('+ amrfinder -p',r['proteins'],'-o',output)
            else:
                with tempfile.TemporaryDirectory(dir=dest) as td:
                    faa=unpack(r['proteins'],Path(td))
                    amr_db=database_root()/'amrfinder'
                    if not amr_db.is_dir() or not any(amr_db.iterdir()):
                        raise FileNotFoundError('AMRFinder database is not installed in '+str(amr_db))
                    execute(['amrfinder','-p',faa,'-o',output,'--threads',threads,'--database',amr_db],False)
    elif name=='metabolic_c':
        # METABOLIC-C legacy pathway; prepare uncompressed reads automatically (large disk use).
        from scripts.validate_inputs import read_sheet
        from scripts.mag_ani import reviewed_mags
        samples=read_sheet(conf['sample_path'],not dry)
        mags=[m['path'] for m in reviewed_mags(conf['mag_review_path'],require_files=not dry)]
        source=conf.get('metabolic_c_script','')
        if not source:raise ValueError('Set metabolic_c_script to METABOLIC-C.pl')
        mag_out=dest/'metabolic_genomes'
        fastq_out=dest/'uncompressed_reads'
        if not dry:
            mag_out.mkdir(parents=True,exist_ok=True)
            fastq_out.mkdir(parents=True,exist_ok=True)
            for m in mags:
                src=Path(m)
                if src.suffix=='.gz':
                    dst=mag_out/(src.stem+'.fa');
                    if not dst.exists():
                        with gzip.open(src,'rb') as r,open(dst,'wb') as w:shutil.copyfileobj(r,w)
                else:
                    dst=mag_out/(src.stem+'.fa')
                    if not dst.exists():dst.symlink_to(src)
        reads=[]
        for r in samples:
            pair=[]
            for k in ('r1','r2'):
                p=Path(r[k]);out=fastq_out/(r['sample']+'_'+k+'.fastq')
                if p.suffix=='.gz':
                    if not dry and not out.exists():
                        with gzip.open(p,'rb') as reader,open(out,'wb') as writer:shutil.copyfileobj(reader,writer)
                    pair.append(str(out))
                else:pair.append(str(p))
            reads.append(','.join(pair))
        sheet=dest/'metabolic_reads.txt'
        if not dry:sheet.write_text('\n'.join(reads)+'\n')
        execute(['perl',str(source),'-in-gn',mag_out,'-r',sheet,'-st','illumina',
                 '-rt','metaG','-t',threads,'-o',dest/'METABOLIC_C'],dry,
                folder=str(Path(source).resolve().parent),
                env=tool_env(conda_root()/'envs'/'metax'/'bin'))
    elif name=='antismash':
        from scripts.mag_ani import reviewed_mags
        mags=[m['path'] for m in reviewed_mags(conf['mag_review_path'],require_files=not dry)]
        for m in mags:
            src=Path(m)
            fasta=src
            if src.name.endswith('.gz'):
                fasta=dest/src.name[:-3]
                if not dry and not fasta.exists():
                    with gzip.open(src,'rb') as reader,open(fasta,'wb') as writer:shutil.copyfileobj(reader,writer)
            out=dest/fasta.stem
            execute(['antismash',str(fasta),'--output-dir',out,'--cpus',threads,'--genefinding-tool','prodigal'],dry)
    else:
        from scripts.integrations import RUNNERS, run_integration
        if name in RUNNERS:
            run_integration(conf, name, dry)
        else:
            raise ValueError('Unsupported module: '+name)
