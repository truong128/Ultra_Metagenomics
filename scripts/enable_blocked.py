"""Install the tools and reference databases that still block scheduler stages.

Does not stop the running nf-core/mag job. GTDB-Tk pplacer stays blocked until
RAM plus swap is at least 96 GB; this host cannot create that swap without sudo.
"""
import json, os, shutil, subprocess, sys, threading, time
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))
from scripts.host_paths import conda_root, database_root

ROOT = database_root()
CONFIG = Path(__file__).resolve().parents[1] / 'config' / 'project.json'
PINS = ROOT / 'pins.json'
STATUS = ROOT / 'ENABLE_STATUS.txt'
BIN = Path.home() / '.local' / 'bin'
ENVS = conda_root() / 'envs'
PIN_LOCK = threading.Lock()


def log(text):
    line = time.strftime('%Y-%m-%d %H:%M:%S') + '  ' + text
    print(line, flush=True)
    STATUS.parent.mkdir(parents=True, exist_ok=True)
    with STATUS.open('a') as handle:
        handle.write(line + '\n')


def run(argv, cwd=None, check=True):
    log('$ ' + ' '.join(str(x) for x in argv))
    result = subprocess.run([str(x) for x in argv], cwd=cwd)
    if check and result.returncode != 0:
        raise RuntimeError(f'exit {result.returncode}: {argv[0]}')
    return result.returncode


def wget(url, dest):
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(1, 6):
        code = run(['wget', '-c', '-O', dest, url], check=False)
        if code == 0 and dest.is_file() and dest.stat().st_size > 100_000:
            log(f'downloaded {dest.name} ({dest.stat().st_size} bytes)')
            return dest
        log(f'wget attempt {attempt} failed for {dest.name}')
        time.sleep(20)
    raise RuntimeError('download failed: ' + url)


def bash(script):
    return run(['bash', '-lc', script])


def pin(values):
    """Merge database paths into project.json and remember them across later rewrites."""
    with PIN_LOCK:
        current = {}
        if PINS.is_file():
            try:
                current = json.loads(PINS.read_text())
            except json.JSONDecodeError:
                current = {}
        current.update(values)
        PINS.write_text(json.dumps(current, indent=2) + '\n')
        for _ in range(10):
            data = json.loads(CONFIG.read_text())
            dbs = data.setdefault('reference_databases', {})
            for key, value in current.items():
                if key == 'metabolic_c_script':
                    data[key] = value
                else:
                    dbs[key] = value
            payload = json.dumps(data, indent=2) + '\n'
            CONFIG.write_text(payload)
            time.sleep(0.3)
            check = json.loads(CONFIG.read_text())
            ok = True
            for key, value in current.items():
                got = check.get(key) if key == 'metabolic_c_script' else check.get('reference_databases', {}).get(key)
                ok = ok and got == value
            if ok:
                return
        raise RuntimeError('project.json did not keep the new database paths')


def repin():
    if PINS.is_file():
        try:
            pin(json.loads(PINS.read_text()))
        except Exception as exc:
            log('repin skipped: ' + str(exc))


def wait_proc(pid):
    log(f'waiting for installer pid {pid}')
    while Path(f'/proc/{pid}').exists():
        time.sleep(30)
    log(f'installer pid {pid} finished')


def link(env, name, dest_name=None):
    src = ENVS / env / 'bin' / name
    if not src.exists():
        log(f'missing binary {src}')
        return False
    BIN.mkdir(parents=True, exist_ok=True)
    target = BIN / (dest_name or name)
    # A symlink here points at the conda binary. Writing through it would
    # replace that binary with a wrapper that calls itself.
    if target.is_symlink() or target.exists():
        target.unlink()
    wrapper = '#!/bin/bash\nexec ' + json.dumps(str(src)) + ' "$@"\n'
    target.write_text(wrapper)
    target.chmod(0o755)
    log(f'linked {target.name}')
    return True


def mamba_env(name, packages):
    env = ENVS / name / 'bin'
    if env.is_dir() and any(env.iterdir()):
        log(name + ' already exists')
        return
    spec = ' '.join(packages)
    bash(f'source {conda_root()}/etc/profile.d/conda.sh && mamba create -y -n {name} -c conda-forge -c bioconda {spec}')


def extract(archive, dest):
    dest.mkdir(parents=True, exist_ok=True)
    log(f'extracting {archive.name}')
    run(['tar', '-xf', archive, '-C', dest])
    return dest


def one_child_dir(path):
    kids = [p for p in path.iterdir() if p.name not in ('.', '..')]
    dirs = [p for p in kids if p.is_dir()]
    if len(dirs) == 1 and len(kids) == 1:
        return dirs[0]
    return path


def step(name, func):
    try:
        log('START ' + name)
        func()
        log('OK ' + name)
    except Exception as exc:
        log(f'FAIL {name}: {exc}')


DOWNLOADS = {
    'checkv': ('https://portal.nersc.gov/CheckV/checkv-db-v1.5.tar.gz', ROOT / 'checkv' / 'checkv-db-v1.5.tar.gz'),
    'eukcc': ('http://ftp.ebi.ac.uk/pub/databases/metagenomics/eukcc/eukcc2_db_ver_1.2.tar.gz', ROOT / 'eukcc' / 'eukcc2_db_ver_1.2.tar.gz'),
    'iphop': ('https://portal.nersc.gov/cfs/m342/iphop/db/iPHoP.latest_rw.tar.gz', ROOT / 'iphop' / 'iPHoP.latest_rw.tar.gz'),
    'metaphlan': ('http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vJan26_CHOCOPhlAnSGB_202605.tar', ROOT / 'metaphlan' / 'mpa_vJan26_CHOCOPhlAnSGB_202605.tar'),
    'kaiju_ref': ('https://kaiju-idx.s3.eu-central-1.amazonaws.com/2026/kaiju_db_refseq_ref_2026-09-17.tgz', ROOT / 'kaiju' / 'kaiju_db_refseq_ref_2026-09-17.tgz'),
    'kaiju_fungi': ('https://kaiju-idx.s3.eu-central-1.amazonaws.com/2026/kaiju_db_fungi_2026-09-18.tgz', ROOT / 'kaiju' / 'kaiju_db_fungi_2026-09-18.tgz'),
    'swissprot': ('https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz', ROOT / 'diamond' / 'uniprot_sprot.fasta.gz'),
}


def start_downloads():
    procs = {}
    for name, (url, dest) in DOWNLOADS.items():
        dest.parent.mkdir(parents=True, exist_ok=True)
        log(f'starting {name}')
        procs[name] = subprocess.Popen(['wget', '-c', '-O', str(dest), url])
    return procs


def install_tools():
    wait_pid = os.environ.get('METAX_WAIT_PID', '').strip()
    if wait_pid and Path('/proc/' + wait_pid).exists():
        wait_proc(int(wait_pid))
    for name in ('skani', 'kraken2', 'bracken', 'checkv', 'dRep', 'minced', 'hmmsearch', 'nonpareil',
                 'exec_annotation', 'seqkit', 'emapper.py', 'download_eggnog_data.py', 'diamond', 'amrfinder'):
        link('metax', name)
    for name in ('tiara',):
        link('metax-tiara', name)
    link('metax-phyloflash', 'phyloFlash.pl')
    link('metax-phyloflash', 'phyloFlash_makedb.pl')
    link('metax-coverm', 'coverm')
    link('metax-whokaryote', 'whokaryote.py')
    extras = [
        ('metax-gtdbtk', ['gtdbtk'], ['gtdbtk']),
        ('metax-dbcan', ['dbcan'], ['run_dbcan']),
        ('metax-fegenie', ['fegenie'], ['FeGenie.py']),
        ('metax-antismash', ['antismash'], ['antismash', 'download-antismash-databases']),
        ('metax-eukcc', ['eukcc'], ['eukcc']),
        ('metax-iphop', ['iphop'], ['iphop']),
        ('metax-humann', ['humann'], ['humann', 'humann_databases']),
    ]
    for env, packages, binaries in extras:
        try:
            mamba_env(env, packages)
        except Exception as exc:
            log(f'FAIL mamba {env}: {exc}')
            continue
        for binary in binaries:
            link(env, binary)


def tool_databases():
    egg = ROOT / 'eggnog'
    if shutil.which('download_eggnog_data.py') and not (egg / 'eggnog.db').is_file():
        egg.mkdir(parents=True, exist_ok=True)
        run(['download_eggnog_data.py', '-y', '--data_dir', egg])
        pin({'eggnog_data': str(egg)})
    dbcan = ROOT / 'dbcan'
    if shutil.which('run_dbcan') and not any(dbcan.glob('**/*')):
        dbcan.mkdir(parents=True, exist_ok=True)
        run(['run_dbcan', 'database', '--db_dir', dbcan])
        pin({'dbcan': str(dbcan)})
    silva = ROOT / 'phyloflash_silva'
    if shutil.which('phyloFlash_makedb.pl') and not any(silva.glob('*')):
        silva.mkdir(parents=True, exist_ok=True)
        run(['phyloFlash_makedb.pl', '--remote'], cwd=silva)
        pin({'phyloflash_silva': str(silva)})
    if shutil.which('amrfinder'):
        run(['amrfinder', '-u'], check=False)
    human = ROOT / 'humann'
    if shutil.which('humann_databases'):
        human.mkdir(parents=True, exist_ok=True)
        run(['humann_databases', '--download', 'chocophlan', 'full', human / 'chocophlan', '--update-config', 'no'])
        run(['humann_databases', '--download', 'uniref', 'uniref90_diamond', human / 'uniref', '--update-config', 'no'])
        pin({'humann_nucleotide': str(human / 'chocophlan'), 'humann_protein': str(human / 'uniref')})
    if shutil.which('download-antismash-databases'):
        run(['download-antismash-databases'], check=False)


def metabolic():
    dest = ROOT / 'metabolic'
    if not (dest / 'METABOLIC-C.pl').is_file():
        dest.mkdir(parents=True, exist_ok=True)
        run(['git', 'clone', '--depth', '1', 'https://github.com/AnantharamanLab/METABOLIC.git', dest])
    os.chdir(dest)
    for archive in ('METABOLIC_hmm_db.tgz', 'METABOLIC_template_and_database.tgz', 'Motif.tgz', 'Accessory_scripts.tgz'):
        path = dest / archive
        if path.is_file():
            run(['tar', '-xf', path, '-C', dest])
    hmms = sorted(dest.rglob('*.hmm'))
    combined = dest / 'cnps_biogeochemical.hmm'
    if hmms and not combined.is_file():
        with combined.open('wb') as out:
            for hmm in hmms:
                out.write(hmm.read_bytes())
                if not hmm.read_bytes().endswith(b'\n'):
                    out.write(b'\n')
    script = dest / 'METABOLIC-C.pl'
    values = {}
    if script.is_file():
        values['metabolic_c_script'] = str(script)
    if combined.is_file():
        values['cnps_hmms'] = str(combined)
    if values:
        pin(values)


def diamond_db(fasta_gz):
    if not fasta_gz or not Path(fasta_gz).is_file():
        raise RuntimeError('SwissProt FASTA was not downloaded')
    if not shutil.which('diamond'):
        raise RuntimeError('diamond is not on PATH')
    out = ROOT / 'diamond' / 'swissprot'
    if not Path(str(out) + '.dmnd').is_file():
        run(['diamond', 'makedb', '--in', fasta_gz, '-d', out])
    dmnd = Path(str(out) + '.dmnd')
    if not dmnd.is_file():
        raise RuntimeError('diamond makedb did not write swissprot.dmnd')
    return dmnd


def finish_archives(paths):
    ready = {}
    if paths.get('checkv'):
        dest = extract(paths['checkv'], ROOT / 'checkv' / 'db')
        ready['checkv'] = str(one_child_dir(dest))
    if paths.get('eukcc'):
        dest = extract(paths['eukcc'], ROOT / 'eukcc' / 'db')
        ready['eukcc'] = str(one_child_dir(dest))
    if paths.get('metaphlan'):
        dest = extract(paths['metaphlan'], ROOT / 'metaphlan' / 'db')
        ready['metaphlan'] = str(dest)
    kaiju_rows = []
    for key, label in (('kaiju_ref', 'refseq_ref'), ('kaiju_fungi', 'fungi')):
        if paths.get(key):
            dest = extract(paths[key], ROOT / 'kaiju' / label)
            kaiju_rows.append((label, str(one_child_dir(dest))))
    if paths.get('iphop'):
        archive = paths['iphop']
        free = shutil.disk_usage(ROOT).free
        need = int(archive.stat().st_size * 1.2)
        if free < need:
            log(f'iphop extract needs about {need} free bytes and the disk has {free}')
        else:
            dest = extract(archive, ROOT / 'iphop' / 'db')
            ready['iphop'] = str(one_child_dir(dest))
            archive.unlink()
    if ready:
        pin({k: v for k, v in ready.items() if k in ('checkv', 'eukcc', 'iphop')})
    dmnd = None
    try:
        dmnd = diamond_db(paths.get('swissprot'))
    except Exception as exc:
        log('FAIL diamond: ' + str(exc))
    if kaiju_rows and ready.get('metaphlan') and dmnd:
        sheet = ROOT / 'taxprofiler' / 'databases.csv'
        sheet.parent.mkdir(parents=True, exist_ok=True)
        lines = ['tool,db_name,db_params,db_type,db_path']
        for label, path in kaiju_rows:
            lines.append(f'kaiju,{label},,short,{path}')
        lines.append(f'metaphlan,mpa_vJan26,,short,{ready["metaphlan"]}')
        lines.append(f'diamond,swissprot,,short,{dmnd}')
        sheet.write_text('\n'.join(lines) + '\n')
        pin({'taxprofiler_sheet': str(sheet)})


def skani_from_gtdb():
    root = ROOT / 'gtdb'
    while True:
        sketches = []
        if root.is_dir():
            for path in root.rglob('*'):
                if path.is_dir() and path.name.lower() in {'sketches', 'skani'}:
                    sketches.append(path)
                elif path.is_file() and 'skani' in path.name.lower():
                    sketches.append(path.parent)
        if sketches:
            pin({'skani_representatives': str(sketches[0])})
            log('skani database ' + str(sketches[0]))
            return
        try:
            out = subprocess.check_output(['ps', '-eo', 'cmd'], text=True, errors='replace')
        except subprocess.CalledProcessError:
            out = ''
        busy = 'gtdbtk_data.tar.gz' in out
        metadata = list(root.rglob('metadata')) if root.is_dir() else []
        if metadata and not busy:
            log('GTDB package is extracted and has no standalone skani sketch directory')
            return
        if not busy and not (root / 'gtdbtk_data.tar.gz').is_file():
            log('GTDB package is not on disk')
            return
        log('waiting for the GTDB package before looking for skani sketches')
        time.sleep(120)


def collect_finished(procs):
    ready = {}
    pending = dict(procs)
    while pending:
        time.sleep(20)
        for name, proc in list(pending.items()):
            if proc.poll() is None:
                continue
            dest = DOWNLOADS[name][1]
            if proc.returncode == 0 and dest.is_file() and dest.stat().st_size > 100_000:
                ready[name] = dest
                log(f'downloaded {name} ({dest.stat().st_size} bytes)')
            else:
                log(f'FAIL download {name}: exit {proc.returncode}')
            del pending[name]
        if pending:
            log('still downloading: ' + ', '.join(sorted(pending)))
    return ready


def main():
    os.chdir(ROOT)
    log('enabling blocked stages')
    procs = start_downloads()
    holder = {}
    def collect():
        holder['paths'] = collect_finished(procs)
    thread = threading.Thread(target=collect, daemon=True)
    thread.start()
    step('tools', install_tools)
    step('tool databases', tool_databases)
    step('METABOLIC and marker HMMs', metabolic)
    thread.join()
    step('extract and taxprofiler sheet', lambda: finish_archives(holder.get('paths', {})))
    step('skani sketches from GTDB', skani_from_gtdb)
    repin()
    log('enable pass finished')


if __name__ == '__main__':
    main()
