"""Download the full PlusPFP Kraken2 index and classify reads only after the other stages finish.

The 16 GB capped index remains reference_databases.kraken2. This queue writes reference_databases.kraken2_full
and does not start classification while run-all or Nextflow is still running.
"""
import json, os, shutil, subprocess, sys, time
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
if str(BASE) not in sys.path:
    sys.path.insert(0, str(BASE))
from scripts.host_paths import database_root
URL = 'https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_20261022.tar.gz'
DB_ROOT = database_root() / 'kraken'
ARCHIVE = DB_ROOT / 'k2_pluspfp_20261022.tar.gz'
EXTRACT = DB_ROOT / 'pluspfp_full'
NOTE = DB_ROOT / 'full_pluspfp_QUEUE.txt'
DONE = {'SUCCEEDED', 'CACHED', 'HELD', 'NOT_APPLICABLE'}


def note(text):
    DB_ROOT.mkdir(parents=True, exist_ok=True)
    line = time.strftime('%Y-%m-%d %H:%M:%S') + '  ' + text
    print(line, flush=True)
    NOTE.write_text(line + '\n')


def pipeline_busy():
    try:
        out = subprocess.check_output(['ps', '-eo', 'cmd'], text=True, errors='replace')
    except (OSError, subprocess.CalledProcessError):
        return True
    for line in out.splitlines():
        if 'queue_kraken_full.py' in line:
            continue
        if 'metax.py run-all' in line or 'metax.py --config' in line:
            return True
        if 'nextflow' in line:
            return True
    return False


def others_done():
    status_path = BASE / 'results' / 'execution' / 'status.json'
    if not status_path.is_file():
        return False
    try:
        stages = json.loads(status_path.read_text()).get('stages', {})
    except json.JSONDecodeError:
        return False
    from scripts.autorun import STAGES
    needed = [s.id for s in STAGES if s.id not in ('taxonomy_full', 'report')]
    return all(stages.get(name, {}).get('status') in DONE for name in needed)


def download():
    DB_ROOT.mkdir(parents=True, exist_ok=True)
    note('downloading full PlusPFP archive; classification is not started')
    for attempt in range(1, 6):
        result = subprocess.run(['wget', '-c', '-O', str(ARCHIVE), URL])
        size = ARCHIVE.stat().st_size if ARCHIVE.is_file() else 0
        if result.returncode == 0 and size > 50_000_000_000:
            note(f'archive downloaded ({size} bytes); waiting for the other stages')
            return
        note(f'wget attempt {attempt} incomplete (exit {result.returncode}, {size} bytes); retrying')
        time.sleep(30)
    raise RuntimeError(f'Full PlusPFP archive is incomplete: {ARCHIVE}')


def publish(db):
    config = BASE / 'config' / 'project.json'
    for _ in range(8):
        data = json.loads(config.read_text())
        data.setdefault('reference_databases', {})['kraken2_full'] = str(db)
        payload = json.dumps(data, indent=2) + '\n'
        config.write_text(payload)
        time.sleep(0.3)
        check = json.loads(config.read_text()).get('reference_databases', {}).get('kraken2_full')
        if check == str(db):
            return
    raise RuntimeError('Could not keep reference_databases.kraken2_full in project.json')


def extract():
    found = next(EXTRACT.rglob('hash.k2d'), None) if EXTRACT.exists() else None
    if found:
        return found.parent
    EXTRACT.mkdir(parents=True, exist_ok=True)
    note('other stages are finished; extracting the 270 GB index')
    subprocess.run(['tar', '-xzf', str(ARCHIVE), '-C', str(EXTRACT)], check=True)
    found = next(EXTRACT.rglob('hash.k2d'), None)
    if not found:
        raise RuntimeError(f'hash.k2d not found under {EXTRACT}')
    publish(found.parent)
    if ARCHIVE.is_file():
        ARCHIVE.unlink()
        note('removed the archive after hash.k2d was extracted')
    return found.parent


def run_full_kraken(conf, log):
    """Classify each sample with the full index. Existing reports are kept."""
    import fcntl, shlex
    db = Path((conf.get('reference_databases') or {}).get('kraken2_full') or '')
    if not (db / 'hash.k2d').is_file():
        raise FileNotFoundError('Full PlusPFP index is not extracted (hash.k2d missing)')
    from scripts.validate_inputs import read_sheet
    rows = read_sheet(conf['sample_path'], False)
    dest = conf['ext'] / 'taxonomy_full'
    dest.mkdir(parents=True, exist_ok=True)
    with (dest / '.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        base = Path(conf['_base'])
        pending = [row for row in rows if not (dest / row['sample'] / 'report.txt').is_file()]
        if not pending:
            log.write_text('Full PlusPFP reports already present\n')
            return
        kraken = shutil.which('kraken2')
        if not kraken:
            raise FileNotFoundError('kraken2 is not on PATH')
        with log.open('a') as out:
            for row in pending:
                sample = dest / row['sample']
                sample.mkdir(parents=True, exist_ok=True)
                def resolve(name, item=row):
                    path = Path(item[name])
                    return path if path.is_absolute() else (base / path)
                r1, r2 = resolve('r1'), resolve('r2')
                cmd = [kraken, '--db', str(db), '--paired', '--memory-mapping', '--threads', '8',
                       '--report', str(sample / 'report.txt'), '--output', str(sample / 'kraken.output'),
                       str(r1), str(r2)]
                if str(r1).endswith('.gz'):
                    cmd.insert(4, '--gzip-compressed')
                wrapper = []
                if shutil.which('nice'):
                    wrapper += ['nice', '-n', '19']
                if shutil.which('ionice'):
                    wrapper += ['ionice', '-c', '3']
                out.write('$ ' + shlex.join(wrapper + cmd) + '\n')
                out.flush()
                subprocess.run(wrapper + cmd, check=True, stdout=out, stderr=subprocess.STDOUT)
                if not (sample / 'report.txt').is_file():
                    raise RuntimeError(f'Kraken2 did not write {sample / "report.txt"}')
        note('full PlusPFP classification finished')


def main():
    import sys
    os.chdir(BASE)
    if str(BASE) not in sys.path:
        sys.path.insert(0, str(BASE))
    download()
    while pipeline_busy() or not others_done():
        why = 'pipeline still running' if pipeline_busy() else 'waiting until every other stage is SUCCEEDED, CACHED, HELD, or NOT_APPLICABLE'
        note(why)
        time.sleep(120)
    import sys
    sys.path.insert(0, str(BASE))
    from metax import load
    conf = load(BASE / 'config' / 'project.json')
    db = extract()
    conf.setdefault('reference_databases', {})['kraken2_full'] = str(db)
    log = BASE / 'results' / 'execution' / 'taxonomy_full.log'
    log.parent.mkdir(parents=True, exist_ok=True)
    run_full_kraken(conf, log)


if __name__ == '__main__':
    main()
