"""Adapters for evidence that nf-core/mag does not already produce.

Each function adds one kind of evidence. None of them re-runs MEGAHIT, Prodigal,
read QC, or the production binners. Duplicate full engines (SqueezeMeta,
metagenome-atlas) are separate benchmark directories and stay out of the
production tables.
"""
import csv
import gzip
import shutil
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]


def _run(argv, dry, folder=None, env=None):
    from scripts.advanced_modules import execute
    execute(argv, dry, folder=folder, env=env)


def _must_db(conf, key, dry):
    from scripts.advanced_modules import must
    return must(conf['reference_databases'], key, dry)


def _products(conf, dry):
    from scripts.advanced_modules import products
    return products(conf, dry)


def _unpack(src, dest_dir):
    src = Path(src)
    dest = Path(dest_dir) / src.name.replace('.gz', '')
    if str(src).endswith('.gz'):
        with gzip.open(src, 'rb') as inp, open(dest, 'wb') as out:
            shutil.copyfileobj(inp, out)
        return dest
    return src


def _reads(conf, dry):
    from scripts.validate_inputs import read_sheet
    rows = read_sheet(conf['sample_path'], not dry)
    resolved = []
    for row in rows:
        item = dict(row)
        for key in ('r1', 'r2'):
            path = Path(row[key]).expanduser()
            item[key] = str(path if path.is_absolute() else (BASE / path).resolve())
        resolved.append(item)
    return resolved


def _note(path, text, dry):
    print(text)
    if not dry:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text + '\n')


def run_whokaryote(conf, dry):
    """Gene-structure contig domain calls. Model S does not re-run Tiara."""
    import tempfile
    data = _products(conf, dry)
    dest = conf['ext'] / 'whokaryote'
    if not dry:
        dest.mkdir(parents=True, exist_ok=True)
    threads = str(conf['threads'])
    for sample in data['samples']:
        out = dest / sample['sample']
        if dry:
            print('+ whokaryote.py --contigs', sample['assembly'], '--gff', sample['gff'], '--model S --outdir', out)
            continue
        out.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=dest) as tmp:
            fasta = _unpack(sample['assembly'], tmp)
            gff = _unpack(sample['gff'], tmp)
            _run(['whokaryote.py', '--contigs', fasta, '--gff', gff, '--outdir', out,
                  '--model', 'S', '--minsize', '5000'], False)
        hits = list(out.glob('whokaryote_predictions_*.tsv'))
        if not hits:
            raise FileNotFoundError(f'Whokaryote produced no prediction table for {sample["sample"]}')


def _nonpareil_prefix_ready(prefix):
    from scripts.autorun import _nonpareil_ready
    return _nonpareil_ready(prefix)


def run_nonpareil(conf, dry):
    dest = conf['ext'] / 'nonpareil'
    if not dry:
        dest.mkdir(parents=True, exist_ok=True)
    for row in _reads(conf, dry):
        prefix = dest / row['sample']
        _run(['nonpareil', '-s', row['r1'], '-T', 'kmer', '-f', 'fastq',
              '-t', str(min(8, int(conf['threads']))), '-o', prefix], dry)
        if not dry and not _nonpareil_prefix_ready(prefix):
            raise FileNotFoundError(f'Nonpareil output missing for {row["sample"]}')


def run_kofam(conf, dry):
    import tempfile
    db = Path(_must_db(conf, 'kofam', dry))
    data = _products(conf, dry)
    dest = conf['ext'] / 'kofam'
    if not dry:
        dest.mkdir(parents=True, exist_ok=True)
    profiles = db / 'profiles'
    ko_list = db / 'ko_list'
    if not dry and (not profiles.is_dir() or not ko_list.is_file()):
        raise FileNotFoundError('KOfam directory must contain profiles/ and ko_list')
    for sample in data['samples']:
        out = dest / (sample['sample'] + '.tsv')
        if dry:
            print('+ exec_annotation', sample['proteins'], '-o', out)
            continue
        with tempfile.TemporaryDirectory(dir=dest) as tmp:
            faa = _unpack(sample['proteins'], tmp)
            _run(['exec_annotation', '-o', out, '--cpu', str(conf['threads']),
                  '-k', ko_list, '--profile', profiles, '-f', 'detail-tsv', faa], False)


def run_catalogue(conf, dry):
    """Nonredundant protein catalogue. Does not replace the per-sample Prodigal calls."""
    import tempfile
    data = _products(conf, dry)
    dest = conf['ext'] / 'catalogue'
    if dry:
        print('+ mmseqs easy-linclust <core proteins>', dest)
        return
    dest.mkdir(parents=True, exist_ok=True)
    merged = dest / 'proteins.faa'
    with merged.open('wb') as out:
        for sample in data['samples']:
            src = Path(sample['proteins'])
            opener = gzip.open if str(src).endswith('.gz') else open
            with opener(src, 'rb') as inp:
                shutil.copyfileobj(inp, out)
    tmp = dest / 'tmp'
    tmp.mkdir(exist_ok=True)
    prefix = dest / 'catalogue'
    _run(['mmseqs', 'easy-linclust', merged, prefix, tmp,
          '--min-seq-id', '0.95', '-c', '0.90', '--cov-mode', '1',
          '--threads', str(conf['threads'])], False)
    if not Path(str(prefix) + '_cluster.tsv').is_file():
        raise FileNotFoundError('MMseqs2 cluster table was not written')
    shutil.rmtree(tmp, ignore_errors=True)


def _bams(conf):
    root = Path(conf['asset_path'] and json_root(conf))
    return sorted(root.glob('**/*.bam'))


def json_root(conf):
    import json
    return json.loads(conf['asset_path'].read_text())['nfcore_root']


def run_coverm(conf, dry):
    """Abundance per sample. MEGAHIT contig names repeat across samples, so one CoverM folder would reject them."""
    from scripts.mag_ani import reviewed_mags
    from scripts.validate_inputs import read_sheet
    dest = conf['ext'] / 'coverm'
    mags = reviewed_mags(conf['mag_review_path'], require_files=not dry)
    root = Path(json_root(conf))
    if not dry:
        dest.mkdir(parents=True, exist_ok=True)
    pieces = []
    for sample in read_sheet(conf['sample_path'], not dry):
        sid = sample['sample']
        chosen = [m for m in mags if sid in m['id']]
        if not chosen:
            continue
        folder = dest / (sid + '_mags')
        bams = [p for p in root.glob('**/*.bam') if p.name.endswith(sid + '-' + sid + '.bam')]
        if not dry:
            if not bams:
                raise FileNotFoundError('No self-mapping BAM for ' + sid)
            folder.mkdir(parents=True, exist_ok=True)
            for mag in chosen:
                src = Path(mag['path'])
                out_fa = folder / (src.name[:-3] if src.name.endswith('.gz') else src.name)
                if out_fa.exists():
                    continue
                if src.name.endswith('.gz'):
                    with gzip.open(src, 'rb') as reader, open(out_fa, 'wb') as writer:
                        shutil.copyfileobj(reader, writer)
                else:
                    out_fa.symlink_to(src)
        piece = dest / (sid + '.tsv')
        _run(['coverm', 'genome', '-d', folder, '-x', 'fa', '-b', *[str(b) for b in bams],
              '-o', piece, '-t', str(conf['threads']), '--min-covered-fraction', '0'], dry)
        pieces.append(piece)
    if dry or not pieces:
        return
    combined = dest / 'mag_abundance.tsv'
    header = None
    with combined.open('w') as out:
        for piece in pieces:
            lines = piece.read_text().splitlines()
            if not lines:
                continue
            if header is None:
                header = lines[0]
                out.write(header + '\n')
            for line in lines[1:]:
                out.write(line + '\n')


def run_drep(conf, dry):
    from scripts.mag_ani import reviewed_mags
    dest = conf['ext'] / 'drep'
    mags = [m['path'] for m in reviewed_mags(conf['mag_review_path'], require_files=not dry)]
    if len(mags) < 2:
        _note(dest / 'NOT_ENOUGH_MAGS.txt',
              'dRep needs at least two QC-threshold MAGs. Dereplication was not run. This is not a species count.',
              dry)
        return
    _run(['dRep', 'dereplicate', dest, '-g', *mags, '-p', str(conf['threads']),
          '-comp', '50', '-con', '10'], dry)


def _eukaryotic_bins(conf):
    root = Path(json_root(conf))
    tables = list(root.glob('**/*tiara*')) + list(root.glob('**/*domain*classification*'))
    euk = []
    for table in tables:
        if not table.is_file() or table.suffix.lower() not in {'.tsv', '.txt', '.csv'}:
            continue
        text = table.read_text(errors='replace')
        for line in text.splitlines():
            cells = [c.strip() for c in line.replace(',', '\t').split('\t')]
            if any(c.lower() in {'eukarya', 'eukaryote', 'eukaryotic'} or 'eukarya' in c.lower() for c in cells[1:]):
                euk.append(cells[0])
    return tables, euk


def run_eukcc(conf, dry):
    db = _must_db(conf, 'eukcc', dry)
    dest = conf['ext'] / 'eukcc'
    if dry:
        print('+ eukcc single --db', db)
        return
    tables, euk_ids = _eukaryotic_bins(conf)
    if not tables:
        raise FileNotFoundError('No Tiara/domain-classification tables under the nf-core output, so EukCC had no eukaryotic bin list.')
    if not euk_ids:
        _note(dest / 'no_eukaryotic_bins.txt',
              'Bin domain classification tables contained no eukaryotic bin. EukCC was not run. This is not a claim that eukaryotes are absent from the reads.',
              False)
        return
    from scripts.mag_ani import reviewed_mags
    mags = reviewed_mags(conf['mag_review_path'])
    chosen = [m for m in mags if any(token and token in Path(m['path']).name for token in euk_ids)]
    if not chosen:
        _note(dest / 'no_eukaryotic_bins.txt',
              'Eukaryotic domain calls did not match any QC-threshold MAG filename. EukCC was not run on prokaryotic MAGs.',
              False)
        return
    dest.mkdir(parents=True, exist_ok=True)
    for mag in chosen:
        out = dest / mag['id']
        _run(['eukcc', 'single', '--db', db, '--out', out, mag['path']], False)


def run_iphop(conf, dry):
    db = _must_db(conf, 'iphop', dry)
    data = _products(conf, dry)
    explicit = conf.get('viral_fasta', '')
    if explicit:
        viral = [Path(explicit)]
    else:
        root = Path(data['nfcore_root'])
        viral = sorted(root.glob('VirusIdentification/geNomad/**/*virus*.fna*'))
    if not viral:
        raise FileNotFoundError('No geNomad viral FASTA. iPHoP was not run and no host was inferred.')
    dest = conf['ext'] / 'iphop'
    if not dry:
        dest.mkdir(parents=True, exist_ok=True)
    _run(['iphop', 'predict', '--fa_file', viral[0], '--db_dir', db, '--out_dir', dest,
          '-t', str(min(8, int(conf['threads'])))], dry)


def run_crispr(conf, dry):
    from scripts.mag_ani import reviewed_mags
    dest = conf['ext'] / 'crispr'
    mags = reviewed_mags(conf['mag_review_path'], require_files=not dry)
    if not dry:
        dest.mkdir(parents=True, exist_ok=True)
    for mag in mags:
        src = Path(mag['path'])
        fasta = src
        if src.name.endswith('.gz'):
            fasta = dest / src.name[:-3]
            if not dry and not fasta.exists():
                with gzip.open(src, 'rb') as reader, open(fasta, 'wb') as writer:
                    shutil.copyfileobj(reader, writer)
        out = dest / (mag['id'] + '.gff')
        _run(['minced', '-gff', str(fasta), out], dry)
        if not dry and not out.is_file():
            out.write_text('##gff-version 3\n# no CRISPR arrays found\n')


def run_taxprofiler(conf, dry):
    db_sheet = _must_db(conf, 'taxprofiler_sheet', dry)
    dest = conf['ext'] / 'taxprofiler'
    sheet = dest / 'samplesheet.csv'
    if not dry:
        dest.mkdir(parents=True, exist_ok=True)
        with sheet.open('w', newline='') as handle:
            writer = csv.writer(handle)
            writer.writerow(['sample', 'run_accession', 'instrument_platform', 'fastq_1', 'fastq_2', 'fasta'])
            for row in _reads(conf, dry):
                writer.writerow([row['sample'], row['sample'], 'ILLUMINA', row['r1'], row['r2'], ''])
    version = str(conf.get('taxprofiler_version', '1.2.4'))
    cmd = ['nextflow', 'run', 'nf-core/taxprofiler', '-r', version,
           '-profile', conf['profile'], '--input', sheet, '--databases', db_sheet,
           '--outdir', dest / 'results', '-resume']
    if conf.get('resource_mode') == 'linux64':
        cmd += ['-c', str(BASE / 'profiles' / 'taxprofiler_linux.config')]
    import os
    env = os.environ.copy()
    # taxprofiler 1.2.4 references manifest inside its config. Nextflow 26's new
    # parser rejects that. The legacy parser accepts the published config.
    env['NXF_SYNTAX_PARSER'] = 'v1'
    _run(cmd, dry, env=env)


def run_humann(conf, dry):
    """Read-level pathway abundance. This is encoded potential, and the default HUMAnN nucleotide reference is human-microbiome biased."""
    db = _must_db(conf, 'humann_nucleotide', dry)
    dest = conf['ext'] / 'humann'
    if not dry:
        dest.mkdir(parents=True, exist_ok=True)
    for row in _reads(conf, dry):
        out = dest / row['sample']
        joined = dest / (row['sample'] + '.fastq.gz')
        if not dry and not joined.exists():
            with joined.open('wb') as handle:
                for key in ('r1', 'r2'):
                    with open(row[key], 'rb') as src:
                        shutil.copyfileobj(src, handle)
        cmd = ['humann', '--input', joined, '--output', out, '--threads', str(conf['threads']),
               '--nucleotide-database', db]
        protein = (conf.get('reference_databases') or {}).get('humann_protein') or ''
        if protein:
            cmd += ['--protein-database', protein]
        from scripts.host_paths import database_root
        mpa = database_root() / 'metaphlan' / 'db'
        cmd += ['--metaphlan-options', '--db_dir ' + str(mpa) + ' -x mpa_vJan26_CHOCOPhlAnSGB_202605']
        _run(cmd, dry)


def run_markers_hmm(conf, dry):
    """HMMER screen of the curated biogeochemical HMM file. A hit is not a complete pathway."""
    import tempfile
    hmm = _must_db(conf, 'cnps_hmms', dry)
    data = _products(conf, dry)
    dest = conf['ext'] / 'markers_hmm'
    if not dry:
        dest.mkdir(parents=True, exist_ok=True)
    for sample in data['samples']:
        out = dest / (sample['sample'] + '.tbl')
        if dry:
            print('+ hmmsearch --tblout', out, hmm, sample['proteins'])
            continue
        with tempfile.TemporaryDirectory(dir=dest) as tmp:
            faa = _unpack(sample['proteins'], tmp)
            _run(['hmmsearch', '--cpu', str(conf['threads']), '--tblout', out, hmm, faa], False)


def run_benchmark_squeezemeta(conf, dry):
    dest = conf['ext'] / 'benchmark_squeezemeta'
    sheet = dest / 'squeezemeta_samples.tsv'
    if not dry:
        dest.mkdir(parents=True, exist_ok=True)
        with sheet.open('w') as handle:
            for row in _reads(conf, False):
                handle.write(f"{row['sample']}\t{row['r1']}\t{row['r2']}\n")
    _run(['SqueezeMeta.pl', '-m', 'sequential', '-p', dest / 'project', '-s', sheet,
          '-t', str(conf['threads'])], dry)


def run_benchmark_atlas(conf, dry):
    dest = conf['ext'] / 'benchmark_atlas'
    if not dry:
        dest.mkdir(parents=True, exist_ok=True)
    _run(['atlas', 'run', '--working-dir', dest, '-j', str(conf['threads'])], dry)


RUNNERS = {
    'whokaryote': run_whokaryote,
    'nonpareil': run_nonpareil,
    'kofam': run_kofam,
    'catalogue': run_catalogue,
    'coverm': run_coverm,
    'drep': run_drep,
    'eukcc': run_eukcc,
    'iphop': run_iphop,
    'crispr': run_crispr,
    'taxprofiler': run_taxprofiler,
    'humann': run_humann,
    'markers_hmm': run_markers_hmm,
    'benchmark_squeezemeta': run_benchmark_squeezemeta,
    'benchmark_atlas': run_benchmark_atlas,
}


def run_integration(conf, name, dry=False):
    if name not in RUNNERS:
        raise ValueError('Unsupported integration: ' + name)
    RUNNERS[name](conf, dry)
