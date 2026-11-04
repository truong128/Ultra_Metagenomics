#!/usr/bin/env python3
"""Reviewed MAG -> conservative skani ANI candidate screen; NEVER auto-assign GTDB species.

Only standard-library parsing; skani and reference sketches are external prerequisites.
"""
import csv
import json
import re
import shutil
import subprocess
from pathlib import Path

FIELDS = ['mag_id', 'fasta', 'approved', 'completeness_pct', 'contamination_pct',
          'gunc_status', 'reviewer_notes']


def mag_id(path):
    name = Path(path).name
    return re.sub(r'\.(fna|fasta|fa)(\.gz)?$', '', name, flags=re.I)


def make_review(data, destination):
    destination = Path(destination)
    if destination.exists():
        raise FileExistsError(f'Review file already exists, will not overwrite: {destination}')
    mags = data.get('mags') or []
    if not mags:
        raise ValueError('No candidate MAG FASTAs in assets.json; run bridge and inspect binners first')
    ids = [mag_id(x) for x in mags]
    if len(ids) != len(set(ids)):
        raise ValueError('MAG basenames not unique; curate stable IDs before review')
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open('w', newline='') as out:
        writer = csv.DictWriter(out, FIELDS, delimiter='\t')
        writer.writeheader()
        for path, name in zip(mags, ids):
            writer.writerow({'mag_id': name, 'fasta': str(Path(path).resolve()), 'approved': 'false',
                             'completeness_pct': '', 'contamination_pct': '',
                             'gunc_status': 'not_run', 'reviewer_notes': 'Fill CheckM2 QC before approval'})
    return destination


def reviewed_mags(review_file, require_files=True):
    review_file = Path(review_file)
    if not review_file.is_file():
        raise FileNotFoundError(f'Run "python3 metax.py mag-review" then curate {review_file}')
    with review_file.open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not set(FIELDS).issubset(reader.fieldnames or []):
            raise ValueError('MAG review TSV missing required columns')
        rows = list(reader)
    eligible = []
    ids = set()
    for row in rows:
        ident = row['mag_id'].strip()
        if not ident or ident in ids:
            raise ValueError('Missing or duplicate mag_id in review sheet')
        ids.add(ident)
        if row['approved'].strip().lower() not in ('true', 'yes', '1'):
            continue
        try:
            comp = float(row['completeness_pct'])
            contamination = float(row['contamination_pct'])
        except ValueError as exc:
            raise ValueError(f'{ident}: approved=true requires numeric CheckM2 completeness/contamination') from exc
        if not (50 <= comp <= 100 and 0 <= contamination <= 10):
            raise ValueError(f'{ident}: fails default MAG gate (completeness >=50%, contamination <=10%)')
        gunc = row['gunc_status'].strip().lower()
        if gunc not in ('pass', 'not_run'):
            raise ValueError(f'{ident}: GUNC status must be pass or not_run to enter ANI screen')
        p = Path(row['fasta'])
        if require_files and not p.is_file():
            raise FileNotFoundError(f'{ident}: MAG FASTA missing {p}')
        eligible.append({'id': ident, 'path': str(p.resolve()), 'gunc': gunc,
                         'completeness': comp, 'contamination': contamination})
    if not eligible:
        raise ValueError('No approved MAGs. Review CheckM2 and mark approved=true explicitly.')
    return eligible


def skani_command(mags, database, output, threads):
    return ['skani', 'search', *[m['path'] for m in mags], '-d', str(database),
            '-o', str(output), '-t', str(threads)]


def parse_search_file(search_file, eligible, report_file, min_ani=95, min_af=50):
    """skani outputs ANI and AF as percentages, e.g. ANI=98.7, AF=71.2.

    These are candidate matches; no automatic GTDB species-radius checks.
    """
    eligible_by_path = {str(Path(m['path']).resolve()): m for m in eligible}
    eligible_by_name = {}
    for m in eligible:
        name = Path(m['path']).name
        if name in eligible_by_name:
            raise ValueError('Duplicate MAG FASTA basenames; cannot disambiguate skani queries')
        eligible_by_name[name] = m
    hits = {m['id']: [] for m in eligible}
    with Path(search_file).open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        required = {'Ref_file', 'Query_file', 'ANI', 'Align_fraction_query'}
        ref_af_col = 'Align_fraction_ref' if 'Align_fraction_ref' in (reader.fieldnames or []) else 'Align_fraction_reference'
        if not reader.fieldnames or not required.issubset(reader.fieldnames) or ref_af_col not in reader.fieldnames:
            raise ValueError('Unsupported skani output; expected tab-separated header with ANI/AF columns')
        for item in reader:
            raw = item['Query_file']
            obj = eligible_by_path.get(str(Path(raw).resolve())) or eligible_by_name.get(Path(raw).name)
            if obj is None:
                raise ValueError('skani output contains a query outside approved MAG inputs: ' + raw)
            vals = [float(item[k]) for k in ('ANI',ref_af_col,'Align_fraction_query')]
            if any(v < 0 or v > 100 for v in vals):
                raise ValueError('ANI and AF must be percentage values from 0 to 100')
            hits[obj['id']].append((vals, item))
    report_file = Path(report_file)
    report_file.parent.mkdir(parents=True, exist_ok=True)
    out_fields = ['mag_id','status','reference_file','ANI_pct','AF_reference_pct','AF_query_pct',
                  'completeness_pct','contamination_pct','gunc_status','note']
    with report_file.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=out_fields, delimiter='\t')
        writer.writeheader()
        for mag in eligible:
            best = sorted(hits[mag['id']], key=lambda v:(v[0][0], min(v[0][1:])), reverse=True)
            accepted = [v for v in best if v[0][0] >= min_ani and v[0][1] >= min_af and v[0][2] >= min_af]
            chosen = accepted[0] if accepted else (best[0] if best else None)
            if chosen:
                (ani, ref_af, query_af), hit = chosen
                status = 'ANI_CANDIDATE_REVIEW' if accepted else 'NO_CONFIDENT_ANI_ASSIGNMENT'
                ref = hit['Ref_file']
            else:
                ani = ref_af = query_af = ''
                ref = ''
                status = 'NO_DATABASE_HIT'
            writer.writerow({'mag_id':mag['id'],'status':status,'reference_file':ref,
                             'ANI_pct':ani,'AF_reference_pct':ref_af,'AF_query_pct':query_af,
                             'completeness_pct':mag['completeness'],'contamination_pct':mag['contamination'],
                             'gunc_status':mag['gunc'],
                             'note':'Candidate ANI match only; not an official GTDB species assignment. GTDB uses species-specific circumscription radii and QC.' if accepted else 'Retain as unresolved; do not force species rank.'})
    return report_file


def run_ani(conf, dry=False):
    db = conf.get('reference_databases', {}).get('skani_representatives', '')
    if not db:
        raise ValueError('Set reference_databases.skani_representatives to an independent skani sketch DB (not an arbitrary GTDB-Tk internal DB)')
    if not dry and not Path(db).is_dir():
        raise FileNotFoundError('Missing skani database directory: '+db)
    selected = reviewed_mags(conf['mag_review_path'], not dry)
    dest = conf['root'] / 'taxonomy_ani'
    raw = dest / 'skani_search.tsv'
    annotated = dest / 'ani_candidates.tsv'
    cmd = skani_command(selected, db, raw, conf['threads'])
    import shlex
    print('+', shlex.join(map(str, cmd)))
    if dry:
        print('Next: parse skani_search.tsv into reviewed ani_candidates.tsv (no GTDB species claims)')
        return
    if not shutil.which('skani'):
        raise FileNotFoundError('skani not found. Install a matching skani release.')
    dest.mkdir(parents=True, exist_ok=True)
    subprocess.run(cmd, check=True)
    parse_search_file(raw, selected, annotated,
                      conf.get('ani',{}).get('min_ani_pct',95.0),
                      conf.get('ani',{}).get('min_af_pct',50.0))
    print('Wrote:', annotated)


def hpc_gtdb_command(conf, dry=False, out_name='taxonomy_gtdb_hpc', cpus=None, pplacer_cpus=1):
    """GTDB-Tk classify_wf. Never uses --full_tree or pplacer_useram.

    On a 64 GB host this still needs swap. It is not launched from the mac64 core.
    """
    selected = reviewed_mags(conf['mag_review_path'], not dry)
    out = conf['root'] / out_name
    batch = out / 'approved_mags_batch.tsv'
    out.mkdir(parents=True, exist_ok=True)
    workers = int(cpus if cpus is not None else min(8, int(conf['threads'])))
    with batch.open('w') as f:
        for m in selected:
            f.write(f"{m['path']}\t{m['id']}\n")
    return ['gtdbtk','classify_wf','--batchfile',str(batch),'--out_dir',str(out/'gtdbtk'),
            '--cpus',str(max(1, workers)),'--pplacer_cpus',str(max(1, int(pplacer_cpus))),
            '--scratch_dir',str(out/'scratch')]
