"""Automated, conservative CheckM2 threshold selection; NOT human-curated MAG validation.

Refuses ambiguous mappings, duplicate QC rows, missing QC, or invalid numbers.
"""
import csv, json, math, re
from pathlib import Path
from scripts.mag_ani import FIELDS, mag_id


def norm(s):
    return re.sub(r'\.(fa|fna|fasta)(\.gz)?$', '', str(s).strip(), flags=re.I).split('/')[-1]


def choose_col(columns, terms, kind):
    if kind=='MAG identifier':
        canonical=[x for x in columns if x.lower().strip().replace(' ','_') in ('bin_id','bin','name')]
        if len(canonical)==1:return canonical[0]
    hits=[x for x in columns if all(t in x.lower().replace(' ','_') for t in terms)]
    if len(hits)==1:return hits[0]
    # nf-core/mag 5.5.0 publishes CheckM2's own header: Completeness, Contamination.
    exact=[x for x in columns if x.lower().replace(' ','_')==terms[0]]
    if not hits and len(exact)==1:return exact[0]
    raise ValueError(f'Cannot unambiguously identify {kind}: {columns}, candidates {hits or exact}')


def auto_review(assets, out, min_completeness=50, max_contamination=10):
    """Accept a deterministic exact-ID join to nf-core/mag CheckM2 summary; never fuzzy match."""
    assets = json.loads(Path(assets).read_text()) if not isinstance(assets,dict) else assets
    out=Path(out)
    if out.exists():return {'status':'EXISTS_KEEP_USER_FILE','path':str(out)}
    mags=assets.get('mags',[])
    if not mags:return {'status':'NO_MAGS','path':None}
    qc=Path(assets['nfcore_root'])/'GenomeBinning'/'QC'/'checkm2_summary.tsv'
    if not qc.is_file():raise FileNotFoundError(f'Cannot auto-select MAGs without CheckM2 summary: {qc}')
    with qc.open() as stream:
        reader=csv.DictReader(stream,delimiter='\t')
        fields=reader.fieldnames or []
        # Field name variations from nf-core/mag may have tool suffixes, never infer positional indices.
        key=choose_col(fields,['bin'],'MAG identifier')
        comp=choose_col(fields,['completeness','checkm2'],'CheckM2 completeness')
        contam=choose_col(fields,['contamination','checkm2'],'CheckM2 contamination')
        rows=list(reader)
    idx={}
    for row in rows:
        k=norm(row[key])
        if k in idx:raise ValueError(f'Duplicate CheckM2 MAG ID: {k}')
        idx[k]=row
    rows_out=[]
    for path in mags:
        k=mag_id(path)
        record=idx.get(k)
        if record is None:raise ValueError(f'Cannot exactly map MAG {k} to CheckM2 summary. Fix identifiers, do not make up QC.')
        try:
            c=float(record[comp]);t=float(record[contam])
        except (ValueError,TypeError) as exc:raise ValueError(f'Missing numeric CheckM2 QC for {k}') from exc
        if not all(map(math.isfinite,[c,t])) or not(0<=c<=100 and 0<=t<=100):raise ValueError(f'Invalid CheckM2 QC for {k}')
        passed=c>=min_completeness and t<=max_contamination
        rows_out.append({'mag_id':k,'fasta':str(Path(path).resolve()),'approved':str(passed).lower(),
                        'completeness_pct':c,'contamination_pct':t,'gunc_status':'not_run',
                        'reviewer_notes':'AUTO_THRESHOLD_ONLY; NOT taxonomic/structural/manual validation; check chimerism and provenance'})
    out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,delimiter='\t',fieldnames=FIELDS);writer.writeheader();writer.writerows(rows_out)
    return {'status':'GENERATED_AUTO_THRESHOLD','path':str(out),'eligible':sum(x['approved']=='true' for x in rows_out),'candidates':len(rows_out)}
