"""Conservative descriptive gene annotation summary, not metabolic activity inference."""
import csv,collections,json,html,re
from pathlib import Path


def summarize_eggnog(conf):
    base=conf['ext']/'eggnog';dest=conf['root']/'summaries';dest.mkdir(parents=True,exist_ok=True)
    assets=json.loads(conf['asset_path'].read_text())
    totals=[];abundance=[]
    for sample in assets['samples']:
        sid=sample['sample']
        files=sorted((base/sid).glob('*.emapper.annotations'))
        if len(files)!=1:raise ValueError(f'Expected exactly 1 eggNOG annotation table for {sid}: {files}')
        genes=set();kohits=collections.defaultdict(set);has_ko=0
        with files[0].open(encoding='utf-8') as f:
            hdr=None
            for line in f:
                if line.startswith('##') or not line.strip():continue
                if line.startswith('#query'):
                    hdr=line.rstrip('\n').lstrip('#').split('\t')
                    if 'KEGG_ko' not in hdr:raise ValueError(f'KEGG_ko column missing in {files[0]}')
                    if 'query' not in hdr:raise ValueError(f'query ID column missing in {files[0]}')
                    continue
                if hdr is None:raise ValueError('eggNOG header not found: '+str(files[0]))
                values=line.rstrip('\n').split('\t')
                if len(values)!=len(hdr):raise ValueError('Malformed eggNOG row: '+str(files[0]))
                row=dict(zip(hdr,values));g=row['query'];genes.add(g)
                kos={x.removeprefix('ko:') for x in row['KEGG_ko'].split(',') if x not in ('','-')}
                kos={k for k in kos if re.fullmatch(r'K\d{5}',k)}
                if kos:has_ko+=1
                for ko in kos:kohits[ko].add(g)
        totals.append([sid,len(genes),has_ko,len(kohits)])
        for ko,identifiers in sorted(kohits.items()):abundance.append([sid,ko,len(identifiers)])
    for name,header,rows in [
        ('eggnog_ko_gene_counts.tsv',['sample','ko','unique_annotated_gene_count'],abundance),
        ('eggnog_summary.tsv',['sample','unique_annotated_genes','genes_with_kos','distinct_kos'],totals)]:
        with (dest/name).open('w',newline='') as output:
            w=csv.writer(output,delimiter='\t');w.writerow(header);w.writerows(rows)
    (dest/'INTERPRETATION.txt').write_text('KO counts are counts of annotated ORFs, not abundance, enzyme activity, metabolic flux, or pathway completeness. Multiple KO assignments can occur per gene.\n')
    return dest
