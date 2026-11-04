#!/usr/bin/env python3
import argparse,json,html
from pathlib import Path

def main():
 p=argparse.ArgumentParser();p.add_argument('--assets',required=True);p.add_argument('--out',required=True);a=p.parse_args()
 data=json.loads(Path(a.assets).read_text());out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
 rows=[]
 for x in data['samples']:
  s=x['assembly_stats'];t=x['proteins_stats'];rows.append({'sample':x['sample'],'contigs':s['sequences'],'total_bases':s['bases'],'N50':s['n50'],'protein_orfs':t['sequences'],'negative_control':x['is_negative_control']})
 result={'sample_metrics':rows,'mags_curated_fasta_count':len(data['mags']),'warnings':data['warnings'],'important':'Counts are technical descriptors. No claim about species richness, flux, activity, or ecological interaction.'}
 (out/'audit.json').write_text(json.dumps(result,indent=2)+'\n')
 lines=['<html><meta charset="utf-8"><title>MangroveMetaX audit</title><style>body{font-family:Arial,sans-serif;margin:4rem;max-width:950px}table{border-collapse:collapse}td,th{padding:10px;border:1px solid #ccc}th{background:#edf5f5}</style><h1>Metagenomics quantitative audit</h1><p>All values are calculated from real bridged FASTA files. No biodiversity or metabolic activity is inferred.</p><table><tr>']
 labels=['sample','contigs','total_bases','N50','protein_orfs','negative_control']
 lines+=['<th>'+html.escape(v)+'</th>' for v in labels];lines+=['</tr>']
 for row in rows:lines+=['<tr>']+['<td>'+html.escape(str(row[k]))+'</td>' for k in labels]+['</tr>']
 lines+=['</table><h2>Warnings</h2><p>'+html.escape('; '.join(result['warnings']) or 'None')+'</p></html>']
 (out/'audit.html').write_text('\n'.join(lines));print(out/'audit.json');print(out/'audit.html')
if __name__=='__main__':main()
