#!/usr/bin/env python3
"""Claim validator: orthogonal evidence, data provenance and language constraints."""
import argparse,csv,json,sys
from pathlib import Path
FIELD=['claim_id','sample','process','gene_id','contig_id','mag_id','taxonomy_rank','taxonomy_label','method','database_version','sequence_identity','alignment_coverage','bit_score','e_value','hmm_score','context_verified','phylogeny_verified','orthogonal_validation','evidence_source','claim_type','reviewer','notes']
TYPES={'dna_potential','rna_expression','protein_expression','biochemical_activity','experimental_flux'}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input');ap.add_argument('--template',action='store_true');ap.add_argument('--output',default='results/claim_qa.json');args=ap.parse_args()
 if args.template:
  target=Path(args.output).with_suffix('.tsv');target.parent.mkdir(parents=True,exist_ok=True)
  with open(target,'w') as f:f.write('\t'.join(FIELD)+'\n')
  print(target);return
 if not args.input:raise ValueError('Provide --input evidence ledger TSV or --template')
 with open(args.input,newline='') as f:
  reader=csv.DictReader(f,delimiter='\t');fields=set(reader.fieldnames or [])
  if set(FIELD)-fields:raise ValueError('Missing columns: '+str(sorted(set(FIELD)-fields)))
  rows=list(reader)
 qa=[]
 for i,r in enumerate(rows,2):
  reasons=[]
  if r['claim_type'] not in TYPES:reasons.append('unknown_claim_type')
  if not r['evidence_source']:reasons.append('no_evidence_source')
  if not r['gene_id'] and not r['mag_id']:reasons.append('no_sequence_or_MAG_id')
  if r['claim_type'] in ('rna_expression','protein_expression','biochemical_activity','experimental_flux') and not r['orthogonal_validation']:reasons.append('activity_requires_RNA_protein_assay_or_flux_data')
  if r['taxonomy_rank'].lower()=='species' and not r['phylogeny_verified']:reasons.append('species_requires_extra_validation')
  if r['context_verified'].lower() not in ('true','false'):reasons.append('context_verified_not_boolean')
  qa.append({'row':i,'claim_id':r['claim_id'],'status':'PASS' if not reasons else 'REVIEW','reasons':reasons})
 output=Path(args.output);output.parent.mkdir(parents=True,exist_ok=True)
 output.write_text(json.dumps({'claims':len(rows),'passing':sum(r['status']=='PASS' for r in qa),'results':qa},indent=2)+'\n')
 print(output)
 if any(r['status']=='REVIEW' for r in qa):sys.exit(3)
if __name__=='__main__':
 try:main()
 except (OSError,ValueError) as e:sys.exit('LEDGER ERROR: '+str(e))
