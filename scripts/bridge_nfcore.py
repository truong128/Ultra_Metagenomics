#!/usr/bin/env python3
"""Strict output adapter for nf-core/mag v5.5.0; never invents or ambiguously chooses filenames."""
import argparse,csv,gzip,json,shutil,sys
from pathlib import Path

def fasta_stats(p):
    op=gzip.open if str(p).endswith('.gz') else open
    n=0;total=0;longest=0;length=0;seq=[]
    with op(p,'rt') as f:
        for line in f:
            if line.startswith('>'):
                if n:seq.append(length)
                n+=1;length=0
            else:length+=len(line.strip())
        if n:seq.append(length)
    total=sum(seq)
    for l in seq:longest=max(longest,l)
    seen=0;n50=0
    for v in sorted(seq,reverse=True):
        seen+=v
        if seen>=total/2:
            n50=v
            break
    return {'sequences':n,'bases':total,'max_length':longest,'n50':n50}

def single(path,what):
    if not path.is_file():raise FileNotFoundError(f'Missing {what}: {path}')
    return path.resolve()
def one_of(candidates,what):
    """Use the one published path that exists. Two existing candidates is an error, not a guess."""
    found=[p for p in candidates if p.is_file()]
    if len(found)==1:return found[0].resolve()
    if len(found)>1:raise FileNotFoundError('Ambiguous '+what+': '+', '.join(str(p) for p in found))
    raise FileNotFoundError('Missing '+what+': '+' or '.join(str(p) for p in candidates))
def build(root,out,rows,assembler='MEGAHIT',globs='GenomeBinning/DASTool/bins/*.fa'):
    root=Path(root).resolve();out=Path(out).resolve()
    if not root.is_dir():raise FileNotFoundError(root)
    items=[]
    for row in rows:
        sid=row['sample']
        # Does not automatically handle coassembly: sample->group mapping requires explicit workflow configuration.
        # nf-core/mag 5.5.0 publishes MEGAHIT-<sample>.contigs.fa.gz and a per-sample Prodigal directory.
        assembly_dir=root/'Assembly'/assembler
        prodigal=root/'Annotation'/'Prodigal'
        f=one_of([
            assembly_dir/f'{sid}.contigs.fa.gz',
            assembly_dir/f'{assembler}-{sid}.contigs.fa.gz',
        ],f'{assembler} assembly {sid}')
        aa=one_of([
            prodigal/f'{assembler}-{sid}.faa.gz',
            prodigal/assembler/sid/f'{assembler}-{sid}_prodigal.faa.gz',
        ],f'Prodigal proteins {sid}')
        gff=one_of([
            prodigal/f'{assembler}-{sid}.gff.gz',
            prodigal/assembler/sid/f'{assembler}-{sid}_prodigal.gff.gz',
        ],f'Prodigal GFF {sid}')
        stats=fasta_stats(f);prot=fasta_stats(aa)
        items.append({'sample':sid,'assembly':str(f),'proteins':str(aa),'gff':str(gff),'assembly_stats':stats,'proteins_stats':prot,'is_negative_control':row['is_negative_control'].lower()=='true'})
    patterns=[globs]
    if globs.endswith('*.fa'):patterns.append(globs+'.gz')
    bins=sorted({p.resolve() for pattern in patterns for p in root.glob(pattern)})
    multiqc=root/'multiqc'/'multiqc_report.html'
    report=root/'GenomeBinning'/'bin_summary.tsv'
    result={'source':'nf-core/mag','assembler':assembler,'nfcore_root':str(root),'samples':items,'mags':[str(x.resolve()) for x in bins], 'bin_summary':str(report.resolve()) if report.exists() else None,'multiqc':str(multiqc.resolve()) if multiqc.exists() else None,'warnings':[]}
    if bins:result['warnings'].append('DAS Tool refined bins are CANDIDATES, not automatically high-quality. Filter by CheckM2/GUNC/rRNA/tRNA and appropriate MIMAG criteria before publication.')
    if not bins:result['warnings'].append('No curated DAS Tool bins found. This is NOT the same as no MAGs: review alternate binner outputs and bin quality.')
    if not multiqc.exists():result['warnings'].append('MultiQC missing: core run may be incomplete.')
    out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,indent=2)+'\n');return result

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',required=True);p.add_argument('--out',required=True);p.add_argument('--samples',required=True);p.add_argument('--assembler',default='MEGAHIT');p.add_argument('--mag-glob',default='GenomeBinning/DASTool/bins/*.fa');a=p.parse_args()
    with open(a.samples) as f:rows=list(csv.DictReader(f,delimiter='\t'))
    data=build(a.root,a.out,rows,a.assembler,a.mag_glob)
    print(f"Bridged {len(data['samples'])} samples, {len(data['mags'])} curated MAGs")
    for w in data['warnings']:print('WARNING:',w,file=sys.stderr)
if __name__=='__main__':
    try:main()
    except (OSError,ValueError) as e:sys.exit('BRIDGE ERROR: '+str(e))
