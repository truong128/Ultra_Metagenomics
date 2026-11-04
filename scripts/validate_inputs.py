#!/usr/bin/env python3
"""Safe stdlib paired-FASTQ validator; reports explicitly sampled vs full scan."""
import argparse,csv,gzip,hashlib,json,re,sys
from pathlib import Path

def opener(p): return gzip.open(p,'rt') if str(p).endswith('.gz') else open(p,'rt')
def norm(s):
    s=s.split()[0][1:] if s.startswith('@') else s.split()[0]
    return re.sub(r'([/.][12])$','',s)
def validate_fastq_pair(a,b,limit=10000):
    counts=0;bases=0
    with opener(a) as f1, opener(b) as f2:
        while True:
            lines1=[f1.readline() for _ in range(4)]
            lines2=[f2.readline() for _ in range(4)]
            if not any(lines1) and not any(lines2): break
            if any(not x for x in lines1+lines2): raise ValueError(f'Truncated or desynchronized FASTQ near pair {counts+1}')
            for ls in [lines1,lines2]:
                if not ls[0].startswith('@') or not ls[2].startswith('+'): raise ValueError('Invalid FASTQ header or +')
                if len(ls[1].rstrip('\n\r'))!=len(ls[3].rstrip('\n\r')): raise ValueError('Sequence/quality length mismatch')
            if norm(lines1[0])!=norm(lines2[0]): raise ValueError(f'Pair header mismatch {lines1[0].strip()} {lines2[0].strip()}')
            counts+=1;bases+=len(lines1[1].strip())+len(lines2[1].strip())
            if limit and counts>=limit: break
    if not counts: raise ValueError('Zero FASTQ pairs')
    return {'pairs_examined':counts,'bases_examined':bases,'full_scan':limit==0}
def read_sheet(p,check_files=True):
    with open(p,newline='') as f:
        reader=csv.DictReader(f,delimiter='\t')
        mandatory={'sample','group','r1','r2','kind','is_negative_control'}
        if not mandatory.issubset(reader.fieldnames or []):raise ValueError(f'Missing fields {sorted(mandatory-set(reader.fieldnames or []))}')
        rows=list(reader)
    if not rows:raise ValueError('Empty samplesheet')
    known=set()
    for row in rows:
        name=row['sample']
        if not re.fullmatch('[A-Za-z][A-Za-z0-9_-]*',name):raise ValueError('Unsafe sample ID '+repr(name))
        if name in known:raise ValueError('Duplicate sample ID '+name)
        known.add(name)
        if row['is_negative_control'].lower() not in ('true','false'):raise ValueError('is_negative_control must be true/false')
        if check_files:
            for typ in ('r1','r2'):
                file=Path(row[typ]).expanduser()
                if not file.is_file():raise FileNotFoundError(f"{name} {typ}: {file}")
    return rows

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--samples',required=True);ap.add_argument('--out');ap.add_argument('--limit',type=int,default=10000);ap.add_argument('--header-only',action='store_true');a=ap.parse_args()
    rows=read_sheet(a.samples,not a.header_only)
    summary={'samples':len(rows),'negative_controls':sum(x['is_negative_control']=='true' for x in rows),'records':[]}
    for r in rows:
        report={'sample':r['sample'],'kind':r['kind'],'negative_control':r['is_negative_control']=='true'}
        if not a.header_only:report.update(validate_fastq_pair(r['r1'],r['r2'],a.limit))
        summary['records'].append(report)
    result=json.dumps(summary,indent=2,ensure_ascii=False)
    if a.out:
        Path(a.out).parent.mkdir(parents=True,exist_ok=True);Path(a.out).write_text(result+'\n')
    print(result)
if __name__=='__main__':
    try:main()
    except (ValueError,OSError) as e:sys.exit('INPUT ERROR: '+str(e))
