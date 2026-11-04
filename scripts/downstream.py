"""Final status reporting and nonduplicative Bracken postprocessing."""
import csv,html,json,os,shlex,subprocess
from pathlib import Path


def run_bracken(conf,log):
    reports=sorted((conf['ext']/'taxonomy').glob('*/report.txt'))
    if not reports:raise FileNotFoundError('No Kraken2 reports for Bracken')
    db=conf['reference_databases']['kraken2']
    dest=conf['ext']/'bracken';dest.mkdir(parents=True,exist_ok=True)
    length=str(conf['reference_databases'].get('bracken_read_length',150))
    with log.open('w') as out:
        for report in reports:
            output=dest/(report.parent.name+'.tsv')
            cmd=['bracken','-d',db,'-i',str(report),'-o',str(output),'-r',length,'-l','S']
            out.write('$ '+shlex.join(cmd)+'\n');out.flush()
            subprocess.run(cmd,check=True,stdout=out,stderr=subprocess.STDOUT)
            if not output.exists():raise RuntimeError(f'Bracken output not produced: {output}')


def write_report(conf,records,signature):
    dest=conf['root']/'execution';dest.mkdir(parents=True,exist_ok=True)
    status={k:v['status'] for k,v in records.items()}
    counts={x:sum(z==x for z in status.values()) for x in sorted(set(status.values()))}
    data={'pipeline':'MangroveMetaX v3.1','mode':conf.get('resource_mode'),'input_signature':signature,
          'status_counts':counts,'stages':records,
          'disclaimer':'Results reflect computational evidence only; DNA does not establish metabolic activity, species richness or causal microbial interactions.',
          'interpretation':'SUCCEEDED = stage exited and expected outputs checked where available; CACHED = matched configuration; HELD = duplicate engine intentionally not mixed into production; BLOCKED = unmet tool/DB/dependency; NOT_APPLICABLE = no QC eligible MAG; FAILED = error.',
          'full_coverage':all(x in ('SUCCEEDED','CACHED','HELD') for x in status.values())}
    registry=Path(__file__).resolve().parents[1]/'docs'/'STEP_REGISTRY.tsv'
    if registry.exists():
        with registry.open() as handle: registry_rows=list(csv.DictReader(handle,delimiter='\t'))
        review_categories={'DESIGN_PROTOCOL','EXPERIMENTAL','CURATED_DOWNSTREAM','OPTIONAL_EXTENSION','PUBLICATION_GATE','POST_AND_CURATED','OPTIONAL_HPC','BENCHMARK_ONLY'}
        data['research_step_registry_total']=len(registry_rows)
        data['registry_not_fully_automated_count']=sum(r.get('STATUS') in review_categories for r in registry_rows)
        data['coverage_warning']='A finished execution DAG does not mean every scientific protocol, experimental, or curation step in the registry was done.'
    audit=conf['root']/'audit'/'audit.json'
    if audit.exists():
        data['assembly_audit']=json.loads(audit.read_text())
    (dest/'report.json').write_text(json.dumps(data,indent=2)+'\n')
    cols=['Stage','Status','Reason','Log']
    content=['<!DOCTYPE html><html><head><meta charset="utf-8"><title>MangroveMetaX Audit</title>',
             '<style>body{font-family:Arial,sans-serif;margin:35px auto;max-width:1200px;padding:15px}td,th{border:1px solid #ccc;padding:7px;text-align:left}table{border-collapse:collapse;width:100%}th{background:#e9eef4}small{color:#4c5869}.failed{color:#aa2433}.ok{color:#21633f}</style>',
             '</head><body><h1>MangroveMetaX v3.0 execution audit</h1>',
             f'<p>Mode: <b>{html.escape(str(data["mode"]))}</b> | Full coverage: <b>{data["full_coverage"]}</b></p>',
             '<p>Every missing reference/tool is listed; unexecuted work is never represented as completed.</p>',
             '<table><thead><tr>'+''.join('<th>'+x+'</th>' for x in cols)+'</tr></thead><tbody>']
    for stage,d in records.items():
        cells=[stage,d['status'],d.get('reason',''),d.get('log','')]
        content.append('<tr>'+''.join('<td>'+html.escape(str(x))+'</td>' for x in cells)+'</tr>')
    content.append('</tbody></table>')
    if 'assembly_audit' in data:
        content.append('<h2>Quantitative audit</h2><pre>'+html.escape(json.dumps(data['assembly_audit'],indent=2))+'</pre>')
    content.append('<h2>Scope limitation</h2><p>'+html.escape(data.get('coverage_warning',''))+'</p>')
    content.append('<h2>Interpretation boundary</h2><p>'+html.escape(data['disclaimer'])+'</p></body></html>')
    (dest/'report.html').write_text('\n'.join(content))
