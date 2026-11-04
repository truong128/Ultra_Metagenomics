"""Dependency-free descriptive SVG plots based only on measured assembly statistics."""
import json,html
from pathlib import Path

def write_svg(conf):
    audit=conf['root']/'audit'/'audit.json'
    if not audit.is_file():raise FileNotFoundError(audit)
    rows=json.loads(audit.read_text())['sample_metrics']
    if not rows:raise ValueError('No samples in audit')
    dest=conf['root']/'figures';dest.mkdir(parents=True,exist_ok=True)
    metrics=[('contigs','Assembled contig counts'),('protein_orfs','Predicted ORF counts'),('N50','Assembly N50 (bp)')]
    for key,caption in metrics:
        values=[int(row[key]) for row in rows];maximum=max(values) or 1
        width=760;height=105+54*len(rows)
        parts=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="{width}" height="{height}">',
          '<rect width="100%" height="100%" fill="#ffffff"/>',f'<text x="16" y="35" font-family="Arial" font-weight="bold" font-size="19" fill="#172431">{html.escape(caption)}</text>',
          '<text x="16" y="58" font-family="Arial" font-size="11" fill="#555">Measured from FASTA, not species richness or biological activity</text>']
        for i,(row,val) in enumerate(zip(rows,values)):
            y=80+54*i;length=420*val/maximum
            parts.extend([f'<text x="16" y="{y+17}" font-family="Arial" font-size="12" fill="#172431">{html.escape(str(row["sample"]))}</text>',
                          f'<rect x="215" y="{y}" width="{length:.2f}" height="22" rx="4" fill="#227d79"/>',
                          f'<text x="{235+length:.2f}" y="{y+16}" font-family="Arial" font-size="12" fill="#172431">{val:,}</text>'])
        parts.append('</svg>')
        (dest/(key+'.svg')).write_text('\n'.join(parts)+'\n')
    (dest/'FIGURE_NOTES.txt').write_text('Descriptive technical metrics only. Differential ecology, causal interaction, abundance normalization and environmental activity are not inferred.\n')
    return dest
