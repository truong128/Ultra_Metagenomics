import json,sys,tempfile,unittest
from pathlib import Path
BASE=Path(__file__).resolve().parents[1];sys.path.insert(0,str(BASE))
from scripts.functional_summary import summarize_eggnog
from scripts.scientific_figures import write_svg

class ReportsTest(unittest.TestCase):
 def test_ko_annotation_counts_and_svg(self):
  with tempfile.TemporaryDirectory() as tmp:
   r=Path(tmp);(r/'extensions/eggnog/A').mkdir(parents=True);(r/'audit').mkdir()
   (r/'assets.json').write_text(json.dumps({'samples':[{'sample':'A'}]}))
   (r/'extensions/eggnog/A/A.emapper.annotations').write_text('## emapper\n#query\tKEGG_ko\nG1\tko:K00001,ko:K00002\nG2\tko:K00001\nG3\t-\n')
   (r/'audit/audit.json').write_text(json.dumps({'sample_metrics':[{'sample':'A','contigs':2,'protein_orfs':3,'N50':1500}]}))
   conf={'root':r,'ext':r/'extensions','asset_path':r/'assets.json'}
   summarize_eggnog(conf);write_svg(conf)
   s=(r/'summaries/eggnog_ko_gene_counts.tsv').read_text()
   self.assertIn('A\tK00001\t2',s);self.assertIn('A\tK00002\t1',s)
   self.assertEqual(len(list((r/'figures').glob('*.svg'))),3)
if __name__=='__main__':unittest.main()
