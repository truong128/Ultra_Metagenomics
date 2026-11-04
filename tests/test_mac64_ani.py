import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

BASE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(BASE))
from metax import load,build_core
from scripts.mag_ani import make_review,reviewed_mags,parse_search_file,skani_command


class TestMac64(unittest.TestCase):
    def test_mac64_limits_and_gtdb_disabled(self):
        cfg=load(BASE/'config/project_mac64.json')
        cmd=build_core(cfg)
        self.assertIn('--skip_gtdbtk',cmd)
        self.assertIn('--skip_semibin',cmd)
        self.assertIn('--max_memory',cmd)
        self.assertIn('40.GB',cmd)
        self.assertIn('-c',cmd)
        self.assertTrue(str(BASE/'profiles/mac64.config') in cmd)
        self.assertIn('--refine_bins_dastool',cmd)
        self.assertIn('--run_checkm2',cmd)
        self.assertIn('--skip_spades',cmd)
        self.assertNotIn('--gtdb_db',cmd)

    def test_unsafe_core_gtdb_refused(self):
        cfg=load(BASE/'config/project_mac64.json')
        cfg['core_enable_gtdbtk']=True
        with self.assertRaisesRegex(ValueError,'prohibited'):build_core(cfg)

    def test_review_gate(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t); mags=root/'mag.fa';mags.write_text('>c1\nATGC\n')
            review=root/'review.tsv'
            make_review({'mags':[str(mags)]},review)
            with self.assertRaisesRegex(ValueError,'No approved'):reviewed_mags(review)
            txt=review.read_text().replace('false\t\t\tnot_run','true\t90\t3\tpass')
            review.write_text(txt)
            mags_approved=reviewed_mags(review)
            self.assertEqual(len(mags_approved),1)
            self.assertEqual(mags_approved[0]['completeness'],90)
            with self.assertRaises(FileExistsError):make_review({'mags':[str(mags)]},review)

    def test_review_refuses_low_quality(self):
        with tempfile.TemporaryDirectory() as t:
            r=Path(t)/'review.tsv'; fasta=Path(t)/'m.fa';fasta.write_text('>c\nATGC\n')
            r.write_text('mag_id\tfasta\tapproved\tcompleteness_pct\tcontamination_pct\tgunc_status\treviewer_notes\n'
                         f'M01\t{fasta}\ttrue\t40\t2\tpass\t\n')
            with self.assertRaisesRegex(ValueError,'fails default MAG gate'):reviewed_mags(r)

    def test_ani_candidates_and_unclassified(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t)
            ms=[{'id':'A','path':str(root/'A.fa'),'gunc':'pass','completeness':95,'contamination':1},
                {'id':'B','path':str(root/'B.fa'),'gunc':'not_run','completeness':88,'contamination':2},
                {'id':'C','path':str(root/'C.fa'),'gunc':'pass','completeness':82,'contamination':1}]
            # skani versions use either Align_fraction_ref or Align_fraction_reference
            for col in ('Align_fraction_ref','Align_fraction_reference'):
                inp=root/f'skani_{col}.tsv';report=root/f'report_{col}.tsv'
                with inp.open('w') as f:
                    f.write(f'Ref_file\tQuery_file\tANI\t{col}\tAlign_fraction_query\tRef_name\tQuery_name\n')
                    f.write(f'ref1.fa\t{root}/A.fa\t97.1\t75\t82\tref1\tA\n')
                    f.write(f'ref2.fa\t{root}/B.fa\t98.0\t12\t98\tref2\tB\n')
                    f.write(f'ref3.fa\t{root}/B.fa\t94.0\t91\t90\tref3\tB\n')
                parse_search_file(inp,ms,report)
                with report.open() as f: rows={r['mag_id']:r for r in csv.DictReader(f,delimiter='\t')}
                self.assertEqual(rows['A']['status'],'ANI_CANDIDATE_REVIEW')
                self.assertEqual(rows['B']['status'],'NO_CONFIDENT_ANI_ASSIGNMENT')
                self.assertEqual(rows['C']['status'],'NO_DATABASE_HIT')
                self.assertIn('not an official GTDB',rows['A']['note'])

    def test_skani_command_without_full_gtdb(self):
        cmd=skani_command([{'id':'x','path':'/tmp/x.fa'}],'/tmp/sketch','/tmp/out.tsv',8)
        self.assertEqual(cmd[:2],['skani','search'])
        self.assertEqual(cmd[cmd.index('-d')+1],'/tmp/sketch')
        self.assertNotIn('gtdbtk',cmd)

if __name__=='__main__':unittest.main()
