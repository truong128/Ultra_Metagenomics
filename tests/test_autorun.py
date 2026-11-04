import csv,json,sys,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
BASE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(BASE))
from scripts.auto_mag_qc import auto_review
from scripts.autorun import STAGES,run_all,missing_reqs
from metax import load

class TestAutoRun(unittest.TestCase):
 def test_auto_checkm2_threshold_and_exact_join(self):
  with tempfile.TemporaryDirectory() as td:
   base=Path(td)
   (base/'GenomeBinning/QC').mkdir(parents=True)
   (base/'GenomeBinning/QC/checkm2_summary.tsv').write_text('Bin Id\tCompleteness (CheckM2)\tContamination (CheckM2)\nMAG-A.fa\t93.2\t2.4\nMAG-B.fa\t38\t3\n')
   for k in ('MAG-A','MAG-B'):(base/(k+'.fa')).write_text('>x\nACGT\n')
   assets={'nfcore_root':str(base),'mags':[str(base/'MAG-A.fa'),str(base/'MAG-B.fa')]}
   sheet=base/'approved.tsv'
   result=auto_review(assets,sheet)
   self.assertEqual(result['eligible'],1)
   with sheet.open() as handle:x=list(csv.DictReader(handle,delimiter='\t'))
   self.assertEqual([v['approved'] for v in x],['true','false'])
   self.assertIn('AUTO_THRESHOLD_ONLY',x[0]['reviewer_notes'])
   self.assertEqual(auto_review(assets,sheet)['status'],'EXISTS_KEEP_USER_FILE')
 def test_auto_qc_refuses_ambiguous_columns(self):
  with tempfile.TemporaryDirectory() as td:
   base=Path(td);(base/'GenomeBinning/QC').mkdir(parents=True)
   (base/'GenomeBinning/QC/checkm2_summary.tsv').write_text('Bin Id\tCompleteness (CheckM2)\tCompleteness second checkm2\tContamination (CheckM2)\nA.fa\t90\t90\t1\n')
   (base/'A.fa').write_text('>x\nAC\n')
   with self.assertRaises(ValueError):auto_review({'nfcore_root':str(base),'mags':[str(base/'A.fa')]},base/'out.tsv')
 def test_scheduler_plan_and_missing_dbs(self):
  conf=load(BASE/'config/project.linux64.json')
  plan=run_all(conf,dry_run=True)
  self.assertEqual(len(plan),len(STAGES))
  self.assertTrue(any('Reference DB missing' in v for v in missing_reqs(next(s for s in STAGES if s.id=='taxonomy'),conf)))
  self.assertEqual(plan[0]['stage'],'validate')
 def test_scheduler_full_execution_and_resume_mock(self):
  with tempfile.TemporaryDirectory() as td:
   conf=load(BASE/'config/project.linux64.json')
   conf['root']=Path(td);conf['raw']=Path(td)/'nfcore_mag';conf['asset_path']=Path(td)/'assets.json';conf['ext']=Path(td)/'extensions';conf['mag_review_path']=Path(td)/'approved.tsv'
   from scripts import autorun
   with patch.object(autorun,'execute_stage',return_value='SUCCEEDED') as fake, patch.object(autorun,'path_status',return_value=True),patch.object(autorun,'missing_reqs',return_value=[]),patch.object(autorun,'mag_state',return_value='READY'):
    result=run_all(conf)
    self.assertTrue(all(v['status'] in ('SUCCEEDED','HELD') for v in result.values()))
    held=sum(1 for s in STAGES if s.id.startswith('benchmark_'))
    self.assertEqual(fake.call_count,len(STAGES)-1-held)
    second=run_all(conf)
    self.assertTrue(all(v['status'] in ('CACHED','SUCCEEDED','HELD') for v in second.values()))
    self.assertEqual(fake.call_count,len(STAGES)-1-held)
   self.assertTrue((Path(td)/'execution/report.json').exists())
 def test_failed_inputs_produce_truthful_audit(self):
  with tempfile.TemporaryDirectory() as td:
   sheet=Path(td)/'samples.tsv'
   sheet.write_text('sample\tgroup\tr1\tr2\tkind\tis_negative_control\nMissing\tG\t/no/such_R1.fastq.gz\t/no/such_R2.fastq.gz\twater\tfalse\n')
   cfg=json.loads((BASE/'config/project.linux64.json').read_text())
   cfg['samples_tsv']=str(sheet)
   proj=Path(td)/'project.json'
   proj.write_text(json.dumps(cfg))
   conf=load(proj)
   conf['root']=Path(td);conf['raw']=Path(td)/'nfcore_mag';conf['asset_path']=Path(td)/'assets.json';conf['ext']=Path(td)/'extensions';conf['mag_review_path']=Path(td)/'approved.tsv'
   with self.assertRaises(RuntimeError):run_all(conf)
   report=json.loads((Path(td)/'execution'/'report.json').read_text())
   self.assertFalse(report['full_coverage'])
   self.assertEqual(report['stages']['validate']['status'],'FAILED')
   self.assertEqual(report['stages']['core']['status'],'BLOCKED')
   self.assertTrue((Path(td)/'execution'/'report.html').is_file())

if __name__=='__main__':unittest.main()
