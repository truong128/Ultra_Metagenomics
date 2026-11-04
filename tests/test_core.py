import gzip,json,sys,tempfile,unittest
from pathlib import Path
BASE=Path(__file__).resolve().parents[1];sys.path.insert(0,str(BASE))
from scripts.validate_inputs import validate_fastq_pair,read_sheet
from scripts.bridge_nfcore import build,fasta_stats
from scripts.advanced_modules import products
from metax import load,build_core

def write_fastq(p,reads):
    with gzip.open(p,'wt') as f:
        for name,seq in reads:f.write(f'@{name}\n{seq}\n+\n'+('I'*len(seq))+'\n')
def fa(p,items):
    with gzip.open(p,'wt') as f:
        for name,seq in items:f.write(f'>{name}\n{seq}\n')
class TestMmx(unittest.TestCase):
 def test_fastq_integrity(self):
    with tempfile.TemporaryDirectory() as tmp:
        a=Path(tmp)/'r1.gz';b=Path(tmp)/'r2.gz'
        write_fastq(a,[('id1/1','ATGC'),('id2/1','TTTT')]);write_fastq(b,[('id1/2','CGAT'),('id2/2','AAAA')])
        x=validate_fastq_pair(a,b,0);self.assertEqual(x['pairs_examined'],2);self.assertTrue(x['full_scan'])
        write_fastq(b,[('id1/2','CGAT'),('different/2','AAAA')])
        with self.assertRaises(ValueError):validate_fastq_pair(a,b,0)
 def test_bridge_counts_and_n50(self):
    with tempfile.TemporaryDirectory() as tmp:
        root=Path(tmp)/'nf';dest=Path(tmp)/'assets.json'
        (root/'Assembly/MEGAHIT').mkdir(parents=True)
        (root/'Annotation/Prodigal').mkdir(parents=True)
        fa(root/'Assembly/MEGAHIT/sample1.contigs.fa.gz',[('ctg1','A'*100),('ctg2','T'*80),('ctg3','G'*20)])
        fa(root/'Annotation/Prodigal/MEGAHIT-sample1.faa.gz',[('p1','M'*10),('p2','G'*12)])
        with gzip.open(root/'Annotation/Prodigal/MEGAHIT-sample1.gff.gz','wt') as f:f.write('##gff-version 3\n')
        x=build(root,dest,[{'sample':'sample1','is_negative_control':'false'}])
        self.assertEqual(x['samples'][0]['assembly_stats']['sequences'],3)
        self.assertEqual(x['samples'][0]['assembly_stats']['bases'],200)
        self.assertEqual(x['samples'][0]['assembly_stats']['n50'],100)
        self.assertEqual(x['samples'][0]['proteins_stats']['sequences'],2)
        self.assertTrue(dest.exists());self.assertEqual(len(x['mags']),0)
        self.assertTrue(x['warnings'])
 def test_core_no_duplicate_assembly(self):
    conf=load(BASE/'config/project.linux64.json');cmd=build_core(conf)
    text=' '.join(str(x) for x in cmd)
    self.assertIn('nf-core/mag',text)
    self.assertIn('--skip_spades',text)
    self.assertIn('--refine_bins_dastool',text)
    self.assertIn('--run_checkm2',text)
    self.assertNotIn('--skip_maxbin2',text) # DAS Tool requires sufficient bin algorithms
    self.assertNotIn('SqueezeMeta',text)
 def test_missing_precomputed_inputs_fail(self):
    with tempfile.TemporaryDirectory() as tmp:
        with self.assertRaises(FileNotFoundError):build(tmp,Path(tmp)/'x.json',[{'sample':'bad','is_negative_control':'false'}])
if __name__=='__main__':unittest.main()
