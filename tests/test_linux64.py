import sys
import unittest
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE))
from metax import load, build_core


class TestLinux64(unittest.TestCase):
    def test_full_binning_and_gtdb_outside_nextflow(self):
        cfg = load(BASE / 'config/project.linux64.json')
        self.assertEqual(cfg['resource_mode'], 'linux64')
        cmd = build_core(cfg)
        self.assertIn('--skip_spades', cmd)
        self.assertIn('--skip_gtdbtk', cmd)
        self.assertNotIn('--gtdb_db', cmd)
        self.assertNotIn('--skip_semibin', cmd)
        self.assertNotIn('--skip_concoct', cmd)
        self.assertNotIn('--skip_comebin', cmd)
        self.assertNotIn('--skip_metabinner', cmd)
        self.assertIn('--semibin_environment', cmd)
        self.assertIn('global', cmd)
        self.assertIn('--run_busco', cmd)
        self.assertIn('--run_checkm2', cmd)
        self.assertIn('--run_gunc', cmd)
        self.assertIn('--run_virus_identification', cmd)
        self.assertIn('--metaeuk_mmseqs_db', cmd)
        self.assertIn('UniRef50', cmd)
        self.assertIn('--refine_bins_dastool', cmd)
        self.assertIn('48.GB', cmd)
        self.assertIn(str(BASE / 'profiles' / 'linux64.config'), cmd)
        self.assertNotIn('genomad', cfg['_core_omitted'])
        self.assertNotIn('metaeuk', cfg['_core_omitted'])
        self.assertNotIn('SqueezeMeta', ' '.join(cmd))

    def test_optional_database_is_attached_when_present(self):
        cfg = load(BASE / 'config/project.linux64.json')
        cfg['reference_databases'] = dict(cfg['reference_databases'])
        cfg['reference_databases']['genomad'] = str(BASE / 'config')
        cmd = build_core(cfg)
        self.assertIn('--run_virus_identification', cmd)
        self.assertIn('--genomad_db', cmd)
        self.assertIn('--genomad_splits', cmd)
        self.assertNotIn('genomad', cfg['_core_omitted'])


if __name__ == '__main__':
    unittest.main()
