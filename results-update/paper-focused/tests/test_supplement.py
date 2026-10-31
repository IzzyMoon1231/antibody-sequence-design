"""Public-package integration tests; reduced budgets are test fixtures."""
import unittest
from pathlib import Path
import copy
import tempfile
import pandas as pd
import workflow

ROOT = Path(__file__).resolve().parents[1]

class SupplementTests(unittest.TestCase):
    def test_opt_in_and_original_tables(self):
        import supplemental as s
        d=workflow.validate_inputs(workflow.read_inputs(ROOT))
        self.assertIsNone(s.run(d))
        r=s.run(d, enabled=True, n_perm=99, n_boot=20)
        self.assertEqual(len(r['comparisons']),4)
        self.assertEqual(len(r['associations']),4)
        self.assertEqual(len(r['families']),22)
        self.assertEqual(r['receipt']['joined_rows'],400)
        self.assertNotIn('significant',r['families'])
        self.assertAlmostEqual(r['associations'].query("test_id=='9NH7:baseline:rho'").rho.iloc[0],.191613,places=5)

    def test_changed_missing_and_stale(self):
        import supplemental as s
        d=workflow.validate_inputs(workflow.read_inputs(ROOT))
        bad=copy.deepcopy(d); bad['main'].loc[0,'oasis']=0
        with self.assertRaises(ValueError): s.run(bad,enabled=True,n_perm=9,n_boot=2)
        with tempfile.TemporaryDirectory(dir=ROOT) as td:
            import shutil
            t=Path(td)
            shutil.copytree(ROOT/'inputs',t/'inputs'); shutil.copyfile(ROOT/'selection.csv',t/'selection.csv')
            f=t/'inputs/cdr_rmsd_results.csv'; raw=pd.read_csv(f)
            raw.loc[raw.cdr=='H3','rmsd_angstrom']+=1
            raw.to_csv(f,index=False)
            changed=workflow.validate_inputs(workflow.read_inputs(t))
            r=s.run(changed,enabled=True,n_perm=9,n_boot=2)
            self.assertTrue(r['associations'].rho.isna().all())
            self.assertTrue(r['families'].candidate_adjusted_p.isna().all())
            self.assertIn('actual predicted-coordinate',r['receipt']['qualification'])
            self.assertGreater(r['comparisons'].ensemble_median.max(),1)
            raw.iloc[1:].to_csv(f,index=False)
            with self.assertRaises(ValueError): workflow.validate_inputs(workflow.read_inputs(t))

    def test_export_and_notebook_controls(self):
        import supplemental as s
        import build_notebook
        nb=build_notebook.notebook()
        text='\n'.join(c.source for c in nb.cells)
        self.assertIn('RUN_EXPLORATORY_STATS = False',text)
        self.assertLess(text.index('Explain the result yourself'),text.index('RUN_EXPLORATORY_STATS'))
        d=workflow.validate_inputs(workflow.read_inputs(ROOT))
        r=s.run(d,enabled=True,n_perm=9,n_boot=2)
        with tempfile.TemporaryDirectory(dir=ROOT) as td:
            with self.assertRaises(ValueError): s.export(Path(td),d,r,review_plan='review only')
            p=s.export(Path(td),d,r,confirm=True,review_plan='review only')
            self.assertTrue((p/'receipt.json').is_file())
            q=s.export(Path(td),d,r,confirm=True,review_plan='review only')
            self.assertNotEqual(p,q)
            r['comparisons'].loc[0,'p']=.5
            with self.assertRaises(ValueError): s.export(Path(td),d,r,confirm=True,review_plan='review only')

if __name__=='__main__': unittest.main()
