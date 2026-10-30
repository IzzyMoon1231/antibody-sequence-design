"""Real retained-input tests; mutated fixtures are controls, not observations."""
import importlib.util
import unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]

class AnalysisTests(unittest.TestCase):
    def test_real_table_matches_checked_arithmetic(self):
        self.assertIsNotNone(importlib.util.find_spec('analysis'), 'Missing portable analysis helper')
        import analysis
        import pandas as pd
        result = analysis.load(ROOT)
        table = analysis.summarize(result)[0].set_index('structure')
        expected = pd.read_csv(ROOT/'tests/expected_table.csv').set_index('structure')
        for col in expected.select_dtypes('number').columns:
            for case in expected.index:
                self.assertAlmostEqual(table.loc[case,col], expected.loc[case,col], places=12)

    def test_render_real_and_changed_inputs(self):
        import analysis
        import pandas as pd
        import tempfile, shutil
        self.assertTrue(hasattr(analysis, 'publish'), 'Missing regenerating reader and chart')
        with tempfile.TemporaryDirectory(dir=ROOT) as tmp:
            dest=Path(tmp)
            shutil.copytree(ROOT/'inputs',dest/'inputs')
            shutil.copyfile(ROOT/'selection.csv',dest/'selection.csv')
            data=analysis.load(dest)
            _, ablation=analysis.summarize(data)
            expected=pd.read_csv(ROOT/'tests/expected_ablation.csv')
            for col in expected.select_dtypes('number').columns:
                import numpy as np
                np.testing.assert_allclose(ablation[col],expected[col],atol=1e-12,equal_nan=True)
            first=analysis.publish(dest)
            before=(dest/'figures/figure1.png').read_bytes()
            before_text=(dest/'READ_RESULTS.md').read_text()
            for st in analysis.CASES:
                p=dest/f'inputs/{st}_ensemble.xlsx'
                import openpyxl
                wb=openpyxl.load_workbook(p); ws=wb['Overview']
                headers=[c.value for c in ws[1]]; col=headers.index('Heavy OASis Identity')+1
                for row in range(2,ws.max_row+1): ws.cell(row,col).value=0.1
                wb.save(p)
            rmsd=pd.read_csv(dest/'inputs/cdr_rmsd_results.csv')
            rmsd.loc[(rmsd.method=='ensemble')&(rmsd.cdr=='H3'),'rmsd_angstrom']=0.4
            rmsd.to_csv(dest/'inputs/cdr_rmsd_results.csv',index=False)
            a=pd.read_csv(dest/'inputs/ablation_scores.csv'); a['id']=a.id.str.replace('baseline','TEMP').str.replace('ensemble','baseline').str.replace('TEMP','ensemble'); a.to_csv(dest/'inputs/ablation_scores.csv',index=False)
            second=analysis.publish(dest)
            self.assertNotEqual(before,(dest/'figures/figure1.png').read_bytes())
            self.assertNotEqual(before_text,(dest/'READ_RESULTS.md').read_text())
            self.assertIn('lower heavy-domain',second['readout'])
            self.assertIn('lower H3',second['readout'])
            self.assertNotIn('clear adverse',second['readout'])
            self.assertIn('Changed-input analysis',(dest/'READ_RESULTS.md').read_text())
            self.assertIn('lower heavy-domain',(dest/'results.html').read_text())
            self.assertTrue(all(x<0 for x in second['bar_values']))
            self.assertEqual(first['bar_values'],[-x for x in second['bar_values']])
            (dest/'inputs/ablation_scores.csv').unlink()
            with self.assertRaisesRegex(ValueError,'Missing input: inputs/ablation_scores.csv'):
                analysis.publish(dest)

    def test_missing_inputs_invalidate_readers(self):
        import analysis, tempfile, shutil
        with tempfile.TemporaryDirectory(dir=ROOT) as tmp:
            dest=Path(tmp); shutil.copytree(ROOT/'inputs',dest/'inputs'); shutil.copyfile(ROOT/'selection.csv',dest/'selection.csv')
            analysis.publish(dest)
            (dest/'inputs/ablation_scores.csv').unlink()
            with self.assertRaisesRegex(ValueError,'Missing input'):
                analysis.publish(dest)
            self.assertFalse((dest/'figures/figure1.png').exists())
            self.assertFalse((dest/'tables/table1.csv').exists(), 'Stale tables must not survive failed update')
            self.assertIn('Update failed',(dest/'READ_RESULTS.md').read_text())

    def test_review_loop_coverage_and_regeneration(self):
        import analysis, tempfile, shutil, pandas as pd
        with tempfile.TemporaryDirectory(dir=ROOT) as tmp:
            dest=Path(tmp); shutil.copytree(ROOT/'inputs',dest/'inputs'); shutil.copyfile(ROOT/'selection.csv',dest/'selection.csv')
            result=analysis.publish(dest)
            self.assertIn('loops', result, 'Complete supporting loop table missing')
            loops=pd.read_csv(dest/'tables/loops.csv')
            self.assertEqual(len(loops),9)
            raw=pd.read_csv(dest/'inputs/cdr_rmsd_results.csv')
            for r in loops.itertuples():
                for method in ('baseline','ensemble'):
                    selected=raw[(raw.structure==r.structure)&(raw.cdr==r.cdr)&(raw.method==method)]
                    self.assertEqual(len(selected),100)
                    self.assertAlmostEqual(getattr(r,method+'_median_A'),selected.rmsd_angstrom.median(),places=12)
            for file in ('READ_RESULTS.md','results.html'):
                text=(dest/file).read_text()
                for phrase in ('H2, L2, L3','H1','not overall scFv preservation','H1/H2 + L1/L2/L3','not an isolated heavy-loop intervention','Paper changes','mean loop pLDDT'):
                    self.assertIn(phrase,text)
            self.assertIn('Heavy-domain', ' '.join(result['table'].columns))
            raw.loc[(raw.structure=='9NFU')&(raw.cdr=='L2')&(raw.method=='ensemble'),'rmsd_angstrom']=0.1
            raw.to_csv(dest/'inputs/cdr_rmsd_results.csv',index=False)
            changed=analysis.publish(dest)
            self.assertNotEqual(result['loop_readout'],changed['loop_readout'])
            self.assertNotIn('H2, L2, L3',changed['loop_readout'])
            self.assertIn('H1, L2',changed['loop_readout'])
            self.assertIn(changed['loop_readout'],(dest/'READ_RESULTS.md').read_text())

    def test_missing_supporting_loop_invalidates_output(self):
        import analysis, tempfile, shutil, pandas as pd
        with tempfile.TemporaryDirectory(dir=ROOT) as tmp:
            dest=Path(tmp); shutil.copytree(ROOT/'inputs',dest/'inputs'); shutil.copyfile(ROOT/'selection.csv',dest/'selection.csv')
            analysis.publish(dest)
            p=dest/'inputs/cdr_rmsd_results.csv'; raw=pd.read_csv(p)
            raw=raw.drop(raw[(raw.structure=='9NFU')&(raw.cdr=='L2')&(raw.method=='ensemble')].index[0])
            raw.to_csv(p,index=False)
            with self.assertRaisesRegex(ValueError,'Complete loops'):
                analysis.publish(dest)
            self.assertFalse((dest/'tables/loops.csv').exists())
            self.assertIn('Update failed',(dest/'READ_RESULTS.md').read_text())

    def test_notebook_builder_exists(self):
        self.assertTrue((ROOT/'build_notebook.py').exists(), 'Need reproducible clean-kernel notebook builder')

    def test_unsupported_structure_error_is_specific(self):
        import analysis, tempfile, shutil, pandas as pd
        with tempfile.TemporaryDirectory(dir=ROOT) as tmp:
            dest=Path(tmp); shutil.copytree(ROOT/'inputs',dest/'inputs'); shutil.copyfile(ROOT/'selection.csv',dest/'selection.csv')
            p=dest/'inputs/cdr_rmsd_results.csv'; d=pd.read_csv(p).drop(columns='cdr'); d.to_csv(p,index=False)
            with self.assertRaisesRegex(ValueError,'cdr_rmsd_results.csv: missing columns: cdr'):
                analysis.publish(dest)

    def test_case_interpretation_tracks_tradeoff(self):
        import analysis
        d=analysis.load(ROOT); t,a=analysis.summarize(d)
        answer=analysis.describe(d,t,a)[0]
        self.assertIn('does not support a blanket no-structural-cost claim',answer)
        t['h3_difference_medians_A']=-0.1
        answer=analysis.describe(d,t,a)[0]
        self.assertNotIn('does not support a blanket no-structural-cost claim',answer)

    def test_scientific_gates(self):
        import analysis, tempfile, shutil, pandas as pd
        changes=[
            ('atoms','Unsupported H3 endpoint'),
            ('missing_key','identical case/method/index'),
            ('nonfinite','Invalid rmsd_angstrom'),
            ('confidence','low_confidence must match'),
            ('duplicate_ablation','unique case_condition'),
            ('unpaired_ablation','100 complete finite reported seed pairs'),
            ('denominator','declared peptide denominator'),
            ('fasta','VH/FASTA mismatch'),
        ]
        for kind,message in changes:
            with self.subTest(kind=kind), tempfile.TemporaryDirectory(dir=ROOT) as tmp:
                dest=Path(tmp); shutil.copytree(ROOT/'inputs',dest/'inputs'); shutil.copyfile(ROOT/'selection.csv',dest/'selection.csv')
                p=dest/'inputs/cdr_rmsd_results.csv'; d=pd.read_csv(p)
                idx=d[(d.method=='ensemble')&(d.cdr=='H3')].index[0]
                if kind=='atoms': d.loc[idx,'n_atoms_per_residue']=1
                if kind=='missing_key': d=d.drop(idx)
                if kind=='nonfinite': d.loc[idx,'rmsd_angstrom']=float('nan')
                if kind=='confidence': d.loc[idx,'low_confidence']=1-d.loc[idx,'low_confidence']
                d.to_csv(p,index=False)
                p=dest/'inputs/ablation_scores.csv'; a=pd.read_csv(p)
                if kind=='duplicate_ablation': a=pd.concat([a,a.iloc[:1]])
                if kind=='unpaired_ablation': a=a.iloc[1:]
                if kind=='denominator': a.loc[0,'Heavy OASis Identity']=0.12345
                a.to_csv(p,index=False)
                if kind=='fasta': (dest/'inputs/9NH7_baseline.fa').write_text('>wrong\nAAAA\n')
                with self.assertRaisesRegex(ValueError,message):
                    analysis.publish(dest)

if __name__ == '__main__': unittest.main()
