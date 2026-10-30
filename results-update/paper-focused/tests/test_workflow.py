"""Student workflow controls. Mutations/subsets are SYNTHETIC, not research."""
from pathlib import Path
import unittest
import nbformat
ROOT = Path(__file__).resolve().parents[1]

class WorkflowTests(unittest.TestCase):
    def test_notebook_has_visible_read_compute_plot_export_steps(self):
        nb = nbformat.read(ROOT/'results.ipynb', as_version=4)
        sources = '\n'.join(c.source for c in nb.cells if c.cell_type == 'code')
        self.assertNotIn('publish(ROOT)', sources, 'Setup must not publish all outputs')
        for step in ('read_inputs(', 'validate_inputs(', '.groupby(', 'plot_panel(', 'plot_ablation(', 'EXPORT = False', 'export_bundle('):
            self.assertIn(step, sources)
        self.assertNotIn('Image(', sources)

    def test_changed_inspected_frames_cannot_inherit_original_provenance(self):
        import workflow
        loaded = workflow.read_inputs(ROOT)
        loaded['books'][('9NH7','ensemble')]['Overview']['Heavy OASis Identity'] = .1
        with self.assertRaisesRegex(ValueError, 'changed after reading'):
            workflow.validate_inputs(loaded)

    def test_fresh_cells_real_and_synthetic_without_saved_figures(self):
        import build_notebook, workflow, shutil, tempfile, pandas as pd, numpy as np
        from nbclient.exceptions import CellExecutionError
        import hashlib
        def snapshot(root):
            return {str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest()
                    for p in root.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
        with tempfile.TemporaryDirectory(prefix='.wf-',dir=ROOT) as tmp:
            dest=Path(tmp)
            for name in ('analysis.py','workflow.py','selection.csv'):
                shutil.copyfile(ROOT/name,dest/name)
            shutil.copytree(ROOT/'inputs',dest/'inputs')
            # No figures/tables copied: fresh kernel cannot view old PNGs.
            before=snapshot(dest)
            nb=build_notebook.execute(build_notebook.notebook(),dest)
            self.assertEqual(before,snapshot(dest), 'Run All with EXPORT=False wrote a file')
            images=lambda n:[o['data']['image/png'] for c in n.cells if c.cell_type=='code' for o in c.outputs if 'image/png' in o.get('data',{})]
            self.assertEqual(len(images(nb)),6)
            # SYNTHETIC subset/reversal control: not new research observations.
            loaded=workflow.read_inputs(dest)
            for row in loaded['selection'].itertuples():
                n=7 if row.method=='baseline' else 5
                sheets=loaded['books'][(row.structure,row.method)]
                with pd.ExcelWriter(dest/row.workbook) as writer:
                    for name,frame in sheets.items():
                        indices=frame.Antibody.str.extract(r'(\d+)(?:_VH)?$')[0].astype(int) - (row.method=='baseline')
                        frame=frame[indices<n].copy()
                        if name=='Overview' and row.method=='ensemble': frame['Heavy OASis Identity']=.1
                        frame.to_excel(writer,sheet_name=name,index=False)
                records=loaded['sequences'][(row.structure,row.method)][:n]
                (dest/row.fasta).write_text(''.join('>'+header+'\n'+seq+'\n' for header,seq in records))
            raw=loaded['structural']; idx=raw.design_id.str.extract(r'(\d+)$')[0].astype(int)
            raw=raw[raw.method.isin(['baseline','ensemble']) & (idx < raw.method.map({'baseline':7,'ensemble':5}).fillna(0))].copy()
            raw.loc[(raw.method=='ensemble')&(raw.cdr=='H3'),'rmsd_angstrom']=.4
            raw.to_csv(dest/'own-structure.csv',index=False)
            a=loaded['ablation']; a=a[a.id.str.extract(r'_s(\d+)$')[0].astype(int)<6].copy()
            a.id=a.id.str.replace('baseline','TEMP').str.replace('ensemble','baseline').str.replace('TEMP','ensemble')
            a.to_csv(dest/'own-ablation.csv',index=False)
            chosen=build_notebook.notebook()
            chosen.cells[2].source=chosen.cells[2].source.replace("'inputs/cdr_rmsd_results.csv'","'own-structure.csv'").replace("'inputs/ablation_scores.csv'","'own-ablation.csv'")
            probe=nbformat.v4.new_code_cell("import numpy as np\nassert not data['original_inputs']\nassert table.baseline_n.eq(7).all() and table.ensemble_n.eq(5).all()\nassert table.biophi_heavy_difference_medians.lt(0).all()\nassert 'lower heavy-domain' in readout\nfor i, st in enumerate(analysis.CASES):\n    for j, method in enumerate(analysis.METHODS):\n        expected = main.query('structure == @st and method == @method').oasis.to_numpy()\n        np.testing.assert_allclose(globals()['panel_'+st+'_similarity'].axes[0].collections[2*j].get_offsets()[:,0], expected)\n    for j, condition in enumerate(('cdrh3','others')):\n        expected = pairs.loc[(st,condition)].difference_pp.to_numpy()\n        np.testing.assert_allclose(ablation_figure.axes[i].collections[2*j].get_offsets()[:,0], expected)\nassert ablation_table.n_pairs.eq(6).all()\nassert ablation_table.query(\"domain == 'H'\").mean_difference.lt(0).all()")
            chosen.cells.append(probe)
            before=snapshot(dest)
            changed=build_notebook.execute(chosen,dest)
            self.assertEqual(before,snapshot(dest))
            self.assertNotEqual(images(nb)[0],images(changed)[0])
            self.assertNotEqual(images(nb)[5],images(changed)[5])
            optional=build_notebook.notebook()
            optional.cells[2].source=chosen.cells[2].source.replace("ABLATION = 'own-ablation.csv'", "ABLATION = None")
            before=snapshot(dest)
            optional_run=build_notebook.execute(optional,dest)
            self.assertEqual(len(images(optional_run)),5)
            self.assertEqual(before,snapshot(dest))
            # Exercise the actual last cell with explicit synthetic reasoning notes.
            for name in ('README.md','SOURCES.md','build_notebook.py'):
                shutil.copyfile(ROOT/name,dest/name)
            exporting=build_notebook.notebook()
            exporting.cells[2].source=chosen.cells[2].source
            for cell in exporting.cells:
                if cell.cell_type=='code' and cell.source.startswith('student_notes ='):
                    cell.source=cell.source.replace(": ''", ": 'SYNTHETIC workflow control, not a research interpretation'")
                if cell.cell_type=='code' and cell.source.startswith('EXPORT = False'):
                    cell.source=cell.source.replace('EXPORT = False','EXPORT = True')
            build_notebook.execute(exporting,dest)
            exported=list((dest/'exports').iterdir())
            self.assertEqual(len(exported),1)
            self.assertTrue((exported[0]/'student-notes.json').is_file())
            self.assertIn('lower heavy-domain',(exported[0]/'READ_RESULTS.md').read_text())
            self.assertIn('lower heavy-domain',(exported[0]/'results.html').read_text())
            self.assertNotIn('Pinned retained-input analysis',(exported[0]/'READ_RESULTS.md').read_text())
            for invalid in ('missing','negative'):
                failing=build_notebook.notebook()
                if invalid=='missing':
                    failing.cells[2].source=failing.cells[2].source.replace("'inputs/cdr_rmsd_results.csv'","'missing.csv'")
                else:
                    raw.loc[raw.index[0],'rmsd_angstrom']=-1
                    raw.to_csv(dest/'invalid.csv',index=False)
                    failing.cells[2].source=failing.cells[2].source.replace("'inputs/cdr_rmsd_results.csv'","'invalid.csv'")
                before=snapshot(dest)
                with self.assertRaises(CellExecutionError): build_notebook.execute(failing,dest)
                self.assertEqual(before,snapshot(dest))
                self.assertEqual(images(failing),[], 'Failure used stale images')

    def test_missing_preview_columns_fail_plainly_at_read(self):
        import workflow, tempfile, shutil, openpyxl
        with tempfile.TemporaryDirectory(prefix='.schema-',dir=ROOT) as tmp:
            dest=Path(tmp)
            shutil.copytree(ROOT/'inputs',dest/'inputs'); shutil.copyfile(ROOT/'selection.csv',dest/'selection.csv')
            p=dest/'inputs/9NH7_baseline.xlsx'
            wb=openpyxl.load_workbook(p); sheet=wb['Overview']
            for cell in sheet[1]:
                if cell.value=='Heavy OASis Identity': cell.value='wrong endpoint'
            wb.save(p)
            with self.assertRaisesRegex(ValueError,'9NH7_baseline.xlsx.*missing columns.*Heavy'):
                workflow.read_inputs(dest)

    def test_wrong_scoring_threshold_is_not_compatible(self):
        import workflow, tempfile, shutil, openpyxl
        with tempfile.TemporaryDirectory(prefix='.threshold-',dir=ROOT) as tmp:
            dest=Path(tmp)
            shutil.copytree(ROOT/'inputs',dest/'inputs'); shutil.copyfile(ROOT/'selection.csv',dest/'selection.csv')
            p=dest/'inputs/9NH7_baseline.xlsx'
            wb=openpyxl.load_workbook(p); sheet=wb['Overview']
            column=[c.value for c in sheet[1]].index('Threshold')+1
            sheet.cell(2,column).value='strict'; wb.save(p)
            with self.assertRaisesRegex(ValueError,'relaxed'):
                workflow.validate_inputs(workflow.read_inputs(dest))

    def test_local_paths_and_export_preservation(self):
        import workflow, tempfile, shutil, hashlib, pandas as pd
        with tempfile.TemporaryDirectory(prefix='.ex-',dir=ROOT) as tmp:
            dest=Path(tmp)
            shutil.copytree(ROOT/'inputs',dest/'inputs'); shutil.copyfile(ROOT/'selection.csv',dest/'selection.csv')
            for bad in ('../selection.csv','https://example.org/a.csv','/etc/passwd'):
                with self.subTest(bad=bad), self.assertRaises(ValueError): workflow.read_inputs(dest,selection=bad)
            (dest/'escape.csv').symlink_to(ROOT/'selection.csv')
            with self.assertRaises(ValueError): workflow.read_inputs(dest,selection='escape.csv')
            (dest/'escape.csv').unlink()
            data=workflow.validate_inputs(workflow.read_inputs(dest))
            self.assertTrue(data['original_inputs'])
            with self.assertRaisesRegex(ValueError,'opt-in'): workflow.export_bundle(dest,data)
            self.assertFalse((dest/'exports').exists())
            first=workflow.export_bundle(dest,data,confirm=True)
            frozen={str(p.relative_to(first)):hashlib.sha256(p.read_bytes()).hexdigest() for p in first.rglob('*') if p.is_file()}
            second=workflow.export_bundle(dest,data,confirm=True)
            self.assertNotEqual(first,second)
            self.assertEqual(frozen,{str(p.relative_to(first)):hashlib.sha256(p.read_bytes()).hexdigest() for p in first.rglob('*') if p.is_file()})
            (dest/'inputs/ablation_scores.csv').unlink()
            with self.assertRaisesRegex(ValueError,'Missing input'): workflow.read_inputs(dest)
            self.assertEqual(frozen,{str(p.relative_to(first)):hashlib.sha256(p.read_bytes()).hexdigest() for p in first.rglob('*') if p.is_file()})
            optional=workflow.validate_inputs(workflow.read_inputs(dest,ablation=None))
            self.assertTrue(optional['ablation'].empty)
            third=workflow.export_bundle(dest,optional,confirm=True)
            self.assertIn('not selected',(third/'figures/figure2.svg').read_text())

if __name__ == '__main__': unittest.main()
