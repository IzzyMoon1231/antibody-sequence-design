"""Artist-level observations; mutations below are synthetic controls."""
import unittest
from pathlib import Path
import numpy as np
import analysis
ROOT=Path(__file__).resolve().parents[1]

class FigureTests(unittest.TestCase):
    def test_main_draws_and_all_loop_quartiles(self):
        self.assertTrue(hasattr(analysis,'plot_main'), 'Main comparison figure missing')
        data=analysis.load(ROOT)
        fig=analysis.plot_main(data)
        for i,st in enumerate(analysis.CASES):
            ax=fig.axes[2*i]
            for j,m in enumerate(analysis.METHODS):
                expected=data['main'].query('structure==@st and method==@m').oasis.to_numpy()
                np.testing.assert_allclose(ax.collections[j*2].get_offsets()[:,0],expected)
            ax=fig.axes[2*i+1]
            for j,m in enumerate(analysis.METHODS):
                loops=['H1','H2','H3']+(['L1','L2','L3'] if st=='9NFU' else [])
                expected=[data['loops'].query('structure==@st and method==@m and cdr==@c').rmsd_angstrom.quantile([.25,.5,.75]).to_numpy() for c in loops]
                np.testing.assert_allclose(ax.collections[j*2+1].get_offsets()[:,0],np.array(expected)[:,1])
                segments=ax.collections[j*2].get_segments()
                np.testing.assert_allclose([[s[0,0],s[1,0]] for s in segments],np.array(expected)[:,[0,2]])
        data['main'].loc[data['main'].method=='ensemble','oasis']=.1
        changed=analysis.plot_main(data)
        np.testing.assert_allclose(changed.axes[0].collections[2].get_offsets()[:,0],.1)
        import matplotlib.pyplot as plt
        plt.close('all')

    def test_ablation_points_are_reported_pairs(self):
        data=analysis.load(ROOT)
        import inspect
        fig=(analysis.plot_ablation(analysis.summarize(data)[1]) if len(inspect.signature(analysis.plot_ablation).parameters)==1 else analysis.plot_ablation(analysis.summarize(data)[1], data['ablation']))
        for ax,st in zip(fig.axes,analysis.CASES):
            self.assertEqual(len(ax.patches),0, 'Replace mean-only bars with actual pairs')
            for j,cond in enumerate(('cdrh3','others')):
                w=data['ablation'].query('structure==@st and condition==@cond').pivot(index='seed',columns='method',values='Heavy OASis Identity')
                expected=100*(w.ensemble-w.baseline).to_numpy()
                np.testing.assert_allclose(ax.collections[j*2].get_offsets()[:,0],expected)
                self.assertAlmostEqual(ax.lines[j+1].get_xdata()[0],expected.mean())
        reversed_data=data['ablation'].copy()
        reversed_data['method']=reversed_data.method.map({'baseline':'ensemble','ensemble':'baseline'})
        rev=analysis.plot_ablation(analysis.summarize(data)[1], reversed_data)
        np.testing.assert_allclose(rev.axes[0].collections[0].get_offsets()[:,0],-fig.axes[0].collections[0].get_offsets()[:,0])
        import matplotlib.pyplot as plt
        plt.close('all')

    def test_reading_surfaces_prioritize_main_and_separate_checklist(self):
        import tempfile, shutil
        with tempfile.TemporaryDirectory(dir=ROOT) as tmp:
            dest=Path(tmp)
            shutil.copytree(ROOT/'inputs',dest/'inputs')
            shutil.copyfile(ROOT/'selection.csv',dest/'selection.csv')
            result=analysis.publish(dest)
            for name in ('figure1','figure2','experiment_map'):
                for ext in ('png','svg'):
                    self.assertTrue((dest/f'figures/{name}.{ext}').exists(), name)
            for name in ('READ_RESULTS.md','results.html'):
                text=(dest/name).read_text()
                self.assertLess(text.index('Figure 1'),text.index('Table 1'))
                self.assertIn('Figure 2',text)
                self.assertIn('AUTHOR_CHECKLIST.md',text)
                self.assertNotIn('## Paper changes',text)
                self.assertIn('not confidence intervals',text)
            self.assertIn('index-paired Wilcoxon',(dest/'AUTHOR_CHECKLIST.md').read_text())
            self.assertIn('Main comparison',(dest/'figures/experiment_map.svg').read_text())
            (dest/'inputs/ablation_scores.csv').unlink()
            with self.assertRaises(ValueError): analysis.publish(dest)
            for name in ('figure1','figure2','experiment_map'):
                self.assertTrue((dest/f'figures/{name}.png').exists(), 'Prior export preserved')

    def test_plain_opening_and_notebook_figure_hierarchy(self):
        d=analysis.load(ROOT); t,a=analysis.summarize(d)
        answer=analysis.describe(d,t,a)[0]
        self.assertTrue(answer.startswith('Ensemble sequences score higher'), answer)
        t['biophi_heavy_difference_medians']=-.2
        self.assertTrue(analysis.describe(d,t,a)[0].startswith('Ensemble sequences score lower'))
        source=(ROOT/'build_notebook.py').read_text()
        self.assertIn('analysis.plot_experiment_map(data)',source)
        self.assertIn('analysis.plot_ablation(ablation_table, ablation)',source)
        self.assertNotIn('Bars show',source)
        self.assertNotIn('## Paper changes',source)

if __name__=='__main__': unittest.main()
