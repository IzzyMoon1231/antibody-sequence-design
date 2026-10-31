"""Bounded novice orientation and experiment-map regression controls."""
import unittest
from pathlib import Path
import analysis
import build_notebook
ROOT = Path(__file__).resolve().parents[1]

class OrientationTests(unittest.TestCase):
    def test_before_start_and_validation_reference(self):
        nb = build_notebook.notebook()
        opening = nb.cells[0].source
        self.assertIn('Before you start', opening)
        orientation = opening.split('Before you start', 1)[1]
        self.assertTrue(200 <= len(orientation.split()) <= 300)
        for term in ('ProteinMPNN', 'AbLang', 'biological replications', 'BioPhi', 'OASis', 'heavy', 'light', 'VHH', 'scFv', 'CDR', 'H3', 'RMSD', 'Å', 'pLDDT', 'ablation'):
            self.assertIn(term, orientation)
        validation = next(c.source for c in nb.cells if c.source.startswith('## 3'))
        self.assertIn('README.md#what-is-checked', validation)
        self.assertIn('stop', validation)
        self.assertNotIn('three backbone atoms', validation)
        self.assertIn('## What is checked', (ROOT/'README.md').read_text())

    def test_map_semantics_counts_layout_and_placement(self):
        import matplotlib.pyplot as plt
        from matplotlib.patches import FancyArrowPatch
        data = analysis.load(ROOT)
        for main_n, ablation_n in ((400,800), (13,26), (13,0)):
            selected = dict(data, main=data['main'].iloc[:main_n], ablation=data['ablation'].iloc[:ablation_n])
            fig = analysis.plot_experiment_map(selected)
            texts = [t for ax in fig.axes for t in ax.texts]
            text = '\n'.join(t.get_text() for t in texts)
            for word in ('Question:', 'Compare:', 'Measure:', 'Contrast:', f'{main_n} draws', 'unpaired', 'No structural'):
                self.assertIn(word, text)
            self.assertIn(f'{ablation_n} draws' if ablation_n else 'not selected', text)
            self.assertIn('reported seed' if ablation_n else 'No ablation comparison', text)
            self.assertEqual(len(fig.axes), 2)
            self.assertLessEqual(fig.get_size_inches()[1], 5.1)
            self.assertGreaterEqual(min(t.get_fontsize() for t in texts), 13)
            self.assertFalse(any(isinstance(p,FancyArrowPatch) for ax in fig.axes for p in ax.patches))
            fig.canvas.draw()
            bounds = fig.bbox
            for t in texts:
                box = t.get_window_extent(fig.canvas.get_renderer())
                self.assertTrue(bounds.contains(box.x0, box.y0) and bounds.contains(box.x1, box.y1), t.get_text())
            plt.close(fig)
        nb = build_notebook.notebook()
        map_cells = [i for i,c in enumerate(nb.cells) if 'analysis.plot_experiment_map(data)' in c.source]
        self.assertEqual(len(map_cells),1)
        validation = next(i for i,c in enumerate(nb.cells) if 'data = validate_inputs(loaded)' in c.source)
        calculation = next(i for i,c in enumerate(nb.cells) if c.source.startswith('## 4'))
        self.assertTrue(validation < map_cells[0] < calculation)
        self.assertIn('map after validation', next(c.source for c in nb.cells if c.source.startswith('## 7')))

    def test_readers_share_orientation_and_current_map(self):
        import tempfile
        with tempfile.TemporaryDirectory(dir=ROOT) as tmp:
            analysis.publish(Path(tmp), data=analysis.load(ROOT))
            for name in ('READ_RESULTS.md','results.html'):
                text=(Path(tmp)/name).read_text()
                self.assertIn('Before you start',text)
                self.assertLess(text.index('Before you start'),text.index('Main comparison'))
                self.assertIn('ProteinMPNN',text)
                self.assertNotIn('Arrows show',text)

    def test_change_to_rerun_route_is_explicit(self):
        guide=(ROOT/'README.md').read_text()
        self.assertIn('## What changed? What do I rerun?', guide)
        section=guide.split('## What changed? What do I rerun?',1)[1].split('\n## ',1)[0]
        for phrase in ('Restart Kernel and Run All', 'sequences', 'rescoring', 'refolding', 'stale', 'ablation', 'formatting', 'Mann', 'correlations', 'multiplicity', 'descriptive'):
            self.assertIn(phrase, section)

if __name__ == '__main__': unittest.main()
