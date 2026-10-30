"""Presentation-only contracts; mutations are synthetic controls."""
import unittest
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import analysis
ROOT = Path(__file__).resolve().parents[1]

class PanelTests(unittest.TestCase):
    def test_four_visible_cells_and_local_interpretations(self):
        import build_notebook
        nb = build_notebook.notebook()
        panels = [(i,c) for i,c in enumerate(nb.cells) if c.cell_type=='code' and 'analysis.plot_panel(' in c.source]
        expected = [(s,m) for s in analysis.CASES for m in ('similarity','structure')]
        self.assertEqual(len(panels), 4)
        for (i,c),(st,metric) in zip(panels,expected):
            self.assertIn(repr(st),c.source)
            self.assertIn(repr(metric),c.source)
            self.assertIn('?',nb.cells[i-1].source)
            self.assertIn('analysis.panel_readout(',nb.cells[i+1].source)
        source = '\n'.join(c.source for c in nb.cells if c.cell_type=='code')
        self.assertNotIn('publish(',source)
        self.assertNotIn('Image(',source)
        self.assertNotIn('plot_main(',source)
        self.assertIn('EXPORT = False',source)

    def test_panel_captions_and_exports_follow_changed_data(self):
        import tempfile
        d=analysis.load(ROOT)
        initial=analysis.panel_readout(d,'9NH7','similarity')
        d['main'].loc[d['main'].method=='ensemble','oasis']=.1
        self.assertNotEqual(initial,analysis.panel_readout(d,'9NH7','similarity'))
        self.assertIn('lower',analysis.panel_readout(d,'9NH7','similarity'))
        self.assertIn('H3',analysis.panel_readout(d,'9NFU','structure'))
        d['original_inputs']=False
        with tempfile.TemporaryDirectory(dir=ROOT) as tmp:
            dest=Path(tmp)
            analysis.publish(dest,data=d)
            for st in analysis.CASES:
                for metric in ('similarity','structure'):
                    name=f'{st}_{metric}'
                    for ext in ('png','svg'):
                        self.assertTrue((dest/f'figures/{name}.{ext}').is_file())
                    self.assertIn(f'figures/{name}.png',(dest/'READ_RESULTS.md').read_text())
                    self.assertIn(analysis.panel_readout(d,st,metric),(dest/'READ_RESULTS.md').read_text())
                    self.assertIn(name,(dest/'results.html').read_text())
            for name in ('READ_RESULTS.md','results.html'):
                text=(dest/name).read_text()
                self.assertIn('figure1.svg',text)
                self.assertNotIn('Pinned retained-input analysis',text)

    def tearDown(self):
        plt.close('all')

    def test_compact_layout_and_exact_panel_artists(self):
        data = analysis.load(ROOT)
        combined = analysis.plot_main(data)
        self.assertEqual(combined.axes[0].get_subplotspec().get_gridspec().get_geometry(), (2, 2))
        for i, case in enumerate(analysis.CASES):
            for row, metric in enumerate(('similarity', 'structure')):
                one = analysis.plot_panel(data, case, metric).axes[0]
                both = combined.axes[row*2+i]
                self.assertIn(case, both.get_title(loc='left'))
                self.assertTrue(both.get_title(loc='left').startswith(chr(65+row*2+i)))
                self.assertEqual(one.get_xlim(), both.get_xlim())
                self.assertEqual(one.get_xlabel(), both.get_xlabel())
                self.assertEqual([t.get_text() for t in one.get_yticklabels()], [t.get_text() for t in both.get_yticklabels()])
                for a,b in zip(one.collections, both.collections, strict=True):
                    np.testing.assert_array_equal(a.get_offsets(), b.get_offsets())
                    np.testing.assert_array_equal(a.get_facecolors(), b.get_facecolors())
                    np.testing.assert_array_equal(a.get_edgecolors(), b.get_edgecolors())
                    if hasattr(a, 'get_segments'):
                        np.testing.assert_array_equal(a.get_segments(), b.get_segments())
                for a,b in zip(one.lines, both.lines, strict=True):
                    np.testing.assert_array_equal(a.get_xydata(), b.get_xydata())
                    self.assertEqual(a.get_color(), b.get_color())
        for row in (0, 1):
            self.assertEqual(combined.axes[2*row].get_xlim(), combined.axes[2*row+1].get_xlim())

    def test_changed_frames_expand_shared_scales_without_clipping(self):
        data = analysis.load(ROOT)
        data['loops'].loc[(data['loops'].structure=='9NFU') & (data['loops'].cdr=='L3'), 'rmsd_angstrom'] = 80
        data['main'].loc[data['main'].method=='ensemble','oasis'] = .1
        figures = [analysis.plot_panel(data, st, 'structure') for st in analysis.CASES]
        self.assertEqual(figures[0].axes[0].get_xlim(), figures[1].axes[0].get_xlim())
        self.assertGreater(figures[0].axes[0].get_xlim()[1],80)
        panel = analysis.plot_panel(data, '9NH7', 'similarity')
        np.testing.assert_allclose(panel.axes[0].collections[2].get_offsets()[:,0], .1)
        with self.assertRaises(ValueError): analysis.plot_panel(data,'other','structure')
        with self.assertRaises(ValueError): analysis.plot_panel(data,'9NH7','other')
