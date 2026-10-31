"""Synthetic arithmetic fixtures, not research observations."""
import unittest
import numpy as np
import stats_helper as s

class Tests(unittest.TestCase):
    def test_unpaired_ties_and_direction(self):
        r=s.compare([1,2],[2,3],n_perm=99,n_boot=20,seed=42)
        # Pairwise E-B: 1,0,2,1; 3 wins, one tie.
        self.assertEqual(r['U_ensemble'],3.5)
        self.assertEqual(r['prob_greater_plus_half_ties'],.875)
        self.assertEqual(r['rank_biserial'],.75)
        self.assertEqual(r['hodges_lehmann'],1)
        self.assertEqual(r['difference_medians'],1)
        self.assertAlmostEqual(r['p'],2/3)
        self.assertEqual(s.compare([1,1],[2,2],n_boot=20)['p'],1/3)
        self.assertEqual(s.compare([1,1],[1,1],n_boot=20)['p'],1)
        with self.assertRaises(ValueError): s.compare([1,np.nan],[2,3])
        with self.assertRaises(ValueError): s.compare([1],[2,3])

    def test_sequence_join_and_association(self):
        import pandas as pd
        scores=pd.DataFrame({'id':['a','b','c'],'sequence':['AA','AB','AC'],'score':[1,2,3]})
        folds=pd.DataFrame({'id':['c','a','b'],'sequence':['AC','AA','AB'],'rmsd':[1,3,2]})
        r=s.associate(scores,folds,n_perm=99)
        self.assertAlmostEqual(r['rho'],-1)
        self.assertAlmostEqual(r['p'],1/3)
        with self.assertRaises(ValueError): s.associate(scores,folds.iloc[:2])
        with self.assertRaises(ValueError): s.associate(scores,pd.concat([folds,folds.iloc[:1]]))
        bad=folds.copy(); bad.loc[0,'sequence']='ZZ'
        with self.assertRaises(ValueError): s.associate(scores,bad)
        bad=scores.copy(); bad['score']=1
        self.assertEqual(s.associate(bad,folds)['status'],'withheld: constant variable')

    def test_family_and_reported_pairs(self):
        import pandas as pd
        self.assertTrue(np.allclose(s.adjust({'a':.01,'b':.04,'c':.03},family=['a','b','c'],method='BH').values.tolist(),[.03,.04,.04]))
        self.assertTrue(np.allclose(s.adjust({'a':.01,'b':.04,'c':.03},family=['a','b','c'],method='Holm').values.tolist(),[.03,.06,.06]))
        with self.assertRaises(ValueError): s.adjust({'a':.01},family=['a','b'],method='BH')
        pairs=pd.DataFrame({'seed':['0','1','0','1'],'method':['baseline','baseline','ensemble','ensemble'],'score':[1,2,2,4]})
        r=s.paired_summary(pairs,pairing_source='synthetic same-seed fixture')
        self.assertEqual(r['mean_paired_difference'],1.5)
        self.assertIsNone(r['p'])
        with self.assertRaises(ValueError): s.paired_summary(pairs.iloc[:3],pairing_source='fixture')
        with self.assertRaises(ValueError): s.paired_summary(pairs,pairing_source='')


if __name__=="__main__": unittest.main()
