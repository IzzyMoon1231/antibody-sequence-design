"""Offline exploratory statistics. No model/scorer imports or filesystem writes.
Draw-weighted conditional inference requires independent/exchangeable sampling.
All effect directions are ensemble minus baseline, not automatically beneficial.
"""
import itertools
import math
import numpy as np
from scipy.stats import rankdata


def vector(x, minimum=2):
    a=np.asarray(x,dtype=float)
    if a.ndim!=1 or len(a)<minimum or not np.isfinite(a).all():
        raise ValueError('Require finite 1-D data with at least %d observations; no silent dropping'%minimum)
    return a


def adjust(pvalues, *, family, method):
    """Require complete named family; never silently correct only available tests."""
    import pandas as pd
    if not family or len(set(family))!=len(family) or set(pvalues)!=set(family):
        raise ValueError('Provide each member of the declared family exactly once')
    p=vector([pvalues[k] for k in family],1)
    if ((p<0)|(p>1)).any(): raise ValueError('p outside [0,1]')
    order=np.argsort(p); m=len(p); ranks=np.arange(1,m+1)
    if method=='BH': q=np.minimum.accumulate((p[order]*m/ranks)[::-1])[::-1]
    elif method=='Holm': q=np.maximum.accumulate(p[order]*(m-ranks+1))
    else: raise ValueError('Explicit supported adjustment: BH or Holm')
    out=np.empty(m); out[order]=np.minimum(q,1)
    return pd.Series(out,index=family)


def paired_summary(frame, *, pairing_source):
    """Descriptive reported seed-pair contrasts only; no execution certification."""
    if not pairing_source or not {'seed','method','score'}<=set(frame.columns):
        raise ValueError('Document the source-defined pairing and provide seed/method/score')
    if frame[['seed','method','score']].isna().any().any() or frame.duplicated(['seed','method']).any():
        raise ValueError('Missing or duplicate reported pair')
    w=frame.pivot(index='seed',columns='method',values='score')
    if set(w.columns)!={'baseline','ensemble'} or w.isna().any().any():
        raise ValueError('Require exactly complete baseline/ensemble pairs')
    d=vector(w.ensemble-w.baseline,1)
    return dict(n_pairs=len(d),mean_paired_difference=float(d.mean()),median_paired_difference=float(np.median(d)),q25_paired_difference=float(np.quantile(d,.25)),q75_paired_difference=float(np.quantile(d,.75)),positive=int((d>0).sum()),negative=int((d<0).sum()),zero=int((d==0).sum()),pairing_source=pairing_source,p=None,status='descriptive reported seed pairs; executed pairing not certified; no test adopted')


def associate(scores, folds, *, n_perm=100000, seed=42):
    """Inputs require real full-sequence evidence, not an ID-only assertion.
    Caller binds score VH to full FASTA and fold sequence to actual PDB/receipt.
    Each call is ONE case/method/endpoint, never pooled across contexts.
    """
    for frame,cols in [(scores,['id','sequence','score']),(folds,['id','sequence','rmsd'])]:
        if not set(cols)<=set(frame.columns) or frame[cols].isna().any().any():
            raise ValueError('Missing join evidence or measurement')
        if frame.id.duplicated().any() or frame.sequence.eq('').any():
            raise ValueError('Duplicate IDs or empty sequences')
    joined=scores.merge(folds,on='id',how='outer',validate='one_to_one',suffixes=('_score','_fold'),indicator=True)
    if not joined._merge.eq('both').all() or not joined.sequence_score.eq(joined.sequence_fold).all():
        raise ValueError('Incomplete or sequence-mismatched score/fold join')
    x,y=vector(joined.score,3),vector(joined.rmsd,3)
    result=dict(n=len(x),unique_sequences=int(joined.sequence_score.nunique()),unique_scores=len(np.unique(x)),status='exploratory; sequence-joined, independent draw assumption')
    if np.ptp(x)==0 or np.ptp(y)==0:
        return dict(result,rho=None,p=None,status='withheld: constant variable')
    if not isinstance(n_perm,int) or n_perm<1: raise ValueError('Positive permutation count required')
    rx,ry=rankdata(x),rankdata(y); rx-=rx.mean(); ry-=ry.mean()
    denom=np.linalg.norm(rx)*np.linalg.norm(ry)
    rho=float(rx@ry/denom); rng=np.random.default_rng(seed)
    exact=len(x)<=8
    permutations=itertools.permutations(ry) if exact else (rng.permutation(ry) for _ in range(n_perm))
    hits=sum(abs(rx@np.asarray(v)/denom)>=abs(rho)-1e-12 for v in permutations)
    p=hits/math.factorial(len(x)) if exact else (hits+1)/(n_perm+1)
    return dict(result,rho=rho,p=p,p_method='exhaustive ranks permutation' if exact else 'Monte Carlo ranks permutation, plus-one',seed=seed,n_permutations=math.factorial(len(x)) if exact else n_perm)


def compare(baseline, ensemble, *, n_perm=100000, n_boot=10000, seed=42):
    b,e=vector(baseline),vector(ensemble)
    if not isinstance(n_perm,int) or n_perm<1 or not isinstance(n_boot,int) or n_boot<2:
        raise ValueError('Positive integer permutations and at least two bootstrap draws required')
    ne,nb=len(e),len(b)
    if ne*nb>1000000: raise ValueError('Pairwise allocation exceeds one million; review implementation before scaling')
    d=e[:,None]-b[None,:]
    u=float((d>0).sum()+.5*(d==0).sum())
    ranks=rankdata(np.r_[e,b]); center=ne*nb/2; observed=abs(u-center)
    rng=np.random.default_rng(seed)
    count=math.comb(ne+nb,ne)
    if ne==nb and (e.min()>b.max() or b.min()>e.max()):
        p=2/count; p_method='exact complete separation, equal arm sizes'
    elif count<=50000:
        hits=sum(abs(ranks[list(ix)].sum()-ne*(ne+1)/2-center)>=observed-1e-12 for ix in itertools.combinations(range(ne+nb),ne))
        p=hits/count; p_method='exhaustive label permutation (ties retained)'
    else:
        hits=0
        for _ in range(n_perm):
            up=rng.permutation(ranks)[:ne].sum()-ne*(ne+1)/2
            hits+=abs(up-center)>=observed-1e-12
        p=(hits+1)/(n_perm+1); p_method='Monte Carlo label permutation, plus-one'
    boot=np.empty(n_boot)
    for i in range(n_boot):
        boot[i]=np.median(rng.choice(e,ne)[:,None]-rng.choice(b,nb)[None,:])
    lo,hi=np.quantile(boot,[.025,.975])
    return dict(n_baseline=nb,n_ensemble=ne,baseline_median=float(np.median(b)),ensemble_median=float(np.median(e)),difference_medians=float(np.median(e)-np.median(b)),hodges_lehmann=float(np.median(d)),hl_ci95_low=float(lo),hl_ci95_high=float(hi),U_ensemble=u,prob_greater_plus_half_ties=u/(ne*nb),rank_biserial=2*u/(ne*nb)-1,p=float(p),p_method=p_method,n_permutations=n_perm,n_bootstrap=n_boot,seed=seed,status='exploratory; conditional on independent draws')
