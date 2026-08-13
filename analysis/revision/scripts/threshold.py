import pandas as pd, numpy as np
from scipy import stats
from sklearn.mixture import GaussianMixture
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

np.random.seed(0)
df = pd.read_csv('./data/Tree_ID_info.csv')
df['site'] = np.where(df['plot']=='BGS','BGS','EMS')
sppmap={'rm':'A.rubrum','bg':'N.sylvatica','ro':'Q.rubra','hem':'T.canadensis'}
df['sp']=df['species'].map(sppmap)
ert=pd.read_csv('./data/ERT_application_results.csv')
m=df.merge(ert,on='tree',how='inner').reset_index(drop=True)
metrics=['mean','median','sd','cv','gini','entropy','cma','radialgradient']

def build_pc1(data, scheme):
    X=data[metrics].astype(float).copy()
    if scheme=='pooled':
        Z=(X-X.mean())/X.std(ddof=1)
    elif scheme=='species':
        Z=X.copy()
        for met in metrics:
            Z[met]=data.groupby('sp')[met].transform(lambda x:(x-x.mean())/x.std(ddof=1))
    Zc=Z.values - Z.values.mean(0)
    U,S,Vt=np.linalg.svd(Zc,full_matrices=False)
    ve=(S**2)/np.sum(S**2)
    pc1=(U*S)[:,0]
    # orient: higher PC1 = lower resistivity => negative corr with mean resistivity
    if np.corrcoef(pc1,data['mean'])[0,1]>0: pc1=-pc1; Vt[0]=-Vt[0]
    load=dict(zip(metrics,Vt[0]))
    return pc1, ve[0], ve[1], load

for scheme in ['species','pooled']:
    pc1,ve1,ve2,load=build_pc1(m,scheme)
    m[f'pc1_{scheme}']=pc1
    print(f"\n=== {scheme} normalization: PC1 {ve1*100:.1f}%, PC2 {ve2*100:.1f}% ===")
    print("  loadings:", {k:round(v,2) for k,v in load.items()})

print("\ncorr(species-PC1, pooled-PC1) =", round(np.corrcoef(m['pc1_species'],m['pc1_pooled'])[0,1],3))

# ---- MODALITY: is there a natural break? ----
print("\n"+"="*66,"\nMODALITY of the anomaly axis (natural break?)\n","="*66)
def modality(x,label):
    x=np.asarray(x).reshape(-1,1)
    g1=GaussianMixture(1,random_state=0).fit(x); g2=GaussianMixture(2,random_state=0,n_init=5).fit(x)
    print(f"{label}: BIC 1-comp={g1.bic(x):.1f}  2-comp={g2.bic(x):.1f}  -> {'2-comp (bimodal) favored' if g2.bic(x)+2<g1.bic(x) else 'unimodal favored'}")
    # skew/kurt
    print(f"   skew={stats.skew(x.ravel()):.2f} kurtosis={stats.kurtosis(x.ravel()):.2f}")
    return g2
g2_sp=modality(m['pc1_species'],'species-PC1')
g2_pl=modality(m['pc1_pooled'],'pooled-PC1')

# ---- Threshold rules on species-PC1 (primary axis) ----
pc=m['pc1_species'].values
mu,med,sd=pc.mean(),np.median(pc),pc.std(ddof=1)
def otsu(x,nbins=256):
    hist,edges=np.histogram(x,bins=nbins); centers=(edges[:-1]+edges[1:])/2
    w=hist/hist.sum(); cum=np.cumsum(w); cummean=np.cumsum(w*centers); gmean=cummean[-1]
    best_t,best_v=None,-1
    for i in range(1,nbins):
        w0=cum[i]; w1=1-w0
        if w0==0 or w1==0: continue
        m0=cummean[i]/w0; m1=(gmean-cummean[i])/w1
        v=w0*w1*(m0-m1)**2
        if v>best_v: best_v,best_t=v,centers[i]
    return best_t
def gmm_cross(gmm):
    # crossover between the two component densities within overlap
    means=gmm.means_.ravel(); order=np.argsort(means)
    lo,hi=means[order]
    xs=np.linspace(lo,hi,2000)
    logp=gmm._estimate_weighted_log_prob(xs.reshape(-1,1))  # weighted comp log-prob
    diff=logp[:,order[0]]-logp[:,order[1]]
    sign=np.sign(diff); idx=np.where(np.diff(sign)!=0)[0]
    return xs[idx[0]] if len(idx) else np.nan
t_otsu=otsu(pc); t_gmm=gmm_cross(g2_sp)
rules={'mean (published)':mu,'median':med,'GMM 2-comp crossover':t_gmm,"Otsu break":t_otsu,
       'mean-0.5SD':mu-0.5*sd,'mean+0.5SD':mu+0.5*sd,'mean+1SD':mu+1*sd}
m['struct']=m['percent_damaged']>1
print("\n"+"="*66,"\nTHRESHOLD SENSITIVITY on species-PC1\n","="*66)
print(f"{'rule':<22}{'thr':>7}{'nAnom':>7}{'nIncip':>7}{'BGSinc%':>9}{'EMSinc%':>9}  species incipient %")
rows=[]
for name,t in rules.items():
    anom=pc>t
    incip=(~m['struct'])&anom
    bgs=incip[m.site=='BGS'].sum()/(m.site=='BGS').sum()*100
    ems=incip[m.site=='EMS'].sum()/(m.site=='EMS').sum()*100
    sppc={s: round(incip[m.sp==s].sum()/(m.sp==s).sum()*100) for s in sorted(m.sp.unique())}
    print(f"{name:<22}{t:>7.2f}{anom.sum():>7}{incip.sum():>7}{bgs:>9.0f}{ems:>9.0f}  {sppc}")
    rows.append((name,t,anom.sum(),incip.sum(),bgs,ems))

# how many trees flip anomaly label across the full rule set
labels=np.vstack([(pc>t).astype(int) for t in rules.values()])
flip=(labels.min(0)!=labels.max(0)).sum()
print(f"\nTrees whose anomaly label is NOT stable across all rules: {flip}/{len(m)}")
print(f"Always-anomalous: {(labels.min(0)==1).sum()}, never-anomalous: {(labels.max(0)==0).sum()}")

# ---- Within vs between: site contrast under each normalization at its own mean ----
print("\n"+"="*66,"\nWITHIN (species) vs BETWEEN (pooled) normalization: site incipient contrast\n","="*66)
for scheme in ['species','pooled']:
    a=m[f'pc1_{scheme}'].values; t=a.mean()
    incip=(~m['struct'])&(a>t)
    bgs=incip[m.site=='BGS'].mean()*100; ems=incip[m.site=='EMS'].mean()*100
    # specialists: do their PC1 carry between-site info? report per-species mean PC1
    sp_means={s: round(a[m.sp==s].mean(),2) for s in sorted(m.sp.unique())}
    print(f"{scheme:>8}: BGS incip {bgs:.0f}% vs EMS {ems:.0f}%  | per-species mean PC1 {sp_means}")

# ---- External anchor demonstration via validation regression ----
print("\n"+"="*66,"\nEXTERNAL ANCHOR demo (validation moisture <-> mean resistivity)\n","="*66)
v=pd.read_csv('./data/hemlock/validation_summary.csv')
v['mean_res']=1000/v['Conductance']
sl=stats.linregress(v['mean_res'],v['moisture'])
print(f"validation: moisture = {sl.intercept:.1f} + {sl.slope:.3f}*mean_res  (r={sl.rrvalue if hasattr(sl,'rrvalue') else sl.rvalue:.2f})")
# invert: resistivity at a moisture anchor
for Mcut in [90,100,110]:
    Rcut=(Mcut-sl.intercept)/sl.slope
    n_below=(m['mean']<Rcut).sum()
    print(f"  moisture>{Mcut}% <=> mean_res<{Rcut:.0f} Ohm-m -> {n_below}/{len(m)} main trees flagged (species-blind, absolute)")

# ---------- FIGURE ----------
fig,axes=plt.subplots(1,2,figsize=(12,4.6))
for ax,scheme,title in [(axes[0],'pc1_species','Within-species normalization (published)'),
                        (axes[1],'pc1_pooled','Pooled normalization (preserves absolute differences)')]:
    a=m[scheme].values
    for s,mk in zip(sorted(m.sp.unique()),['o','s','D','^']):
        ax.hist(a[m.sp==s],bins=np.linspace(a.min(),a.max(),18),alpha=0.0)  # keep range
    # kde
    xs=np.linspace(a.min()-0.3,a.max()+0.3,300)
    kde=stats.gaussian_kde(a)
    ax.plot(xs,kde(xs),color='k',lw=1.5)
    ax.hist(a,bins=16,density=True,color='0.85',edgecolor='white')
    cols={'mean (published)':'#c0392b','median':'#2980b9','GMM 2-comp crossover':'#27ae60',"Otsu break":'#8e44ad'}
    if scheme=='pc1_species':
        for name,c in cols.items():
            t=rules[name]
            if np.isfinite(t): ax.axvline(t,color=c,ls='--',lw=1.6,label=f"{name}={t:.2f}")
    else:
        ax.axvline(a.mean(),color='#c0392b',ls='--',lw=1.6,label=f"mean={a.mean():.2f}")
    ax.set_title(title,fontsize=10); ax.set_xlabel('ERT PC1 (higher = lower resistivity / more heterogeneous)')
    ax.set_ylabel('density'); ax.legend(fontsize=7,frameon=False)
plt.tight_layout()
plt.savefig('analysis/revision/output/CJFR-threshold-rules.png',dpi=170)
print("\nsaved figure")
