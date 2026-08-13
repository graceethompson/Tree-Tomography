import pandas as pd, numpy as np
from scipy import stats
import statsmodels.formula.api as smf
import statsmodels.api as sm
from statsmodels.stats.anova import anova_lm

pd.set_option('display.width',160); pd.set_option('display.max_columns',30)

df = pd.read_csv('./data/Tree_ID_info.csv')
df['site'] = np.where(df['plot']=='BGS','BGS','EMS')
df['decay_bin'] = (df['percent_damaged']>0).astype(int)
sppmap={'rm':'A.rubrum','bg':'N.sylvatica','ro':'Q.rubra','hem':'T.canadensis'}
df['sp']=df['species'].map(sppmap)

print("="*70,"\nN by species x site\n","="*70)
print(pd.crosstab(df['sp'],df['site'],margins=True))
print("\nDecay present counts by species x site:")
print(pd.crosstab([df['sp'],df['site']],df['decay_bin']))

# ---------- A. STEM SIZE CONFOUND (Reviewer #8) ----------
print("\n"+"="*70,"\nA. DBH by site & species (Reviewer #8 stem-size confound)\n","="*70)
print(df.groupby('site')['dbh'].agg(['n' if False else 'count','mean','std','min','max']))
print(df.groupby('sp')['dbh'].agg(['count','mean','std','min','max']))
t=stats.ttest_ind(df[df.site=='BGS'].dbh, df[df.site=='EMS'].dbh)
print(f"\nDBH BGS vs EMS: BGS mean={df[df.site=='BGS'].dbh.mean():.1f}, EMS mean={df[df.site=='EMS'].dbh.mean():.1f}, t-test p={t.pvalue:.4f}")
# does dbh correlate with percent damaged?
r=stats.pearsonr(df.dbh,df.percent_damaged); rs=stats.spearmanr(df.dbh,df.percent_damaged)
print(f"dbh vs percent_damaged: Pearson r={r[0]:.3f} p={r[1]:.3f}; Spearman rho={rs[0]:.3f} p={rs[1]:.3f}")
# among decayed only
d=df[df.decay_bin==1]
r2=stats.pearsonr(d.dbh,d.percent_damaged)
print(f"among decayed (n={len(d)}): dbh vs pct_damaged Pearson r={r2[0]:.3f} p={r2[1]:.3f}")
# dbh vs decay presence
print("dbh mean by decay presence:", df.groupby('decay_bin')['dbh'].mean().round(1).to_dict())
glm_dbh=smf.glm('decay_bin ~ dbh', df, family=sm.families.Binomial()).fit()
print(f"decay~dbh logistic: beta={glm_dbh.params['dbh']:.3f} p={glm_dbh.pvalues['dbh']:.3f}")

# ---------- B. ANOVA one-way vs two-way (Reviewer #6) ----------
print("\n"+"="*70,"\nB. ANOVA on percent_damaged (Reviewer #6)\n","="*70)
m1=smf.ols('percent_damaged ~ site',df).fit()
print("One-way (site):"); print(anova_lm(m1).round(4))
m2=smf.ols('percent_damaged ~ C(sp) + site',df).fit()
print("\nTwo-way additive (species+site):"); print(anova_lm(m2).round(4))
m3=smf.ols('percent_damaged ~ C(sp)*site',df).fit()
print("\nTwo-way with interaction:"); print(anova_lm(m3).round(4))

# ---------- C. Binomial GLM presence (site vs species) ----------
print("\n"+"="*70,"\nC. Decay presence GLM (site vs species)\n","="*70)
gnull=smf.glm('decay_bin ~ 1',df,family=sm.families.Binomial()).fit()
gsite=smf.glm('decay_bin ~ site',df,family=sm.families.Binomial()).fit()
gspp=smf.glm('decay_bin ~ C(sp)',df,family=sm.families.Binomial()).fit()
gfull=smf.glm('decay_bin ~ C(sp)+site',df,family=sm.families.Binomial()).fit()
def lrt(a,b):
    stat=2*(b.llf-a.llf); dfd=b.df_model-a.df_model; return stat, dfd, stats.chi2.sf(stat,dfd)
print("Null->Site LRT:", [round(x,4) for x in lrt(gnull,gsite)])
print("Site->Full (add species):", [round(x,4) for x in lrt(gsite,gfull)])
print("Null->Species LRT:", [round(x,4) for x in lrt(gnull,gspp)])
print("Species->Full (add site):", [round(x,4) for x in lrt(gspp,gfull)])
print("Full model site coef:", round(gfull.params.get('site[T.EMS]',np.nan),3), "p=",round(gfull.pvalues.get('site[T.EMS]',np.nan),4))
print("Prevalence by site:", df.groupby('site')['decay_bin'].mean().round(3).to_dict())

# ---------- D. Threshold sensitivity: reconstruct PC1 & classification ----------
print("\n"+"="*70,"\nD. Reconstruct ERT PC1 (species z-scores) + classification thresholds\n","="*70)
ert=pd.read_csv('./data/ERT_application_results.csv')
m=df.merge(ert,on='tree',how='inner')
print(f"Merged main trees with ERT: {len(m)}")
metrics=['mean','median','sd','cv','gini','entropy','cma','radialgradient']
# species-specific z-scores
z=m.copy()
for met in metrics:
    z[met+'_z']=m.groupby('sp')[met].transform(lambda x:(x-x.mean())/x.std(ddof=1))
Z=z[[met+'_z' for met in metrics]].values
Z=Z-Z.mean(0)
# PCA
U,S,Vt=np.linalg.svd(Z,full_matrices=False)
scores=U*S
ve=(S**2)/np.sum(S**2)
# orient PC1 so higher = lower resistivity (neg loading on mean). loading sign:
load_mean=Vt[0][metrics.index('mean')]
pc1=scores[:,0]*(-1 if load_mean>0 else 1)
z['pc1']=pc1
print(f"PC1 var explained={ve[0]*100:.1f}%, PC2={ve[1]*100:.1f}%")
thr=pc1.mean()
z['struct_loss']=z['percent_damaged']>1
z['anom']=z['pc1']>thr
# classification
def classify(r):
    if not r['struct_loss'] and not r['anom']: return 'I_NoDecay'
    if not r['struct_loss'] and r['anom']: return 'II_Incipient'
    if r['struct_loss'] and r['anom']: return 'III_Active'
    return 'IV_Cavity'
z['cls']=z.apply(classify,axis=1)
print("Classification counts (threshold=mean, PC1>mean):")
print(z['cls'].value_counts())
print("By site:"); print(pd.crosstab(z['site'],z['cls'],normalize='index').round(2))
# threshold sensitivity: shift threshold by +/- fractions of SD
print("\nThreshold sensitivity for 'incipient' assignment among structurally-sound trees:")
sd1=pc1.std()
for delta,lab in [(-0.5,'mean-0.5SD'),(0,'mean'),(0.5,'mean+0.5SD'),(1,'mean+1SD')]:
    t2=pc1.mean()+delta*sd1
    anom=pc1>t2
    inc=((~z['struct_loss'])&anom).sum()
    # site split of incipient
    bgs_inc=((~z['struct_loss'])&anom&(z['site']=='BGS')).mean()
    print(f"  {lab}: n_incipient={inc}, n_anomalous_total={anom.sum()} / {len(z)}")
# how many trees flip class if threshold moves +/-0.5SD
base=z['cls'].copy()
anom2=pc1>(pc1.mean()+0.5*sd1)
z['cls2']=z.apply(lambda r:classify(r),axis=1)
flips=(base!= pd.Series([classify({'struct_loss':s,'anom':a}) for s,a in zip(z['struct_loss'],anom2)],index=z.index)).sum()
print(f"Trees changing class when threshold -> mean+0.5SD: {flips} of {len(z)}")

# 1% SoT threshold sensitivity
print("\n1% SoT threshold: trees by damage band:")
bands=pd.cut(df.percent_damaged,[-.1,0,1,5,100],labels=['0','(0-1]','(1-5]','>5'])
print(bands.value_counts().sort_index())

# ---------- E. Validation: PC1 vs mean resistivity vs moisture (Reviewer #3) ----------
print("\n"+"="*70,"\nE. Hemlock validation: PC1 vs mean resistivity vs moisture (Reviewer #3)\n","="*70)
v=pd.read_csv('./data/hemlock/validation_summary.csv')
# metrics available: Conductance? uses Median,SD,CV,Gini,Entropy,CMA,RadialGradient + Conductance. mean resistivity = 1000/Conductance
v['mean_res']=1000/v['Conductance']
vmet=['mean_res','Median','SD','CV','Gini','Entropy','CMA','RadialGradient']
Vz=v[vmet].apply(lambda x:(x-x.mean())/x.std(ddof=1)).values
Vz=Vz-Vz.mean(0)
Uu,Ss,Vtt=np.linalg.svd(Vz,full_matrices=False)
vpc1=(Uu*Ss)[:,0]
# orient: higher = lower resistivity
if np.corrcoef(vpc1,v['mean_res'])[0,1]>0: vpc1=-vpc1
rpc=stats.pearsonr(vpc1,v['moisture']); rmean=stats.pearsonr(v['mean_res'],v['moisture'])
print(f"n={len(v)} hemlocks")
print(f"PC1 vs moisture: r={rpc[0]:.3f} p={rpc[1]:.3f}")
print(f"mean resistivity vs moisture: r={rmean[0]:.3f} p={rmean[1]:.3f}")
print(f"Conductance vs moisture: r={stats.pearsonr(v['Conductance'],v['moisture'])[0]:.3f}")
for met in ['Median','CV','Gini','Entropy','CMA','RadialGradient']:
    rr=stats.pearsonr(v[met],v['moisture']); print(f"  {met} vs moisture: r={rr[0]:.3f} p={rr[1]:.3f}")

# ---------- F. Carbon: 47% vs 26% and mean damage among structural-loss ----------
print("\n"+"="*70,"\nF. Carbon framing numbers\n","="*70)
print(f"% trees any anomaly (incipient+active+cavity) = {(z['cls']!='I_NoDecay').mean()*100:.0f}%")
print(f"% trees structural loss (>1% SoT) = {(df.percent_damaged>1).mean()*100:.0f}%")
sl=df[df.percent_damaged>1]
print(f"Among structural-loss trees: mean pct damage overall={sl.percent_damaged.mean():.1f}%")
print("  by site:", sl.groupby('site')['percent_damaged'].mean().round(1).to_dict())
print(f"Mean pct damage all trees by site:", df.groupby('site')['percent_damaged'].mean().round(2).to_dict())
