import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
df=pd.read_csv('./data/Tree_ID_info.csv')
df['site']=np.where(df['plot']=='BGS','BGS','EMS')
sm={'rm':'A.rubrum','bg':'N.sylvatica','ro':'Q.rubra','hem':'T.canadensis'}; df['sp']=df['species'].map(sm)
ert=pd.read_csv('./data/ERT_application_results.csv')
m=df.merge(ert,on='tree').reset_index(drop=True)
mets=['mean','median','sd','cv','gini','entropy','cma','radialgradient']
Z=m[mets].copy()
for c in mets: Z[c]=m.groupby('sp')[c].transform(lambda x:(x-x.mean())/x.std(ddof=1))
Zc=Z.values-Z.values.mean(0); U,S,Vt=np.linalg.svd(Zc,full_matrices=False); pc=(U*S)[:,0]
if np.corrcoef(pc,m['mean'])[0,1]>0: pc=(-pc); Vt[0]=-Vt[0]
m['pc1']=pc; m['res_z']=Z['mean']; m['struct']=m['percent_damaged']>1
mu=pc.mean()

# ---- PC1 loadings: does PC1 "contain" CMA? ----
print("PC1 loadings (does PC1 include CMA?):")
for met,l in sorted(zip(mets,Vt[0]),key=lambda t:-abs(t[1])): print(f"   {met:14s}{l:+.2f}")
print(f"   -> CMA correlation with PC1: r={np.corrcoef(m['cma'],pc)[0,1]:+.2f}   (CMA loads on PC2, not PC1)")
print(f"   -> mean resistivity vs PC1: r={np.corrcoef(m['mean'],pc)[0,1]:+.2f}")

def four(anom, struct):
    return np.where(~struct&~anom,'I',np.where(~struct&anom,'II',np.where(struct&anom,'III','IV')))
# Scheme A: original PC1 (anomaly = wet/heterogeneous = PC1>mean)
A=four(pc>mu, m['struct'].values)
# Scheme B: CMA (anomaly = central moisture accumulation > 0.33, the paper's own 'even' point)
B=four(m['cma'].values>0.33, m['struct'].values)
# Scheme C: absolute mean resistivity (anomaly = wetter than sample median)
C=four(m['mean'].values<np.median(m['mean']), m['struct'].values)
# Scheme D: species-normalized resistivity only (anomaly = res_z<0)
D=four(m['res_z'].values<0, m['struct'].values)
m['A']=A; m['B']=B; m['C']=C; m['D']=D

from collections import Counter
for name,arr in [('A PC1 (original)',A),('B CMA>0.33',B),('C abs-resistivity',C),('D species-resistivity',D)]:
    print(f"\n{name}: {dict(Counter(arr))}")

def agree(x,y): return (x==y).mean()*100
print("\n--- AGREEMENT with original PC1 scheme (A) ---")
for name,arr in [('B CMA',B),('C abs-res',C),('D sp-res',D)]:
    print(f"  A vs {name}: {agree(A,arr):.0f}% identical  ({(A!=arr).sum()} of 57 trees change category)")
print(f"  A vs B, restricted to the moisture-axis calls only (I/II among sound; III/IV among damaged):")
sound=~m['struct'].values; dmgd=m['struct'].values
print(f"     I/II (sound):  {agree(A[sound],B[sound]):.0f}% ({(A[sound]!=B[sound]).sum()}/{sound.sum()} differ)")
print(f"     III/IV (damaged): {agree(A[dmgd],B[dmgd]):.0f}% ({(A[dmgd]!=B[dmgd]).sum()}/{dmgd.sum()} differ)")

# confusion A vs B
print("\nConfusion A (rows) vs B (cols):")
print(pd.crosstab(pd.Series(A,name='PC1'),pd.Series(B,name='CMA')))

# ---------------- FIGURE: phase diagram under 3 binnings ----------------
catcol={'I':'#3b6fb0','II':'#5aa02c','III':'#d98a1f','IV':'#b83232'}
def panel(ax,x,anom_line,xlabel,cats,title,xinvert=False,logx=False):
    y=np.sqrt(np.clip(m['percent_damaged'],0,None))
    ax.axhline(np.sqrt(1),color='#888',ls='--',lw=1)
    ax.axvline(anom_line,color='#888',ls='--',lw=1)
    for cat,c in catcol.items():
        idx=cats==cat
        ax.scatter(x[idx], y[idx]+np.random.normal(0,0.03,idx.sum()), c=c, s=42, edgecolor='white',lw=.5,label=cat)
    yt=[0,1,5,10,20,30]; ax.set_yticks([np.sqrt(t) for t in yt]); ax.set_yticklabels(yt)
    ax.set_xlabel(xlabel,fontsize=9); ax.set_title(title,fontsize=9.5)
    if xinvert: ax.invert_xaxis()
    if logx: ax.set_xscale('log')
np.random.seed(1)
fig,axes=plt.subplots(1,3,figsize=(16,5.2))
panel(axes[0],pc,mu,'ERT PC1 (species-normalized composite)',A,'A. Original — bin by PC1')
panel(axes[1],m['cma'].values,0.33,'CMA — central moisture accumulation',B,'B. Bin by CMA (>0.33)')
panel(axes[2],m['mean'].values,np.median(m['mean']),'mean resistivity (Ω·m, absolute)',C,'C. Bin by absolute resistivity',xinvert=True)
axes[0].set_ylabel('SoT structural loss %  (√ scale)')
axes[0].legend(title='Category',fontsize=8,frameon=False,loc='upper left')
fig.suptitle('Does the ERT axis choice change the binning? Same trees, three anomaly criteria',fontsize=12,y=1.01)
plt.tight_layout(); plt.savefig('analysis/revision/output/CJFR-binning-compare.png',dpi=150,bbox_inches='tight')

# second fig: PC1 diagram colored by CMA value (continuous) to show CMA is a SEPARATE axis
fig2,ax=plt.subplots(figsize=(7.6,5.6))
y=np.sqrt(np.clip(m['percent_damaged'],0,None))
ax.axhline(np.sqrt(1),color='#888',ls='--',lw=1); ax.axvline(mu,color='#888',ls='--',lw=1)
sc=ax.scatter(pc,y+np.random.normal(0,0.03,len(y)),c=m['cma'],cmap='RdYlBu',s=60,edgecolor='k',lw=.4,vmin=0,vmax=0.5)
yt=[0,1,5,10,20,30]; ax.set_yticks([np.sqrt(t) for t in yt]); ax.set_yticklabels(yt)
ax.set_xlabel('ERT PC1'); ax.set_ylabel('SoT structural loss % (√)')
plt.colorbar(sc,label='CMA (blue=wet core, red=dry core)')
ax.set_title('PC1 phase diagram colored by CMA:\nCMA does NOT line up with the PC1 axis (r=+0.05) — it is orthogonal (PC2)')
plt.tight_layout(); plt.savefig('analysis/revision/output/CJFR-pc1-cma.png',dpi=150,bbox_inches='tight')
print("\nsaved figures")
