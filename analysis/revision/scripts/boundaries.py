import pandas as pd, numpy as np
from scipy import stats
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

np.random.seed(0)
df=pd.read_csv('./data/Tree_ID_info.csv')
df['site']=np.where(df['plot']=='BGS','BGS','EMS')
sppmap={'rm':'A.rubrum','bg':'N.sylvatica','ro':'Q.rubra','hem':'T.canadensis'}; df['sp']=df['species'].map(sppmap)
ert=pd.read_csv('./data/ERT_application_results.csv')
m=df.merge(ert,on='tree',how='inner').reset_index(drop=True)
metrics=['mean','median','sd','cv','gini','entropy','cma','radialgradient']
Z=m[metrics].copy()
for met in metrics: Z[met]=m.groupby('sp')[met].transform(lambda x:(x-x.mean())/x.std(ddof=1))
Zc=Z.values-Z.values.mean(0); U,S,Vt=np.linalg.svd(Zc,full_matrices=False)
pc=(U*S)[:,0]
if np.corrcoef(pc,m['mean'])[0,1]>0: pc=-pc
m['pc']=pc; dmg=m['percent_damaged'].values
mu,sd=pc.mean(),pc.std(ddof=1)

def classify(struct,anom):
    return np.where(~struct&~anom,'I', np.where(~struct&anom,'II', np.where(struct&anom,'III','IV')))

# ---- the two axes have different shape ----
print("SoT structural axis (percent_damaged): distribution")
print("  values >0:", sorted(dmg[dmg>0]))
bands=pd.cut(dmg,[-.1,0,1,3,5,10,100],labels=['0','(0,1]','(1,3]','(3,5]','(5,10]','>10'])
print("  bands:", bands.value_counts().sort_index().to_dict())
print(f"ERT moisture axis: unimodal continuum, SD={sd:.2f}, all trees have a value\n")

base=classify(dmg>1, pc>mu)
from collections import Counter
print("Baseline categories:", dict(Counter(base)))

# ======== BOUNDARY-BY-BOUNDARY STABILITY ========
# Each boundary crosses ONE axis. Vary that axis' threshold over a plausible range; count flips.
print("\n"+"="*80,"\nBOUNDARY STABILITY: which category edges move when thresholds move?\n","="*80)

# --- SoT-axis boundaries: I<->IV (anom=0) and II<->III (anom=1). Vary tau_SoT 1..5% ---
print("\n[SoT-axis boundaries]  I<->IV and II<->III  (vary structural-loss threshold 1%..5%)")
for tau in [0.5,1,2,3,5,10]:
    c=classify(dmg>tau, pc>mu)
    ch=(c!=base).sum()
    print(f"  tau_SoT={tau:>4}% : categories={dict(Counter(c))}  trees changed vs baseline={ch}")
n_soT_band=((dmg>1)&(dmg<=5)).sum()
print(f"  -> trees in the 1%-5% 'ambiguous structural' band: {n_soT_band}  (gap: 0 trees in (0,1])")

# --- ERT-axis boundaries: I<->II (struct=0) and III<->IV (struct=1). Vary tau_ERT mean±1SD ---
print("\n[ERT-axis boundaries]  I<->II and III<->IV  (vary moisture-anomaly threshold mean±1SD)")
for k in [-1,-0.5,-0.25,0,0.25,0.5,1]:
    tau=mu+k*sd
    c=classify(dmg>1, pc>tau)
    ch=(c!=base).sum()
    print(f"  tau_ERT=mean{k:+.2f}SD : categories={dict(Counter(c))}  trees changed vs baseline={ch}")

# Decompose ERT instability by boundary, over mean±0.5SD
lo,hi=mu-0.5*sd, mu+0.5*sd
sound=dmg<=1; dmgd=dmg>1
inband=(pc>lo)&(pc<hi)
print(f"\nERT 'ambiguous' band = mean±0.5SD = ({lo:.2f},{hi:.2f})")
print(f"  trees in band overall: {inband.sum()}")
print(f"    of which structurally SOUND  (I<->II at risk): {(inband&sound).sum()}")
print(f"    of which structurally DAMAGED (III<->IV at risk): {(inband&dmgd).sum()}")
print(f"  structurally-damaged trees total: {dmgd.sum()}; their PC1 range: {pc[dmgd].min():.2f}..{pc[dmgd].max():.2f}")
print(f"  cavities (IV) at baseline: {(base=='IV').sum()}; active (III): {(base=='III').sum()}")

# explicit per-boundary flip counts over full ERT sweep and SoT sweep
def flips_between(a,b,base_c,alt_c):
    return int(((base_c==a)&(alt_c==b)|(base_c==b)&(alt_c==a)).sum())
print("\nPer-boundary flips (baseline vs a perturbed threshold):")
altERT=classify(dmg>1, pc>(mu+0.5*sd))
altSoT=classify(dmg>5, pc>mu)
print(f"  I<->II  (ERT axis, +0.5SD): {flips_between('I','II',base,altERT)}")
print(f"  III<->IV(ERT axis, +0.5SD): {flips_between('III','IV',base,altERT)}")
print(f"  II<->III(SoT axis, 1->5%) : {flips_between('II','III',base,altSoT)}")
print(f"  I<->IV  (SoT axis, 1->5%) : {flips_between('I','IV',base,altSoT)}")

# ======== FIGURE: phase diagram with threshold BANDS ========
fig,ax=plt.subplots(figsize=(8.2,6.4))
ysqrt=np.sqrt(np.clip(dmg,0,None))
# bands
ax.axvspan(lo,hi,color='#f39c12',alpha=0.13,zorder=0)              # ERT ambiguous band (wide)
ax.axhspan(np.sqrt(1),np.sqrt(5),color='#7f8c8d',alpha=0.18,zorder=0) # SoT ambiguous band (thin, empty)
ax.axvline(mu,color='#c0392b',ls='--',lw=1.4,zorder=2)
ax.axhline(np.sqrt(1),color='#333',ls='--',lw=1.4,zorder=2)
cols={'A.rubrum':'#2c7fb8','T.canadensis':'#41b6c4','N.sylvatica':'#7fcdbb','Q.rubra':'#d95f0e'}
mk={'A.rubrum':'^','T.canadensis':'o','N.sylvatica':'s','Q.rubra':'D'}
for s in cols:
    idx=m.sp==s
    ax.scatter(pc[idx],ysqrt[idx]+np.random.normal(0,0.03,idx.sum()),marker=mk[s],c=cols[s],s=44,edgecolor='white',lw=.5,label=s,zorder=3)
ax.set_xlabel('ERT PC1  (moisture axis — CONTINUUM, no natural break)',fontsize=10)
ax.set_ylabel('SoT structural loss %  (√ scale — BIMODAL, gap at 0→6%)',fontsize=10)
yt=[0,1,5,10,20,30]; ax.set_yticks([np.sqrt(t) for t in yt]); ax.set_yticklabels(yt)
# quadrant labels
ax.text(pc.min()+0.2,np.sqrt(28),'IV: Cavity',fontsize=10,color='#7f0000',weight='bold')
ax.text(pc.max()-1.4,np.sqrt(28),'III: Active',fontsize=10,color='#7f4f00',weight='bold')
ax.text(pc.min()+0.2,0.05,'I: No decay',fontsize=10,color='#26456e',weight='bold')
ax.text(pc.max()-1.6,0.05,'II: Incipient',fontsize=10,color='#4d7f00',weight='bold')
ax.text(mu, np.sqrt(34),'moisture threshold\n(wide band = many trees flip)',fontsize=7.3,ha='center',color='#c0392b')
ax.text(pc.max()-0.1, np.sqrt(3), 'structural threshold\n(band empty = stable)',fontsize=7.3,ha='right',color='#333')
ax.legend(fontsize=8,frameon=False,loc='center right')
ax.set_title('Boundary stability: the structural (SoT) split is crisp; the moisture (ERT) split is fuzzy',fontsize=10.2)
plt.tight_layout(); plt.savefig('analysis/revision/output/CJFR-boundary-fig.png',dpi=170)
print("\nsaved figure")
