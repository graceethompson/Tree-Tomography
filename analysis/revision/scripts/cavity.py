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
Zc=Z.values-Z.values.mean(0); U,S,Vt=np.linalg.svd(Zc,full_matrices=False)
pc=(U*S)[:,0]
if np.corrcoef(pc,m['mean'])[0,1]>0: pc=-pc
m['pc1']=pc; mu=pc.mean()
dm=m[m['percent_damaged']>1].copy()
dm['base']=np.where(dm['pc1']>mu,'III active','IV cavity')

# Physically-motivated cavity signals among damaged trees:
# (a) crack detected  (b) dry core: low CMA (few low-resistivity pixels in inner third)
print("Damaged trees — is there a physical 'dry cavity' signature?")
print(f"CMA (central moisture accumulation; LOW = dry core): III active mean={dm[dm.base=='III active'].cma.mean():.2f}, "
      f"IV cavity mean={dm[dm.base=='IV cavity'].cma.mean():.2f}")
# alternative dry-core cavity rule
dm['dry_core']=dm['cma']<0.15
dm['crack']=dm['crack_detected']=='yes'
print("\nProposed physical cavity rule = structural loss AND (crack OR dry core CMA<0.15):")
dm['cav_phys']=dm['crack']|dm['dry_core']
print(dm[['tree','sp','percent_damaged','crack','cma','dry_core','pc1','base','cav_phys']].sort_values('pc1').to_string(index=False))
print(f"\nComposite (PC1<mean) cavities: {(dm.base=='IV cavity').sum()}  | physical (crack/dry-core) cavities: {dm.cav_phys.sum()}")
print("Agreement between the two cavity definitions:", (( dm.base=='IV cavity')==dm.cav_phys).sum(),"/",len(dm))

# gap check on PC1 among damaged trees
s=np.sort(dm.pc1.values); g=np.diff(s)
print(f"\nGap at the III/IV boundary: cavities at PC1<= {s[2]:.2f}; next active at {s[3]:.2f}; "
      f"threshold mean={mu:.2f} sits in a {s[3]-s[2]:.2f}-wide gap")

# ---------- figure ----------
fig,ax=plt.subplots(figsize=(9,5.2))
cols={'A.rubrum':'#2c7fb8','T.canadensis':'#41b6c4','N.sylvatica':'#7fcdbb','Q.rubra':'#d95f0e'}
for _,r in dm.iterrows():
    ax.scatter(r.pc1,r.cma,s=90+r.percent_damaged*10,c=cols[r.sp],
               edgecolor='k' if r.crack else 'white',lw=2 if r.crack else .6,zorder=3,alpha=.9)
    ax.annotate(f"{int(r.tree)}",(r.pc1,r.cma),fontsize=6.5,xytext=(4,4),textcoords='offset points')
ax.axvline(mu,color='#c0392b',ls='--',lw=1.3); ax.text(mu+0.05,0.02,'III/IV threshold\n(PC1 mean)',fontsize=7,color='#c0392b')
ax.axhspan(0,0.15,color='#8e6b3f',alpha=0.10); ax.text(4.6,0.075,'dry core\n(cavity-like)',fontsize=7,ha='right',color='#6b4f2f')
ax.set_xlabel('ERT PC1  (low = drier/normal → cavity side; high = wet → active)')
ax.set_ylabel('CMA — central moisture accumulation\n(low = dry core)')
ax.set_title('The 15 structurally-damaged trees: cavity vs active\n(marker size = % structural loss; black outline = crack detected)',fontsize=10)
from matplotlib.lines import Line2D
leg=[Line2D([0],[0],marker='o',color='w',markerfacecolor=cols[s],markersize=8,label=s) for s in cols]
leg.append(Line2D([0],[0],marker='o',color='w',markerfacecolor='grey',markeredgecolor='k',markeredgewidth=2,markersize=8,label='crack detected'))
ax.legend(handles=leg,fontsize=7.5,frameon=False,loc='upper left')
plt.tight_layout(); plt.savefig('analysis/revision/output/CJFR-cavity-fig.png',dpi=170)
print("saved figure")
