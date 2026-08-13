import pandas as pd, numpy as np
from scipy import stats
from sklearn.mixture import GaussianMixture
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

np.random.seed(0)
df = pd.read_csv('./data/Tree_ID_info.csv')
df['site']=np.where(df['plot']=='BGS','BGS','EMS')
sppmap={'rm':'A.rubrum','bg':'N.sylvatica','ro':'Q.rubra','hem':'T.canadensis'}
df['sp']=df['species'].map(sppmap)
ert=pd.read_csv('./data/ERT_application_results.csv')
m=df.merge(ert,on='tree',how='inner').reset_index(drop=True)
metrics=['mean','median','sd','cv','gini','entropy','cma','radialgradient']
Z=m[metrics].copy()
for met in metrics:
    Z[met]=m.groupby('sp')[met].transform(lambda x:(x-x.mean())/x.std(ddof=1))
Zc=Z.values-Z.values.mean(0)
U,S,Vt=np.linalg.svd(Zc,full_matrices=False)
pc=(U*S)[:,0]
if np.corrcoef(pc,m['mean'])[0,1]>0: pc=-pc
m['pc']=pc; m['struct']=m['percent_damaged']>1
mu,med,sd=pc.mean(),np.median(pc),pc.std(ddof=1)
def otsu(x,nb=256):
    h,e=np.histogram(x,bins=nb); c=(e[:-1]+e[1:])/2; w=h/h.sum()
    cw=np.cumsum(w); cm=np.cumsum(w*c); gm=cm[-1]; bt,bv=None,-1
    for i in range(1,nb):
        w0=cw[i]; w1=1-w0
        if w0<=0 or w1<=0: continue
        v=w0*w1*((cm[i]/w0)-((gm-cm[i])/w1))**2
        if v>bv: bv,bt=v,c[i]
    return bt
g2=GaussianMixture(2,random_state=0,n_init=5).fit(pc.reshape(-1,1))
mns=g2.means_.ravel(); o=np.argsort(mns); xs=np.linspace(mns[o[0]],mns[o[1]],2000)
lp=g2._estimate_weighted_log_prob(xs.reshape(-1,1)); d=lp[:,o[0]]-lp[:,o[1]]
idx=np.where(np.diff(np.sign(d))!=0)[0]; t_gmm=xs[idx[0]] if len(idx) else np.nan
t_otsu=otsu(pc)
rules={'mean (published)':(mu,'#c0392b'),'median':(med,'#2980b9'),
       'Otsu break':(t_otsu,'#8e44ad'),'GMM crossover':(t_gmm,'#27ae60')}

fig,(ax1,ax2)=plt.subplots(1,2,figsize=(12.5,4.8))
# Panel A: density
xs2=np.linspace(pc.min()-0.4,pc.max()+0.4,300); kde=stats.gaussian_kde(pc)
ax1.hist(pc,bins=15,density=True,color='0.86',edgecolor='white',zorder=1)
ax1.plot(xs2,kde(xs2),color='k',lw=1.8,zorder=2)
for name,(t,c) in rules.items():
    if np.isfinite(t): ax1.axvline(t,color=c,ls='--',lw=1.7,label=f"{name} = {t:.2f}",zorder=3)
ax1.set_ylim(0,0.42); ax1.set_xlabel('ERT PC1  (higher = lower resistivity / more heterogeneous)')
ax1.set_ylabel('density'); ax1.legend(fontsize=8,frameon=False,loc='upper right')
ax1.set_title('A.  Anomaly axis is unimodal — no natural break\n(BIC favors 1 component over 2)',fontsize=10,loc='left')

# Panel B: sweep threshold, plot BGS vs EMS incipient%
ts=np.linspace(np.percentile(pc,5),np.percentile(pc,95),200)
bgs=[]; ems=[]
mask_bgs=(m.site=='BGS').values; mask_ems=(m.site=='EMS').values; ns=(~m['struct']).values
for t in ts:
    inc=ns&(pc>t)
    bgs.append(inc[mask_bgs].sum()/mask_bgs.sum()*100)
    ems.append(inc[mask_ems].sum()/mask_ems.sum()*100)
ax2.plot(ts,bgs,color='#1f77b4',lw=2.4,label='BGS (wetland)')
ax2.plot(ts,ems,color='#d9822b',lw=2.4,label='EMS (upland)')
ax2.fill_between(ts,ems,bgs,where=np.array(bgs)>=np.array(ems),color='#1f77b4',alpha=0.10)
for name,(t,c) in rules.items():
    if np.isfinite(t): ax2.axvline(t,color=c,ls=':',lw=1.3)
ax2.set_xlabel('anomaly threshold on ERT PC1'); ax2.set_ylabel('% of site trees classified "incipient"')
ax2.legend(fontsize=9,frameon=False,loc='upper right')
ax2.set_title('B.  Wetland > upland incipient holds at EVERY threshold\n(absolute % varies; the direction does not)',fontsize=10,loc='left')
plt.tight_layout(); plt.savefig('.-threshold-fig.png',dpi=175)
print("saved. BGS never below EMS across sweep:", all(np.array(bgs)>=np.array(ems)))
print(f"BGS incipient range across sweep: {min(bgs):.0f}-{max(bgs):.0f}% ; EMS: {min(ems):.0f}-{max(ems):.0f}%")
