import pandas as pd, numpy as np
from scipy import stats
from sklearn.mixture import GaussianMixture
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

np.random.seed(0)
df=pd.read_csv('./data/Tree_ID_info.csv')
df['site']=np.where(df['plot']=='BGS','BGS','EMS')
sppmap={'rm':'A.rubrum','bg':'N.sylvatica','ro':'Q.rubra','hem':'T.canadensis'}
df['sp']=df['species'].map(sppmap)
ert=pd.read_csv('./data/ERT_application_results.csv')
m=df.merge(ert,on='tree',how='inner').reset_index(drop=True)
m['res']=m['mean']            # mean resistivity, Ohm-m (absolute)
m['struct']=m['percent_damaged']>1
ns=(~m['struct']).values
mb=(m.site=='BGS').values; me=(m.site=='EMS').values

# ---- validation calibration: moisture ~ resistivity (n=12) ----
v=pd.read_csv('./data/hemlock/validation_summary.csv')
v['res']=1000/v['Conductance']
sl=stats.linregress(v['res'],v['moisture'])
n=len(v); xbar=v['res'].mean(); Sxx=((v['res']-xbar)**2).sum()
resid=v['moisture']-(sl.intercept+sl.slope*v['res']); rstd=np.sqrt((resid**2).sum()/(n-2))
def pred_moist(x): return sl.intercept+sl.slope*np.asarray(x,float)
def pred_se(x):    return rstd*np.sqrt(1+1/n+(np.asarray(x,float)-xbar)**2/Sxx)   # prediction SE
def res_for_moist(Mcut):  # invert: resistivity giving predicted moisture=Mcut, with band
    x=(Mcut-sl.intercept)/sl.slope
    # approximate inverse band via slope
    band=pred_se(x)/abs(sl.slope)
    return x, x-1.96*band, x+1.96*band
m['pmoist']=pred_moist(m['res'])
r100=res_for_moist(100)[0]; print(f"Calibration: moisture(%) = {sl.intercept:.1f} {sl.slope:+.4f}*res   r={sl.rvalue:.2f}, n={n}, resid s={rstd:.1f}%")
print(f"Predicted moisture range across 57 trees: {m['pmoist'].min():.0f}-{m['pmoist'].max():.0f}%")

print("\n--- Absolute mean resistivity (Ohm-m) by species & site ---")
print(m.groupby('sp')['res'].agg(['count','mean','min','max']).round(0))
print(m.groupby('site')['res'].agg(['mean','median']).round(0).to_dict())

# ---- modality on absolute resistivity ----
def modality(x,label):
    x=np.asarray(x,float).reshape(-1,1)
    b1=GaussianMixture(1,random_state=0).fit(x).bic(x)
    b2=GaussianMixture(2,random_state=0,n_init=5).fit(x).bic(x)
    print(f"{label}: BIC 1={b1:.1f} 2={b2:.1f} -> {'bimodal' if b2+2<b1 else 'unimodal'}; skew={stats.skew(x.ravel()):.2f}")
modality(m['res'],'abs resistivity')
modality(m['pmoist'],'pred moisture')

# ---- candidate ABSOLUTE thresholds (anomaly = LOW resistivity / HIGH moisture) ----
def otsu(x,nb=256,invert=False):
    h,e=np.histogram(x,bins=nb); c=(e[:-1]+e[1:])/2; w=h/h.sum()
    cw=np.cumsum(w); cm=np.cumsum(w*c); gm=cm[-1]; bt,bv=None,-1
    for i in range(1,nb):
        w0=cw[i]; w1=1-w0
        if w0<=0 or w1<=0: continue
        val=w0*w1*((cm[i]/w0)-((gm-cm[i])/w1))**2
        if val>bv: bv,bt=val,c[i]
    return bt
res=m['res'].values
t_otsu=otsu(res)
g2=GaussianMixture(2,random_state=0,n_init=5).fit(res.reshape(-1,1))
mns=np.sort(g2.means_.ravel()); xs=np.linspace(mns[0],mns[1],2000)
o=np.argsort(g2.means_.ravel()); lp=g2._estimate_weighted_log_prob(xs.reshape(-1,1))
dd=lp[:,o[0]]-lp[:,o[1]]; ix=np.where(np.diff(np.sign(dd))!=0)[0]; t_gmm=xs[ix[0]] if len(ix) else np.nan

# absolute anomaly = res < threshold
rules=[
 ('Sample median resistivity', np.median(res)),
 ('Sample mean resistivity', res.mean()),
 ('Otsu break (data-driven)', t_otsu),
 ('GMM crossover (data-driven)', t_gmm),
 ('Lit: silver-fir wetwood ~200', 200.0),
 ('Lit: silver-fir low-ER ~166', 166.0),
 ('Validation: moisture>90% (<505)', res_for_moist(90)[0]),
 ('Validation: moisture>100% (<402)', res_for_moist(100)[0]),
 ('Validation: moisture>110% (<299)', res_for_moist(110)[0]),
]
print("\n"+"="*92)
print(f"{'ABSOLUTE threshold rule':<34}{'res<':>7}{'~moist>':>9}{'nAnom':>7}{'nIncip':>7}{'BGSinc%':>9}{'EMSinc%':>9}")
print("="*92)
tbl=[]
for name,t in rules:
    anom=res<t
    inc=ns&anom
    bgs=inc[mb].sum()/mb.sum()*100; ems=inc[me].sum()/me.sum()*100
    mo=pred_moist(t)
    print(f"{name:<34}{t:>7.0f}{mo:>9.0f}{anom.sum():>7}{inc.sum():>7}{bgs:>9.0f}{ems:>9.0f}")
    tbl.append(dict(rule=name,res_thresh=round(t),moist_equiv=round(mo),n_anom=int(anom.sum()),
                    n_incip=int(inc.sum()),BGS_inc=round(bgs),EMS_inc=round(ems)))
pd.DataFrame(tbl).to_csv('analysis/revision/output/CJFR-absolute-thresholds.csv',index=False)

# stability of anomaly label across absolute rule set
lab=np.vstack([(res<t).astype(int) for _,t in rules])
print(f"\nAlways-anomalous (all abs rules): {(lab.min(0)==1).sum()}, never: {(lab.max(0)==0).sum()}, unstable: {(lab.min(0)!=lab.max(0)).sum()}")

# does site DIRECTION hold on absolute axis?
grid=np.linspace(np.percentile(res,5),np.percentile(res,95),200)
bgs=[]; ems=[]
for t in grid:
    inc=ns&(res<t); bgs.append(inc[mb].sum()/mb.sum()*100); ems.append(inc[me].sum()/me.sum()*100)
bgs=np.array(bgs); ems=np.array(ems)
print(f"\nSite DIRECTION on absolute axis: BGS>=EMS at {np.mean(bgs>=ems)*100:.0f}% of thresholds; "
      f"BGS>EMS at {np.mean(bgs>ems)*100:.0f}%; crossovers where EMS>BGS at {np.mean(ems>bgs)*100:.0f}%")
print(f"BGS incip range {bgs.min():.0f}-{bgs.max():.0f}%, EMS {ems.min():.0f}-{ems.max():.0f}%")

# ---------------- FIGURE (3 panels) ----------------
fig,axes=plt.subplots(1,3,figsize=(16.5,4.9))
spp_order=['A.rubrum','T.canadensis','N.sylvatica','Q.rubra']
cols={'A.rubrum':'#2c7fb8','T.canadensis':'#41b6c4','N.sylvatica':'#7fcdbb','Q.rubra':'#d95f0e'}

# A: absolute resistivity by species (the confound made visible)
ax=axes[0]
for i,s in enumerate(spp_order):
    d=m.loc[m.sp==s,'res']
    ax.scatter(d, np.random.normal(i,0.08,len(d)), color=cols[s], s=34, alpha=.8, edgecolor='white',lw=.5)
    ax.plot([d.mean(),d.mean()],[i-.28,i+.28],color=cols[s],lw=2.5)
ax.set_yticks(range(4)); ax.set_yticklabels([s+('*' if s in('N.sylvatica','Q.rubra') else '') for s in spp_order])
ax.set_xlabel('mean resistivity (Ω·m)  — lower = wetter'); ax.invert_xaxis()
ax.set_title('A.  Absolute axis is species-confounded\n(* = single-site specialist)',fontsize=10,loc='left')

# B: predicted-moisture distribution by site with candidate thresholds
ax=axes[1]
for s,msk,c in [('BGS',mb,'#1f77b4'),('EMS',me,'#d9822b')]:
    xs2=np.linspace(m['pmoist'].min()-5,m['pmoist'].max()+5,300)
    ax.plot(xs2,stats.gaussian_kde(m.loc[msk,'pmoist'])(xs2),color=c,lw=2,label=s)
    ax.scatter(m.loc[msk,'pmoist'],np.full(msk.sum(),-0.001)+ (0.0004 if s=='EMS' else 0),color=c,s=14,alpha=.6)
for Mc,c in [(90,'#999'),(100,'#555'),(110,'#111')]:
    ax.axvline(Mc,color=c,ls=':',lw=1.2); ax.text(Mc,ax.get_ylim()[1]*0.92,f'{Mc}%',fontsize=7,ha='center')
ax.set_xlabel('predicted moisture content (%)  [from validation calibration]'); ax.set_ylabel('density')
ax.legend(fontsize=9,frameon=False); ax.set_title('B.  Absolute (predicted moisture) axis by site',fontsize=10,loc='left')

# C: site direction sweep on absolute resistivity, anchors marked
ax=axes[2]
ax.plot(grid,bgs,color='#1f77b4',lw=2.4,label='BGS (wetland)')
ax.plot(grid,ems,color='#d9822b',lw=2.4,label='EMS (upland)')
ax.fill_between(grid,ems,bgs,where=bgs>=ems,color='#1f77b4',alpha=.10)
ax.fill_between(grid,ems,bgs,where=ems>bgs,color='#d9822b',alpha=.18)
anchor={'median':np.median(res),'mean':res.mean(),'Otsu':t_otsu,'moist>100%':r100,'wetwood~200':200}
for nm,t in anchor.items():
    ax.axvline(t,color='#444',ls=':',lw=1.0); ax.text(t,ax.get_ylim()[1]*0.98,nm,rotation=90,va='top',ha='right',fontsize=6.5,color='#444')
ax.set_xlabel('anomaly threshold: resistivity below … (Ω·m)'); ax.set_ylabel('% of site trees "incipient"')
ax.legend(fontsize=9,frameon=False,loc='upper left'); ax.invert_xaxis()
ax.set_title('C.  Site direction across absolute thresholds\n(shaded orange = where upland exceeds wetland)',fontsize=10,loc='left')
plt.tight_layout(); plt.savefig('analysis/revision/output/CJFR-absolute-fig.png',dpi=170)
print("\nsaved figure + CSV")
