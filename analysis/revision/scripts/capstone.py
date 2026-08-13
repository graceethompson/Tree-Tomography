import pandas as pd, numpy as np
df=pd.read_csv('./data/Tree_ID_info.csv'); df['site']=np.where(df['plot']=='BGS','BGS','EMS')
sm={'rm':'A.rubrum','bg':'N.sylvatica','ro':'Q.rubra','hem':'T.canadensis'}; df['sp']=df['species'].map(sm)
ert=pd.read_csv('./data/ERT_application_results.csv'); m=df.merge(ert,on='tree').reset_index(drop=True)
mets=['mean','median','sd','cv','gini','entropy','cma','radialgradient']
Z=m[mets].copy()
for c in mets: Z[c]=m.groupby('sp')[c].transform(lambda x:(x-x.mean())/x.std(ddof=1))
Zc=Z.values-Z.values.mean(0); U,S,Vt=np.linalg.svd(Zc,full_matrices=False); pc=(U*S)[:,0]
if np.corrcoef(pc,m['mean'])[0,1]>0: pc=-pc
m['pc1']=pc; m['struct']=m['percent_damaged']>1
mb=(m.site=='BGS').values&(~m['struct'].values); me=(m.site=='EMS').values&(~m['struct'].values)
nb=(~m['struct'].values&(m.site=='BGS').values).sum(); ne=(~m['struct'].values&(m.site=='EMS').values).sum()
# species-standardized wetness (higher=wetter)
xdev=-m.groupby('sp')['mean'].transform(lambda x:(x-x.median())/x.std(ddof=1)).values
# pooled PC1
Zp=(m[mets]-m[mets].mean())/m[mets].std(ddof=1); Zp=Zp.values-Zp.values.mean(0)
Up,Sp,Vtp=np.linalg.svd(Zp,full_matrices=False); pcp=(Up*Sp)[:,0]
if np.corrcoef(pcp,m['mean'])[0,1]>0: pcp=-pcp

sd=pc.std(); sdx=xdev.std(); sdp=pcp.std()
# each axis: (values array where higher=more anomalous/wet, list of thresholds spanning -1..+1 SD around its center)
axes={
 'PC1 (within-sp)': (pc, [pc.mean()+k*sd for k in (-1,-0.5,-0.25,0,0.25,0.5,1)]),
 'PC1 (pooled)':    (pcp,[pcp.mean()+k*sdp for k in (-1,-0.5,-0.25,0,0.25,0.5,1)]),
 'species-median resistivity': (xdev,[0+k*sdx for k in (-1,-0.5,-0.25,0,0.25,0.5,1)]),
 'absolute resistivity': (-m['mean'].values, sorted(-np.quantile(m['mean'].values,[0.35,0.42,0.5,0.58,0.65]))),
 'CMA': (m['cma'].values, [0.25,0.30,0.33,0.36,0.40]),
}
print("SITE DIRECTION ROBUSTNESS: incipient among structurally-sound trees, BGS vs EMS")
print(f"{'axis':28s}{'#thr':>5}{'BGS>EMS holds':>15}{'BGS% range':>16}{'EMS% range':>16}")
total=0; held=0
for name,(v,ths) in axes.items():
    bs=[]; es=[]; ok=0
    for t in ths:
        anom=v>t
        b=(anom&mb).sum()/nb*100; e=(anom&me).sum()/ne*100
        bs.append(b); es.append(e); ok+= (b>=e)
    total+=len(ths); held+=ok
    print(f"{name:28s}{len(ths):>5}{ok:>8}/{len(ths):<6}{min(bs):>6.0f}-{max(bs):<9.0f}{min(es):>6.0f}-{max(es):<9.0f}")
print(f"\nACROSS ALL {total} axis×threshold combinations: BGS>=EMS incipient in {held}/{total} = {held/total*100:.0f}%")

# structural prevalence stability by SoT threshold
print("\nStructural-loss prevalence by site, across SoT thresholds:")
for tau in [1,2,3,5]:
    s=m['percent_damaged']>tau
    b=s[m.site=='BGS'].mean()*100; e=s[m.site=='EMS'].mean()*100
    print(f"  SoT>{tau}%: BGS {b:.0f}% vs EMS {e:.0f}%  (EMS>BGS: {e>b})")

# species severity ordering stability (mean % damage among decayed, by species) — bootstrap-free descriptive
print("\nSpecies severity (mean % damage among decayed trees) — ordering:")
dd=m[m.percent_damaged>0]
print("  ", dd.groupby('sp')['percent_damaged'].mean().round(1).sort_values(ascending=False).to_dict())
print("  decayed n per species:", dd.groupby('sp').size().to_dict())
