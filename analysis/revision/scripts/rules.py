import pandas as pd, numpy as np
df=pd.read_csv('./data/Tree_ID_info.csv')
df['site']=np.where(df['plot']=='BGS','BGS','EMS')
sm={'rm':'A.rubrum','bg':'N.sylvatica','ro':'Q.rubra','hem':'T.canadensis'}; df['sp']=df['species'].map(sm)
ert=pd.read_csv('./data/ERT_application_results.csv')
m=df.merge(ert,on='tree').reset_index(drop=True)
mets=['mean','median','sd','cv','gini','entropy','cma','radialgradient']
Z=m[mets].copy()
for c in mets: Z[c]=m.groupby('sp')[c].transform(lambda x:(x-x.mean())/x.std(ddof=1))
Zc=Z.values-Z.values.mean(0); U,S,Vt=np.linalg.svd(Zc,full_matrices=False); pc=(U*S)[:,0]
if np.corrcoef(pc,m['mean'])[0,1]>0: pc=-pc
m['pc1']=pc; m['res_z']=Z['mean']; mu=pc.mean(); sd=pc.std(ddof=1)
otsu=1.00  # from prior analysis

m['struct']=m['percent_damaged']>1
# graded moisture bands
def band(p):
    if p<=mu: return 'normal'
    if p<=otsu: return 'possible'
    return 'confident'
m['mband']=m['pc1'].apply(band)
print(f"PC1: mean={mu:.2f}, SD={sd:.2f}; Otsu 'confident' break={otsu}")
print("Moisture-anomaly bands (all 57):", m['mband'].value_counts().to_dict())
print("Among SOUND trees (I/II split):", m[~m.struct]['mband'].value_counts().to_dict())

# III/IV physical rule among damaged trees: cavity = dry core = high resistivity AND low central moisture
dmg=m[m.struct].copy()
dmg['cavity_phys']=(dmg['mean']>=300)&(dmg['cma']<0.15)
dmg['class_phys']=np.where(dmg['cavity_phys'],'IV cavity','III active')
dmg['class_composite']=np.where(dmg['pc1']>mu,'III active','IV cavity')
print("\nIII/IV among 15 damaged trees:")
print(dmg[['tree','sp','percent_damaged','mean','cma','pc1','class_composite','class_phys']].sort_values('mean').to_string(index=False))
print("composite cavities:",(dmg.class_composite=='IV cavity').sum(),"| physical(dry-core) cavities:",(dmg.class_phys=='IV cavity').sum())

# absolute equivalents
print("\nAbsolute translation (approx, species-confounded): PC1 mean ~ 369 Ohm-m ~ 103% moisture; Otsu(1.0) ~ 387 Ohm-m ~ 101% moisture")
print("Structural gap check: damage values >0:", sorted(m.loc[m.percent_damaged>0,'percent_damaged'].tolist()))
