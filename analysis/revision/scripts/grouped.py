import pandas as pd, numpy as np, glob
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from PIL import Image
df=pd.read_csv('./data/Tree_ID_info.csv'); df['site']=np.where(df['plot']=='BGS','BGS','EMS')
sm={'rm':'A.rubrum','bg':'N.sylvatica','ro':'Q.rubra','hem':'T.canadensis'}; df['sp']=df['species'].map(sm)
ert=pd.read_csv('./data/ERT_application_results.csv'); m=df.merge(ert,on='tree').reset_index(drop=True)
mets=['mean','median','sd','cv','gini','entropy','cma','radialgradient']
Z=m[mets].copy()
for c in mets: Z[c]=m.groupby('sp')[c].transform(lambda x:(x-x.mean())/x.std(ddof=1))
Zc=Z.values-Z.values.mean(0); U,S,Vt=np.linalg.svd(Zc,full_matrices=False); pc=(U*S)[:,0]
if np.corrcoef(pc,m['mean'])[0,1]>0: pc=-pc
m['pc1']=pc; mu=pc.mean(); m['struct']=m['percent_damaged']>1
m['res_dev']=m.groupby('sp')['mean'].transform(lambda x:(x-x.median())/x.std(ddof=1))
def four(a): return np.where(~m.struct&~a,'I',np.where(~m.struct&a,'II',np.where(m.struct&a,'III','IV')))
schemes={
 'A_PC1':('Bin by PC1 (original)', pc>mu),
 'B_CMA':('Bin by CMA > 0.33', m['cma'].values>0.33),
 'C_absres':('Bin by absolute resistivity < sample median', m['mean'].values<np.median(m['mean'])),
 'D_spmed':('Bin by within-species median split (wetter than species median)', (m['res_dev']<0).values),
}
ERT='./images/main_ERT'
def find(t):
    g=[x for x in sorted(glob.glob(f'{ERT}/{t}_*.jpg')) if 'labtest' not in x]; return g[0]
catinfo=[('I','No decay','#3b6fb0'),('II','Incipient','#4e9a2c'),('III','Active','#d98a1f'),('IV','Cavity','#b83232')]
NC=8
def montage(key):
    title,anom=schemes[key]; cats=four(anom); m['_c']=cats
    # build row plan
    plan=[]  # ('hdr',text,color) or ('img',[tids])
    for code,nm,col in catinfo:
        ids=m.loc[m._c==code,'tree'].tolist()
        plan.append(('hdr',f"{code} — {nm}   (n={len(ids)})",col))
        for i in range(0,max(len(ids),1),NC): plan.append(('img',ids[i:i+NC]))
    units=[0.34 if r[0]=='hdr' else 1.0 for r in plan]; tot=sum(units)
    fig=plt.figure(figsize=(NC*1.7, tot*1.55))
    ycur=1.0
    for (kind,*rest),u in zip(plan,units):
        h=u/tot; ytop=ycur; ycur-=h
        if kind=='hdr':
            txt,col=rest
            fig.patches.append(plt.Rectangle((0.01,ycur+0.01),0.98,h-0.02,transform=fig.transFigure,color=col,alpha=0.16,zorder=0))
            fig.text(0.02,ycur+h/2,txt,fontsize=12,fontweight='bold',color=col,va='center')
        else:
            ids=rest[0]
            for k,t in enumerate(ids):
                x=0.01+k*(0.98/NC); w=0.98/NC*0.94; hh=h*0.80
                ax=fig.add_axes([x,ycur+0.02,w,hh]); ax.imshow(Image.open(find(t)).crop((60,30,470,500))); ax.axis('off')
                row=m[m.tree==t].iloc[0]
                ax.set_title(f"#{t} {row['sp'][:4]}\nSoT{row['percent_damaged']} ρ{row['mean']:.0f}",fontsize=6.3)
    fig.suptitle(title+'   —   ERT scans grouped by resulting category',fontsize=13,y=1.0)
    out=f'analysis/revision/output/CJFR-grouped-{key}.png'; plt.savefig(out,dpi=85,bbox_inches='tight'); plt.close(); print('saved',out)
for k in schemes: montage(k)
