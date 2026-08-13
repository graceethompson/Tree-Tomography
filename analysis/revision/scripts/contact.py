import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from PIL import Image
import pandas as pd, numpy as np, glob
df=pd.read_csv('./data/Tree_ID_info.csv')
sm={'rm':'A.rubrum','bg':'N.sylvatica','ro':'Q.rubra','hem':'T.canadensis'}; df['sp']=df['species'].map(sm)
ERT='./images/main_ERT'
def find(tid):
    g=[x for x in sorted(glob.glob(f'{ERT}/{tid}_*.jpg')) if 'labtest' not in x]
    return g[0] if g else None
rows=df.sort_values('tree').reset_index(drop=True)
def sheet(sub,fname,title):
    n=len(sub); ncol=6; nrow=(n+ncol-1)//ncol
    fig,axes=plt.subplots(nrow,ncol,figsize=(ncol*2.4,nrow*2.6))
    axes=np.atleast_2d(axes)
    for i,(_,r) in enumerate(sub.iterrows()):
        ax=axes.flat[i]; p=find(r['tree'])
        if p: ax.imshow(Image.open(p).crop((60,30,470,500)))
        ax.set_title(f"#{r['tree']} {r['sp'][:4]}\nSoT {r['percent_damaged']}%",fontsize=8)
        ax.axis('off')
    for j in range(n,nrow*ncol): axes.flat[j].axis('off')
    fig.suptitle(title,fontsize=12,y=1.0)
    plt.tight_layout(); plt.savefig(fname,dpi=95,bbox_inches='tight'); print('saved',fname,n,'trees')
sheet(rows.iloc[:30],'analysis/revision/output/contact1.png','ERT tomograms 1/2  (blue=wet, red=dry)')
sheet(rows.iloc[30:],'analysis/revision/output/contact2.png','ERT tomograms 2/2  (blue=wet, red=dry)')
