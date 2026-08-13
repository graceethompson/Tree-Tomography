import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
from PIL import Image
import glob, os

ERT='./images/main_ERT'; SOT='./images/main_SoT'
def find(d,tid):
    g=sorted(glob.glob(f'{d}/{tid}_*.jpg'))
    g=[x for x in g if '_2.jpg' not in x and '_jon' not in x and '30ii' not in x]
    return g[0]

# tree: (id, category-with-color, species, %dmg, mean-res, CMA, PC1, one-line visual read)
trees=[
 (544,'I — No decay','#26456e','Q. rubra',0,672,0.33,-2.15,'dry (red) throughout, no void → sound'),
 (418,'II — Incipient','#4d7f00','A. rubrum',0,135,0.41,2.44,'wet (blue) interior but NO structural loss → incipient/wetwood'),
 (217,'III — Active','#7f4f00','T. canadensis',11,128,0.24,5.01,'wet (blue) interior WITH structural loss → active decay'),
 (880,'IV — Cavity','#7f0000','Q. rubra',30,408,0.10,-0.42,'DRY (red) core + large SoT void → open/desiccated cavity'),
 (420,'IV — Cavity','#7f0000','A. rubrum',16,337,0.04,-1.08,'dry (red) central mass + crack → cavity'),
 (893,'III vs IV (borderline)','#555','T. canadensis',12,454,0.40,-0.30,'composite labels "cavity", but ERT shows a WET blue core → really Active'),
 (380,'III — Active (not cavity)','#7f4f00','A. rubrum',34,81,0.05,3.95,'cracked & 34% loss, but ERT uniformly WET → Active, not a cavity'),
]

n=len(trees); ncol=2; nrow=(n+ncol-1)//ncol
fig=plt.figure(figsize=(15, 4.7*nrow))
outer=gridspec.GridSpec(nrow,ncol,figure=fig,hspace=0.30,wspace=0.12)
for i,(tid,cat,col,sp,dmg,res,cma,pc1,read) in enumerate(trees):
    r,c=divmod(i,ncol)
    inner=gridspec.GridSpecFromSubplotSpec(1,2,subplot_spec=outer[r,c],wspace=0.04)
    for j,(path,tag) in enumerate([(find(SOT,tid),'SoT (structure)'),(find(ERT,tid),'ERT (moisture)')]):
        ax=fig.add_subplot(inner[j]); ax.imshow(Image.open(path)); ax.axis('off')
        ax.set_title(tag,fontsize=8,color='#333')
    axt=fig.add_subplot(outer[r,c]); axt.axis('off')
    axt.set_title(f"Tree {tid}  ·  {cat}\n{sp}  ·  SoT loss {dmg}%  ·  mean ρ {res} Ω·m  ·  CMA {cma}  ·  PC1 {pc1:+.2f}\n“{read}”",
                  fontsize=9.2,color=col,fontweight='bold',pad=26,linespacing=1.5)
fig.suptitle('What the categories look like: paired SoT (structural loss) + ERT (moisture) tomograms\n'
             'Cavity = DRY red core + void   ·   Active = WET blue interior + loss   ·   Incipient = WET interior, no loss   ·   No-decay = DRY, no loss',
             fontsize=12,y=0.995)
plt.savefig('analysis/revision/output/CJFR-scan-montage.png',dpi=115,bbox_inches='tight')
print('saved montage with',n,'trees')
