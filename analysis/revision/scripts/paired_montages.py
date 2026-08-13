# Paired SoT+ERT scan montages, grouped by classification scheme.
# Run from repo root. Outputs to analysis/revision/output/.
import pandas as pd, numpy as np, glob
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from PIL import Image

df = pd.read_csv('data/Tree_ID_info.csv'); df['site'] = np.where(df['plot'] == 'BGS', 'BGS', 'EMS')
sm = {'rm': 'A.rubrum', 'bg': 'N.sylvatica', 'ro': 'Q.rubra', 'hem': 'T.canadensis'}
df['sp'] = df['species'].map(sm)
ert = pd.read_csv('data/ERT_application_results.csv')
m = df.merge(ert, on='tree').reset_index(drop=True)
mets = ['mean', 'median', 'sd', 'cv', 'gini', 'entropy', 'cma', 'radialgradient']
Z = m[mets].copy()
for c in mets:
    Z[c] = m.groupby('sp')[c].transform(lambda x: (x - x.mean()) / x.std(ddof=1))
Zc = Z.values - Z.values.mean(0)
U, S, Vt = np.linalg.svd(Zc, full_matrices=False)
pc = (U * S)[:, 0]
if np.corrcoef(pc, m['mean'])[0, 1] > 0:
    pc = -pc
m['pc1'] = pc; m['struct'] = m['percent_damaged'] > 1
m['xdev'] = -m.groupby('sp')['mean'].transform(lambda x: (x - x.median()) / x.std(ddof=1))

ERT = 'images/main_ERT'; SOT = 'images/main_SoT'
def find(d, t):
    g = [x for x in sorted(glob.glob(f'{d}/{t}_*.jpg'))
         if all(s not in x for s in ('labtest', '_2.jpg', '_jon', '30ii', 'dummy', 'test'))]
    return g[0]

def four(anom, struct):
    return np.where(~struct & ~anom, 'I', np.where(~struct & anom, 'II',
           np.where(struct & anom, 'III', 'IV')))

catinfo = [('I', 'No decay', '#3b6fb0'), ('II', 'Incipient', '#4e9a2c'),
           ('III', 'Active', '#d98a1f'), ('IV', 'Cavity', '#b83232')]
schemes = {
    'PC1': ('Bin by PC1 (published scheme)', four(pc > pc.mean(), m['struct'].values)),
    'spmed': ('Bin by within-species median split', four((m['xdev'] > 0).values, m['struct'].values)),
    '6cell': ('6-cell: species-median +/-0.5 SD dead band (defaults to row baseline)',
              np.where(~m['struct'] & (m['xdev'] > 0.5), 'II',
              np.where(~m['struct'] & ~(m['xdev'] > 0.5), 'I',
              np.where(m['struct'] & (m['xdev'] < -0.5), 'IV', 'III')))),
}
PER = 4
def montage(key):
    title, cats = schemes[key]; m['_c'] = cats
    plan = []
    for code, nm, col in catinfo:
        ids = m.loc[m._c == code, 'tree'].tolist()
        plan.append(('hdr', f"{code} — {nm}  (n={len(ids)})", col))
        for i in range(0, max(len(ids), 1), PER):
            plan.append(('img', ids[i:i + PER]))
    units = [0.32 if r[0] == 'hdr' else 1.0 for r in plan]; tot = sum(units)
    fig = plt.figure(figsize=(PER * 3.5, tot * 1.7)); ycur = 1.0
    for (kind, *rest), u in zip(plan, units):
        h = u / tot; ycur -= h
        if kind == 'hdr':
            txt, col = rest
            fig.patches.append(plt.Rectangle((0.01, ycur + 0.01), 0.98, h - 0.02,
                               transform=fig.transFigure, color=col, alpha=0.16))
            fig.text(0.02, ycur + h / 2, txt, fontsize=12, fontweight='bold', color=col, va='center')
        else:
            for k, t in enumerate(rest[0]):
                base = 0.01 + k * (0.98 / PER); cellw = 0.98 / PER
                row = m[m.tree == t].iloc[0]
                for jj, d in enumerate([SOT, ERT]):
                    ax = fig.add_axes([base + jj * (cellw * 0.48), ycur + 0.03, cellw * 0.45, h * 0.74])
                    ax.imshow(Image.open(find(d, t)).crop((60, 30, 470, 500))); ax.axis('off')
                fig.text(base + cellw * 0.48, ycur + h * 0.86,
                         f"#{t} {row['sp'][:4]} SoT{row['percent_damaged']} ρ{row['mean']:.0f}",
                         fontsize=6.6, ha='center')
    fig.suptitle(title + '   —   SoT (left) + ERT (right) per tree', fontsize=13, y=1.0)
    out = f'analysis/revision/output/CJFR-paired-{key}.png'
    plt.savefig(out, dpi=80, bbox_inches='tight'); plt.close(); print('saved', out)

for k in schemes:
    montage(k)
