# Per-tree category assignments under each classification scheme.
# Writes analysis/revision/tables/scheme_assignments.csv. Run from repo root.
import pandas as pd, numpy as np

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
m['pc1'] = pc.round(3)
m['struct'] = m['percent_damaged'] > 1
m['xdev'] = (-m.groupby('sp')['mean'].transform(lambda x: (x - x.median()) / x.std(ddof=1))).round(3)

def four(anom, struct):
    return np.where(~struct & ~anom, 'I', np.where(~struct & anom, 'II',
           np.where(struct & anom, 'III', 'IV')))

m['cat_pc1'] = four(pc > pc.mean(), m['struct'].values)
m['cat_species_median'] = four((m['xdev'] > 0).values, m['struct'].values)
m['cat_absolute'] = four(m['mean'].values < np.median(m['mean']), m['struct'].values)
m['cat_6cell'] = np.where(~m['struct'] & (m['xdev'] > 0.5), 'II',
                 np.where(~m['struct'] & ~(m['xdev'] > 0.5), 'I',
                 np.where(m['struct'] & (m['xdev'] < -0.5), 'IV', 'III')))

out = m[['tree', 'sp', 'site', 'dbh', 'percent_damaged', 'mean', 'cma', 'pc1', 'xdev',
         'cat_pc1', 'cat_species_median', 'cat_absolute', 'cat_6cell']]
out.to_csv('analysis/revision/output/scheme_assignments.csv', index=False)
print(f"wrote {len(out)} rows; 3-way identical (pc1/median/absolute): "
      f"{((out.cat_pc1==out.cat_species_median)&(out.cat_species_median==out.cat_absolute)).sum()}/57")
