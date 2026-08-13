# Statistical tests supporting the revision:
#  1. Fisher exact tests on the incipient site contrast (by threshold; generalists-only)
#  2. Calibration-anchor uncertainty: CI on the mean calibration line vs single-tree PI
# Run from repo root.
import pandas as pd, numpy as np
from scipy import stats

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
sd = pc.std(ddof=1); sound = ~m['struct']

print("1. Fisher exact: incipient site contrast by ERT threshold")
for k, lab in [(-0.5, 'mean-0.5SD'), (0, 'mean (published)'), (0.5, 'mean+0.5SD'), (1, 'mean+1SD')]:
    t = pc.mean() + k * sd; inc = sound & (pc > t)
    a = int(inc[m.site == 'BGS'].sum()); b = int((sound & (m.site == 'BGS') & ~inc).sum())
    c = int(inc[m.site == 'EMS'].sum()); d = int((sound & (m.site == 'EMS') & ~inc).sum())
    OR, p = stats.fisher_exact([[a, b], [c, d]])
    print(f"  {lab:18s}: BGS {a}/{a+b} vs EMS {c}/{c+d}   OR={OR:.2f}  p={p:.3f}")

inc0 = sound & (pc > pc.mean())
a = int(inc0[m.site == 'BGS'].sum()); n1 = int((m.site == 'BGS').sum())
c = int(inc0[m.site == 'EMS'].sum()); n2 = int((m.site == 'EMS').sum())
OR, p = stats.fisher_exact([[a, n1 - a], [c, n2 - c]])
print(f"  As % of ALL trees (31% vs 11%): {a}/{n1} vs {c}/{n2}  OR={OR:.2f}  p={p:.3f}")

g = m[m.sp.isin(['A.rubrum', 'T.canadensis'])]
incg = (~g['struct']) & (g['pc1'] > pc.mean())
a = int(incg[g.site == 'BGS'].sum()); n1 = int((g.site == 'BGS').sum())
c = int(incg[g.site == 'EMS'].sum()); n2 = int((g.site == 'EMS').sum())
OR, p = stats.fisher_exact([[a, n1 - a], [c, n2 - c]])
print(f"  GENERALISTS ONLY: BGS {a}/{n1} ({a/n1*100:.0f}%) vs EMS {c}/{n2} ({c/n2*100:.0f}%)  OR={OR:.2f}  p={p:.3f}")

print("\n2. Absolute anchor uncertainty (hemlock calibration, n=12)")
v = pd.read_csv('data/hemlock/validation_summary.csv'); v['res'] = 1000 / v['Conductance']
sl = stats.linregress(v['res'], v['moisture'])
n = len(v); xbar = v['res'].mean(); Sxx = ((v['res'] - xbar) ** 2).sum()
rstd = np.sqrt((((v['moisture'] - (sl.intercept + sl.slope * v['res']))) ** 2).sum() / (n - 2))
tcrit = stats.t.ppf(0.975, n - 2)
for M in [90, 100, 110]:
    x = (M - sl.intercept) / sl.slope
    ci = rstd * np.sqrt(1 / n + (x - xbar) ** 2 / Sxx) / abs(sl.slope) * tcrit
    pi = rstd * np.sqrt(1 + 1 / n + (x - xbar) ** 2 / Sxx) / abs(sl.slope) * tcrit
    print(f"  moisture>{M}% -> res<{x:.0f} Ohm-m; 95% CI (stand-level anchor): {x-ci:.0f}..{x+ci:.0f}; "
          f"95% PI (single tree): {x-pi:.0f}..{x+pi:.0f}")
