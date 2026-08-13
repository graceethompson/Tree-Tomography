import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

fig,ax=plt.subplots(figsize=(13.5,9)); ax.set_xlim(0,100); ax.set_ylim(0,100); ax.axis('off')
def box(x,y,w,h,text,fc,ec='#333',fs=9.5,tc='#111',bold=False):
    ax.add_patch(FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle='round,pad=0.6,rounding_size=1.5',
                 fc=fc,ec=ec,lw=1.4))
    ax.text(x,y,text,ha='center',va='center',fontsize=fs,color=tc,fontweight='bold' if bold else 'normal',linespacing=1.35)
def arrow(x1,y1,x2,y2,label=None,lx=0,ly=0,col='#444'):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',mutation_scale=15,lw=1.5,color=col))
    if label: ax.text((x1+x2)/2+lx,(y1+y2)/2+ly,label,fontsize=8.5,color=col,ha='center',fontweight='bold')

# input
box(50,95,54,7,'Each tree:  SoT % structural damage  +  ERT resistivity metrics (8, incl. mean ρ and CMA)','#eef2f7',fs=9)
arrow(50,91.5,50,88)
# split 1 SoT
box(50,84,40,8,'AXIS 1 — STRUCTURAL LOSS  (SoT)\nIs SoT % damage  >  1% ?','#d6e4f0',fs=9.5,bold=True)
ax.text(50,77.5,'crisp / bimodal: 0% or ≥6%, empty gap between  ·  stable for any cut in 1–5% (2 trees change)',ha='center',fontsize=7.3,color='#555',style='italic')

# NO branch -> moisture anomaly (I/II)
arrow(31,84,16,70,'NO  (no structural loss)',lx=-3,ly=6,col='#26456e')
box(16,66,30,9,'AXIS 2a — MOISTURE ANOMALY\nspecies-normalized ERT PC1\n(continuum — graded)','#e8f0e2',fs=8.8,bold=True)
# graded outcomes
arrow(16,61.5,10,50)
box(10,45,17,8,'I — No decay\nPC1 ≤ mean (0)','#dbe4ef',fs=8.5,tc='#26456e',bold=True)
arrow(16,61.5,27,50)
box(30,45,20,10,'II — Incipient\n"possible": 0 < PC1 ≤ 1.0\n"confident": PC1 > 1.0','#e5efd6',fs=8.2,tc='#4d7f00',bold=True)
ax.text(20,37,'≈ 30 normal / 7 possible / 5 confident (of 39 sound trees)\nincipient ≠ decay: could be early decay OR bacterial wetwood',
        ha='center',fontsize=7,color='#555',style='italic')

# YES branch -> core wetness (III/IV)
arrow(69,84,84,70,'YES  (structural loss present)',lx=6,ly=6,col='#7f4f00')
box(84,66,30,9,'AXIS 2b — CORE WETNESS\nis the damaged core WET or DRY ?','#f3e6d6',fs=8.8,bold=True)
arrow(84,61.5,72,50)
box(69,44,26,11,'III — Active decay\nWET core:\nmean ρ < 300 Ω·m  OR  CMA ≥ 0.15','#f0e2cf',fs=8,tc='#7f4f00',bold=True)
arrow(84,61.5,95,50)
box(92,44,15,11,'IV — Cavity\nDRY core:\nmean ρ ≥ 300\nAND CMA < 0.15','#f0d6d6',fs=7.8,tc='#7f0000',bold=True)
ax.text(84,35.5,'physically-grounded refinement of the composite cut\n(reclassifies 893: wet core → Active, not Cavity)',
        ha='center',fontsize=7,color='#555',style='italic')

# sensitivity footer band
box(50,22,92,12,
 'SENSITIVITY / STATUS\n'
 '• SoT 1% threshold: ROBUST — reclassifies only 2/57 trees across 1–5%; report equipment detection floor.\n'
 '• PC1 anomaly cut (I↔II): OPERATIONAL CONVENTION on a continuum — ±0.5 SD (±0.95) reclassifies 11–14 trees; report graded bands + full sweep. Direction (BGS>EMS) invariant to threshold & normalization.\n'
 '• Core-wetness cut (III↔IV): PHYSICALLY grounded but thresholds (300 Ω·m, CMA 0.15) are data-suggested, not validated.',
 '#fbfbf6',ec='#aaa',fs=8,tc='#222')
box(50,10,92,7,
 'ABSOLUTE / FUTURE:  PC1 mean ≈ 369 Ω·m ≈ 103% predicted moisture; Otsu break ≈ 387 Ω·m ≈ 101%.  Absolute cut premature — n=12 calibration 95% band ≈ 5–800 Ω·m.\n'
 'True incipient-decay vs wetwood, and per-species absolute thresholds, require destructive discs/cross-sections.',
 '#f0f0ea',ec='#aaa',fs=7.6,tc='#333')
ax.set_title('Decay classification — decision rules & thresholds',fontsize=13,fontweight='bold',y=1.0)
plt.tight_layout(); plt.savefig('analysis/revision/output/CJFR-decision-tree.png',dpi=150,bbox_inches='tight')
print('saved flowchart')
