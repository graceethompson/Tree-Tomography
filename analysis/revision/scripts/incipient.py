import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
df=pd.read_csv('./data/Tree_ID_info.csv'); df['site']=np.where(df['plot']=='BGS','BGS','EMS')
sm={'rm':'A.rubrum','bg':'N.sylvatica','ro':'Q.rubra','hem':'T.canadensis'}; df['sp']=df['species'].map(sm)
ert=pd.read_csv('./data/ERT_application_results.csv'); m=df.merge(ert,on='tree').reset_index(drop=True)
# species-standardized resistivity deviation; ORIENT so higher = WETTER (right)
m['xdev']=-m.groupby('sp')['mean'].transform(lambda x:(x-x.median())/x.std(ddof=1))
m['struct']=m['percent_damaged']>1
sound=~m['struct']; dmgd=m['struct']; mb=(m.site=='BGS').values; me=(m.site=='EMS').values

# ---- WHERE DOES INCIPIENT START? sensitivity of the sound->incipient cut ----
print("Incipient onset = sound tree confidently WETTER than its species median.")
print("(xdev = species-standardized; wetter = higher; +0.5 SD = confidently wet)")
print(f"{'cut (SD wetter than median)':32s}{'nII':>5}{'BGS%':>7}{'EMS%':>7}")
for cut,lab in [(0,'median (any wetter)'),(0.25,'+0.25 SD'),(0.5,'+0.5 SD'),(1.0,'+1.0 SD')]:
    inc=sound&(m['xdev']>cut)
    b=inc[mb].sum()/sound[mb].sum()*100; e=inc[me].sum()/sound[me].sum()*100
    print(f"  {lab:30s}{inc.sum():>5}{b:>7.0f}{e:>7.0f}")

# ---- internal 'positive control': moisture of structurally-confirmed decay (active) trees ----
# active = damaged & wet (xdev>0). Their xdev shows 'decay-associated moisture' range.
active=m[dmgd & (m['xdev']>0)]
print(f"\nInternal reference — structurally-damaged WET (active) trees, n={len(active)}:")
print(f"  their xdev (wetness): median={active['xdev'].median():.2f}, range {active['xdev'].min():.2f}..{active['xdev'].max():.2f}")
print(f"  -> a sound tree wetter than ~+0.5 SD overlaps the active-decay moisture range (supports +0.5 SD onset)")

# ---- 6-cell / dead-band classification ----
def classify(cut):
    wet=m['xdev']>cut; dry=m['xdev']<-cut
    cat=np.where(sound & wet,'II',
        np.where(sound & ~wet,'I',
        np.where(dmgd & dry,'IV','III')))
    return cat
for cut in [0.5,1.0]:
    from collections import Counter
    print(f"\n6-cell scheme, dead band ±{cut} SD:  {dict(Counter(classify(cut)))}")

# ---- FIGURE: 6-cell dead-band phase diagram (two dead-band widths) ----
catcol={'I':'#3b6fb0','II':'#4e9a2c','III':'#d98a1f','IV':'#b83232'}
np.random.seed(2)
fig,axes=plt.subplots(1,2,figsize=(15,6),sharey=True)
for ax,cut in zip(axes,[0.5,1.0]):
    cat=classify(cut); y=np.sqrt(np.clip(m['percent_damaged'],0,None))
    ax.axvspan(-cut,cut,color='0.85',alpha=0.6,zorder=0)
    ax.axvline(0,color='#666',ls='-',lw=1); ax.axhline(np.sqrt(1),color='#666',ls='-',lw=1)
    for c,col in catcol.items():
        idx=cat==c
        ax.scatter(m['xdev'][idx], y[idx]+np.random.normal(0,0.03,idx.sum()), c=col, s=48, edgecolor='white',lw=.5,zorder=3,label=c)
    yt=[0,1,5,10,20,30]; ax.set_yticks([np.sqrt(t) for t in yt]); ax.set_yticklabels(yt)
    ax.set_xlabel('moisture anomaly  (species-standardized;  drier ←   0 = species median   → wetter)',fontsize=9)
    # cell labels
    xl=ax.get_xlim()
    ax.text(xl[0]+0.2,np.sqrt(31),'IV Cavity\n(damaged, confidently dry)',fontsize=8,color=catcol['IV'],weight='bold',va='top')
    ax.text(xl[1]-0.2,np.sqrt(31),'III Active\n(damaged, wet)',fontsize=8,color=catcol['III'],weight='bold',ha='right',va='top')
    ax.text(xl[0]+0.2,0.05,'I No decay\n(sound, not-wet)',fontsize=8,color=catcol['I'],weight='bold')
    ax.text(xl[1]-0.2,0.05,'II Incipient\n(sound, confidently wet)',fontsize=8,color=catcol['II'],weight='bold',ha='right')
    ax.text(0,np.sqrt(34),f'dead band ±{cut} SD\n(transition → defaults to row baseline)',fontsize=7.5,ha='center',color='#555')
    ax.set_title(f'Dead band ±{cut} SD',fontsize=10)
axes[0].set_ylabel('SoT structural loss %  (√ scale)')
axes[0].legend(title='Category',fontsize=8.5,frameon=False,loc='center left')
fig.suptitle('6-cell scheme: 4 confident corners + a "normal" transition column that defaults to the row baseline\n'
             '(sound→No-decay, damaged→Active). Only "confidently wet sound" (II) and "confidently dry damaged" (IV) are off-baseline.',
             fontsize=11,y=1.02)
plt.tight_layout(); plt.savefig('analysis/revision/output/CJFR-6cell.png',dpi=140,bbox_inches='tight')
print("\nsaved 6-cell figure")
