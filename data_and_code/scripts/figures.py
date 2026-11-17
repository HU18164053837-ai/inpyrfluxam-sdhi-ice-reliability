"""Journal layouts from unchanged saved results; never call prediction models."""
from pathlib import Path
import sys
BASE=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(BASE/'tools/pylibs'),str(BASE/'phase1c/tools/pylibs')]
import pandas as pd,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
R=Path(__file__).resolve().parents[2];I=R/'data_and_code/inputs'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
C={'inpyrfluxam':'#0072B2','isopyrazam':'#D55E00'}
sp=['Oncorhynchus mykiss','Pimephales promelas','Lepomis macrochirus','Cyprinus carpio','Danio rerio','Poecilia reticulata','Oryzias latipes']
short={s:s.split()[0][0]+'. '+s.split()[1] for s in sp}
p=pd.read_csv(I/'phase1e/all90_predictions_recomputed.csv');u=pd.read_csv(I/'phase1g_b/endpoint_uncertainty_assessment.csv')
t=pd.read_csv(I/'phase1g_b/target_species_summary.csv');t=t[(t.scenario=='ALL_AVAILABLE') & (t.validation_tier=='EXTERNAL_CANDIDATE')].copy()
d=pd.read_csv(I/'phase1g_b/compound_direction_summary.csv');d=d[(d.scenario=='ALL_AVAILABLE') & (d.validation_tier=='EXTERNAL_CANDIDATE')]
m=pd.read_csv(I/'phase1g_b/compound_equal_target_summary.csv');m=m[(m.scenario=='ALL_AVAILABLE') & (m.validation_tier=='EXTERNAL_CANDIDATE')]
def label(ax,s):ax.text(-.12,1.045,s,transform=ax.transAxes,fontweight='bold',fontsize=12,va='bottom')
def save(fig,n):
 for ext in ['pdf','svg','png']:fig.savefig(R/'figures'/f'{n}.{ext}',dpi=360,bbox_inches='tight',metadata={'Creator':'Research figure'} if ext=='pdf' else None)
 plt.close(fig)
fig,(a,b)=plt.subplots(1,2,figsize=(8.3,4.0),gridspec_kw={'width_ratios':[1.05,1.4]})
for ci,c in enumerate(C):
 z=u[u.compound==c].set_index('species')
 for j,s in enumerate(sp):
  if s not in z.index:continue
  r=z.loc[s];y=j+(ci-.5)*.18;x=r.observed_LC50_ug_L
  if pd.notna(r.CI95_lower_ug_L):a.plot([r.CI95_lower_ug_L,r.CI95_upper_ug_L],[y,y],color=C[c],lw=1)
  a.scatter(x,y,color=C[c],s=27,marker='o' if ci==0 else 's',zorder=3)
a.set_xscale('log');a.set_yticks(range(7),[short[s] for s in sp],fontstyle='italic');a.invert_yaxis();a.set_xlabel('Measured 96-h LC₅₀ (μg/L)');a.axvline(100,c='.55',ls=':',lw=.8)
a.legend(handles=[Line2D([],[],color=C[c],marker='o' if i==0 else 's',ls='',label=c.capitalize()) for i,c in enumerate(C)],loc='lower left',frameon=False,fontsize=7)
comp=list(C)+['fluxapyroxad','benzovindiflupyr','bixafen','sedaxane','penthiopyrad','pydiflumetofen']
end=p[['compound','target','observed_ug_L']].drop_duplicates();matrix=np.full((8,7),np.nan)
for j,c in enumerate(comp):
 for k,s in enumerate(sp):
  z=end[(end.compound==c)&(end.target==s)]
  if len(z):matrix[j,k]=z.iloc[0].observed_ug_L
b.imshow(np.where(np.isnan(matrix),np.nan,1),cmap=matplotlib.colors.ListedColormap(['#edf0f2']),aspect='auto')
for j in range(8):
 for k in range(7):
  if np.isfinite(matrix[j,k]):b.text(k,j,f'{matrix[j,k]:g}',ha='center',va='center',fontsize=7,color=C.get(comp[j],'#333333'))
b.set_xticks(range(7),[short[s] for s in sp],rotation=65,ha='right',fontstyle='italic',fontsize=7)
counts=p.groupby('compound').size();b.set_yticks(range(8),[f'{c.capitalize()} ({counts[c]})' for c in comp]);b.axhline(1.5,c='.5',ls='--',lw=.8);b.set_xlabel('Measured LC₅₀ (μg/L); directions in parentheses')
label(a,'A');label(b,'B');fig.tight_layout(w_pad=2.7);save(fig,'Fig1')
fig,ax=plt.subplots(2,2,figsize=(7.8,7.0))
for j,c in enumerate(C):
 z=p[p.compound==c].copy();a=ax[0,j];x=np.geomspace(10,2000,150);a.fill_between(x,x/5,x*5,color='#eeeeee');a.plot(x,x,c='.4',lw=1)
 for good,marker in [(True,'o'),(False,'x')]:
  q=z[z.guidance_pass==good];a.scatter(q.observed_ug_L,q.predicted_ug_L,c=C[c],s=24,marker=marker,alpha=.8)
 a.set(xscale='log',yscale='log',xlim=(10,2000),ylim=(5,5000),xlabel='Observed LC₅₀ (μg/L)',ylabel='Predicted LC₅₀ (μg/L)');a.text(.04,.96,c.capitalize(),transform=a.transAxes,va='top',fontweight='bold')
 a=ax[1,j];z=z.sort_values(['target','surrogate']).reset_index(drop=True)
 for k,r in z.iterrows():
  obs=r.observed_ug_L;a.plot(np.log10([r.PI95_lower/obs,r.PI95_upper/obs]),[k,k],c='.78',lw=.6)
  a.plot(np.log10([r.CI95_lower/obs,r.CI95_upper/obs]),[k,k],c=C[c],lw=1.3)
  a.scatter(r.signed_log10_error,k,c=C[c],s=12,marker='o' if r.guidance_pass else 'x',zorder=3)
 a.axvline(0,c='.35',lw=.8);a.set(xlabel='Signed log₁₀ error',ylabel='Saved direction index',xlim=(-4,4),ylim=(len(z),-1))
 a.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
for a,l in zip(ax.flat,'ABCD'):label(a,l)
fig.tight_layout(h_pad=2,w_pad=2);save(fig,'Fig2')
fig,ax=plt.subplots(2,2,figsize=(8.0,7.4),gridspec_kw={'height_ratios':[1.4,1]})
t=t.sort_values(['compound','target']).reset_index(drop=True)
for k,r in t.iterrows():
 ax[0,0].plot([r.identity_MAE_log10,r.webice_MAE_log10],[k,k],c='.65',lw=1)
 ax[0,0].scatter(r.identity_MAE_log10,k,marker='s',s=23,facecolor='white',edgecolor=C[r.compound],zorder=3)
 ax[0,0].scatter(r.webice_MAE_log10,k,s=23,c=C[r.compound],zorder=3)
 delta=100*(r.webice_class_agreement_rate-r.identity_class_agreement_rate)
 ax[0,1].barh(k,delta,color=C[r.compound],height=.62)
for a in ax[0]:a.set_yticks(range(len(t)),[short[s] for s in t.target],fontstyle='italic',fontsize=7);a.invert_yaxis();a.axhline(6.5,c='.6',ls=':',lw=.8)
ax[0,0].set_xlabel('Target mean absolute log₁₀ error');fig.legend(handles=[Line2D([],[],c='k',marker='o',ls='',label='Web-ICE'),Line2D([],[],c='k',marker='s',mfc='white',ls='',label='Identity transfer')],frameon=False,fontsize=7,loc='upper center',bbox_to_anchor=(.5,1.025),ncol=2)
ax[0,1].set_xlabel('Agreement change (percentage points)');ax[0,1].axvline(0,c='.4',lw=.8);ax[0,1].set_xlim(-70,70)
xs=[];labs=[]
for j,c in enumerate(C):
 r=d[d.compound==c].iloc[0];n=int(r.n_directions)
 for k,method in enumerate(['webice','identity']):
  x=j*3+k;xs.append(x);labs.append(c.capitalize()+'\n'+('Web-ICE' if k==0 else 'Identity'))
  mismatch=int(r[method+'_class_mismatch_n']);under=int(r[method+'_hazard_category_underestimated_n']);match=n-mismatch
  bottom=0
  for v,color in [(match,'#c7d9e8'),(under,'#d89065'),(mismatch-under,'#dddddd')]:
   ax[1,0].bar(x,v,bottom=bottom,color=color,width=.72)
   if v:ax[1,0].text(x,bottom+v/2,str(v),ha='center',va='center',fontsize=8)
   bottom+=v
 r=m[m.compound==c].iloc[0]
 for k,method in enumerate(['webice','identity']):
  ax[1,1].scatter(j*3+k,100*r['equal_target_'+method+'_hazard_category_underestimated_rate'],marker='o' if k==0 else 's',c=C[c],s=35)
ax[1,0].set_xticks(xs,['Web-ICE','Identity']*2,fontsize=7);ax[1,0].set_ylabel('Dependent prediction directions');ax[1,0].legend(handles=[matplotlib.patches.Patch(color=co,label=la) for co,la in [('#c7d9e8','Class match'),('#d89065','Hazard underestimated'),('#dddddd','Hazard overestimated')]],frameon=False,fontsize=7,loc='upper right')
ax[1,0].set_ylim(0,46)
ax[1,1].set_xticks(xs,['Web-ICE','Identity']*2,fontsize=7);ax[1,1].set_ylabel('Target-equal hazard underestimation (%)');ax[1,1].set_ylim(0,60)
for a in ax[1]:
 for j,c in enumerate(C):a.text(j*3+.5,-.16,c.capitalize(),transform=a.get_xaxis_transform(),ha='center',va='top',fontsize=8)
for a,l in zip(ax.flat,'ABCD'):label(a,l)
fig.tight_layout(h_pad=2.2,w_pad=2.1);save(fig,'Fig3')
print('Created 3 figures in PDF/SVG/PNG from frozen tables.')
