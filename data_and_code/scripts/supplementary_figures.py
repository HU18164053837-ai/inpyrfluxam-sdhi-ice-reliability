from pathlib import Path
import csv,sys,math
ROOT=Path(__file__).resolve().parents[2];BASE=ROOT.parent
sys.path[:0]=[str(BASE/'tools/pylibs'),str(BASE/'phase1c/tools/pylibs')]
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt,numpy as np
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42})
blue='#0072B2';orange='#D55E00';gray='#777777'
def read(p):
 with Path(p).open(encoding='utf-8-sig') as f:return list(csv.DictReader(f))
def f(r,k):return float(r[k])
def out(fig,n):
 name='FigS2' if n==3 else 'FigS1'
 for ext in ['pdf','svg','png']:fig.savefig(ROOT/'figures'/f'{name}.{ext}',dpi=360,bbox_inches='tight')
 plt.close(fig)
comp=['inpyrfluxam','isopyrazam','fluxapyroxad','benzovindiflupyr','bixafen','sedaxane','penthiopyrad','pydiflumetofen']
strategy=read(ROOT/'data_and_code/inputs/phase1e/strategy_by_compound.csv');names=['A_quality','B_domain','C_geomean','C_median','D_rainbow','D_zebrafish'];labels=['A Minimum MSE','B Minimum leverage','C Log-mean','C Log-median','D Trout','D Zebrafish']
fig=plt.figure(figsize=(10,8.4),layout='constrained');gs=fig.add_gridspec(2,2,height_ratios=[1,1.2]);ax=[fig.add_subplot(gs[0,0]),fig.add_subplot(gs[0,1]),fig.add_subplot(gs[1,:])]
for j,c in enumerate(['inpyrfluxam','isopyrazam']):
 for i,s in enumerate(names):
  z=next(x for x in strategy if x['scenario']=='PRIMARY' and x['compound']==c and x['strategy']==s);ax[j].scatter(f(z,'MAE_log10'),i,c=blue);ax[j].text(f(z,'MAE_log10')+.015,i,f"{z['n_directions']}/{z['n_possible_targets']} targets",va='center',fontsize=10)
 ax[j].set(yticks=range(6),yticklabels=labels,xlabel='Mean absolute log$_{10}$ error',xlim=(0,1.1));ax[j].invert_yaxis();ax[j].text(.04,1.03,c,transform=ax[j].transAxes)
di=read(ROOT/'data_and_code/inputs/phase1e/trout_fathead_direction_and_reverse.csv')
for i,c in enumerate(comp):
 z=next(x for x in di if x['compound']==c and x['surrogate']=='Oncorhynchus mykiss');rev=next(x for x in di if x['compound']==c and x['surrogate']=='Pimephales promelas');ax[2].scatter(f(z,'signed_log10_error'),i,c=blue,marker='o');ax[2].scatter(f(rev,'signed_log10_error'),i,c=orange,marker='s')
ax[2].axhline(1.5,c='#bbbbbb',lw=.7,ls=':');ax[2].axvline(0,c=gray,lw=.8);ax[2].set(yticks=range(8),yticklabels=comp,xlabel='Signed log$_{10}$ error',xlim=(-1,1));ax[2].invert_yaxis();ax[2].scatter([],[],c=blue,label='Trout to fathead');ax[2].scatter([],[],c=orange,marker='s',label='Fathead to trout');ax[2].legend(fontsize=10,loc='upper center',bbox_to_anchor=(.5,-.2))
for i,a in enumerate(ax):a.text(-.15,1.12,'ABC'[i],transform=a.transAxes,fontweight='bold')
out(fig,3)
fig=plt.figure(figsize=(10,7.6),layout='constrained');gs=fig.add_gridspec(2,2,height_ratios=[1,1.2]);ax=[fig.add_subplot(gs[0,:]),fig.add_subplot(gs[1,0]),fig.add_subplot(gs[1,1])];rr=read(ROOT/'supplementary/tables/core_pocket_alignment_phase1c.csv');order=['Oncorhynchus_mykiss','Pimephales_promelas','Danio_rerio','Oryzias_latipes'];sites=list(dict.fromkeys((x['subunit'],x['human_position']) for x in rr));arr=np.ones((4,len(sites)))
for i,s in enumerate(order):
 for j,site in enumerate(sites):
  z=next(x for x in rr if x['species']==s and (x['subunit'],x['human_position'])==site);arr[i,j]=z['fish_aa']==z['human_aa'];ax[0].text(j,i,z['fish_aa'],ha='center',va='center',fontsize=11)
ax[0].imshow(arr,cmap=matplotlib.colors.ListedColormap(['#e6ab5b','#d8e7f3']),vmin=0,vmax=1,aspect='auto',zorder=-1);ax[0].set(yticks=range(4),yticklabels=['Trout*','Fathead','Zebrafish','Medaka'],xticks=range(len(sites)),xticklabels=[s+p for s,p in sites]);plt.setp(ax[0].get_xticklabels(),rotation=60,ha='right');ax[0].set_xlabel('Human template chain and position')
red=read(ROOT/'supplementary/tables/redocking_summary.csv')
for i,(p,score) in enumerate([('9KC4','vina'),('9KC4','vinardo'),('9QDM','vina'),('9QDM','vinardo')]):
 z=[x for x in red if x['pdb']==p and x['scoring']==score];ax[1].scatter([i-.1,i,i+.1],[f(x,'top_RMSD_A') for x in z],c=blue if score=='vina' else orange)
ax[1].axhline(2,c=gray,ls='--',lw=.8);ax[1].set(xticks=range(4),xticklabels=['Human\nVina','Human\nVinardo','Yeast\nVina','Yeast\nVinardo'],ylabel='Native-ligand RMSD (angstrom)')
pr=read(ROOT/'supplementary/tables/local_proxy_summary.csv')
for score,col,delta in [('vina',blue,-.08),('vinardo',orange,.08)]:
 z=[x for x in pr if x['scoring']==score];ax[2].errorbar([['I','L','M'].index(x['proxy'])+delta for x in z],[f(x,'mean') for x in z],yerr=[[f(x,'mean')-f(x,'min') for x in z],[f(x,'max')-f(x,'mean') for x in z]],fmt='o',c=col,label=score,capsize=3)
ax[2].set(xticks=[0,1,2],xticklabels=['I','L','M'],xlabel='Local human-scaffold substitution',ylabel='Docking score (not free energy)');ax[2].legend(fontsize=10)
for i,a in enumerate(ax):a.text(-.12,1.05,'ABC'[i],transform=a.transAxes,fontweight='bold')
out(fig,4)
