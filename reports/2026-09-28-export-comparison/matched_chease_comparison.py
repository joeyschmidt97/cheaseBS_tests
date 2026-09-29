import sys,pathlib,json,csv,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,'C:/Users/joesc/git')
from TPED.projects.discharge_tools.src.cheasebs_runner import _load_run,_shear_profile,_at_radius
from TPED.projects.discharge_tools.src.filetypes.gfile_data import GFileData
root=pathlib.Path('C:/Users/joesc/git/cheaseBS_tests/data');out=pathlib.Path(__file__).parent/'chease_matched';out.mkdir(exist_ok=True)
rows=[];fig,axes=plt.subplots(2,2,figsize=(11,7))
for row,shot in enumerate(['129015','129038']):
 for col,rho in enumerate([.9,.95]):
  ax=axes[row,col]
  for era in ['new','old']:
   vals=[]
   for alpha in [.8,.9,1.,1.1]:
    path=(root/'lleppin_NSTX_chease_tests'/('NSTX'+shot)/'variations'/('omt'+str(alpha).replace('.','p')+'_omne1p0') if era=='new' else root/'129015_129038'/'pm0.2_ncscal4_warmup0'/shot/f'omt_{alpha:.3f}')
    d=_load_run(str(path));g=GFileData(d['files']['gfile_after']).gfile_to_xarray();val=float(_at_radius(_shear_profile(g),rho));vals.append(val)
    rows.append(dict(shot=shot,era=era,rho_tor=rho,omt=alpha,omne=1.,shear=val,path=str(path),**d['scalars']))
   span=100*np.ptp(vals)/abs(vals[2]);label=('new: jparallel / rhot' if era=='new' else 'old: I-star / rhop')+f' (span {span:.2f}%)'
   ax.plot([.8,.9,1.,1.1],vals,'o-' if era=='new' else 's--',label=label)
  ax.set_title(f'NSTX {shot}: rho_tor={rho}');ax.set_xlabel('omt factor (omne=1)');ax.set_ylabel('magnetic shear');ax.grid(alpha=.25);ax.legend(fontsize=8)
fig.suptitle('NSTX matched scan factors: new exports versus September 25 runs\nBoth NCSCAL=4; other solver settings and provenance differ',fontsize=12)
fig.tight_layout(rect=(0,.03,1,.92));fig.text(.015,.012,'Span = (maximum - minimum) / |unity shear|. Same TPED parser and magnetic-shear calculation. Not an isolated replay-method experiment.',fontsize=8)
fig.savefig(out/'nstx_matched_old_new.png',dpi=160)
(out/'nstx_matched_old_new.json').write_text(json.dumps(rows,indent=2))
print(out)
