"""Regenerate exported CHEASE-BS comparisons using the TPED plotting CLI."""
import os, sys, json, pathlib, subprocess, csv, re, argparse, importlib.util
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
parser=argparse.ArgumentParser()
parser.add_argument('--data-root',type=pathlib.Path,required=True)
parser.add_argument('--tped',type=pathlib.Path,required=True)
parser.add_argument('--outdir',type=pathlib.Path,required=True)
parser.add_argument('--resume',action='store_true')
args=parser.parse_args()
sys.path.insert(0,str(args.tped.parent))
from TPED.projects.discharge_tools.src.cheasebs_runner import _load_run,_shear_profile,_at_radius
from TPED.projects.discharge_tools.src.filetypes.gfile_data import GFileData
cli=args.tped/'projects/discharge_tools/scripts/compare_cheasebs_runs.py'
out=args.outdir;out.mkdir(parents=True,exist_ok=True)
configs=sorted(args.data_root.glob('lleppin_*/*/*/scan_config.json'))
allrows=[];campaigns=[]; commands=[]; failures=[]
if args.resume and (out/'case_values.json').exists():
 allrows=json.loads((out/'case_values.json').read_text())
 campaigns=json.loads((out/'campaign_summary.json').read_text())
finished={(c['shot'],c['campaign']) for c in campaigns}
for k,config in enumerate(configs,1):
 campaign=config.parent; baseline=campaign.parent; shot=baseline.name
 paths=sorted(campaign.glob('omt*/metadata.json'))
 if not paths: continue
 dest=out/shot/campaign.name;dest.mkdir(parents=True,exist_ok=True)
 print(f'[{k}/{len(configs)}] {shot}/{campaign.name}: {len(paths)} cases',flush=True)
 command=[str(cli),'--baseline',str(baseline),'--variations',str(campaign),'--outdir',str(dest),
          '--pedestal','0.8,0.98','--slices','0.85','0.9','0.95','0.975','--quiet',
          '--title',f'{shot} / {campaign.name}']
 commands.append([sys.executable]+command)
 if (shot,campaign.name) in finished and (dest/"q_shear.png").exists():
  print("  already complete",flush=True);continue
 with (dest/'plot.log').open('w',encoding='utf-8') as log:
  completed=subprocess.run([sys.executable]+command,stdout=log,stderr=subprocess.STDOUT,
                            env=dict(os.environ,PYTHONUTF8='1'),text=True)
 if completed.returncode:
  failures.append(str(campaign));print('PLOT FAILED',campaign,flush=True)
 rows=[];curves=[]; source=None
 for meta_path in paths:
  meta=json.loads(meta_path.read_text());run=meta_path.parent
  d=_load_run(str(run),str(baseline));files=d['files'];sc=d['scalars']
  if source is None: source=GFileData(files['gfile_before']).gfile_to_xarray()
  if source is None: raise ValueError('Source EQDSK could not be parsed: '+str(files['gfile_before']))
  g=GFileData(files['gfile_after']).gfile_to_xarray();shear=_shear_profile(g)
  mods=meta.get('modification',{})
  row=dict(shot=shot,campaign=campaign.name,case=run.name,path=str(run),
           ncscal=d['cfg'].get('ncscal'),coordinate=d['cfg'].get('coordinate'),
           replay=d['cfg'].get('replay_representation'),
           rhotTopPed=mods.get('omt',{}).get('rhotTopPed'),
           omt=mods.get('omt',{}).get('alpha'),omne=mods.get('omne',{}).get('alpha'),
           **sc)
  for rho in [.85,.9,.95,.975]:
   row[f'q_{rho}']=float(np.interp(rho,g.rho_tor,np.abs(g.q.values)))
   row[f's_{rho}']=float(_at_radius(shear,rho))
  row['q0']=float(abs(g.q.values[0]))
  rows.append(row);curves.append((run.name,g,shear,sc))
 allrows+=rows
 summary=dict(shot=shot,campaign=campaign.name,cases=len(rows),converged=sum(r['converged'] is True for r in rows),
              iterations=[min(r['iterations'] for r in rows if r['iterations'] is not None),max(r['iterations'] for r in rows if r['iterations'] is not None)],
              ncscal=rows[0]['ncscal'],rhotTopPed=rows[0]['rhotTopPed'])
 for rho in [.85,.9,.95,.975]:
  vals=[r[f's_{rho}'] for r in rows];unity=next((r[f's_{rho}'] for r in rows if r['omt']==1 and r['omne']==1),None)
  summary[f's_range_{rho}']=[min(vals),max(vals)]
  summary[f's_span_unity_pct_{rho}']=100*(max(vals)-min(vals))/abs(unity) if unity else None
 campaigns.append(summary)
 # Compact view reuses exactly the q parser and shear helper used by the CLI.
 fig,axes=plt.subplots(2,2,figsize=(13,9))
 sx=source.rho_tor.values;sy=np.abs(source.q.values);ss=_shear_profile(source)
 for i,(name,g,sh,sc) in enumerate(curves):
  color=plt.get_cmap('turbo')(i/max(len(curves)-1,1));style='-' if sc['converged'] else '--'
  label=name+(' [not converged]' if sc['converged'] is False else '')
  for a in axes[:,0]: a.plot(g.rho_tor,np.abs(g.q.values),color=color,ls=style,lw=1.2,label=label)
  for a in axes[:,1]: a.plot(*sh,color=color,ls=style,lw=1.2,label=label)
 for a in axes[:,0]:a.plot(sx,sy,'k:',lw=2,label='source EFIT');a.set_ylabel('q')
 for a in axes[:,1]:a.plot(*ss,'k:',lw=2,label='source EFIT');a.set_ylabel('magnetic shear')
 for a in axes[0]:a.set_xlim(0,1)
 for a in axes[1]:
  a.set_xlim(.8,.98)
  # Match the CLI zoom autoscaling to the displayed radius window.
  y=[]
  for line in a.lines:
   x,v=line.get_data();x=np.asarray(x);v=np.asarray(v);m=(x>=.8)&(x<=.98)&np.isfinite(v);y.extend(v[m])
  if y:
   span=max(y)-min(y);a.set_ylim(min(y)-.06*span,max(y)+.06*span)
 for a in axes.flat:a.grid(alpha=.2);a.set_xlabel('rho_tor')
 axes[0,0].set_title('Full q profile');axes[0,1].set_title('Full shear profile')
 axes[1,0].set_title('Pedestal q');axes[1,1].set_title('Pedestal shear')
 handles,labels=axes[0,0].get_legend_handles_labels()
 fig.legend(handles,labels,loc='center left',bbox_to_anchor=(.77,.5),fontsize=7)
 fig.suptitle(f'{shot} / {campaign.name}\nNCSCAL={summary["ncscal"]}; rhotTopPed={summary["rhotTopPed"]}; '
              f'{summary["converged"]}/{len(rows)} solver-converged; iterations {summary["iterations"]}',fontsize=12)
 fig.tight_layout(rect=(0,.025,.78,.92));fig.text(.02,.012,'Same TPED q/shear calculation. Dashed colored curves: saved solver flag false. Separatrix-adjacent derivatives need resolution checks.',fontsize=8)
 fig.savefig(dest/'q_shear.png',dpi=150);plt.close(fig)
 (out/'campaign_summary.json').write_text(json.dumps(campaigns,indent=2),encoding='utf-8')
 (out/'case_values.json').write_text(json.dumps(allrows,indent=2),encoding='utf-8')
(out/'commands.json').write_text(json.dumps(commands,indent=2),encoding='utf-8')
with (out/'case_values.csv').open('w',newline='',encoding='utf-8') as f:
 writer=csv.DictWriter(f,fieldnames=list(allrows[0]));writer.writeheader();writer.writerows(allrows)
(out/'failures.json').write_text(json.dumps(failures,indent=2),encoding='utf-8')
print('DONE',len(campaigns),'campaigns',len(allrows),'cases;',len(failures),'plot failures',flush=True)
