"""Finish the report index and paired checks after generate_report.py."""
import sys,pathlib,json,shutil,hashlib,os
import numpy as np
sys.path.insert(0,'C:/Users/joesc/git')
from TPED.projects.discharge_tools.src.cheasebs_runner import _load_run,_shear_profile
from TPED.projects.discharge_tools.src.filetypes.gfile_data import GFileData
out=pathlib.Path('C:/Users/joesc/git/cheaseBS_tests/reports/2026-09-28-export-comparison')
root=out.parents[1]/'data'
rows=json.loads((out/'case_values.json').read_text());campaigns=json.loads((out/'campaign_summary.json').read_text())
assert len(campaigns)==27 and len(rows)==323,(len(campaigns),len(rows))
assert not json.loads((out/'failures.json').read_text())
# Save the generators and matched NSTX analysis as reviewable, reproducible artifacts.
tmp=pathlib.Path(os.environ['TEMP'])
if (tmp/'build_chease_export_report.py').exists(): shutil.copyfile(tmp/'build_chease_export_report.py',out/'generate_report.py')
if (tmp/'matched_chease_comparison.py').exists(): shutil.copyfile(tmp/'matched_chease_comparison.py',out/'matched_chease_comparison.py')
if pathlib.Path(__file__).resolve() != (out/'finish_report.py').resolve(): shutil.copyfile(__file__,out/'finish_report.py')
(out/'chease_matched').mkdir(exist_ok=True)
for p in ((tmp/'chease_matched').iterdir() if (tmp/'chease_matched').exists() else []):
 if p.is_file():shutil.copyfile(p,out/'chease_matched'/p.name)
# Compare only cases whose saved solver flags both pass, on a shared rho_tor grid.
x=np.linspace(.8,.98,181);paired=[]
for shot in ['DIIID153764','DIIID162940','DIIID174082']:
 for a,b in [('chease','chease_ncscal4'),('chease_ncscal4','chease_rhotTopPed0p3_ncscal4')]:
  ar={r['case']:r for r in rows if r['shot']==shot and r['campaign']==a and r['converged'] is True}
  br={r['case']:r for r in rows if r['shot']==shot and r['campaign']==b and r['converged'] is True}
  values=[]
  for name in sorted(ar.keys() & br.keys()):
   ds=[]
   for r in [ar[name],br[name]]:
    d=_load_run(r['path']);ds.append(GFileData(d['files']['gfile_after']).gfile_to_xarray())
   ga,gb=ds;qa=np.interp(x,ga.rho_tor,np.abs(ga.q));qb=np.interp(x,gb.rho_tor,np.abs(gb.q));sa=_shear_profile(ga);sb=_shear_profile(gb)
   values.append(dict(case=name,max_relative_q_pct=float(np.max(np.abs(qb/qa-1))*100),max_absolute_shear_difference=float(np.max(np.abs(np.interp(x,*sb)-np.interp(x,*sa))))))
  paired.append(dict(shot=shot,a=a,b=b,both_converged=len(values),cases=values))
(out/'diiid_matched_checks.json').write_text(json.dumps(paired,indent=2))
# Make a four-scenario shear comparison using the same TPED helper.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
for shot in ['DIIID153764','DIIID162940','DIIID174082']:
 fig,axes=plt.subplots(2,2,figsize=(13,8),sharex=True,sharey=True)
 names=['chease','chease_ncscal4','chease_rhotTopPed0p3','chease_rhotTopPed0p3_ncscal4'];limits=[]
 case_names=sorted({r['case'] for r in rows if r['shot']==shot})
 colors={name:plt.get_cmap('turbo')(i/max(len(case_names)-1,1)) for i,name in enumerate(case_names)}
 for ax,cam in zip(axes.flat,names):
  group=[r for r in rows if r['shot']==shot and r['campaign']==cam]
  for r in group:
   d=_load_run(r['path']);g=GFileData(d['files']['gfile_after']).gfile_to_xarray();sh=_shear_profile(g)
   ax.plot(*sh,color=colors[r['case']],ls='-' if r['converged'] else '--',lw=1.1,label=r['case'])
   m=(sh[0]>=.8)&(sh[0]<=.98)&np.isfinite(sh[1]);limits.extend(sh[1][m])
  ax.set_title(f'NCSCAL={group[0]["ncscal"]}, rhotTopPed={group[0]["rhotTopPed"]}\n{sum(r["converged"] is True for r in group)}/{len(group)} saved solver flags pass',fontsize=10)
  ax.set_xlim(.8,.98);ax.grid(alpha=.2);ax.set_xlabel('rho_tor');ax.set_ylabel('magnetic shear')
 span=max(limits)-min(limits);axes[0,0].set_ylim(min(limits)-.06*span,max(limits)+.06*span)
 handles=[plt.Line2D([],[],color=colors[n],label=n) for n in case_names]
 fig.legend(handles=handles,loc='center left',bbox_to_anchor=(.79,.5),fontsize=8)
 fig.suptitle(f'{shot}: NCSCAL and core-change comparison',fontsize=14)
 fig.tight_layout(rect=(0,.04,.79,.94));fig.text(.02,.015,'Same colors and axes across panels. Dashed colored curves: solver flag false. Missing cases are not counted as failures.',fontsize=9)
 fig.savefig(out/f'{shot}_settings_comparison.png',dpi=150);plt.close(fig)
# Check identical source data for old/new NSTX.
hashes=[]
for shot in ['129015','129038']:
 old=_load_run(str(root/'129015_129038/pm0.2_ncscal4_warmup0'/shot/'omt_1.000'))['files']
 new=_load_run(str(root/'lleppin_NSTX_chease_tests'/('NSTX'+shot)/'variations/omt1p0_omne1p0'))['files']
 for kind,a,b in [('source',old['gfile_before'],new['gfile_before'])]+[(f'profile_{s}',old['profiles_before'][s],new['profiles_before'][s]) for s in 'eiz']:
  ha,hb=[hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest() for p in [a,b]]
  hashes.append(dict(shot=shot,kind=kind,old=a,new=b,old_sha256=ha,new_sha256=hb,equal=ha==hb))
(out/'nstx_source_hashes.json').write_text(json.dumps(hashes,indent=2))
text='''# DIII-D and NSTX exported CHEASE-BS comparison - 2026-09-28

Replotted **323 saved cases in 27 campaigns** with TPED's `compare_cheasebs_runs.py`: full diagnostic figures, pedestal zooms, pressure-chain plots, and compact q/shear views. The 27 scan configs request 354 cases; 31 lack exported case metadata. Missing cases are not assigned a failure status. Source data were not changed and no new equilibria were solved.

## Main observations

- **The new NSTX scans have a clear shear response.** All 14 cases for each shot report convergence: 129015 takes 5-14 iterations, 129038 takes 6-11. This differs from the earlier NCSCAL=4 I-star scans, which usually stopped in two iterations.
- **A matched-factor comparison changes the interpretation.** For omt factors 0.8, 0.9, 1.0, 1.1 at omne=1 and rho_tor=0.9, shear spans 16.98% in new 129015 versus 3.97% previously, and 7.33% in new 129038 versus 2.90% previously. Shear decreases with omt in both new scans, while it increases in both old scans at this radius. Span means (maximum-minimum)/abs(reconstructed unity shear).
- **This is not an isolated test of NCSCAL or replay representation.** Both NSTX sets use NCSCAL=4, and source EQDSKs plus all three reference profiles are byte-identical. But the new runs use jparallel/rhot, max_iter=20, and 0.001 Ip/bootstrap/q tolerances; the older set uses I-star/rhop, max_iter=50, 0.002 Ip and 0.01 bootstrap/q tolerances, with different code/build provenance and scan implementation. Their reconstructed unity states also differ.
- **The DIII-D metadata qualify the slides' convergence wording.** For NCSCAL=4 and rhotTopPed=0.3, saved flags pass in 6/14 cases for 153764 and 14/14 for 162940 and 174082. The eight 153764 failures stop at the 12-iteration cap. In NCSCAL=1/rhotTopPed=0.3, 174082 has only six saved cases, all passing; that is not evidence that all configured cases converged.
- **Fixed-amplitude controls need separate interpretation.** Their outer Ip gate can remain unsatisfied by design. A false solver flag in those controls is not, by itself, proof that bootstrap/q iteration failed. The exports contain final scalars but no full iteration histories, so trajectory/closure claims remain limited.

[Matched NSTX figure](chease_matched/nstx_matched_old_new.png) | [Exact matched values](chease_matched/nstx_matched_old_new.json) | [Identical-source checks](nstx_source_hashes.json)

![Matched NSTX scan](chease_matched/nstx_matched_old_new.png)

## Direct DIII-D comparisons

[153764 settings figure](DIIID153764_settings_comparison.png) | [162940 settings figure](DIIID162940_settings_comparison.png) | [174082 settings figure](DIIID174082_settings_comparison.png)

Paired checks below use only cases with passing saved flags in both campaigns, interpolated onto 181 common rho_tor points over 0.8-0.98. Differences are maxima over those radii and matched cases. A visually similar q profile need not have identical shear; the outer shear comparison is especially sensitive to differentiation and resolution. These are output comparisons, not numerical-error bounds.

| Shot | Comparison | Matched passing cases | Max q difference | Max absolute shear difference |
|---|---|---:|---:|---:|
'''
for pair in paired:
 values=pair['cases'];label='NCSCAL 1 vs 4, top=0.8' if pair['a']=='chease' else 'top=0.8 vs 0.3, NCSCAL=4'
 text+=f"| {pair['shot']} | {label} | {len(values)} | {max(v['max_relative_q_pct'] for v in values):.3f}% | {max(v['max_absolute_shear_difference'] for v in values):.3f} |\n"
text+='''
[Case-by-case paired checks](diiid_matched_checks.json).

## Plot index and coverage

Every campaign links to a compact q/shear figure and the original TPED diagnostic layout. Each also has `plot.log`. Failed saved flags are identified in diagnostic legends and dashed in compact/settings figures. Campaigns remain separate to avoid mixing constraints. Pressure-chain plots show paired omt/omne scans as scatter points, without a misleading one-parameter connecting line or fitted slope.

The full diagnostics retain TPED's delta-panel rho_tor<=0.95 cap. Absolute shear and pressure-chain samples extend farther out; those outer values need grid/derivative checks before physics interpretation. Pressure-chain x axes are solved/source pressure, not solved/requested pressure error. Compact q/shear plots call the same TPED parser and shear helper as the diagnostics.

| Shot | Campaign | NCSCAL | top | Saved / configured | Passing flags | Iterations | Plots |
|---|---|---:|---:|---:|---:|---|---|
'''
for c in campaigns:
 shot=c['shot'];cam=c['campaign'];p=out/shot/cam
 cfgpath=next(root.glob(f'lleppin_*/{shot}/{cam}/scan_config.json'));cfg=json.loads(cfgpath.read_text());expected=len(cfg.get('modifications',{}))
 pngs=list(p.glob('cheasebs_comparison*'));full=next(x for x in pngs if 'rho' not in x.name and 'pchain' not in x.name);zoom=next(x for x in pngs if 'rho' in x.name);chain=next(x for x in pngs if 'pchain' in x.name)
 links=' / '.join(f'[{label}]({path.relative_to(out).as_posix()})' for label,path in [('q/shear',p/'q_shear.png'),('full',full),('zoom',zoom),('pressure',chain)])
 text+=f"| {shot} | {cam} | {c['ncscal']} | {c['rhotTopPed']} | {c['cases']}/{expected} | {c['converged']} | {c['iterations'][0]}-{c['iterations'][1]} | {links} |\n"
text+='''
TPED implementation: [b5170e7](https://github.com/drdrhatch/TPED/commit/b5170e7) on `gene_datatree`.

## Reproduce / use different paths

Run from the TPED checkout (PowerShell example):

```powershell
.\\.venv\\Scripts\\python.exe projects/discharge_tools/scripts/compare_cheasebs_runs.py `
  --baseline C:/Users/joesc/git/cheaseBS_tests/data/lleppin_NSTX_chease_tests/NSTX129015 `
  --variations C:/Users/joesc/git/cheaseBS_tests/data/lleppin_NSTX_chease_tests/NSTX129015/variations `
  --outdir C:/Users/joesc/git/cheaseBS_tests/reports/nstx129015-comparison `
  --pedestal 0.8,0.98 --slices 0.85 0.9 0.95 0.975
```

- `--baseline DIR`: directory containing the source gfile and unmodified `profiles_*` or `REF_profiles_*`. Use the shot directory, not the exported `variations/baseline` folder, which contains metadata and a preview but no source EQDSK.
- `--baseline-gfile FILE`: select an explicit baseline EQDSK when needed; reference profiles come from `--baseline` or the file's parent. This also allows an explicitly chosen reconstructed reference. Label that reference accordingly when presenting the result.
- `--variations DIR ...`: campaign roots or explicit run directories. Existing positional paths still work. Use one shot per explicit baseline.
- `--slices RHO ...`: pressure-chain radii without alpha contouring. Existing `--chapman` and `--alpha` remain supported; an explicit baseline controls Chapman fitting too.
- Metadata resolves relocated output/profile basenames and supplies settings, transform labels and final summaries. No synthetic iteration histories are constructed. Namelist NCSCAL is read locally where not explicitly overridden in metadata.
- DIII-D 174082 exposed a pre-existing GFileData validation bug: the identifier before the dimensions was required to be 3. Its valid source uses 0. Validation and dimension parsing now accept other integer identifiers and both 15/16-character fields; 174082 itself parses as width 16.

All 15 focused regression tests passed, including legacy layouts, both profile naming conventions, explicit overrides, wrong-shot/ambiguous baseline rejection, final summaries without traces, and EFIT header variants.

[Commands used](commands.json) | [Case CSV](case_values.csv) | [Case JSON](case_values.json) | [Campaign summaries](campaign_summary.json) | [Generator](generate_report.py)

Regenerate the campaign plots with the same Python environment:

```powershell
python generate_report.py --data-root C:/Users/joesc/git/cheaseBS_tests/data --tped C:/Users/joesc/git/TPED --outdir . --resume
```

Omit `--resume` to rebuild every campaign. `matched_chease_comparison.py` writes the older/new NSTX comparison into `chease_matched/`. `finish_report.py` rebuilds paired checks, settings figures and this index for this dated dataset.
'''
(out/'README.md').write_text(text,encoding='utf-8')
print('Report finished:',out,'PNG count',len(list(out.rglob('*.png'))))
