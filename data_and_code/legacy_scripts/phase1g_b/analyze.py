"""Bounded Phase1G-B analysis of frozen point predictions. No model fitting/calls."""
from pathlib import Path
import csv,json,math,statistics as st,hashlib,collections,datetime
ROOT=Path(__file__).resolve().parents[1]
P=ROOT/'config/analysis_plan_frozen.json'
assert hashlib.sha256(P.read_bytes()).hexdigest()==(ROOT/'config/analysis_plan_frozen.sha256').read_text().strip()
PLAN=json.loads(P.read_text(encoding='utf8'));B=[100.,1000.,10000.,100000.]
NAMES=['Very highly toxic','Highly toxic','Moderately toxic','Slightly toxic','Practically nontoxic']
def read(p):return list(csv.DictReader(p.open(encoding='utf-8-sig')))
def save(name,rows):
 if not rows:return
 keys=list(dict.fromkeys(k for r in rows for k in r))
 with (ROOT/'results'/name).open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
def rank(v,alternative=False):
 assert math.isfinite(v) and v>0
 if alternative:return 1+sum(v>b for b in B)
 if v<100:return 1
 if v<=1000:return 2
 if v<=10000:return 3
 if v<=100000:return 4
 return 5
def near(v):
 d,b=min((abs(math.log10(v/b)),b) for b in B)
 return dict(nearest_threshold_ug_L=b,threshold_distance_abs_log10=d,nearest_threshold_fold_distance=10**d,within_factor2_threshold=d<=math.log10(2),exact_threshold=v in B)
def num(v):
 try:x=float(v);return x if math.isfinite(x) else None
 except (ValueError,TypeError):return None
def flag(v):return str(v).lower()=='true'
def mean(a):return st.mean(a) if a else ''
def med(a):return st.median(a) if a else ''
def tier(c):return 'EXTERNAL_CANDIDATE' if c in PLAN['primary_compounds'] else 'CONDITIONAL_INTERNAL_SENSITIVITY'
raw=read(ROOT/'inputs/all90_predictions_frozen.csv');trial=read(ROOT/'inputs/toxicity_trial_evidence_frozen.csv')
bytrial={(r['compound'],r['species']):r for r in trial if flag(r['in_prediction_panel'])}
point={}
for r in raw:
 for s,v in [(r['target'],float(r['observed_ug_L'])),(r['surrogate'],float(r['input_ug_L']))]:
  key=(r['compound'],s)
  if key in point:assert math.isclose(point[key],v,rel_tol=1e-12)
  else:point[key]=v
end=[];byend={}
for (c,s),v in sorted(point.items()):
 t=bytrial.get((c,s),{});lo=num(t.get('CI95_lower_ug_L'));hi=num(t.get('CI95_upper_ug_L'))
 if t:assert math.isclose(v,float(t['frozen_LC50_ug_L']),rel_tol=1e-12)
 if lo is not None and hi is not None:
  assert 0<lo<=v<=hi
  cats=list(range(rank(lo),rank(hi)+1));status='CI_SINGLE_CLASS' if len(cats)==1 else 'CI_MULTIPLE_CLASSES'
 else:lo=hi=None;cats=[];status='CI_MISSING'
 notes=t.get('concentration_stability_HSE','')
 issues=[]
 if 'Departed' in notes:issues.append('MEASURED_CONCENTRATION_DEVIATION')
 if c=='inpyrfluxam' and s=='Lepomis macrochirus':issues.append('REPLICATE_MORTALITY_VARIATION')
 alt=num(t.get('alternative_summary_LC50_ug_L'))
 e=dict(compound=c,validation_tier=tier(c),species=s,observed_LC50_ug_L=v,observed_LC50_mg_L=v/1000,point_class_rank=rank(v),point_class=NAMES[rank(v)-1],CI95_lower_ug_L=lo if lo is not None else '',CI95_upper_ug_L=hi if hi is not None else '',CI_status=status,CI_possible_class_ranks=';'.join(map(str,cats)),CI_possible_classes=';'.join(NAMES[k-1] for k in cats),CI_spans_or_touches_cutoffs=';'.join(str(b) for b in B if lo is not None and lo<=b<=hi),specific_quality_limitation=bool(issues),specific_quality_flags=';'.join(issues),quality_notes=notes,original_report_available=False,original_report_id=t.get('original_report_id','NOT_RESOLVED_TRAINING_AGGREGATE'),source=t.get('frozen_source','Web-ICE standardized species-compound aggregate'),CI_source=t.get('alternative_source','NO_EXPERIMENTAL_CI_AVAILABLE') if cats else 'NO_EXPERIMENTAL_CI_AVAILABLE',source_url=t.get('supplementary_source_url','https://www3.epa.gov/webice/documents/Chemicals_for_Aquatics.xlsx'),concentration_basis=t.get('basis_HSE','NOT_ESTABLISHED_FROM_AGGREGATE'),life_stage_and_water_conditions='NOT_FULLY_VERIFIED',alternative_regulatory_point_ug_L=alt if alt is not None else '',alternative_point_class_rank=rank(alt) if alt is not None else '',alternative_point_class_unchanged=(rank(alt)==rank(v)) if alt is not None else '',CI_interpretation='Target assay confidence set only; no individual probabilities, source-input or model uncertainty propagation',**near(v))
 e.update(duration_h=t.get('duration_h','NOT_RESOLVED_FROM_STANDARDIZED_AGGREGATE'),exposure_regime=t.get('regime','NOT_RESOLVED_FROM_STANDARDIZED_AGGREGATE'),test_material=t.get('test_material','TRAINING_STANDARDIZED_COMPOUND_AGGREGATE'),reported_purity_percent=t.get('purity_percent_HSE','NOT_REPORTED'),reported_GLP=t.get('GLP_HSE','NOT_RESOLVED_FROM_AGGREGATE'),regulator_reliability=t.get('ERMA_reliability','NOT_REASSIGNED'))
 end.append(e);byend[(c,s)]=e
save('endpoint_uncertainty_assessment.csv',[e for e in end if e['validation_tier']=='EXTERNAL_CANDIDATE'])
save('internal_endpoint_uncertainty_assessment.csv',[e for e in end if e['validation_tier']!='EXTERNAL_CANDIDATE'])
def assess(e,p):
 k=rank(p);obs=e['point_class_rank'];mis=k!=obs
 cats=[int(x) for x in e['CI_possible_class_ranks'].split(';') if x]
 if not cats:status='POINT_MISMATCH_CI_MISSING' if mis else 'POINT_MATCH_CI_MISSING';outside=None;lower=None
 else:
  outside=k not in cats;lower=k>max(cats)
  status=('POINT_MISMATCH_TARGET_CI_INCOMPATIBLE' if outside else 'POINT_MISMATCH_CI_COMPATIBLE') if mis else 'POINT_MATCH_CI_COMPATIBLE'
 return dict(predicted_class_rank=k,predicted_class=NAMES[k-1],point_class_mismatch=mis,signed_class_shift=k-obs,absolute_class_shift=abs(k-obs),category_hazard_underestimated=k>obs,category_hazard_overestimated=k<obs,target_CI_incompatible=outside if outside is not None else '',hazard_underestimated_for_all_target_CI_classes=lower if lower is not None else '',uncertainty_status=status)
haz=[];base=[];calc=[]
for r in raw:
 c=r['compound'];s=r['surrogate'];t=r['target'];e=byend[(c,t)];se=byend[(c,s)]
 obs=float(r['observed_ug_L']);pred=float(r['predicted_ug_L']);identity=float(r['input_ug_L'])
 a=assess(e,pred);z=assess(e,identity)
 h=dict(compound=c,validation_tier=tier(c),validation_original=r['validation'],surrogate=s,target=t,observed_LC50_ug_L=obs,predicted_LC50_ug_L=pred,observed_class_rank=e['point_class_rank'],observed_class=e['point_class'],**a,target_CI_status=e['CI_status'],target_CI_lower_ug_L=e['CI95_lower_ug_L'],target_CI_upper_ug_L=e['CI95_upper_ug_L'],target_CI_possible_class_ranks=e['CI_possible_class_ranks'],target_specific_quality_limitation=e['specific_quality_limitation'],target_quality_flags=e['specific_quality_flags'],surrogate_CI_status=se['CI_status'],surrogate_specific_quality_limitation=se['specific_quality_limitation'],guidance_pass=flag(r['guidance_pass']),in_domain=flag(r['in_domain']),MSE=float(r['MSE']),R2=float(r['R2']),slope=float(r['slope']),n_training=int(float(r['n_training'])),**{'observed_'+k:v for k,v in near(obs).items()},**{'predicted_'+k:v for k,v in near(pred).items()})
 haz.append(h)
 we=math.log10(pred/obs);ie=math.log10(identity/obs)
 assert math.isclose(we,float(r['signed_log10_error']),abs_tol=1e-12)
 assert math.isclose(10**abs(we),float(r['fold_error']),rel_tol=1e-12)
 b=dict(compound=c,validation_tier=tier(c),surrogate=s,target=t,observed_LC50_ug_L=obs,surrogate_observed_LC50_ug_L=identity,webice_prediction_ug_L=pred,identity_prediction_ug_L=identity,webice_fold_error=10**abs(we),identity_fold_error=10**abs(ie),webice_signed_log10_error=we,identity_signed_log10_error=ie,webice_absolute_log10_error=abs(we),identity_absolute_log10_error=abs(ie),delta_absolute_log10_error_webice_minus_identity=abs(we)-abs(ie),delta_fold_error_webice_minus_identity=10**abs(we)-10**abs(ie),webice_LC50_overestimated=pred>obs,identity_LC50_overestimated=identity>obs,observed_class_rank=rank(obs),webice_class_rank=rank(pred),identity_class_rank=rank(identity),webice_class_mismatch=a['point_class_mismatch'],identity_class_mismatch=z['point_class_mismatch'],delta_class_mismatch_webice_minus_identity=int(a['point_class_mismatch'])-int(z['point_class_mismatch']),webice_hazard_category_underestimated=a['category_hazard_underestimated'],identity_hazard_category_underestimated=z['category_hazard_underestimated'],webice_target_CI_incompatible=a['target_CI_incompatible'],identity_target_CI_incompatible=z['target_CI_incompatible'],webice_uncertainty_status=a['uncertainty_status'],identity_uncertainty_status=z['uncertainty_status'],target_CI_status=e['CI_status'],surrogate_CI_status=se['CI_status'],guidance_pass=h['guidance_pass'],in_domain=h['in_domain'],MSE=h['MSE'],n_training=h['n_training'],target_within_factor2_threshold=e['within_factor2_threshold'],target_specific_quality_limitation=e['specific_quality_limitation'])
 base.append(b);calc.append((h,b))
save('hazard_classification_results.csv',[e for e in haz if e['validation_tier']=='EXTERNAL_CANDIDATE'])
save('baseline_comparison.csv',[e for e in base if e['validation_tier']=='EXTERNAL_CANDIDATE'])
save('internal_hazard_classification_results.csv',[e for e in haz if e['validation_tier']!='EXTERNAL_CANDIDATE'])
save('internal_baseline_comparison.csv',[e for e in base if e['validation_tier']!='EXTERNAL_CANDIDATE'])
scenarios=PLAN['scenarios'];selected={}
for scenario in scenarios:
 name=scenario.split(' ')[0];z=[]
 for h,b in calc:
  keep=True
  if name in ['GUIDANCE_AND_DOMAIN','STRICT_QUALITY']:keep=h['guidance_pass'] and h['in_domain']
  if name=='STRICT_QUALITY':keep=keep and h['MSE']<=.5 and h['n_training']>=10
  if name=='TARGET_AWAY_FACTOR2':keep=not h['observed_within_factor2_threshold']
  if name.startswith('EXCLUDE_BLUEGILL'):
   if h['compound']!='inpyrfluxam':continue
   keep=h['target']!='Lepomis macrochirus' and (name!='EXCLUDE_BLUEGILL_EITHER_ROLE' or h['surrogate']!='Lepomis macrochirus')
  if keep:z.append((h,b))
 selected[name]=z
def summ(pairs):
 h=[h for h,b in pairs];b=[b for h,b in pairs];n=len(b)
 out=dict(n_directions=n,n_targets=len({x['target'] for x in b}),n_surrogates=len({x['surrogate'] for x in b}))
 if not n:return out
 out.update(n_target_CI_available=sum(x['target_CI_status']!='CI_MISSING' for x in b),n_target_CI_single_class=sum(x['target_CI_status']=='CI_SINGLE_CLASS' for x in b),n_target_CI_multiple_classes=sum(x['target_CI_status']=='CI_MULTIPLE_CLASSES' for x in b))
 for m in ['webice','identity']:
  mis=[x for x in b if x[m+'_class_mismatch']];ci=[x for x in mis if x[m+'_target_CI_incompatible'] is True];amb=[x for x in mis if x[m+'_target_CI_incompatible'] is False];missing=[x for x in mis if x[m+'_target_CI_incompatible']=='']
  out.update({m+'_median_FE':med([x[m+'_fold_error'] for x in b]),m+'_MAE_log10':mean([x[m+'_absolute_log10_error'] for x in b]),m+'_median_signed_log10_error':med([x[m+'_signed_log10_error'] for x in b]),m+'_LC50_overestimation_n':sum(x[m+'_LC50_overestimated'] for x in b),m+'_LC50_overestimation_rate':mean([int(x[m+'_LC50_overestimated']) for x in b]),m+'_class_mismatch_n':len(mis),m+'_class_mismatch_rate':len(mis)/n,m+'_class_agreement_rate':1-len(mis)/n,m+'_hazard_category_underestimated_n':sum(x[m+'_hazard_category_underestimated'] for x in b),m+'_hazard_category_underestimated_rate':mean([int(x[m+'_hazard_category_underestimated']) for x in b]),m+'_point_mismatch_target_CI_incompatible_n':len(ci),m+'_point_mismatch_CI_compatible_n':len(amb),m+'_point_mismatch_CI_missing_n':len(missing),m+'_target_CI_incompatible_rate_among_CI_available':len(ci)/out['n_target_CI_available'] if out['n_target_CI_available'] else ''})
 out.update(delta_MAE_webice_minus_identity=out['webice_MAE_log10']-out['identity_MAE_log10'],delta_class_mismatch_rate=out['webice_class_mismatch_rate']-out['identity_class_mismatch_rate'],n_direction_abs_error_webice_better=sum(x['delta_absolute_log10_error_webice_minus_identity']<-1e-12 for x in b),n_direction_abs_error_identity_better=sum(x['delta_absolute_log10_error_webice_minus_identity']>1e-12 for x in b),max_webice_class_shift=max(x['absolute_class_shift'] for x in h),n_webice_multiclass_shifts=sum(x['absolute_class_shift']>=2 for x in h),n_mismatches_with_both_points_away_factor2=sum(x['point_class_mismatch'] and not x['observed_within_factor2_threshold'] and not x['predicted_within_factor2_threshold'] for x in h))
 return out
target=[];surrogate=[];compound=[];macro=[]
for scenario,pairs in selected.items():
 for c in sorted(set(e['compound'] for e in end)):
  if scenario.startswith('EXCLUDE_BLUEGILL') and c!='inpyrfluxam':continue
  q=[x for x in pairs if x[0]['compound']==c];z=dict(scenario=scenario,compound=c,validation_tier=tier(c),**summ(q));compound.append(z)
  targets=sorted(s for cc,s in byend if cc==c)
  ts=[]
  for s in targets:
   sub=[x for x in q if x[0]['target']==s];e=byend[(c,s)];a=dict(scenario=scenario,compound=c,validation_tier=tier(c),target=s,observed_ug_L=e['observed_LC50_ug_L'],observed_class_rank=e['point_class_rank'],target_CI_status=e['CI_status'],target_CI_possible_class_ranks=e['CI_possible_class_ranks'],**summ(sub));target.append(a)
   if sub:ts.append(a)
  for s in targets:
   sub=[x for x in q if x[0]['surrogate']==s];surrogate.append(dict(scenario=scenario,compound=c,validation_tier=tier(c),surrogate=s,**summ(sub)))
  m=dict(scenario=scenario,compound=c,validation_tier=tier(c),n_total_target_endpoints=len(targets),n_covered_targets=len(ts),n_directions=len(q),n_targets_webice_MAE_better=sum(x['delta_MAE_webice_minus_identity']<-1e-12 for x in ts),n_targets_identity_MAE_better=sum(x['delta_MAE_webice_minus_identity']>1e-12 for x in ts),n_targets_MAE_tied=sum(abs(x['delta_MAE_webice_minus_identity'])<=1e-12 for x in ts),n_targets_webice_class_agreement_better=sum(x['delta_class_mismatch_rate']<-1e-12 for x in ts),n_targets_identity_class_agreement_better=sum(x['delta_class_mismatch_rate']>1e-12 for x in ts),n_targets_class_agreement_tied=sum(abs(x['delta_class_mismatch_rate'])<=1e-12 for x in ts),n_targets_with_any_webice_point_mismatch=sum(x['webice_class_mismatch_n']>0 for x in ts),n_targets_with_any_webice_target_CI_incompatible_mismatch=sum(x['webice_point_mismatch_target_CI_incompatible_n']>0 for x in ts),n_targets_all_webice_directions_point_mismatch=sum(x['webice_class_mismatch_n']==x['n_directions'] for x in ts))
  for name in ['webice_MAE_log10','identity_MAE_log10','delta_MAE_webice_minus_identity','webice_class_mismatch_rate','identity_class_mismatch_rate','webice_class_agreement_rate','identity_class_agreement_rate','webice_LC50_overestimation_rate','identity_LC50_overestimation_rate','webice_hazard_category_underestimated_rate','identity_hazard_category_underestimated_rate']:
   m['equal_target_'+name]=mean([x[name] for x in ts]);m['target_min_'+name]=min([x[name] for x in ts]) if ts else '';m['target_max_'+name]=max([x[name] for x in ts]) if ts else ''
  macro.append(m)
save('target_species_summary.csv',target);save('surrogate_species_summary.csv',surrogate);save('compound_direction_summary.csv',compound);save('compound_equal_target_summary.csv',macro)
# Same official cutoffs; assess exact-boundary ownership only, not shifted thresholds.
boundary=[]
for c in sorted(set(x['compound'] for x in haz)):
 q=[x for x in haz if x['compound']==c];b=[x for x in base if x['compound']==c]
 for alt in [False,True]:
  boundary.append(dict(compound=c,validation_tier=tier(c),boundary_rule='ALL_EQUALS_HIGHER_HAZARD' if alt else 'CURRENT_EPA_PRIMARY',n_directions=len(q),webice_mismatch_n=sum(rank(x['predicted_LC50_ug_L'],alt)!=rank(x['observed_LC50_ug_L'],alt) for x in q),identity_mismatch_n=sum(rank(x['identity_prediction_ug_L'],alt)!=rank(x['observed_LC50_ug_L'],alt) for x in b),n_exact_target_cutoff=sum(x['observed_exact_threshold'] for x in q),n_exact_webice_cutoff=sum(x['predicted_exact_threshold'] for x in q)))
save('boundary_ownership_sensitivity.csv',boundary)
scope=[dict(compound=c,validation_tier=tier(c),n_species_compound_endpoints=sum(e['compound']==c for e in end),n_directions=sum(h['compound']==c for h in haz),n_CI_available=sum(e['compound']==c and e['CI_status']!='CI_MISSING' for e in end),n_CI_single_class=sum(e['compound']==c and e['CI_status']=='CI_SINGLE_CLASS' for e in end),n_CI_multiple_classes=sum(e['compound']==c and e['CI_status']=='CI_MULTIPLE_CLASSES' for e in end),n_specific_quality_flags=sum(e['compound']==c and e['specific_quality_limitation'] for e in end)) for c in sorted(set(e['compound'] for e in end))]
save('analysis_units_and_endpoint_coverage.csv',scope)
log=dict(analysis_time_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),plan_sha256=hashlib.sha256(P.read_bytes()).hexdigest(),n_rows=len(haz),primary_hazard_rows=sum(e['validation_tier']=='EXTERNAL_CANDIDATE' for e in haz),internal_hazard_rows=sum(e['validation_tier']!='EXTERNAL_CANDIDATE' for e in haz),n_unique_direction_keys=len({(x['compound'],x['surrogate'],x['target']) for x in haz}),n_unique_endpoints=len(end),external_candidates=2,new_model_fits=0,new_WebICE_queries=0)
(ROOT/'evidence/analysis_execution.json').write_text(json.dumps(log,indent=2),encoding='utf8')
for z in compound:
 if z['scenario']=='ALL_AVAILABLE' and z['validation_tier']=='EXTERNAL_CANDIDATE':print(z['compound'],'directions',z['n_directions'],'point mismatches',z['webice_class_mismatch_n'],'target-CI incompatible',z['webice_point_mismatch_target_CI_incompatible_n'],'median FE ICE/identity',round(z['webice_median_FE'],6),round(z['identity_median_FE'],6))
for z in macro:
 if z['scenario']=='ALL_AVAILABLE' and z['validation_tier']=='EXTERNAL_CANDIDATE':print('TARGET_MACRO',z['compound'],'MAE ICE/identity',round(z['equal_target_webice_MAE_log10'],6),round(z['equal_target_identity_MAE_log10'],6))
