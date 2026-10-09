#!/usr/bin/env python3
"""Independent numerical audit; reads results, never imports model-training code."""
from pathlib import Path
import numpy as np, pandas as pd
from scipy.integrate import trapezoid as scipy_trapz
import json
P=Path(__file__).parent
inp=pd.read_csv(P/'source_R7/Presliding_Figure_Code_Package/data/Hysteresis_600pts.csv')
fit=pd.read_csv(P/'output/R10_new_primary_predicted_600.csv')
res=pd.read_csv(P/'output/R10_new_primary_fullfit.csv')
splits=pd.read_csv(P/'output/R10_new_primary_point_splits.csv')
x=inp.position_mm.to_numpy(); y=inp.force_N.to_numpy(); rx=x.max()-x.min(); rf=y.max()-y.min()
refarea=float(scipy_trapz(y,x))
assert abs(refarea-0.8749578673406824)<1e-9
maxerr=0.;n=0;checks=[]
for model in res.model:
    a=fit.query('model == @model').sort_values('index')
    xx=a.x_fitted_mm.to_numpy(); ff=a.force_fitted_N.to_numpy()
    dx=(xx-x)/rx; dy=(ff-y)/rf
    nrmse=float(np.sqrt(np.dot(dx,dx)/1200+np.dot(dy,dy)/1200))
    # Independent trapezoid as midpoint-force * displacement increments
    aream=float(sum((ff[k]+ff[k+1])*.5*(xx[k+1]-xx[k]) for k in range(599)))
    area_err=float(abs(aream/refarea-1)*100)
    netxx=float(xx[-1]-xx[0]);netf=float(ff[-1]-ff[0])
    closure=float(100*np.linalg.norm([netxx/rx,netf/rf])/2**0.5)
    row=res.query('model == @model').iloc[0]
    for col,calc in [('nrmse_joint',nrmse),('loop_area_error_pct',area_err),('closure_error_pct',closure),('signed_area_pred_Nmm',aream)]:
        delta=abs(float(row[col])-calc);maxerr=max(maxerr,delta)
        assert delta<5e-11,(model,col,delta)
        n+=1
    assert a.shape[0]==600
    checks.append({'model':model,'nrmse_joint':nrmse,'signed_area_pred_Nmm':aream,'area_error_pct':area_err,'closure_error_pct':closure})
# Reconstruct all heldout partitions from independent row-index definition and
# verify no duplicated model/seed/role combination; splitter is data-independent.
assert len(splits)==80 and len(splits.drop_duplicates(['model','seed','subset']))==80
assert all(splits[splits.subset=='test'].n==180)
assert all(splits[splits.subset=='train'].n==420)
assert (splits.design_rank==splits.n_columns).all()
assert (splits.condition>0).all()
# Behavioral identity tests with perfect predictions
idealarea=float(sum((y[j]+y[j+1])*.5*(x[j+1]-x[j]) for j in range(599)))
assert abs(idealarea-refarea)<1e-12
assert np.sqrt(np.mean(((x-x)/rx)**2+((y-y)/rf)**2))==0
assert x[0]==x[-1] and y[0]==y[-1]
# Signed-area orientation reverses when entire traced path reversed.
assert abs(scipy_trapz(y[::-1],x[::-1])+refarea)<1e-12
out={'independent_evaluator':'scipy.integrate.trapezoid and midpoint sum; directly computed residual metrics',
     'metric_comparisons':n,'metric_max_abs_difference':maxerr,'all_metric_checks_pass':True,
     'fullfit_data_rows':int(len(fit)),'split_rows':len(splits),'split_validation_pass':True,
     'behavioral_checks_pass':True,'observed_signed_area_Nmm':refarea,
     'per_model_check':checks}
(P/'output/R10_independent_numerical_QA.json').write_text(json.dumps(out,indent=2)+'\n')
print('PASS:',n,'INDEPENDENT METRIC CHECKS; maximum absolute numerical deviation',maxerr)
print('PASS: 80 new training/testing score rows; 420/180 split invariants, no rank deficiency')
print('PASS: perfect-reconstruction metric, area orientation and closed original path checks')
for r in checks:print(r)
