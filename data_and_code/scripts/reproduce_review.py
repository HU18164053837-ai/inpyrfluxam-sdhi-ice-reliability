"""Portable reproduction of saved-point statistics; no toxicity prediction."""
from pathlib import Path
import csv,shutil,subprocess,sys,json,math,re
R=Path(__file__).resolve().parents[2];I=R/'data_and_code/inputs';D=R/'data_and_code/reproduction_check'
for x in ['scripts','config','inputs','results','evidence']:(D/x).mkdir(parents=True,exist_ok=True)
shutil.copy2(R/'data_and_code/legacy_scripts/phase1g_b/analyze.py',D/'scripts/analyze.py')
for p in (R/'data_and_code/config').glob('analysis_plan*'):shutil.copy2(p,D/'config'/p.name)
for p in (I/'phase1g_b_source').glob('*.csv'):shutil.copy2(p,D/'inputs'/p.name)
subprocess.run([sys.executable,str(D/'scripts/analyze.py')],check=True)
def read(p):
 with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def compare(actual, frozen):
 assert len(actual)==len(frozen), 'Record count differs'
 for i,(a,b) in enumerate(zip(actual,frozen)):
  assert list(a)==list(b), 'Column schema/order differs'
  for field,x in a.items():
   y=b[field]
   if x==y: continue
   # IDs, labels, integers, booleans and missing values stay exact.
   identity=bool(re.search(r'(^|_)(id|name|species|cas|unit|version)($|_)',field,re.I)) or field.lower() in {'identity','surrogate','target','compound','model','direction','record','endpoint'}
   # Categories and frozen integer-valued fields must remain exact.
   strict_category=bool(re.search(r'(^|_)(category|class)$',field,re.I))
   frozen_integer=bool(re.fullmatch(r'[+-]?\d+',y))
   floating=not (identity or strict_category or frozen_integer) and any(c in x+y for c in '.eE')
   try:
    fx,fy=float(x),float(y)
   except ValueError:
    floating=False
   assert floating and math.isfinite(fx) and math.isfinite(fy) and math.isclose(fx,fy,rel_tol=1e-12,abs_tol=1e-12), (i,field,x,y)
count=0
for p in (D/'results').glob('*.csv'):
 compare(read(p),read(I/'phase1g_b'/p.name))
 count+=1
assert count==12, ('Expected 12 tables',count)
print('PASS: 12 frozen statistical tables match; floats rel_tol=1e-12, abs_tol=1e-12; other fields and record order exact.')
