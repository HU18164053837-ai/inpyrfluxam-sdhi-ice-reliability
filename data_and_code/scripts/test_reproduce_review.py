"""Minimal standard-library regression tests; no analysis runs on import."""
from pathlib import Path
import ast,math,re,unittest
source=Path(__file__).with_name('reproduce_review.py').read_text(encoding='utf8')
node=next(x for x in ast.parse(source).body if isinstance(x,ast.FunctionDef) and x.name=='compare')
scope={'math':math,'re':re};exec(compile(ast.Module(body=[node],type_ignores=[]),'comparator','exec'),scope)
compare=scope['compare']
class ComparatorTests(unittest.TestCase):
 def test_identity_prefixed_numeric_metrics(self):
  for field in ['identity_median_FE','identity_MAE','identity_mean_FE']:
   with self.subTest(field=field):compare([{field:'1.9714285714285715'}],[{field:'1.9714285714285713'}])
 def test_strict_identity_and_category(self):
  for field in ['identity','compound','surrogate','target','hazard_category','record_id','species_name']:
   with self.subTest(field=field):
    with self.assertRaises(AssertionError):compare([{field:'1.00000000000001'}],[{field:'1.0'}])
 def test_integer_count_exact(self):
  for actual in ['3','2.00000000000001']:
   with self.assertRaises(AssertionError):compare([{'count':actual}],[{'count':'2'}])
 def test_excess_tolerance_rejected(self):
  with self.assertRaises(AssertionError):compare([{'identity_median_FE':'1.971429'}],[{'identity_median_FE':'1.9714285714285713'}])
 def test_record_order_and_schema(self):
  with self.assertRaises(AssertionError):compare([{'compound':'B'},{'compound':'A'}],[{'compound':'A'},{'compound':'B'}])
  with self.assertRaises(AssertionError):compare([{'other':'1.0'}],[{'value':'1.0'}])
 def test_nonfinite_and_missing_rejected(self):
  for value in ['NaN','inf','']:
   with self.assertRaises(AssertionError):compare([{'metric':value}],[{'metric':'1.0'}])
if __name__=='__main__':unittest.main(verbosity=2)
