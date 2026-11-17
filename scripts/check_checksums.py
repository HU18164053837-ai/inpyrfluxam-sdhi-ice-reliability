from pathlib import Path
import hashlib
root=Path(__file__).resolve().parents[1]
count=0
for line in (root/'provenance/checksums.sha256').read_text().splitlines():
 expected,name=line.split('  ',1)
 assert hashlib.sha256((root/name).read_bytes()).hexdigest()==expected,name
 count+=1
print('PASS:',count,'public files verified')
