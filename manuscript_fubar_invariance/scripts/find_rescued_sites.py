import json
from Bio import SeqIO

def find_proximal_rescued(gene):
    json_path = f"manuscript_fubar_invariance/data/mammalian/{gene}.B-STILL.json"
    fas_path = f"manuscript_fubar_invariance/data/mammalian/{gene}.fas"
    with open(json_path) as f:
        data = json.load(f)
    content = data["MLE"]["content"]["0"]
    records = list(SeqIO.parse(fas_path, "fasta"))
    
    rescued = []
    for i, row in enumerate(content):
        # row[12] is proximal EBF
        # row[9] is alpha=beta=0 EBF
        if row[12] >= 10 and row[9] < 1:
            col = [str(r.seq[i*3:i*3+3]).upper() for r in records if len(r.seq) >= i*3+3]
            col = [c for c in col if "-" not in c and "N" not in c]
            if len(set(col)) > 1:
                rescued.append((i+1, row[12], col))
    return rescued

print("MECP2 rescued sites:")
for site, ebf, col in find_proximal_rescued("MECP2")[:5]:
    variants = {c: col.count(c) for c in set(col)}
    print(f"  Site {site}: EBF_prox={ebf:.2f}, variants={variants}")

print("\nCFTR rescued sites:")
for site, ebf, col in find_proximal_rescued("CFTR")[:5]:
    variants = {c: col.count(c) for c in set(col)}
    print(f"  Site {site}: EBF_prox={ebf:.2f}, variants={variants}")
