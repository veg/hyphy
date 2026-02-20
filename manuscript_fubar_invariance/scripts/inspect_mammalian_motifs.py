import json
import os
from Bio.Nexus import Nexus
from Bio.Seq import Seq

DATA_DIR = "manuscript_fubar_invariance/data/mammalian"

def find_motif(nex_path, motif_aa):
    nex = Nexus.Nexus(nex_path)
    # Use first taxon as reference
    ref_seq = str(nex.matrix[nex.taxlabels[0]])
    ref_aa = str(Seq(ref_seq.replace("-", "N")).translate())
    
    pos = ref_aa.find(motif_aa)
    if pos != -1:
        # pos is 0-indexed AA position
        return range(pos, pos + len(motif_aa))
    return []

def inspect():
    # 1. ACE2: HEMGH motif
    ace2_nex = os.path.join(DATA_DIR, "ACE2.nex")
    ace2_json = os.path.join(DATA_DIR, "ACE2.B-STILL.json")
    hemgh_pos = find_motif(ace2_nex, "HEMGH")
    print(f"\nACE2 HEMGH positions: {list(hemgh_pos)}")
    if hemgh_pos:
        with open(ace2_json) as f:
            data = json.load(f)
        content = data["MLE"]["content"]["0"]
        for p in hemgh_pos:
            row = content[p]
            print(f"  Site {p+1}: alpha={row[0]:.3f}, beta={row[1]:.3f}, EBF_prox={row[12]:.2f}")

    # 2. TP53: R175, R248, R273 (using human numbering as proxy if match found)
    # Human TP53: R175 (CGC), R248 (CGG), R273 (CGT)
    # We'll just look for CGC, CGG, CGT in the alignment and see EBFs
    tp53_nex = os.path.join(DATA_DIR, "TP53.nex")
    tp53_json = os.path.join(DATA_DIR, "TP53.B-STILL.json")
    print(f"\nTP53 Arginine (CpG) positions:")
    with open(tp53_json) as f:
        data = json.load(f)
    content = data["MLE"]["content"]["0"]
    nex = Nexus.Nexus(tp53_nex)
    ref_seq = str(nex.matrix[nex.taxlabels[0]])
    for i in range(len(content)):
        codon = ref_seq[i*3:i*3+3].upper()
        if codon in ["CGC", "CGG", "CGT", "CGA"]:
            row = content[i]
            if row[12] > 50: # Only show high constraint
                print(f"  Site {i+1} ({codon}): alpha={row[0]:.3f}, beta={row[1]:.3f}, EBF_prox={row[12]:.2f}")

    # 3. XBP1: Synonymous constraint (alpha=0)
    xbp1_json = os.path.join(DATA_DIR, "XBP1.B-STILL.json")
    print(f"\nXBP1 Synonymous Constraint (EBF[alpha=0] >= 50):")
    with open(xbp1_json) as f:
        data = json.load(f)
    content = data["MLE"]["content"]["0"]
    for i, row in enumerate(content):
        # Index 10 is EBF[alpha=0]
        if row[10] >= 50:
            print(f"  Site {i+1}: alpha={row[0]:.3f}, beta={row[1]:.3f}, EBF[a=0]={row[10]:.2f}, EBF_prox={row[12]:.2f}")

if __name__ == "__main__":
    inspect()
