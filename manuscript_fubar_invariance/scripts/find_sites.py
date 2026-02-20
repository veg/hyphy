from Bio.Seq import Seq

def get_taxon_seq(nex_path, taxon_idx):
    with open(nex_path) as f:
        lines = f.readlines()
        matrix_start = -1
        for i, line in enumerate(lines):
            if "MATRIX" in line.upper():
                matrix_start = i + 1
                break
        if matrix_start == -1: return ""
        seq_lines = [l.strip() for l in lines[matrix_start:] if l.strip() and not l.strip().startswith(";")]
        if taxon_idx < len(seq_lines):
            line = seq_lines[taxon_idx]
            parts = line.split()
            if len(parts) > 1: return parts[1]
            return parts[0]
        return ""

# TP53 hg38 is at index 66 (0-indexed)
hg38_tp53 = get_taxon_seq("manuscript_fubar_invariance/data/mammalian/TP53.nex", 66)
aa_tp53 = str(Seq(hg38_tp53.replace("-", "N")).translate())
# R175 is in 'VVVRC' (human index 175)
p175 = aa_tp53.find("VVVRC") + 3 # VVVRC -> 171-175. R is 175.
# RPIL starts at R248
p248 = aa_tp53.find("RPIL")
# FEVRRCP: R is 273.
p273 = aa_tp53.find("FEVRR") + 3

print(f"TP53 site 175: {p175 + 1}")
print(f"TP53 site 248: {p248 + 1}")
print(f"TP53 site 273: {p273 + 1}")

# ACE2 hg38
hg38_ace2 = get_taxon_seq("manuscript_fubar_invariance/data/mammalian/ACE2.nex", 66)
aa_ace2 = str(Seq(hg38_ace2.replace("-", "N")).translate())
p_hemgh = aa_ace2.find("HEMGH")
print(f"ACE2 HEMGH (374-378): {p_hemgh + 1}-{p_hemgh + 5}")
