from Bio import SeqIO

def get_hg38_xbp1():
    for record in SeqIO.parse("manuscript_fubar_invariance/data/mammalian/XBP1.fas", "fasta"):
        if record.id.upper() == "HG38":
            return str(record.seq).upper()
    return ""

seq = get_hg38_xbp1()
motif = "GTCCGCAGCAGGTGCAG"
pos = seq.find(motif)
if pos != -1:
    print(f"Motif {motif} found at nt index {pos} (Site {pos//3 + 1})")
    print(f"Surrounding sequence: {seq[pos-10:pos+len(motif)+10]}")
else:
    print("Motif not found")
