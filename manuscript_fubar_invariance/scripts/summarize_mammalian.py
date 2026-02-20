import json
import os
import glob
import pandas as pd

DATA_DIR = "manuscript_fubar_invariance/data/mammalian"

def summarize():
    json_files = glob.glob(os.path.join(DATA_DIR, "*.B-STILL.json"))
    
    summary = []
    for jf in sorted(json_files):
        gene = os.path.basename(jf).replace(".B-STILL.json", "")
        with open(jf) as f:
            data = json.load(f)
        
        content = data["MLE"]["content"]["0"]
        # Index 12 is proximal EBF
        # Index 9 is alpha=beta=0 EBF
        # Index 10 is alpha=0 EBF
        # Index 11 is beta=0 EBF
        
        ebfs = [row[12] for row in content]
        max_ebf = max(ebfs)
        mean_ebf = sum(ebfs) / len(ebfs)
        sig_sites = len([e for e in ebfs if e >= 10])
        
        summary.append({
            "Gene": gene,
            "Sites": len(ebfs),
            "Sig Sites (EBF>=10)": sig_sites,
            "Mean EBF": round(mean_ebf, 2),
            "Max EBF": round(max_ebf, 2)
        })
        
        # Specific checks
        if gene == "ACE2":
            # HEMGH motif is usually around 374-378 in human? 
            # We need to find where it is in the alignment.
            pass
            
    df = pd.DataFrame(summary)
    print(df.to_string(index=False))

if __name__ == "__main__":
    summarize()
