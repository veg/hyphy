import os
import subprocess
import glob

DATA_DIR = "manuscript_fubar_invariance/data/mammalian"
HYPHY_PATH = "./hyphy"
BSTILL_BF = "res/TemplateBatchFiles/SelectionAnalyses/B-STILL.bf"

def run_cmd(cmd):
    print(f"Executing: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error in command: {' '.join(cmd)}")
        print(f"Stderr: {result.stderr}")
    return result.stdout

def process():
    nex_files = glob.glob(os.path.join(DATA_DIR, "*.nex"))
    
    for nex in sorted(nex_files):
        gene = os.path.basename(nex).replace(".nex", "")
        print(f"\n=== Processing {gene} ===")
        
        out_json = os.path.join(DATA_DIR, f"{gene}.B-STILL.json")
        # Run B-STILL
        # Using the same parameters as before
        cmd = [
            HYPHY_PATH, BSTILL_BF,
            "--code", "Universal",
            "--alignment", nex,
            "--output", out_json,
            "--grid", "20",
            "--method", "Variational-Bayes",
            "--posterior", "0.9",
            "--non-zero", "No",
            "--ebf", "10",
            "--radius-threshold", "0.5"
        ]
        run_cmd(cmd)
        print(f"Completed {gene}")

if __name__ == "__main__":
    process()
