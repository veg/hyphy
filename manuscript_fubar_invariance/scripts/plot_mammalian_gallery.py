import json
import os
import matplotlib.pyplot as plt
import numpy as np

DATA_DIR = "manuscript_fubar_invariance/data/mammalian"
GENES = ["ACE2", "TP53", "XBP1", "MECP2", "CFTR", "CALM1"]

def load_ebf(gene):
    path = os.path.join(DATA_DIR, f"{gene}.B-STILL.json")
    with open(path) as f:
        data = json.load(f)
    content = data["MLE"]["content"]["0"]
    # Proximal EBF is index 12
    # alpha=0 EBF is index 10
    return [row[12] for row in content], [row[10] for row in content]

fig, axes = plt.subplots(3, 2, figsize=(15, 12), sharey=False)
axes = axes.flatten()

for i, gene in enumerate(GENES):
    ax = axes[i]
    prox_ebfs, alpha0_ebfs = load_ebf(gene)
    x = np.arange(1, len(prox_ebfs) + 1)
    
    ax.plot(x, prox_ebfs, color='black', alpha=0.7, label='Proximal EBF')
    if gene == "XBP1":
        ax.plot(x, alpha0_ebfs, color='red', alpha=0.5, label='Synonymous EBF')
        ax.legend()
        
    ax.set_title(f"{gene}", fontsize=14)
    ax.set_xlabel("Codon Site")
    ax.set_ylabel("EBF")
    ax.set_yscale('symlog', linthresh=1)
    ax.grid(True, alpha=0.3)
    
    # Highlight specific regions
    if gene == "ACE2":
        # HEMGH around 374-378
        ax.axvspan(374, 378, color='yellow', alpha=0.3, label='HEMGH')
    elif gene == "XBP1":
        # RNA hairpins around 160-170
        ax.axvspan(160, 170, color='red', alpha=0.2, label='RNA Hairpin')
    elif gene == "TP53":
        # Highlight sites 175, 248, 273 if possible, but they are single sites.
        # Just vertical lines
        ax.axvline(175, color='blue', alpha=0.5, linestyle='--', label='R175')
        ax.axvline(248, color='green', alpha=0.5, linestyle='--', label='R248')

plt.tight_layout()
plt.savefig("manuscript_fubar_invariance/mammalian_gallery.png", dpi=300)
print("Saved gallery plot to manuscript_fubar_invariance/mammalian_gallery.png")
