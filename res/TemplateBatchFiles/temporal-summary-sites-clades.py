import sys
import json
import argparse
import operator
import csv
import re
import os
from datetime import datetime

arguments = argparse.ArgumentParser(description='Summarize selection analysis over time.')

arguments.add_argument('-d', '--directory',   help = 'Directory to scan', required = True, type = str)
arguments.add_argument('-n', '--no_header',   help = 'Do not produce CSV headers', required = False, action = 'store_true')
settings = arguments.parse_args()

report_name = re.compile ("sequences\.([^\.]+)\.FEL\.json")
clades = {}

for root, dirs, files in os.walk(settings.directory):
    for analysis_file  in files:
        m = report_name.match (analysis_file)
        if m:
            dir_name = os.path.split (os.path.split (root)[0])[1]
            if not dir_name in clades:
                clades [dir_name] = {}

            clades[dir_name][m.group(1)] = os.path.join (root, analysis_file)
            

writer = csv.writer (sys.stdout)
if not settings.no_header:
    writer.writerow (['gene','clade','site','seqs', 'p', 'alpha','beta','T_int','T_total'])

for clade, files in clades.items():
    for gene, meme in files.items():
        if gene != 'range':
       
            with open (meme, "r") as fh:
                try:
                    meme_json = json.load (fh)
                    T = 0
                    TT = 0
                    for bl,cl in meme_json["tested"]["0"].items():
                        if cl == "test":
                            T += meme_json["branch attributes"]["0"][bl]["Global MG94xREV"]
                        TT += meme_json["branch attributes"]["0"][bl]["Global MG94xREV"]
                    seqs = str (meme_json ["input"]["number of sequences"])
                    for i, row in enumerate (meme_json["MLE"]["content"]["0"]):
                        if row [4] <= 0.05:
                            writer.writerow ([gene, clade, "%d" % (i+1), seqs, "%g" % row[4], "%g" % row[0], "%g" % row[1], "%g" % T, "%g" % TT])
                except:
                    print ("Failed loading %s" % meme, file = sys.stderr)
            
    
                
    
