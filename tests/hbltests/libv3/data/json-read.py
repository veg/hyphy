import json 
import sys
import csv

with open (sys.argv[1], 'r') as fh:
     fubar_results = json.load (fh)
     headers = fubar_results['MLE']['headers']
     rows = fubar_results['MLE']['content']["0"]
     csv_writer = csv.writer (sys.stdout)
     csv_writer.writerow ([k[0] for k in headers])
     for r in rows:
        csv_writer.writerow (r)
