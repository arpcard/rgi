# Compare BLASTP bitscore to Diamond bitscore RGI results
# In[2]:

import sys
import os
import csv, glob
import pandas as pd


def main():

    # In[3]:

    files = glob.glob('*-diamond.txt')
    print(files)
    files = [i.split('-diamond.txt')[0] for i in files]

    print(files)
    # exit("Debug...")
    # In[4]:

    all_missing = list()
    all_added = list()
    
    for fh in files:
        try:
            blast_file = open('%s.txt' % (fh), 'r')
            diamond_file = open('%s-diamond.txt' % (fh), 'r')
            blast_file.close()
            diamond_file.close()
        except FileNotFoundError:
            print("No file: %s" % (fh))
            pass
        else:
            with open('%s.txt' % (fh), 'r') as blast_file, open('%s-diamond.txt' % (fh), 'r') as diamond_file:
                blast_file = csv.reader(blast_file, delimiter='\t')
                h1 = next(blast_file, None)
                diamond_file = csv.reader(diamond_file, delimiter='\t')
                h2 = next(diamond_file, None)
                blast_genes = dict((row[8], [row[2],row[3]]) for row in blast_file)
                diamond_genes = dict((row[8], [row[2],row[3]]) for row in diamond_file)

            missing_genes = [k for k,v in blast_genes.items() if v not in diamond_genes.values()]
            added_genes = [k for k,v in diamond_genes.items() if v not in blast_genes.values()]
            all_missing.extend(missing_genes)
            all_added.extend(added_genes)
    else:
        all_missing = list(set(all_missing))
        all_added = list(set(all_added))

    # In[5]:

    print(set(all_missing).symmetric_difference(all_added))

if __name__ == '__main__':
    main()