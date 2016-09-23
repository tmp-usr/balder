import glob
import os,sys
import pdb
import operator
import math
from itertools import combinations


import numpy as np
from pandas import *

kallisto_output_dir= "."


dirs= glob.glob(os.path.join( kallisto_output_dir, '*' ))

tpms= []

for dirname in dirs:
    
    sample_name= dirname.lstrip("outputs_")
    tpm_file= os.path.join(dirname, "rnaseq.%s.genes.txt" %sample_name)
    
    try:
        df= read_csv(tpm_file, sep="\t", names= ["genes", sample_name])
        indexes= df['genes']
        indexes_1= indexes.apply(str.split, args= ('.'))
        indexes_2= [i[0] for i in indexes_1]
        df['genes']= indexes_2
        df= df.groupby('genes').sum()
        tpms.append(df)
    
    except Exception, e:
        print e


df_tpms= concat(tpms, axis= 1, join= "outer")

final_counts= np.ceil(df_tpms).astype(int)

final_counts.columns = ["X"+item if item.startswith('2') else item for item in final_counts.columns  ]

sample_names= [item.split('__')[0] for item in final_counts.columns]
final_counts.columns = sample_names
sample_categories= set(sample_names)

for x,y in combinations(sample_categories, 2):    
    sub_df= final_counts[[x,y]]
    res= deSeq2(f, sub_df)
    
    print x,y 
    print res



