import os
import pandas as pd
import sys
path=sys.argv[1]
suff_1=sys.argv[2]
suff_2=sys.argv[3]
index=sys.argv[4]
if os.path.isdir('%s/Result/' % path):
    pass
else:
    os.makedirs('%s/Result/' % path)
files1=sorted([i for i in os.listdir(path) if suff_1 in i])
files2=sorted([i for i in os.listdir(path) if suff_2 in i])
x=[]
for i in range(len(files1)):
    x.append('kallisto quant -i %s -o %s/Result/%s %s %s &' % (index,path,files1[i][:-11],folder+files1[i],folder+files2[i])) 
open('kallisto.sh','a+').write('\n'.join(x))
