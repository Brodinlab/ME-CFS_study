import os
import pandas as pd
import sys
path=sys.argv[1]
gtf_file=sys.argv[2]
counts=0
tpm=0
tipe='HSA'
fopen=open(gtf_file,'r').readlines()
def find_id(string,lst):
    try:
        x=len(string)+lst.index(string)+2
        y=lst[x:].index('\";')+x
        res=lst[x:y]
    except ValueError:
        res=float('NaN')
    return res
t2g=[i.rstrip('\n').split('\t')[8] for i in fopen if ('\ttranscript' in i)]
t2g=[i for i in t2g if 'protein_coding' in i]

t2g_2=[]
for i in t2g:
    temp=[find_id('gene_id',i),find_id('gene_version',i),find_id('gene_name',i),find_id('transcript_id',i)+'.'+find_id('transcript_version',i),find_id('gene_biotype',i),find_id('transcript_biotype',i)]
    t2g_2.append(temp)
t2g=pd.DataFrame(t2g_2)
t2g.columns=['gene_id','gene_version','gene_name','transcript_id','gene_biotype','transcript_biotype']
t2g.index=t2g.transcript_id
path='%s/Result' % path
t2g=t2g[t2g.transcript_biotype == 'protein_coding']
t2g=t2g[t2g.gene_biotype == 'protein_coding']
files=os.listdir(path)
files = [i for i in files if ((os.path.isdir(path+i)) & ~(i.startswith('.')))]
# files=[i for i in files if ('boot_' not in i)]
for i in files:
    if ('.ipynb_checkpoints' == i):
        continue
    temp=pd.read_table(path+i+'/abundance.tsv',index_col=0)
    if type(counts) == int:
        counts = temp.est_counts
        tpm=temp.tpm
    else:
        counts=pd.concat([counts,temp.est_counts],1)
        tpm=pd.concat([tpm,temp.tpm],1)
counts.columns=files
tpm.columns=files
counts.to_csv(path+'raw_counts_transcript.txt',sep='\t')
tpm.to_csv(path+'TPM_transcript.txt',sep='\t')
pd.concat([counts,t2g['gene_id']],1).groupby('gene_id').sum().to_csv(path+'raw_counts_gene_ID.txt',sep='\t')
pd.concat([counts,t2g['gene_name']],1).groupby('gene_name').sum().to_csv(path+'raw_counts_gene_name.txt',sep='\t')
pd.concat([tpm,t2g['gene_id']],1).groupby('gene_id').sum().to_csv(path+'TPM_gene_ID.txt',sep='\t')
pd.concat([tpm,t2g['gene_name']],1).groupby('gene_name').sum().to_csv(path+'TPM_gene_name.txt',sep='\t')
