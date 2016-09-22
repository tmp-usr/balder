import re
import numpy as np

level_labels= ["domain","superkingdom","phylum","class","order","family","genus","species","strain"]

class DataPrep(object):
    """
Data preparation involves 
    1: the standardization of sample identifiers in metadata and otu files.
    2: setting up level labels for the taxonomy file
    
    """
    
    def __init__(self,df_metadata, df_otu, df_taxonomy, ssu_type= "16S", ranks= ranks):
        self.df_metadata= df_metadata
        self.df_abundance= df_abundance
        self.df_taxonomy = df_taxonomy
        
        self.ssu_type= ssu_type
        self.ranks= ranks


    def f_metadata(self):
        self.df_metadata.index= map(str.lower, df_metadata.index)


    def f_otu(self, level_labels):
        self.df_abundance.index= map(str.lower, df.index)


    def f_taxonomy(self):
        for rank in self.ranks: self.df_taxonomy[rank]= None
        for i in df_taxonomy.index:
            for rank in self.ranks:
                try:
                    self.df_taxonomy[rank][i]= self.df_taxonomy[1][i].\
                    split(';')[self.ranks.index(rank)]
                except:
                    continue
        self.df_taxonomy= self.df_taxonomy.drop(1, axis= 1)






trash = """

def f_metadata_16S(df):
    df.index= map(str.lower, df.index)
    return df

def f_otu_1"6S(df):
    df.index= map(str.lower, df.index)
    df.index= [s.replace('16s','') for s in df.index]
    return df
    
def f_taxonomy_16S(df):
    for col in cols: df[col]= None
    for i in df.index:
        for col in cols:
            try:
                df[col][i]= df[1][i].split(';')[cols.index(col)]
            except:
                continue
    df= df.drop(1,axis= 1)
    return df


def f_metadata_18S(df):
    df.index= map(str.lower, df.index)
    df.index= [s.replace('18s','') for s in df.index]
    return df


def f_otu_18S(df):
    df.index= map(str.lower, df.index)
    df.index= [s.replace('18s','') for s in df.index]
    
    exp= re.compile('[a-z]*_[0-9]{1,2}[a-z]')
    for index in df.index:
        results= re.findall(exp,index)
        if results != []:
            org= results[0][:-1]
            results.insert(0,org)
            df.ix[org]= np.ceil(df.ix[results].mean())
            df=df.drop(results[1:],axis=0)
    return df
    
    
def f_taxonomy_18S(df):
    cols= ["domain","division", "kingdom","phylum","class","order","family","genus","species","strain"]
    for col in cols: df[col]= None
    for i in df.index:
        for col in cols:
            try:
                df[col][i]= df[1][i].split(';')[cols.index(col)]
            except:
                continue
    return df.drop(1,axis= 1)

"""

