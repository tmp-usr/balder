from pandas import DataFrame, Series

class TaxonMiner(object):
    def __init__(self, annotation_taxonomy):
        self.annotation_taxonomy = annotation_taxonomy

    def get_taxa_richness_by_level(self, otu_data, level):
        df= self.get_level_taxa_by_otu_data(otu_data, level)
        #df.to_csv('level_otus.tsv',sep="\t", index_label="Name")
        groups= DataFrame(Series({k:len(v) for k,v in df.groupby('taxon').groups.iteritems()}), columns= ["richness"])
        groups= groups.sort('richness', ascending= False)
        return df.groupby('taxon').sum().T, groups
    
    def get_level_taxa_by_otu_data(self, otu_data, level):
        new_df=DataFrame(index= otu_data.T.index, columns= otu_data.index.insert(0,'taxon'))
        for index in otu_data.T.index:
            taxon = self.annotation_taxonomy.ix[index][level]
            new_df['taxon'][index]= taxon
        #print index
            for col in otu_data.index:
                new_df[col][index]= otu_data[index][col]
        #new_df= new_df.dropna(axis=0, how='all')
        return new_df
    
    def get_taxon_by_otu_data(self, otu_data, level, taxon):
        if level == "unknown":
            new_df=DataFrame(index= otu_data.T.index, columns= otu_data.index.insert(0,'taxon'))
            for index in otu_data.T.index:
                hier= self.annotation_taxonomy.ix[index]
                if taxon in hier.tolist():
                    new_df['taxon'][index] = taxon
                    for col in otu_data.index:
                        new_df[col][index] = otu_data[index][col]
            
            return new_df.dropna()
        else:
            df_taxon= self.get_level_taxa_by_otu_data(otu_data, level)
            df_taxon= df_taxon[df_taxon['taxon'].str.lower() == taxon.lower()].drop('taxon',1)
            df_taxon= df_taxon.convert_objects(convert_numeric=True).T 
            df_taxon=  df_taxon.loc[~(df_taxon==0).all(axis=1)] 
        
        return df_taxon 
    
    def filter_by_total_otu_count(self, count):
        new_df= self.otu_data.ix[self.otu_data.T.sum() > count]
        filtered_sample_names= set(self.otu_data.index).difference(set(new_df.index))
        filtered_samples= df_otu.ix[filtered_sample_names]
        return new_df, filtered_samples


    def rarefy(self):
        pass


    def rarefy2(self):
        """
            taken from the below link
            http://stackoverflow.com/posts/18967204/revisions
        """
        prng = RandomState(seed) # reproducible results
        noccur = np.sum(M, axis=1) # number of occurrences for each sample
        nvar = M.shape[1] # number of variables
        depth = np.min(noccur) # sampling depth

        Mrarefied = np.empty_like(M)
        for i in range(M.shape[0]): # for each sample
            p = M[i] / float(noccur[i]) # relative frequency / probability
            choice = prng.choice(nvar, depth, p=p)
            Mrarefied[i] = np.bincount(choice, minlength=nvar)

        return Mrarefied




