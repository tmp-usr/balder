class DiversityIndex(object):
    
    def __init__(self):
        pass

    
    def calculate_nmds(self):    
        r_nmds_data= pandas_df_to_r_df(nmds_data)
        r_analysis_type= self.r.initRVariable('analysis_type', analysis_type)
        r_nmds= self.r.fMetaMDS(r_nmds_data, r_analysis_type)
        r_scores= self.r.fMDSScores(r_nmds) 
       
        stress= list(r_nmds.rx('stress')[0])[0]
        scores= r_matrix_to_dataframe(r_scores)
        return scores


    def calculate_indices(self):
        #df= r_matrix_to_dataframe(diversity_data)

        #df_merged_data= df.join(self.df_metadata)
        metadata= metadata.ix[diversity_data.index]
        r_diversity_data= pandas_df_to_r_df(diversity_data)
        indice_data= self.r.fBuildDiversityIndices(r_diversity_data)
        
        #index_data= map(np.array, indice_data)
        
        Indice_Data= dict(zip(("chao1", "shannon", "richness", "pielou"), indice_data))
        for indice, Index in self.r.divIndexes.iteritems():
            indice_data= Indice_Data[indice]
            ind= self.r.divIndexes[indice]
            ind.data= indice_data
            
            ####### try ######## 
            ind_data= list(ind.data)
            names= diversity_data.index
        
            s1= Series(index= names, data= ind_data )
            s2= metadata
            
            df_indice= concat([s2,s1], axis=1)
            df_indice.columns= ['metadata', indice] 

            #### except pass ####
            
            
            #if r:
            #    self.r.fPlotDiversityIndex(ind.type_short, ind.type_long, metadata, ind.data, ind.title)
            #else:
            
            #self.plotDiversityIndex(ind.type_short, ind.type_long,metadata_type, df_indice, ind.data)
            

###### reporting ######

            #if indice == "richness":
            #        fRichness= os.path.join(self.dirTables,  'richness.tsv') 

            #       ind_data= list(ind.data)
            #        names= list(ind.data.names)
                    
            #        s1= Series(index= names, data= ind_data )
            #        s2= metadata
                    
            #        df_indice= concat([s2,s1], axis=1)
            #        df_indice.columns= ['metadata', indice] 

            #        df_indice.sort('metadata').to_csv(fRichness, sep='\t', index_label="Sample")
                    #fOut.write('\n\n'.join((line1, line2)))

#######################




    def write_indices(self):
        pass

    #def plot_indices(self):
    #    pass

