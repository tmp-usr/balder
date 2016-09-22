class Comparor(object):
    def __init__(self, otu_data, index_controls, index_treatments):
        self.otu_data= otu_data   
        self.g1= otu_data.ix[index_controls]
        self.g2= otu_data.ix[index_treatments]

    
    def compare(self, method="mwu"):
        if method == "mwu":
            return mannWhitneyUTest(self.g1.T, self.g2.T)

    
    def mannWhitneyU(self):




    def result(self):
            
        ## we have to perform the below operations over normalized values
        mean1= self.otu_data[self.g1.index].T.mean()
        mean2= self.otu_data[self.g2.index].T.mean()
        std1= self.otu_data[self.g1.index].T.std()
        std2= self.otu_data[self.g2.index].T.std()
        fold_change= mean1.apply(np.log2) / mean2.apply(np.log2)

        df_result = concat([mean1, std1, mean2, std2, fold_change, s_p_values, s_fdr], axis=1 )
        df_result.columns= ["mean1", "std1", "mean2", "std2", "log fold_change", "p_values", "fdr"]
        
        return df_result



            #control= metadata.query('Triclosan == 0.00')['Triclosan'].index
            #treatments= metadata.groupby('Triclosan')

            #treatment_groups= {g:v for g,v in treatments.groups.iteritems() if (len(v) >1) & (g > 0.0)}
            
            
            ### means and fold changes will be calculated from this 
            #df_relative= self.to_relative(df_otu).T 



            #treatment_group_otus= {}
            
            #for treatment, samples in treatment_groups.iteritems():
            #    g1= df_comparison.ix[control]
            #    g2= df_comparison.ix[samples]
            

                
                
                
              
                #treatment_group_otus[str(treatment)]= all_results_below.index

                ### Test for over-representation: Fisher's exact test
                
                #len_results_below= len(results_below)
                #len_results_above= len(results_above)
                

