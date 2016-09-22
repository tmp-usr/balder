import scipy.stats as st
from pandas import Series
from stats import correctBHFDR

import pdb

class Representator(object):
    """
        handles the detection of over- or under-represented taxa.
        employs Fisher's exact test on OTU data.
    """
    def __init__(self, df_annotation_taxonomy, taxon_miner):
        self.df_annotation_taxonomy = df_annotation_taxonomy
        self.taxon_miner= taxon_miner

    def test_fisher(self, results, method, level, taxa):
        """
            applies fisher's exact test to a dataframe of analysis results having columns
            "p-value","fdr" "stat" and OTUs as index.
        """

        all_results_above= results[results['fdr'] > 0.05]
        all_results_below= results[results['fdr'] <= 0.05]
        ### alternative querying in pandas
        # all_results_below= df_result.query("fdr <= 0.05")
        # all_results_above= df_result.query("fdr > 0.05")

        if method == "regression":
            results_decreasing_below= all_results_below[all_results_below['stat'] < 0.0]
            results_increasing_below= all_results_below[all_results_below['stat'] > 0.0]

            results_decreasing_above= all_results_above[all_results_above['stat'] < 0.0]
            results_increasing_above= all_results_above[all_results_above['stat'] > 0.0]

        elif method == "comparison":
     
            results_decreasing_below= all_results_below[all_results_below['mean2'] < all_results_below['mean1'] ]
            results_increasing_below= all_results_below[all_results_below['mean2'] > all_results_below['mean1'] ]

            results_decreasing_above= all_results_above[all_results_above['mean2'] < all_results_above['mean1'] ]
            results_increasing_above= all_results_above[all_results_above['mean2'] > all_results_above['mean1'] ]

        
        fisher_decreasing_increasing=[] 

        for results_below, results_above in zip((results_decreasing_below, results_increasing_below), (results_decreasing_above, results_increasing_above)):

            taxon_indices= []
            taxon_pvalues= []
            for taxon in taxa:
                ## this has nothing to do with abundances.
                ## we use the function to get the number of OTUs in the 
                ## corresponding phylum.
                
                len_all_results_below= len(all_results_below)
                len_all_results_above= len(all_results_above)
                
                taxon_below= self.taxon_miner.get_taxon_by_otu_data(results_below.T, level, taxon).T
                taxon_above= self.taxon_miner.get_taxon_by_otu_data(all_results_above.T, level, taxon).T
                
                len_taxon_below= len(taxon_below)
                len_taxon_above= len(taxon_above)
            
                ### fisher's exact test: columns: phylum, non-phylum
                ###                      rows   : below, above
                len_non_taxon_below= len_all_results_below - len_taxon_below
                len_non_taxon_above= len_all_results_above - len_taxon_above
                
                oddsratio, p_value_taxon = st.fisher_exact([[len_taxon_below, len_non_taxon_below], 
                                                        [len_taxon_above, len_non_taxon_above]])
                taxon_indices.append(taxon)
                taxon_pvalues.append(p_value_taxon)
            
            #fdr_fisher= correctBHFDR(phylum_pvalues) 
            result_fisher= Series(taxon_pvalues, index= taxon_indices, name="p-value")
            #result_fisher= Series(fdr_fisher, index= phylum_indices, name="fdr")
            result_fisher_final= result_fisher[result_fisher <= 0.2].order()
             
            fisher_decreasing_increasing.append(result_fisher_final)
        
        return fisher_decreasing_increasing

########## reporting ####         
#            if i == 0:
#                fisher_corr_file= os.path.join(self.dirTables, "fisher_decreasing_corr.xls")
#            else:
#                fisher_corr_file= os.path.join(self.dirTables, "fisher_increasing_corr.xls")
#                
#            result_fisher_final.to_csv(fisher_corr_file, sep="\t",index_label="Phylum", header= True, float_format="%.5g")
        
#########################

            #i+=1


####### reporting ########
#        corr_file_increasing= os.path.join(self.dirTables, "deseq2_increasing_corr.xlsx")
#        results_increasing_below.sort('fdr').to_excel(corr_file_increasing, index_label="OTU", header= True, float_format="%.5g")

#        corr_file_decreasing= os.path.join(self.dirTables, "deseq2_decreasing_corr.xlsx")
#        results_decreasing_below.sort('fdr').to_excel(corr_file_decreasing, index_label="OTU", header= True, float_format="%.5g")
        
##########################
########## reporting ##########
                
#                deseq2_file_more= os.path.join(self.dirTables, "deseq2_more_%s.xlsx"% treatment)
#                deseq2_file_less= os.path.join(self.dirTables, "deseq2_less_%s.xlsx"% treatment)
                
#                results_more_below.sort('fdr').to_excel(deseq2_file_more, index_label="OTU", header= True, float_format="%.5g")
#                results_less_below.sort('fdr').to_excel(deseq2_file_less, index_label="OTU", header= True, float_format="%.5g")

###############################

