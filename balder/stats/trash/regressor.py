class Regressor(object):
    
    def __init__(self, method="deseq2", df_otu, metadata):
        self.result= self.regress(method, df_otu, metadata)


    def regress(self, method, df_otu, metadata):
        if method == "spearman":
            pass

        elif method == "deseq2":
            df_corr= self.r.correlateDeseq2(df_otu, metadata)
        
        corr_p_values= df_corr['p_value'] 
        ####
        corr_fdr = correctBHFDR(corr_p_values)
        tmp_fdr= [None]*len(corr_fdr)
         
        for i,value in corr_fdr.iteritems():
            tmp_fdr[i]= value 

        ###

        #!!!results deseq2 increasing and decerasing

        return DataFrame({"p-value":corr_p_values, "fdr": corr_fdr, "stat": df_corr['stat'] }, index= df_otu.index).sort('fdr')

