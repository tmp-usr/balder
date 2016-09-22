import os
import shutil

from pandas import Series

from data_processor import DataProcessor
from stats.stats import deseq2Corr, deseq2Test, deseq2_regression
from stats.stats import spearmanCorr, mwuTest 
from stats.representator import Representator
from stats.microbial_ecology.taxon_miner import TaxonMiner
from helper.radapter import RAdapter
from reporter.report import Report

import pdb

class Analysis(object):
    """
        tables:
            + Use one threshold for the diversity plots and the nmds 
            OTU, phylum etc. taxon abundance:
                1. correlations along the gradient
                2. low concentration comparisons
                3. high concentration comparisons
                for 
                    a. bacteria
                    b. eukaryotes

    
        plots:
            1. community response
            2. nmds
            3. diversity
            4. stacked bar
            5. rarefaction
    
    """
    #### constants
    global FDR_CUT_OFF
    FDR_CUT_OFF= 0.2

    def __init__(self, metadata_filepath, otu_data_filepath, annotation_taxonomy_filepath,  
                       metadata_index_col, ssu_type, analysis_type, level):
        
        ### init attributes
        self.dp= DataProcessor(metadata_filepath, otu_data_filepath, annotation_taxonomy_filepath,  
                       metadata_index_col, ssu_type)

        self.r = RAdapter()
        self.reporter= Report(ssu_type, analysis_type, level)

        #### init variables
        self.analysis_type= analysis_type 
        self.ssu_type= ssu_type
        self.level= level
        
        ### analyses
        if analysis_type == "level_deseq2":
            self.build_level_tables()
        elif analysis_type == "level_mwu":
            self.build_level_tables(comparison='mwu', regression="spearman")
        elif analysis_type == "fisher":
            self.build_fisher_tables()
        
        #self.build_plots()

    
    def build_fisher_tables(self):
        
        #### below step is done to extract all e.g phyla at phylum level
        #self.dp.get_abundance_by_level(self.level)

        otu_data= self.dp.otu_data
        metadata= self.dp.metadata.ix[otu_data.index]
        
        #### phyla names to be used in the fisher's test
        
        self.taxon_miner= TaxonMiner(self.dp.annotation_taxonomy)
        tmp, richness= self.taxon_miner.get_taxa_richness_by_level(otu_data, self.level)
        taxa= richness.index
        ####

        #### regression ####
        result_corr= deseq2_regression(self.r.fCorrelateDeseq2, otu_data.T, metadata)
        
        self.representator= Representator(self.dp.annotation_taxonomy, self.taxon_miner) 
        fisher_decreasing_increasing_corr= self.representator.test_fisher(result_corr, "regression", self.level, taxa)
        # continue with a fisher's exact test

        ####### comparative analyses ########
        control_samples= self.dp.get_control_samples()
        treatment_groups= self.dp.get_treatment_groups()

        control_group= otu_data.ix[control_samples]
        
        ### hypothesis tests ###
        fisher_treatments= {}
        for treatment, treatment_samples in treatment_groups.iteritems():
            treatment_group= otu_data.ix[treatment_samples]
            result_comp= deseq2Test(self.r.fTestDeseq2, control_group, treatment_group)
           
            fisher_decreasing_increasing_treatment= self.representator.test_fisher(result_comp, "comparison", self.level, taxa)
            fisher_treatments[treatment]= fisher_decreasing_increasing_treatment
    
        
        self.reporter.write_table(fisher_decreasing_increasing_corr[0], self.level, self.level, False)
        self.reporter.write_table(fisher_decreasing_increasing_corr[1], self.level, self.level, True)

        for treatment, decreasing_increasing_comparison in fisher_treatments.iteritems():
            self.reporter.write_table(decreasing_increasing_comparison[0], self.level, self.level, False, treatment)
            self.reporter.write_table(decreasing_increasing_comparison[1], self.level, self.level, True, treatment)



    def build_level_tables(self, comparison= "deseq2", regression="deseq2"):
        
        #### below step is done to extract all e.g phyla at phylum level
        #self.dp.get_abundance_by_level(self.level)

        otu_data= self.dp.otu_data
        metadata= self.dp.metadata.ix[otu_data.index]
        
         
        #### get level abundances
        self.taxon_miner= TaxonMiner(self.dp.annotation_taxonomy)
        taxon_data= self.taxon_miner.get_level_taxa_by_otu_data(otu_data, self.level)
        taxon_data= taxon_data.groupby('taxon').sum().T
        
        ### writing the abundance ranks ###
        taxon_data_rel= taxon_data.T
        taxon_data_rel= taxon_data_rel/taxon_data_rel.sum()
        taxon_ranks= taxon_data_rel.T.mean().order(ascending=False)
        self.reporter.write_abundance(taxon_ranks)    
         
        #### regression ####
        if regression == "deseq2":
            result_corr= deseq2_regression(self.r.fCorrelateDeseq2, taxon_data.T, metadata)
        elif regression == "spearman":
            result_corr= spearmanCorr(taxon_data_rel, metadata['Triclosan'])

        result_corr= result_corr[result_corr['pvalue'] < FDR_CUT_OFF]
        
        result_corr_decreasing= result_corr[result_corr['stat'] < 0.0]
        result_corr_increasing= result_corr[result_corr['stat'] > 0.0]
        
        
        ####### comparative analyses ########
        control_samples= self.dp.get_control_samples()
        treatment_groups= self.dp.get_treatment_groups()

        if comparison == "deseq2":
            control_group= taxon_data.ix[control_samples]
        
        elif comparison == "mwu":
            control_group= taxon_data_rel.T.ix[control_samples]


        ### hypothesis tests ###
        level_comparisons= {}
        for treatment, treatment_samples in treatment_groups.iteritems():
            
            if comparison == "deseq2":
                treatment_group= taxon_data.ix[treatment_samples]
                result_comp= deseq2Test(self.r.fTestDeseq2, control_group, treatment_group)
            
            elif comparison == "mwu":
                treatment_group= taxon_data_rel.T.ix[treatment_samples]
                result_comp= mwuTest(control_group, treatment_group)

            result_comp= result_comp[result_comp['pvalue'] < FDR_CUT_OFF]

            result_comp_less= result_comp[result_comp['mean2'] < result_comp['mean1'] ]
            result_comp_more= result_comp[result_comp['mean2'] > result_comp['mean1'] ]
            
            level_comparisons[treatment]=[result_comp_less, result_comp_more]
            
        
        #### adding the abundance rank of the taxa to the results
        self.add_ranks_to_results(result_corr_decreasing, taxon_ranks)
        self.add_ranks_to_results(result_corr_increasing, taxon_ranks)
        
        #### writing results
        self.reporter.write_table(result_corr_decreasing, self.level, self.level, False)
        self.reporter.write_table(result_corr_increasing, self.level, self.level, True)

        for treatment, less_more_comparison in level_comparisons.iteritems():
            
            #### adding the abundance rank of the taxa to the results
            self.add_ranks_to_results(less_more_comparison[0], taxon_ranks)
            self.add_ranks_to_results(less_more_comparison[1], taxon_ranks)
            
            #### writing results
            self.reporter.write_table(less_more_comparison[0], self.level, self.level, False, treatment)
            self.reporter.write_table(less_more_comparison[1], self.level, self.level, True, treatment)

    
    def add_ranks_to_results(self, results, otu_ranks):
        rank_list= [] 
        for taxon in results.index:
            rank= otu_ranks.index.tolist().index(taxon)+ 1
            rank_list.append(rank)
        results['rank']=Series(rank_list, index= results.index)
        return results


    def build_plots(self):
        pass




