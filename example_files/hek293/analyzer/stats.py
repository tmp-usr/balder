import math
from collections import namedtuple,OrderedDict

import scipy.stats as st
import numpy as np
from pandas import DataFrame, Series, concat 
import rpy2.robjects as r
from pandas.rpy.common import convert_to_r_dataframe
import statsmodels.sandbox.stats.multicomp as mt

import pdb

from numpy import array, empty                                                                        
#from balder.helper.radapter import pandas_df_to_r_df

msg1= '%s is a single group comparison test. Check the number of groups in the parameter and re-run.'
msg2= '%s is a two group comparison test. Check the number of groups in the parameter and re-run.'
msg3= '%s is a multiple group comparison test. Check the number of groups in the parameter and re-run.'

#### Correction



def correctBonferroni(pValues):
    ''' Bonferroni correction implementation '''    
    corrected = []
    for pValue in pValues:
      correctedValue = pValue * len(pValues)
      corrected.append(correctedValue)
    return corrected



def correctBHFDR(pValues):
    values = [ (pvalue, i) for i, pvalue in enumerate(pValues) ]                                      
    values.sort()
    values.reverse()                                                                                  
    new_values = []
    for i, vals in enumerate(values):                                                                 
        rank = n - i
        pvalue, index = vals                                                                          
        new_values.append((n/rank) * pvalue)                                                          
    for i in xrange(0, int(n)-1):  
        if new_values[i] < new_values[i+1]:                                                           
            new_values[i+1] = new_values[i]                                                           
    for i, vals in enumerate(values):
        pvalue, index = vals
        new_pvalues[index] = new_values[i]                                                                                                                  
    return new_pvalues



def correctBHFDR(pValues):
    """ Benjamini-Hochberg FDR
        p_values: keys, index ; values, p_value 
    """
    p_values= {k:v for k,v in enumerate(pValues)}
    nComparison= len(p_values)
    modifier= nComparison
    
    indexed_p_values= OrderedDict({(j[1],i):j[0] for i, j in enumerate(p_values.items())})
    
    indexed_p_values= OrderedDict(sorted(indexed_p_values.items(), reverse=True))
    
    #print indexed_p_values
    for (value, i),ind in indexed_p_values.items():
        #print i,ind    
        new_value= p_values[ind] * nComparison / float(modifier)
        p_values[ind]= new_value if new_value <=1 else 1 
        modifier -= 1 

    fdr= [0]*len(p_values)
     
    for i,value in p_values.iteritems():
        fdr[i]= value 

    

    #fdr= mt.multipletests(pValues, alpha= 0.2, method="fdr_bh")[1]
    #fdr= pValues 
    return fdr






def deseq2_two_tailed(r_function, g1, g2):
    
    ### prep data for test 
    g1= g1.dropna()
    g2= g2.dropna()

    joint_index= g1.T.index.union(g2.T.index)
    group_data= concat([g1.T,g2.T], axis=1, join_axes= [joint_index])
    
    group_list= []
    for i,g in zip((0,1),(g1, g2)): group_list+= [np.int8(i)]* len(g.T.columns)
    
    r_group_list = r.IntVector(group_list)
    #r_group_list= r(group_list)
    #try:
    r_group_data= convert_to_r_dataframe(group_data)
    #except:
    #    continue

    r_p_values= r_function(r_group_data, r_group_list)
    
    p_values= list(r_p_values)
    return joint_index, [None]*len(joint_index), p_values


def deseq2_regression(r_function, df_otu, df_metadata):
    """
        pandas v > 0.14.0 ruins this function
        do not upgrade
    """
    r_df_otu= convert_to_r_dataframe(df_otu) 
    r_df_metadata= convert_to_r_dataframe(df_metadata)
    r_results= r_function(r_df_otu, r_df_metadata)
    test_stat= list(r_results[0])
    p_value= list(r_results[1])
    fdr= list(r_results[2])
    return DataFrame({ "stat":test_stat,"pvalue":p_value, "fdr":fdr }, index= df_otu.index)

    
def compareDataframes(r_function, function, *groups):
    g1= groups[0]
    g2= groups[1]
    #Results=[]
    #Result=namedtuple('Result','index pvalue teststat mean1 std1 mean2 std2 foldchange')
    joint_index, test_stat, p_values= function(r_function, *groups)
    
    g1_r= (g1.T/g1.T.sum()).T
    g2_r= (g2.T/g2.T.sum()).T
    
    
    ResultDf= DataFrame({ "pvalue":p_values,"stat":test_stat, "mean1":g1_r.mean(), "std1":g1_r.std(), "mean2": g2_r.mean(), "std2":g2_r.std()}, index= joint_index)
    ResultDf['fdr']= Series(correctBHFDR(ResultDf['pvalue']), index=joint_index)
    fold_change= ResultDf["mean2"].apply(np.log10) - ResultDf['mean1'].apply(np.log10)
    
    ResultDf['fold_change']= fold_change
    

    #ResultDf['fdr'] = ResultDf['pvalue']
    #pdb.set_trace()
    ### add fold change
    #return ResultDf.sort_values(by="p_value") # never versions
    
    return ResultDf.sort("pvalue")
    
def compareIndexes(function, g1, g2):
    
    Results=[]
    Result=namedtuple('Result','index pvalue stat mean1 std1 mean2 std2 fold_change')
    
    for index, samples in g1.T.iterrows():
        group1=samples
        mean1= group1.mean()
        std1= group1.std()
        
        group2=g2.T.ix[index]
        std2= group2.std()
        mean2= group2.mean()
        ######
        
        fold_change= np.log10(mean2) - np.log10(mean1)
        
        try:
            tstat, pvalue= function(group1, group2)
            Results.append(Result(index, pvalue, tstat, mean1, std1, mean2, std2, fold_change))
        except:
            tstat, pvalue= 0, 0
    
    ResultDf= DataFrame(Results, columns= Result._fields).set_index('index')
    ResultDf['bonferroni']= Series(correctBonferroni(ResultDf['pvalue']), index=ResultDf.index)
    ResultDf['fdr']= Series(correctBHFDR(ResultDf['pvalue']), index=ResultDf.index)
    #return ResultDf.sort_values(by='fdr')
    return ResultDf.sort('fdr')



def applyCorrelation(function, g1, metadata, r_function= None):
    Results=[]
    Result=namedtuple('Result','index pvalue stat mean std')
    
    for index, samples in g1.iterrows():
        sample1=samples
        mean1= sample1.mean()
        std1= sample1.std()
        try:
            if r_function is not None:
                rho, pvalue= function(r_function, sample1, metadata)
            else:
                rho, pvalue= function(sample1, metadata)
            Results.append(Result(index, pvalue, rho, mean1, std1))
        except:
            rho, pvalue= 0, 0
    
    ResultDf= DataFrame(Results, columns= Result._fields).set_index('index')
    ResultDf['bonferroni']= Series(correctBonferroni(ResultDf['pvalue']), index=ResultDf.index)
    ResultDf['fdr']= Series(correctBHFDR(ResultDf['pvalue']), index=ResultDf.index)
    #return ResultDf.sort_values(by= 'fdr') # newer versions
    return ResultDf.sort('fdr')



def filterResults(Results, pvalue=0.05, mean= 0.0, foldchange=1):
    return Results[(Results['fdr'] <= pvalue) & 
                     (Results['mean1'] >= mean) &
                     (Results['mean2'] >= mean) &
                     (Results['foldchange'] >= foldchange )
                     ]


def tTest(*args ):
    
    if len(args) != 2:
        raise Exception(msg2 % 't-test')
   
    var= "equal"

    if var=='equal':
        return compareIndexes(st.ttest_ind, *args)
    else: # alternative would be var='unequal'
        return compareIndexes(welchTest, *args)

def welchTest(*args):
    ''' calculates welch's t-test for unequal sample sizes and unequal variances'''
    
    if len(args) != 2:
        raise Exception(msg2 % 'Welch\'s')
    
    x1 = np.mean(a)
    x2 = np.mean(b)
    v1 = np.var(a)
    v2 = np.var(b)
    n1 = len(a)
    n2 = len(b)
    t = (x1-x2)/np.sqrt(v1/n1 + v2/n2)
    df  = (v1/n1 + v2/n2)**2
    df /= ((v1/n1)**2/(n1-1.) +(v2/n2)**2/(n2-1.))
    prob = st.t.sf(np.abs(t), df)*2
    return t, prob



def deseq2Test(r_function, *args):
    if len(args) != 2:
        raise Exception(msg2 % 'DeSeq2 Test')
    return compareDataframes( r_function, deseq2_two_tailed, *args)


def mwuTest(*args):
    if len(args) != 2:
        raise Exception(msg2 % 'Mann Whitney U')
    return compareIndexes(st.mannwhitneyu, *args)

##  Tests for equal variances
def bartletTest(*args):
    if len(args) != 2:
        raise Exception(msg2 % 'Bartlet')
   
    return  compareIndexes(st.bartlett, *args)

def leveneTest(*args):
    if len(args) != 2:
        raise Exception(msg2 % 'Levene')
    return compareIndexes(st.levene, *args)

##  Tests for Normality

def shapiroTest(*args):
    if len(args) != 1:
        raise Exception(msg1 % 'Shapiro')

    return  st.shapiro(args)


#### Correlation
def pearsonCorr(x,y):
    return applyCorrelation(st.pearsonr, x , y)

def spearmanCorr(x,y):
    return applyCorrelation(st.spearmanr, x,y)

def deseq2Corr(r_function, group, metadata):
    """
        deprecated
    """
    return applyCorrelation(r_function, deseq2_regression, group, metadata)
      





def ApplyCorrelation(dataSet, metadata, corrType, pCutOff = 0.05, rCutOff = 0.00):
    corrDataContent={}
    i=0
    for feature,v in dataSet.iteritems():
        if corrType == 0:
            result= pearsonCorr(metadata, v)
        else:
            result= spearmanCorr(metadata, v) 
        
        
        if result[1] <= pCutOff and abs(result[0]) >= rCutOff: 
            row=(feature, result[0],result[1])
            corrDataContent[i] = row
            i+=1 
    return corrDataContent
