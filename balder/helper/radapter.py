from collections import namedtuple

import rpy2.robjects as r
import numpy as np
from numpy.random import RandomState

from pandas import DataFrame, Series, read_csv
from pandas.rpy.common import convert_to_r_dataframe


from r2py_variable import r_imports, testDeseq2, correlateDeseq2



################################################################################
def r_matrix_to_dataframe(matrix):
    cols = list(matrix.colnames)
    rows = list(matrix.rownames)
    try:
        return DataFrame(np.array(matrix), index=rows, columns=cols)
    except:
        return DataFrame(np.array(matrix).transpose(), index=rows, columns=cols)

################################################################################
def pandas_df_to_r_df(pandas_df):
    return convert_to_r_dataframe(pandas_df)

################################################################################
def pandas_df_to_named_r_df(pandas_df, r_name):
    eval("%s = convert_to_r_dataframe(pandas_df)" % r_name)
    eval("return %s" % r_name)


class RAdapter(object):
    """
        handles the communication with R
    """
    
    def __init__(self):
        self.init_vars()
        self.init_r_imports()
        #self.initRVariables()
        self.init_r_functions()
        #self.init_div_indexes()
    
    def init_vars(self):
        self.div_indexes= {}


    def init_r_var(self):
        r.globalenv[var]= value
        return r.r[var]
    
    def init_r_imports(self):
        self.imports= r.r(r_imports)


    def init_r_functions(self):
        #self.fPrepareData= r.r(prepareData)
        #self.fMetaMDS= r.r(runMetaMDS)
        #self.fMDSScores= r.r(mdsScores)
        #self.fDrawLogaxis= r.r(drawLogAxis)
        #self.fBuildDiversityIndices= r.r(buildDiversityIndices)
        #self.fPlotDiversityIndex= r.r(plotDiversityIndex)
        #self.fCalCommunityResponse= r.r(calcCommunityResponse)
        #self.fPlotRarefaction= r.r(plotRarefaction)
        #self.fTestEdgeR= r.r(testEdgeR)
        self.fTestDeseq2= r.r(testDeseq2)
        #self.fNormalizeDeseq2= r.r(normalizeDeseq2)
        self.fCorrelateDeseq2= r.r(correlateDeseq2)

    

    def transform_dataframe_by_r_function(self, df):
        r_df= pandas_df_to_r_df(df)  
        if library != "":
            import_=r.r("library(%s)"%library)
        function_=r.r(function)
        r_updated_df= function_(r_df, *args)
        return r_matrix_to_dataframe(r_updated_df)

 
    def normalize_dataframe_by_alogorithm(self, method= "deseq2"):
        r_df= pandas_df_to_r_df(df) 
        groups= [1]*len(df.columns)
        r_groups= r.IntVector(groups)
        r_updated_df= self.fNormalizeDeseq2(r_df, r_groups)
        return r_matrix_to_dataframe(r_updated_df)
    
    
    def apply_r_function_to_dataframe(self,df, function, library="", *args):

        r_df= pandas_df_to_r_df(df)  
        if library != "":
            import_=r.r("library(%s)"%library)
        function_=r.r(function)
        return function_(r_df, *args)


