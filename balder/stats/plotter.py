class Plotter(object):
    
    def __init__(self):
        pass

    def community_response(self, diversity_data, metadata):

        title= self.r.initRVariable('title', self.analysis_type)
        metadata= metadata.ix[diversity_data.index]
            

        
        r_diversity_data= pandas_df_to_r_df(diversity_data)
        
        r_dis_matrix= self.r.fPlotCommunityResponse(r_diversity_data)
        dis_matrix= r_matrix_to_dataframe(r_dis_matrix)  

        ### build a new dataframe with triclosan values
        A={int(b):a for a,b in zip(diversity_data.T ,set(list(dis_matrix['Var1'])))}
        new_df= DataFrame(columns= ["Var1","Var2","value","Triclosan1","Triclosan2"], index= dis_matrix.index)
        
        for index in dis_matrix.index:
            try:
                new_df['Var1'][index]= A[dis_matrix['Var1'][index]] 
                new_df['Var2'][index]= A[dis_matrix['Var2'][index]] 
                new_df['Triclosan1'][index]= self.r.df_metadata.ix[new_df['Var1'][index]]['Triclosan'] 
                new_df['Triclosan2'][index]= self.r.df_metadata.ix[new_df['Var2'][index]]['Triclosan'] 
                new_df['value'].loc[index]= dis_matrix['value'].loc[index]
            except Exception,e:
                print e
                continue
        
        df_comparisons= new_df[new_df['value'] != 0 ]


        #!!! For Alex!!! I think you missed the following condition in the filters "and Triclosan2 > 0.1"
        df_controls= df_comparisons.query("Triclosan1 < 0.1  and Triclosan2 < 0.1")
        df_treatments= df_comparisons.query("Triclosan2 < 0.1  and Triclosan1 > 0.1")


        cline= df_controls['value'].mean()
        cline_std= df_controls['value'].std()

        df_treatments['value']= df_treatments['value'].astype(float)
        

        averages= df_treatments.groupby('Var1').aggregate(np.mean)
        metadata= self.r.df_metadata.ix[averages.index]['Triclosan']



        # adding a negligible value to the controls so that logs can be meaningful
        metadata[metadata == 0.0] = metadata[metadata == 0.0] +0.1
        ### plot the distance matrix
        
        plt.figure()
        
        plt.scatter(np.log10(metadata), 1- averages, s=80, c="k")
        
        ax= plt.gca()
        x_min, x_max= ax.get_xlim()
        ax.set_xlim(x_min, x_max)
        
        
        #plt.plot(, 'k-')

        a=np.log10(metadata)

        a=list(a)
        #print a
        a.insert(0, x_min)
        a.insert(-1, x_max)
        a=np.array(a)

        x= a
        y= a*0+(1-cline)
        y1= a*0+(1-(cline+cline_std))
        y2= a*0+(1-(cline-cline_std))
        
        plt.plot(x, y,  'k-', alpha= 0.6)
        plt.plot(x, y1 , 'k--', alpha=0.4)
        plt.plot(x, y2, 'k--', alpha=0.4)
        
        ax.set_xlabel('Triclosan concentration [nM]')
        ax.set_ylabel('1 - Bray-Curtis distance')
        ax.set_title('Community similarity in comparison to controls')

        labels= ['%.3g' %s for s in [10**x for x in np.arange(-1,3, 0.5)]]
        labels[0]=""
         
        ax.set_xticklabels(labels)
        ax.grid(True, color="gray", alpha= 0.5)

        plt.savefig(os.path.join(self.analysis_dir, 'community_response.png'))
        plt.close()

    
    def diversity_index(self, type_short, type_long, metadata_type, df_indice, index_data):
 ):
               
        #metadata.replace(0.0, 0.1)
        metadata= df_indice['metadata']
   

        # adding a negligible value to the controls so that logs can be meaningful
        metadata[metadata == 0.0] = metadata[metadata == 0.0] +0.1
        
        env_transformed= np.log10(metadata) 
        m, b = np.polyfit(env_transformed, index_data, 1)
        plt.figure()    
        plt.plot(env_transformed, index_data, 'ko')
        plt.xlim([-1.5, 3])

        ax= plt.gca()

        #print metadata_type, type_long, type_short, title
        ax.set_xlabel(metadata_type)
        ax.set_ylabel(type_long)
        ax.set_title( "%s - %s" % (self.analysis_type, type_long))
        labels= ['%.3g' %s for s in [10**x for x in np.arange(-1.5,3,0.5)]]
        labels[0]=""
        labels[1]="0"
        ax.set_xticklabels(labels)

        ax.grid(True, color="gray", alpha= 0.5)
        plt.savefig(os.path.join(self.analysis_dir, '%s.png' %type_short))

        plt.plot(env_transformed, m*env_transformed + b, 'k-')
        
        #plt.show()
        
        plt.savefig(os.path.join(self.analysis_dir, '%s_regression.png' %type_short))
        
        
        annotations= ["%s - %s" % (ind, df_indice.ix[ind]['metadata']) for ind in df_indice.index]
        #annotations= [ "%s - %s" % (sample, con) for (sample, con) in zip(self.r.df_metadata.index, map(str, list(metadata)))]
        
        
        for label, x, y in zip(annotations, env_transformed, index_data):
            plt.annotate(
                label, 
                xy = (x, y), xytext = (random.randint(-50,-1), random.randint(1,50)),
                textcoords = 'offset points', ha = 'right', va = 'bottom',
                bbox = dict(boxstyle = 'round,pad=0.8', fc = 'yellow', alpha = 0.2),
                arrowprops = dict(facecolor= "white", width= 2, alpha= 0.7), fontsize=10)
        #"""
        plt.savefig( os.path.join(self.analysis_dir, '%s_with_sample_names.png'% type_short))
        plt.close()

   

        
    def abundance(self, otu_data, analysis_type):
 
        #if self.ssu_type == "16S":
        ncol=1
        width_factor= 0.85
        #elif self.ssu_type == "18S":
        #    ncol=2
        #    width_factor= 0.65
        
        
        if analysis_type == self.level:
            f= plt.figure(figsize= (15,10))
            ax= data_frame.plot(kind= "bar", ax= f.gca(), stacked= True, ylim= (0,1), grid= False, colormap="Set1")
            ax.legend(fontsize= 11, ncol= ncol, loc= "best", bbox_to_anchor=(1.,1.)) 
            pos1= ax.get_position()
            pos2= [pos1.x0 - 0.05, pos1.y0 +0.075,  pos1.width*width_factor, pos1.height*0.9] 
            ax.set_position(pos2)
            ax.set_title("%s Phylum Relative Abundances" %self.ssu_type)
            ax.set_ylabel("Relative abundance")
            f.savefig( os.path.join(self.analysis_dir, "phylum_abundance.png"))
            plt.close()


    
    def rarefaction(self, otu_data):
        r_analysis_type= self.r.initRVariable('analysis_type', self.analysis_type)
        r_analysis_dir= self.r.initRVariable('analysis_dir', self.analysis_dir)
        r_nmds_data= pandas_df_to_r_df(nmds_data)
        rarefaction= self.r.fPlotRarefaction(r_nmds_data, r_analysis_type, r_analysis_dir)  


    def nmds(self, otu_data, metadata, analysis_type):
        
        r_nmds_data= pandas_df_to_r_df(nmds_data)
        r_analysis_type= self.r.initRVariable('analysis_type', analysis_type)
        r_nmds= self.r.fMetaMDS(r_nmds_data, r_analysis_type)
        r_scores= self.r.fMDSScores(r_nmds) 
       
        stress= list(r_nmds.rx('stress')[0])[0]
        scores= r_matrix_to_dataframe(r_scores)
        
        plt.figure()
        dim1= scores["NMDS1"]
        dim2= scores["NMDS2"]
    
        metadata= metadata.ix[nmds_data.index]
        
        c_data = plt.cm.jet(np.log10(metadata))
        cm = plt.cm.get_cmap('Greys')
       
        # adding a negligible value to the controls so that logs can be meaningful
        metadata[metadata == 0.0] = metadata[metadata == 0.0] +0.1


        sc = plt.scatter(dim1,dim2, c= np.log10(metadata), vmin=min(np.log10(metadata)), vmax=max(np.log10(metadata)),s=100, cmap=cm)
        
        cb= plt.colorbar(sc, shrink= 0.6, aspect= 10, pad= 0.1, norm=matplotlib.colors.LogNorm())

        cb.ax.tick_params(labelsize='x-small',top=False)
        cb.set_label('Triclosan concentration [nM]', labelpad=-70)
        
        labels= np.sort(np.array(list(set(metadata))))

        ind = np.log10(labels)
        ind_labels= map(str, labels)
        ind_labels[0] = "0.0" 
        cb.set_ticks(ind)
        cb.set_ticklabels(ind_labels)

        ax= plt.gca()
        ax.set_xlabel('NMDS1')
        ax.set_ylabel('NMDS2')
        ax.set_title('NMDS Ordination of %s OTU Abundance' %self.analysis_type)
        
        left, width = 1.05, .5
        bottom, height = 0.05, .5
        
        ax.grid(True, color="gray", alpha= 0.5)
        ax.text(left, bottom, 'Stress: %.2g'% stress, fontsize= 11, 
                horizontalalignment='left', verticalalignment='bottom', 
                transform= ax.transAxes)# bbox={ 'facecolor':'white','pad':10})
        
        
        plt.savefig( os.path.join(self.analysis_dir,'nmds.png'))
        annotations= map(str, list(metadata))
        for label, x, y in zip(annotations, dim1, dim2):
            plt.annotate(
                label, 
                xy = (x, y), xytext = (random.randint(-50,-1), random.randint(1,50)),
                textcoords = 'offset points', ha = 'right', va = 'bottom',
                bbox = dict(boxstyle = 'round,pad=0.8', fc = 'yellow', alpha = 0.2),
                arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'), 
                fontsize= 10)
        #"""
        plt.savefig( os.path.join(self.analysis_dir, 'nmds_with_sample_names.png'))
        
        plt.close()



    def venn(self, sets= 2):

        venn_file= os.path.join(self.dirTables, "venn.png")
        
        if self.analysis_type == "16S":
            venn3(map(int, treatment_group_otus.values()), treatment_group_otus.keys())
        else:
            venn2(map(int, treatment_group_otus.values()), treatment_group_otus.keys())
        
        plt.savefig(venn_file)   

