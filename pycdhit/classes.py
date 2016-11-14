# -*- coding: utf-8 -*-

#import argparse
import re
import pprint
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist,squareform
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import fcluster
from scipy.cluster.hierarchy import dendrogram
from scipy.cluster import hierarchy #do we need this?
import matplotlib #do we need this?
import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt
import logging
import string
from fatool import *


class cdhit_read(object):
    """Represents single CDHIT read with name pb and length"""
    
    def __init__(self, name = None, pb = None, length = None):
        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger(__name__)
        self.name = name
        self.pb = pb
        self.length = length
        logger.debug('creating cdhit_read:'+self.name+' '+str(self.pb)+' '+self.length)
    
   
    def to_json(self):
        """ creating json string form self values"""
        logger = logging.getLogger(__name__)
        logger.debug('read.to_json()'+ str({'len':self.length, 'name':self.name, 'pb':self.pb}))
        return {'len':self.length, 'name':self.name, 'pb':self.pb}
        #{sample:XYZ, data:[ {gene:XYZ, reads:[ 0:{f1:4151aa, name:TB, pb:*}, 1:{} ]} {gene:ABC, reads[]} ]}
        
    def __eq__(self, other):
        if self.name != other.name:
            return 0
        if self.pb != other.pb:
            return 0
        if self.length != other.length:
            return 0
        return 1
        
    def __cmp__(self, other):
        if self.name != other.name:
            return 0
        if self.pb != other.pb:
            return 0
        if self.length != other.length:
            return 0
        return 1
    
class cdhit_cluster(object):
    """Represents whole CDHIT cluster with many reads"""
    def __init__(self, name, reads=None, gn=None):
        self.name = name#Cluster
        """name obtains cluster value of cd-hit file"""
        self.gene_name = gn#Gene name
        """gene_name obtains name of sequence from cd-hit file"""
        if reads:
            self.reads = reads#reads list
            """reads list of reads contained by cluster in cd-hit file"""
        else:
            self.reads = []
    
    def to_json(self):
        """returns json representation of cd-hit cluster"""
        logger = logging.getLogger(__name__)
        logger.debug('starting cluster.to_json')
        logger.debug(self.reads)
        return {'gene':self.gene_name, 'cluster':self.name, 'reads':[r.to_json() for r in self.reads]}
        
    def get_single_value(self):
        """returns match score of first read if read exists, else 0"""
        if len(self.reads) > 1:
            return self.reads[1].pb
        else:
            return 0
            
    def get_label(self):
        """returns label of cluster in this case gene_name"""
        #print self.gene_name
        return self.gene_name
    
    def append(self, read):
        """appends read to reads list in cluster"""
        self.reads.append(read)

    
class cdhit_result(object):
    """represents whole cd-hit file with all clusters as list"""
    def __init__(self, name=None, data=None):
        """initialize cdhit result object"""
        self.name = name
        """ name result_set usualy equal to files name"""
        if data:
            self.data = data
        else:
            self.data = []
        
    def to_json(self):
        """returns json representation ot cd-hit result"""
        return {'sample':self.name, 'data':[r.to_json() for r in self.data]}
        
    def to_df(self):
        """transforms cdht result into pandas dataFrame"""
        #  creation of dict {df.label:df.value} for df creation purpouse
        d = {}
        for r in self.data:
            d[r.get_label()] = r.get_single_value()
        #index=name of row in this case eq file name which eq sampel name
        df = pd.DataFrame(data=d, index=[self.name])
        return df
        
        
    def get_thold_labels(self, th, dir):
        """
        return labels thats score mets usere defined condition defined in th and dir
        th is vlaue treshold
        dir is to specify whether labels with values:
        greater than threshold will be returned (dir >= 0) 
        or lower and equal (dir < 0)
        """
        labels = []
        if dir >= 0:
            for r in self.data:
                if r.get_single_value() > th:
                    labels.append(r.get_label())
        else:
            for r in self.data:
                if r.get_single_value() <= th:
                    labels.append(r.get_label())
        return labels
        
    def append(self, cluster):
        self.data.append(cluster)
    
    def load_from_file(self, file, name = None):
        if name:
            self.name = name
        if isinstance(file, str):
            with open(file) as f:
                self.load_content(f.read().splitlines())
        else:
            with file as f:
                self.load_content(f.read().splitlines())
        
    
    def load_content(self, cnt):
        logger = logging.getLogger(__name__)
        pat = re.compile('>(.+)\.\.\.')
        reads = []
        # analizing each line of file
        for l in cnt:
            # if line begins with > new cluster
            if l[0] == '>':
                # if reads from prev cluster add to cluster and append it, get name for new cluster.
                if reads:
                    self.append(cdhit_cluster(name, reads, reads[0].name))
                    name = l[1:].strip()
                    reads = []
                #no reads it is first cluster
                else:
                    logger.debug('appling first name');
                    name = l[1:].strip()
            # reads (does not starts with >)
            else:
                # looking for gene name
                m = pat.search(l)
                if m:
                    gname = m.group(1)
                else:
                    gname = 'unk'
                
                reads.append(cdhit_read(gname, l.split('... ')[1].strip(' at%'), l[2:].split(',')[0].strip()))
                #{len:l[2:].split(',')[0].strip(), 'name':name, 'pb':l.split('... ')[1].strip(' at%')}
        #self.
            
   
class cdhit_set(object):
    def __init__(self, rlist=None):
        #print rlist
        if rlist:
            self.result_list = rlist
        else:
            self.result_list = []
        self.samples = []
        self.labels = []
        self.all_zeros = [] #labels where there are all 0 in columns
        self.all_non_zeros = [] #labels where there are all non 0 in columns
        self.zeros = {}
        self.nzeros = {}
        #clusters
        
    #def get_non_zeros(self):
    
    def append(self, cdr):
        print 'appending result:'
        print cdr
        self.result_list.append(cdr)
        
    
    def analyze(self):
        print self.result_list[0]
        #labels = ''
        labels = set(self.result_list[0].get_thold_labels(0, -1))
        #print labels
        for r in self.result_list[1:]:
            labels = labels & set(r.get_thold_labels(0, -1))
        #print 'labels as list'
        #print list(labels)
        self.all_zeros = list(labels)
        
        labels = None
        labels = set(self.result_list[0].get_thold_labels(0, 1))
        for r in self.result_list[1:]:
            labels = labels & set(r.get_thold_labels(0, 1))
        self.all_non_zeros = list(labels)


    def to_df(self, skip_zeros = None, skip_non_zeros = None):
        #print len(self.result_list)
        for s in self.result_list:
            #self.samples.append(s.name)
            try:
                df  = df.append(s.to_df())
            except NameError:
                df = s.to_df()
                continue
        if skip_zeros:
            #print 'diff:'
            #print len(set(self.all_zeros) - set(df.columns))
            #print len(self.all_zeros)
            #for r in self.all_zeros:
            #    if r not in df.columns:
            #        print r
            df = df.drop(self.all_zeros, 1)
        if skip_non_zeros:
            #print 'non 0'
            #for r in self.all_non_zeros:
            #    if r not in df.columns:
            #        print r
            #print '0 & n0'
            #print set(self.all_zeros)&set(self.all_non_zeros)
            df = df.drop(self.all_non_zeros, 1)
        return df
        
    def make_dendrogram(self, df, fig_size, font_top_size=1):
        '''
        preparing and priniting to file dendrogram basing on cdhit_set data
        df data frame with data from cdht_set
        fig_size figure size
        '''
        df = df[df.columns].astype(float)
        
        # columns pairwise distance
        col_dists = pdist(df.T, metric='euclidean')
        col_clusters = linkage(col_dists, method='complete')
        ''' printing clusters'''
        print '\n------------------\n'
        #for cc in col_clusters:
        #    print cc
            
        #    print '--- eo cc ---'
        #print sch.leaves_list(col_clusters)
        print '\n------------------\n'
        fc = fcluster(col_clusters, 0.7*col_dists.max(), 'distance')
        #print string.join(['ccf','']+map(str, fc),'\t')
        print fc
        #col_dists.max()
            #print r
            #print '\n---\n'
        
        
        ''' ---- end of printing clusters ----'''
        fig = plt.figure(figsize=fig_size)
        
        axd2 = fig.add_axes([0.22,0.8,0.72,0.2]) # x, y, width, height
        col_dendr = dendrogram(col_clusters, distance_sort='ascending', labels=df.columns, leaf_font_size=str(font_top_size))
        '''test'''
        
        #print col_dendr['ivl']
        #print '\n-----list----\n'
        #print sch.leaves_list(col_clusters)
        #print '\n---------\n'
        #print sch.fclusterdata(col_clusters, 0.7*col_dists.max())
        #print col_dendr['leaves']
        print '\n---------\n'
        #print sch.leaves_list(col_clusters)
        #print pd.DataFrame({'leaves':df.columns, 'clusters':fc})
        #print pd.DataFrame({'leaves':df.columns, 'clusters':fc}).sort_index(by=['clusters'])#.to_csv('cluster_info.csv', ';')
        self.clusters = pd.DataFrame({'leaves':df.columns, 'clusters':fc}).sort_index(by=['clusters'])#.to_csv('cluster_info.csv', ';')
        
        
        
        row_clusters = linkage(pdist(df, metric='euclidean'), method='complete')
        axd1 = fig.add_axes([0.01,0.001,0.18,0.8])# x-pos, y-pos, width, height+
        row_dendr = dendrogram(row_clusters, orientation='left',  
                           count_sort='ascending')
                           #labels=samples)
                           #color_threshold=np.inf) # makes dendrogram black
        # no labels on plot axes
        axd1.set_xticks([])
        axd1.set_yticks([])
        axd2.set_xticks([])
        axd2.set_yticks([])
        
        # remove axes spines from dendrogram
        for i,j in zip(axd1.spines.values(), axd2.spines.values()):
            i.set_visible(False)
            j.set_visible(False)
        
        # sort acording to samples dendrogram
        df_rowclust = df.ix[row_dendr['leaves'][::-1]]
        # sort acording to columns dendrogram
        df_rowclust = df_rowclust[df_rowclust.columns[col_dendr['leaves']]]
        
        axm = fig.add_axes([0.22,0.001,0.9,0.8])# x-pos, y-pos, width, height+
        cax = axm.matshow(df_rowclust, interpolation='nearest', cmap=plt.cm.YlGnBu, aspect='auto')
        fig.colorbar(cax)
        
        plt.xticks(range(len(list(df_rowclust.columns))), list(df_rowclust.columns), rotation=90)#, rotation=90#range(len(list(df_rowclust.columns)))
        plt.yticks(range(len(list(df_rowclust.index))), list(df_rowclust.index))
        plt.setp(axm.get_xticklabels()[::1], fontsize=int(font_top_size))
        
        return plt
            
    
    def to_json(self):
        pass
        
        
    
    
        
    

    
