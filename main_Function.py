# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 20:27:04 2020

@author: somranian
"""

import os
import PC2P_ParallelMultiprocess
import PC2P_ParallelRay
import PC2P_Sequential
import networkx as nx
import numpy as np
import pandas as pd


p = os.getcwd()
"""As an example here the PIPS_Corum_Graph is called!"""
path = p + "\\Human\\PIPS\\PIPS_Corum_Graph.txt"
PIPS_Corum = nx.read_weighted_edgelist(path, create_using = nx.Graph(), nodetype = str)
G = PIPS_Corum.copy()

""" To run code sequentially, we need to call Find_CNP from PC2P_Sequential.
    To run code parallel in Windows and Unix, we nee to call Find_CNP from PC2P_ParallelMultiprocess
    To run code parallel in Linux and Mac, we nee to call Find_CNP from PC2P_ParallelRay """

edge_cut = PC2P_Sequential.Find_CNP(G)
#PC2P_ParallelMultiprocess.Find_CNP(PIPS_Corum)
#PC2P_ParallelRay.Find_CNP(PIPS_Corum)

""" To save the result clusters in Graph format"""
G_copy = G.copy()
G_copy.remove_edges_from(edge_cut)

nx.write_edgelist(G_copy, "STRING_CNPPredicted_V4.edgelist.gz", data=False)
#--- If the edges are weighted ----------
nx.write_weighted_edgelist(G_copy, 'STRING_CNPPredicted_V4.weighted.edgelist')


""" We save each predicted cluster in one line """
G_cnp_components = list(nx.connected_components(G_copy))
G_cnp_components.sort(key=len, reverse=True)

with open('G_PredictedClusters.txt', 'w') as f:
    for item in G_cnp_components:
        for node in item:
            f.write("%s " % node)
        f.write("\n")

