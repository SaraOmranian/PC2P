# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 07:51:58 2020

@author: SaRa Omranian
"""

import networkx as nx
import itertools as itert
from operator import itemgetter
from networkx.algorithms.flow import shortest_augmenting_path
import ray
import psutil
num_cpus = psutil.cpu_count(logical=False)
ray.init(num_cpus=num_cpus)


def edgeCutSet_V2(cnp , G):
    """ This function find edges that connect the component to the rest of graph."""
    edge_cnp = cnp.edges()
    edge_G = G.edges(cnp.nodes())
    edgeSet = list(edge_G - edge_cnp)
    return(edgeSet)
    
def _is_separating_set(G, cut):
    """Assumes that the input graph is connected and return TRUE if removing cut from
        the input graph makes it disconnected."""
    if len(cut) == len(G) - 1:
        return True
    H = G.copy()
    H.remove_nodes_from(cut)
    if nx.is_connected(H):
        return False
    return True

def my_cut_nodes_V4(G,mixed_label = False):
    """ This function finds nodes with minimum degree and removes their neighbors of this node from G """
    cut_node = []
    degree = G.degree()
    sorted_nodes = sorted(degree, key=itemgetter(1), reverse=True)
    #-----if we have articulation point----------------------
    articu_points = list(nx.articulation_points(G))
    if articu_points:
        [cut_node.append([v]) for v in articu_points]
        return(cut_node)
    #----------------------------------------------------------
    #-------Find all nodes with delta degree and remove their neighbors
    delta = sorted_nodes[-1][1]
    X = [n for n, d in sorted_nodes[:delta]]
    if _is_separating_set(G, X):
        cut_node.append(X)
    node_deltaDegree = [j[0] for i,j in enumerate(degree) if j[1]==delta]
    for n in node_deltaDegree:
        cut_node.append(list(G.neighbors(n)))
    if not mixed_label:
        cut_node_set = set(tuple(sorted(x)) for x in cut_node)
        cut_node = [ list(x) for x in cut_node_set ]
    else:
        sorted_x = []
        for i in range(len(cut_node)):
            intList=sorted([i for i in cut_node[i] if type(i) is int])
            strList=sorted([i for i in cut_node[i] if type(i) is str])
            sorted_x.append(intList + strList) 
        cut_node =  list(set(tuple(i) for i in (sorted_x)))
    return(cut_node)

""" This part compute the score for each component """

def coherentCutRatio_V1(deg_cnp, deg_G):
    deg_cnp = [di[1] for di in deg_cnp] #change the set to list
    deg_G = [dg[1] for dg in deg_G] #change the degree set to list
    if sum(deg_cnp)==0:
         score = 100
    else:
        score = sum(deg_G)/sum(deg_cnp)
    return(score)  

#compute ratio of coherent cut which needs to be minimum
def coherentCutRatio_V4(deg_cnp, deg_G,cluster_coef):
    deg_cnp = [di[1] for di in deg_cnp] #change the set to list
    deg_G = [dg[1] for dg in deg_G] #change the degree set to list
    if sum(deg_cnp)==0:
        score = 100
    else:
        r = (sum(deg_G)/sum(deg_cnp)) - 1
        if cluster_coef!=0:
            score = r/cluster_coef
    return(score)  

@ray.remote
def second_Neighb(G,v):
    neighbor1 = list(G.neighbors(v))
    neighbor2 = neighbor1[:]
    [neighbor2.extend(list(G.neighbors(n))) for n in neighbor1]
    neighbor2 = list(set(neighbor2))
    return(neighbor2)

@ray.remote
def CNP(G,v ,mixed_label = False):  
    cnp_v = G.copy()
    neighbor1 = list(cnp_v.neighbors(v))
    neighbor1.extend([v])
    induced_N1 = cnp_v.subgraph(neighbor1)
    #----------- Compute Clustering Coefficient
    length_n1 = len(neighbor1)
    length_e1 = len(list(induced_N1.edges()))
    if length_n1==1 or length_n1==0:
        cc=0
    else:
        cc = 2*length_e1/(length_n1*(length_n1-1))
    #-----Calculating clustering coefficient is finished here!-------
    cutRatioN1_V4 = coherentCutRatio_V4(nx.degree(induced_N1),nx.degree(cnp_v,induced_N1.nodes()),cc)
    cutRatioN1_V1 = coherentCutRatio_V1(nx.degree(induced_N1),nx.degree(cnp_v,induced_N1.nodes()))
    cutRatioN1 = (cutRatioN1_V4+cutRatioN1_V1)/2
    neighbor2 = neighbor1[:]
    [neighbor2.extend(list(cnp_v.neighbors(n))) for n in neighbor1]
    neighbor2 = list(set(neighbor2))
    induced_N2 = cnp_v.subgraph(neighbor2).copy()
    induced_N2 = induced_N2.copy()
    Complement_indN2 = nx.complement(induced_N2)
    if(not nx.is_connected(Complement_indN2)):#we find the CNP
        #----------- Compute Clustering Coefficient
        length_n2 = len(neighbor2)
        length_e2 = len(list(induced_N2.edges()))
        if length_n2==1 or length_n2==0:
            cc=0
        else:
            cc = 2*length_e2/(length_n2*(length_n2-1))
        #-----Calculating clustering coefficient is finished here!-------
        cutRatioN2_V4 = coherentCutRatio_V4(nx.degree(induced_N2),nx.degree(cnp_v,induced_N2.nodes()),cc)
        cutRatioN2_V1 = coherentCutRatio_V1(nx.degree(induced_N2),nx.degree(cnp_v,induced_N2.nodes()))
        cutRatioN2 = (cutRatioN2_V4+cutRatioN2_V1)/2
    else:
        #find nodes to disconnect the complement of N2, change elemnt of node cuts to list
#                 cut_node = [list(n) for n in list(nx.all_node_cuts(Complement_indN2,flow_func=shortest_augmenting_path))]
        cut_node = my_cut_nodes_V4(Complement_indN2, mixed_label)
        if (len(cut_node)==0):
            return
        elif len(cut_node)>1: #we may find differnt cut_node_set 
            #if we have more than one choose a set with minimum weight
            minweight = 1000
            for n in cut_node:
                #we calculate the score for each set seperately and then choose one with the minimum score
                #at the end we choose the set with minimum weight
                temp_N2 = induced_N2.copy()
                temp_N2.remove_nodes_from(n)
                #----------- Compute Clustering Coefficient
                length_n2 = len(temp_N2.nodes())
                length_e2 = len(list(temp_N2.edges()))
                if length_n2==1 or length_n2==0:
                    cc=0
                else:
                    cc = 2*length_e2/(length_n2*(length_n2-1))
                #-----Calculating clustering coefficient is finished here!-------
                cutRatioN2_V4 = coherentCutRatio_V4(nx.degree(temp_N2),nx.degree(cnp_v,temp_N2.nodes()),cc)
                cutRatioN2_V1 = coherentCutRatio_V1(nx.degree(temp_N2),nx.degree(cnp_v,temp_N2.nodes()))
                temp_score = (cutRatioN2_V4+cutRatioN2_V1)/2
                if temp_score < minweight:
                    minweight = temp_score
                    minWeightNode = n
            cut_node = minWeightNode
        else:
            #flatten the cut node
            cut_node = [node for sublist in cut_node for node in sublist]
        induced_N2.remove_nodes_from(cut_node)
        #----------- Compute Clustering Coefficient
        length_n2 = len(induced_N2.nodes())
        length_e2 = len(list(induced_N2.edges()))
        if length_n2==1 or length_n2==0:
            cc=0
        else:
            cc = 2*length_e2/(length_n2*(length_n2-1))
        #-----Calculating clustering coefficient is finished here!-------
        cutRatioN2_V4 = coherentCutRatio_V4(nx.degree(induced_N2),nx.degree(cnp_v,induced_N2.nodes()),cc)
        cutRatioN2_V1 = coherentCutRatio_V1(nx.degree(induced_N2),nx.degree(cnp_v,induced_N2.nodes()))
        cutRatioN2 = (cutRatioN2_V4+cutRatioN2_V1)/2
    if cutRatioN1 < cutRatioN2:
        min_ratio = cutRatioN1
        cnp_nodes = induced_N1
    else:
        min_ratio = cutRatioN2
        cnp_nodes = induced_N2
    result = {cnp_nodes:[v,min_ratio]}
    return(result)

def Find_CNPs_V2(G, mixed_label = False):
    #Find all component of G
    G_components = list(nx.connected_components(G))
    G_temp = G.copy()
    edge_cut = []
    rounds = 1
    while len(G_components) != 0:
        componentOfG = G.subgraph(G_components[0]).copy() #we get the first component and find cnp
        if rounds==1:
            if len(componentOfG.nodes()) <= 3:
                del G_components[0]
                continue  
            nodes = list(componentOfG.nodes())
            result_objects_1 = ray.get([CNP.remote(componentOfG,v,mixed_label) for v in nodes])
            results_1 = result_objects_1
            scores_1 = [list(r.values())[0][1] for r in result_objects_1]
            indx_1 = scores_1.index(min(scores_1))
            subgrf_1 = list(results_1[indx_1].keys())[0]
            edge_cut.append(edgeCutSet_V2(subgrf_1,componentOfG)) 
            #finding second neighbors for all cnp nodes
            result_objects2 = ray.get([second_Neighb.remote(componentOfG,v) for v in subgrf_1.nodes()])
            secondNeighb = []
            [secondNeighb.extend(s) for s in result_objects2]
            secondNeighb = set(secondNeighb)
            Nodesto_NextRound = secondNeighb - set(subgrf_1.nodes())
            updated_results = [results_1[i] for i,r in enumerate(results_1) if not(list(r.values())[0][0] in secondNeighb)]
            G_temp.remove_nodes_from(subgrf_1.nodes())
            G_components = list(nx.connected_components(G_temp))
            print("end of round: ", rounds)
            del nodes
            rounds += 1
        else:
            if len(componentOfG.nodes()) <= 3:
                for i,r in enumerate(updated_results):
                    if ( list(r.values())[0][0] in list(componentOfG.nodes()) ):
                        del updated_results[i]
                for i,n in enumerate(list(Nodesto_NextRound)):
                    if ( n in list(componentOfG.nodes()) ):
                        Nodesto_NextRound.remove(n)
                del G_components[0]
                continue   
    #I get the intersect incase we have multiple components. 
#    At the end I have to add the nodes which are not present in the component to nodesto_NextRound
            nodes = Nodesto_NextRound.intersection(set(componentOfG.nodes()))
            if not nodes:
                nodes = set(componentOfG.nodes())
                nodes_diff = Nodesto_NextRound 
            else:
                nodes_diff = Nodesto_NextRound - set(componentOfG.nodes())
            result_objects = ray.get([CNP.remote(componentOfG,v,mixed_label) for v in nodes])
            results = result_objects
            results.extend(updated_results)
            scores = [list(r.values())[0][1] for r in results]
#            print('scores: ',scores)
            indx = scores.index(min(scores))
            subgrf = list(results[indx].keys())[0]
            edge_cut.append(edgeCutSet_V2(subgrf,G_temp))
            #finding second neighbors for all cnp nodes
            result_objects2 = ray.get([second_Neighb.remote(G_temp,v) for v in subgrf.nodes()])
            results2 = result_objects2
            secondNeighb = []
            [secondNeighb.extend(s) for s in results2]
            secondNeighb = set(secondNeighb)
            Nodesto_NextRound = secondNeighb - set(subgrf.nodes())
            if nodes_diff:
                Nodesto_NextRound = Nodesto_NextRound|nodes_diff
            updated_results = [results[i] for i,r in enumerate(results) if not (list(r.values())[0][0] in secondNeighb)]
            G_temp.remove_nodes_from(subgrf.nodes())
            G_components = list(nx.connected_components(G_temp))
            print("end of round: ", rounds)
            rounds += 1
    edge_cut = [edge for sublist in edge_cut for edge in sublist]
    if not mixed_label:
        edge_cut = list(set(tuple(sorted(x)) for x in edge_cut))
    else:
        sorted_x = []
        for i in range(len(edge_cut)):
            intList=sorted([i for i in edge_cut[i] if type(i) is int])
            strList=sorted([i for i in edge_cut[i] if type(i) is str])
            sorted_x.append(intList + strList) 
        edge_cut =  list(set(tuple(i) for i in (sorted_x)))
    return(edge_cut)
