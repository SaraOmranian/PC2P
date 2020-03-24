# -*- coding: utf-8 -*-
"""
Created on Sat Mar  7 17:13:57 2020

@author: somranian
"""

""" Nepusz, T., Yu, H. & Paccanaro, A., 2012. Detecting overlapping protein complexes in protein-protein
    interaction networks. Nature Methods, Volume 9, pp. 471-472."""
from mwmatching import maxWeightMatching


def canonical_protein_name(name):
    """Returns the canonical name of a protein by performing a few simple
    transformations on the name."""
    return name.strip().upper()

def is_numeric(x):
    """Returns whether the given string can be interpreted as a number."""
    try:
        float(x)
        return True
    except:
        return False

def matching_score(set1, set2):
    """Calculates the matching score between two sets (e.g., a cluster and a complex)
    using the approach of Bader et al, 2001"""
    return len(set1.intersection(set2))**2 / (float(len(set1)) * len(set2))

def maximum_matching_ratio(reference, predicted, score_threshold=0.2):
    scores = {}
    count = 0
    n = len(reference)
    for id1, c1 in enumerate(reference):
        for id2, c2 in enumerate(predicted):
            score = matching_score(c1, c2)
            if score == 1:
                count = count + 1
            if score <= score_threshold:
                continue

            scores[id1, id2+n] = score
    print('score is:' ,scores)
    print('fully matched V1',count,'out of',n)
    inpt = [(v1, v2, w) for (v1, v2), w in scores.items()]
    mates = maxWeightMatching(inpt)
    score = sum(scores[i, mate] for i, mate in enumerate(mates) if i < mate)
    return score / n

def clusteringwise_sensitivity(reference, predicted):
    num, den = 0., 0.
    for complx in reference:
        den += len(complx)
        num += max(len(complx.intersection(cluster)) for cluster in predicted)
    if den == 0.:
        return 0.
    return num / den

def positive_predictive_value(reference, predicted):
    num, den = 0., 0.
    for cluster in predicted:
        isects = [len(cluster.intersection(compl)) for compl in reference]
        isects.append(0.)
        num += max(isects)
        den += sum(isects)
    if den == 0.:
        return 0.
    return num / den

def accuracy(reference, predicted):
    return (clusteringwise_sensitivity(reference, predicted) * \
            positive_predictive_value(reference, predicted)) ** 0.5


def fraction_matched(reference, predicted, score_threshold=0.25):
    result = 0

    for id1, c1 in enumerate(reference):
        for id2, c2 in enumerate(predicted):
            score = matching_score(c1, c2)
            if score > score_threshold:
                result += 1
                break

    return result / len(reference)

def clusteringwise_separation(reference, predicted):
    intersections = {}
    marginal_sums = [0.] * len(predicted), [0.] * len(reference)
    for i, cluster in enumerate(predicted):
        for j, compl in enumerate(reference):
            isect = len(cluster.intersection(compl))
            if isect > 0:
                intersections[i, j] = isect
            marginal_sums[0][i] += isect
            marginal_sums[1][j] += isect

    separations_complex = [0.] * len(reference)
    separations_cluster = [0.] * len(predicted)
    for i, cluster in enumerate(predicted):
        s = marginal_sums[0][i]
        for j, compl in enumerate(reference):
            isect = intersections.get((i, j), 0)
            if isect == 0:
                continue
            val = float(isect * isect) / (s * marginal_sums[1][j])
            separations_complex[j] += val
            separations_cluster[i] += val

    avg_sep_complex = sum(separations_complex) / len(separations_complex)
    avg_sep_cluster = sum(separations_cluster) / len(separations_cluster)
    return (avg_sep_complex * avg_sep_cluster) ** 0.5

def read_network(fname,l):
    known_proteins = list()
    for line in open(fname):
        parts = [canonical_protein_name(part) for part in line.strip().split() if not is_numeric(part)]
        if(len(parts)>=l):
            known_proteins.append(set(parts))
    return known_proteins
###########################################################################

def jaccard(set1,set2):
    return ((len(set1.intersection(set2)))/(float(len(set1))+len(set2)-len(set1.intersection(set2))))

def precision_Jaccard(reference,predicted, threshold=0.5):
    counter = 0
    for pred in predicted:
        for ref in reference:
            score = jaccard(ref,pred)
            if score >= threshold:
                counter += 1
                break
    return counter / float(len(predicted))

def recall_Jaccard(reference,predicted, threshold=0.5):
    counter = 0
    for ref in reference:
        for pred in predicted:
            score = jaccard(ref,pred)
            if score >= threshold:
                counter += 1
                break
    return counter / float(len(reference))

def F_measure_Jaccard(reference,predicted, threshold = 0.5):
    p = precision_Jaccard(reference,predicted,threshold)
    r = recall_Jaccard(reference,predicted,threshold)
    return( (2*p*r)/(p + r))      
    
