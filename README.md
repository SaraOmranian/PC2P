# PC2P
Protein-Protein interaction Network - Clustering Algorithm 

Prediction of protein complexes from protein-protein interaction (PPI) networks rely on algorithms for network community detection that seek to identify dense subgraphs in protein interaction networks. However, recently assembled gold standards in yeast and human indicate that protein complexes also form sparse subgraphs. The contribution of our study is five-fold: (1) we formalize the concept of a protein complex as a biclique spanned subgraph, that addresses the density issue of existing approaches; (2) we propose a parameter-free approximation algorithm, termed Protein Complexes from Coherent Partition (PC2P), that solves the network partitioning into biclique spanned subgraphs by removing the smallest number of edges in a given PPI network; (3) we demonstrate that the resulting clusterings are of high modularity, thus reflecting the local structures in PPI networks, and (4) we show that PC2P outperforms seven seminal approaches with respect to a composite score combining five performance measures.

[1] Collins, S. R. et al., 2007. Toward a Comprehensive Atlas of the Physical Interactome of Saccharomyces cerevisiae. Molecular & Cellular Proteomics, Volume 6, p. 439–450.
[2] Gavin, A. et al., 2006. Proteome survey reveals modularity of the yeast cell machinery.. Nature.
[3] Krogan, N. J. et al., 2006. Global landscape of protein complexes in the yeast Saccharomyces cerevisiae. Nature, Volume 440, pp. 637-643.
[4] McDowall, M. D., Scott, M. S. & Barton, G. J., 2008. PIPs: human protein–protein interaction prediction database. Nucleic Acids Research, 11, Volume 37, pp. D651-D656.
[5] Nepusz, T., Yu, H. & Paccanaro, A., 2012. Detecting overlapping protein complexes in protein-protein interaction networks. Nature Methods, Volume 9, pp. 471-472.
[6] Szklarczyk, D. et al., 2014. STRING v10: protein–protein interaction networks, integrated over the tree of life. Nucleic Acids Research, 10, Volume 43, pp. D447-D452.
