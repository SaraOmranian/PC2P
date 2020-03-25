# PC2P 

### Protein Complexes from Coherent Partition


PC2P is a network clustering algorithm that can be used to find protein complexes from protein-protein interaction networks(PPINs).

PC2P is written in 2 formats: sequential and parallel.

The parallel version is implemented with 2 packages in python: Ray and Multiprocess. The Ray implementation can be run on Mac and Linux operating systems, however, for Windows machine, you can use Multiprocess implementation.

To run the program you need to call one of these implementations and pass your network to the function, please see main_Function.py for more details.

There is an implementation with Jupiter, in which I put some examples. In these examples, you can find the original network and the predicted clusters after applying PC2P along with the total number of removed edges to obtain these clusters. 

You can find PPINs for Human and Yeast in their respective folders with their gold standards, which we used to compare the predicted clusters with reference complexes.

If you want to do further analysis and check how well the predicted clusters represent the reference complexes, please run PredictedClusters_Analysis.py in the Analysis folder for the predicted clusters from PC2P.

## References

[1] Collins, S. R. et al., 2007. Toward a Comprehensive Atlas of the Physical Interactome of Saccharomyces cerevisiae. Molecular & Cellular Proteomics, Volume 6, p. 439–450.

[2] Gavin, A. et al., 2006. Proteome survey reveals modularity of the yeast cell machinery. Nature.

[3] Krogan, N. J. et al., 2006. Global landscape of protein complexes in the yeast Saccharomyces cerevisiae. Nature, Volume 440, pp. 637-643.

[4] McDowall, M. D., Scott, M. S. & Barton, G. J., 2008. PIPs: human protein-protein interaction prediction database. Nucleic Acids Research, 11, Volume 37, pp. D651-D656.

[5] Nepusz, T., Yu, H. & Paccanaro, A., 2012. Detecting overlapping protein complexes in protein-protein interaction networks. Nature Methods, Volume 9, pp. 471-472.

[6] Szklarczyk, D. et al., 2014. STRING v10: protein-protein interaction networks, integrated over the tree of life. Nucleic Acids Research, 10, Volume 43, pp. D447-D452.





[1] Collins, S. R. et al., 2007. Toward a Comprehensive Atlas of the Physical Interactome of Saccharomyces cerevisiae. Molecular & Cellular Proteomics, Volume 6, p. 439–450.

[2] Gavin, A. et al., 2006. Proteome survey reveals modularity of the yeast cell machinery. Nature.

[3] Krogan, N. J. et al., 2006. Global landscape of protein complexes in the yeast Saccharomyces cerevisiae. Nature, Volume 440, pp. 637-643.

[4] McDowall, M. D., Scott, M. S. & Barton, G. J., 2008. PIPs: human protein–protein interaction prediction database. Nucleic Acids Research, 11, Volume 37, pp. D651-D656.

[5] Nepusz, T., Yu, H. & Paccanaro, A., 2012. Detecting overlapping protein complexes in protein-protein interaction networks. Nature Methods, Volume 9, pp. 471-472.

[6] Szklarczyk, D. et al., 2014. STRING v10: protein–protein interaction networks, integrated over the tree of life. Nucleic Acids Research, 10, Volume 43, pp. D447-D452.

