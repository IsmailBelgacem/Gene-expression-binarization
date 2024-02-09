# Gene-expression-binarization-R-programing
Gene expression binarization is crucial for synthesizing Boolean gene regulatory network models (GRNs). Boolean formulas describe the evolution of the dynamics of Boolean GRNs, and the problem of the gene expression binarization is essential to deduce these formulas. The binarization problem is complicated because it depends on the gene's functional roles, and there are always zones of uncertainty. Moreover, the gene expression experiments data generally have only some instantaneous experiments. In this project, we developed a novel gene expression data binarization method. This method is based on the gene expression regulations. Our proposed method is applicable for instantaneous gene expression data sets, even if we have only one measurement at the steady state. Moreover, if the measured expression value for a gene is not provided, then based on the Boolean regulation rules, it can also be binarized using the binary values of the gene neighbors. The algorithm for our binarization method has been tested and validated as well. Using our binarization method, we proved that the genes are correctly binarized based on the ODE simulations of artificial examples of gene regulatory networks or well-known examples of Boolean biological networks. For the distribution to the community, we provide the code of this binarization method implemented using the R language. We also tested here our algorithm on real RNA-seq gene expression data.
