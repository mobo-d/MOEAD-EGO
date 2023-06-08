# MOEA/D-EGO

Open-source implementations for methods presented in the following papers: 

* **Q. Zhang, W. Liu, E. Tsang, and B. Virginas. Expensive multiobjective optimization by MOEA/D with Gaussian process model. IEEE Transactions on Evolutionary Computation, 2009. [[PDF](https://ieeexplore.ieee.org/abstract/document/5353656)]** <br/>

* **L. Zhao and Q. Zhang. Exact Formulas  for the Computation of Expected  Tchebycheff Improvement. Proceedings of the IEEE Congress on Evolutionary  Computation, 2023.** <br/>


Our implementation of MOEA/D-EGO differs slightly from the original in two aspects:

* The FuzzyCM is removed. FuzzyCM is an approximation  method  for GP modeling. It could  reduce the computational time when training the GP models. However, MOEA/D-EGO without FuzzyCM could perform better  in terms of solution quality.  Interested readers can find more related discussions in Section VII-E of [MOEA/D-EGO](https://ieeexplore.ieee.org/abstract/document/5353656). 

* An adaptive adjustment strategy for $z^*$ is used. More related discussions can be found in the [Supplementary File of DirHV-EGO](https://ieeexplore.ieee.org/document/10093980). 

## Acknowledgements
* This implementation is based on [PlatEMO](https://github.com/BIMK/PlatEMO).
* For GP modeling, we employe the [DACE toolbox](https://www.omicron.dk/dace.html).
* For Design of Experiment methods, we utilize existing implementations from the [SURROGATES Toolbox](https://sites.google.com/site/felipeacviana/surrogates-toolbox).
