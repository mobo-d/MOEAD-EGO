# MOEA/D-EGO

Open-source implementations for methods presented in the following papers: 

* **Q. Zhang, W. Liu, E. Tsang, and B. Virginas. Expensive multiobjective optimization by MOEA/D with Gaussian process model. IEEE Transactions on Evolutionary Computation, 2009. [[PDF](https://ieeexplore.ieee.org/abstract/document/5353656)]** <br/>

* **L. Zhao and Q. Zhang. Exact Formulas  for the Computation of Expected  Tchebycheff Improvement. Proceedings of the IEEE Congress on Evolutionary  Computation, 2023.** <br/>

The Java Code of MOEA/D-EGO (written by Wudong Liu) is avaliable at [this website](https://sites.google.com/view/moead/resources). Our implementation differs slightly from the vanilla MOEA/D-EGO in two aspects:

* The FuzzyCM is removed. FuzzyCM is an approximation  method  for GP modeling. It could  reduce the computational time when training the GP models. However, MOEA/D-EGO without FuzzyCM could perform better  in terms of solution quality.  Interested readers can find more related discussions in Section VII-E of [MOEA/D-EGO](https://ieeexplore.ieee.org/abstract/document/5353656). 

* An adaptive adjustment strategy for $z^*$ is used. More related discussions can be found in the [Supplementary File of DirHV-EGO](https://ieeexplore.ieee.org/document/10093980).  

## Usage 

Matlab >= 2018a

### Quick Start

* The `run_MOEAD_EGO.m` provides the basic script to run experiments on ZDT and DTLZ.

### Advanced usage

* Download [PlatEMO](https://github.com/BIMK/PlatEMO) (version 4.6) and read PlatEMO's [User Manual](https://github.com/BIMK/PlatEMO/blob/master/PlatEMO/manual.pdf) to familiarize yourself with how to use this platform.
* Copy the folders within "**Algorithms**" into the directory at **"PlatEMO/Algorithms/"**. Next, add all of the subfolders contained within the "PlatEMO" directory to the MATLAB search path.
* In the MATLAB command window, type **`platemo()`** to run PlatEMO using the GUI.
* Select the label "**expensive**" and choose the algorithm **"MOEA-D-EGO"**.
  * Default setting of `batch size`: 5.
* Select a problem and set appropriate parameters.
  * e.g., ZDT1, N=200, M=2, D=8, maxFE=200.
  * e.g., Inverted DTLZ2,  N=210, M=3, D=6, maxFE=300.


If you have any questions or feedback, please feel free to contact  liazhao5-c@my.cityu.edu.hk and qingfu.zhang@cityu.edu.hk.

## Acknowledgements
* This implementation is based on [PlatEMO](https://github.com/BIMK/PlatEMO).
* For GP modeling, we employe the [DACE toolbox](https://www.omicron.dk/dace.html).
 
