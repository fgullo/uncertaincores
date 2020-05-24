# Core Decomposition of Uncertain Graphs
<strong>Copyright 2013-2020, Francesco Bonchi (ISI Foundation, Italy), Francesco Gullo (UniCredit, Italy), Andreas Kaltenbrunner (NTENT, Spain), Yana Volkovich (Xandr, USA)</strong>
<p>

This is a Java software package that implements algorithms described in the paper "Core Decomposition of Uncertain Graphs", KDD 2014.<br>
This package is free for research, academic and non-profit making purposes only. If you use this piece of software for your work and got something published please include the citation reported below. The software may not be sold or redistributed without prior approval. One may make copies of the software for their use provided that the copies, are not sold or distributed, are used under the same terms and conditions. As unestablished research software, this code is provided on an "as is" basis without warranty of any kind, either expressed or implied. The downloading, or executing any part of this software constitutes an implicit agreement to these terms. These terms and conditions are subject to change at any time without prior notice.


<strong>TERMS OF USAGE:</strong>
The following paper should be cited in any research product whose findings are based on the code here distributed:

- F. Bonchi, F. Gullo, A. Kaltenbrunner, Y. Volkovich.<br>
[Core Decomposition of Uncertain Graphs](https://doi.org/10.1145/2623330.2623655).<br>
Proceedings of the 20th ACM SIGKDD International Conference on Knowledge Discovery and Data Mining (KDD 2014), pp. 1316-1325. New York, NY, USA - August 24-27, 2014.
<p>
  
## Running

### `KEtaCores_DynamicProgramming_BigDecimals_IncrementalInitialEtaDegrees.java`<br>

This class implements Algorithm 2 in the paper. It can be run as follows:

```
java KEtaCores_DynamicProgramming_BigDecimals_IncrementalInitialEtaDegrees $DATASET $ETA $RUNS
```
where `$DATASET` is the dataset, `$ETA` is the `\eta \in [0,1]` threshold of the resulting (k,\eta)-core decomposition, and `$RUNS` is the number of runs the algorithm is executed for.

**IMPORTANT**: the dataset file should be in the same folder as the class.

<br>

### `KEtaCores_DynamicProgramming_BigDecimals_IncrementalInitialEtaDegrees_LB.java`

This class implements an imporved (yet faster) version of Algorithm 2 in the paper. 
It is based on the use of a lower-bound on \eta-degrees, which is described in detail in Section 4 in the paper. 
It can be run as follows:

```
java KEtaCores_DynamicProgramming_BigDecimals_IncrementalInitialEtaDegrees_LB $DATASET $ETA $RUNS
```
where `$DATASET` is the dataset, `$ETA` is the `\eta \in [0,1]` threshold of the resulting (k,\eta)-core decomposition, and `$RUNS` is the number of runs the algorithm is executed for.

The class needs a file containing incomplete-beta-function values for the dataset at hand. 
Such a file should be placed into the `beta_function_values` folder, and should be named `$DATASET_beta_matlab.txt_output`.
A Matlab script to compute such a file is in `matlab/beta_function_values_v3.m` (but it can clearly computed by using any other script or language).


**IMPORTANT**: the dataset file and the `beta_function_values` folder should be in the same folder as the class.
