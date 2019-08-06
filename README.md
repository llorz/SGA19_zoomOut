# ZoomOut: Spectral Upsampling for Efficient Shape Correspondence
This is an example code for the paper "ZoomOut: Spectral Upsampling for Efficient Shape Correspondence" by Simone Melzi, Jing Ren, Emanuele Rodol\`a, Abhishek Sharma, Peter Wonka and Maks Ovsjanikov.

In this paper, we propose a new regularizer, the complex resolvent Laplacian commutativity, for functional map pipeline. This term is
theoretically justified and can improve the quality of the computed functional maps and the corresponding recovered point-wise maps
before and after refinement.

<p align="center">
  <img align="center"  src="/figs/eg_cat_wolf.png">
</p>


Main Function
--------------------
```
[C12, Mask] = compute_fMap_complRes(S1,S2,B1,B2,Ev1,Ev2,fct1,fct2,mask_type)

% Input:
%   S1: the source mesh with the new basis B1, and the corresponding eigenvalues Ev1 (k1 Eigen-functions)
%   S2: the target mesh with the new basis B2, and the corresponding eigenvalues Ev2 (k2 Eigen-functions)
%   fct1: the descriptors of shape S1
%   fct2: the descriptors of shape S2
%   mask_type: 'standard', 'slant', or 'complRes' 
% Output:
%   C12: a functional map from S1 -> S2 (k2-by-k1 matrix)
%   Mask: the penalty mask used to regularize C12
```
- Compute a functional map with different penalty mask term:
  - "**standard**": the standard Laplacian commutativity term (formulated as a mask)
  - "**slant**": a heuristic slanted diagonal penalty mask proposed in ["Partial Functinal Correspondences"](https://arxiv.org/abs/1506.05274)
  - "**complRes**": our proposed complex resolvent Laplacian mask


Comments
-------------------------
- The script ```example.m``` shows how to use the above function (to reproduce the Fig.2 in our paper). 
- You can find our paper [here](https://www.dropbox.com/s/ctvor2e25eaaev6/2019SGP_Structured_Regularization_fMap.pdf?dl=0)
- The dataset used in the paper can be found [here](https://github.com/llorz/SGA18_orientation_BCICP_dataset)
- Please let us know (jing.ren@kaust.edu.sa, melzismn@gmail.com) if you have any question regarding the algorithms/paper ฅ^•ﻌ•^ฅ or you find any bugs in the code \_(°ω°｣∠)\_
