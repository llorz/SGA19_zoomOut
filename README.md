# ZoomOut: Spectral Upsampling for Efficient Shape Correspondence
This is an example code for the paper "ZoomOut: Spectral Upsampling for Efficient Shape Correspondence" by Simone Melzi, Jing Ren, Emanuele Rodolà, Abhishek Sharma, Peter Wonka, and Maks Ovsjanikov.

In this paper, we propose a very simple and efficient map refinement method called "zoomOut" for refining maps or correspondences by iterative upsampling in the spectral domain.

<p align="center">
  <img align="center"  src="/figs/eg_cat_wolf.png">
</p>


Main Functions
--------------------
```
T12_refined = zoomOut_refine(B1, B2, T12, para)

% Input:
%   B1: The LB basis of the source shape S1
%   B2: The LB basis of the target shape S2
%   T12: the initial point-wise map from S1 to S2
%   para: a structure stores the following parameters
%       k_init: the initial dimension of functional map
%       k_step: the step size of chaning the functional map
%       k_final: the final dimension (zoomOut) of the functional map
% Output:
%   T12_refined: the refined point-wise map
```

```
T12_refined = zoomOut_refine_fast(S1, S2, T12, para, seed)

% Input:
%   S1: The source shape S1 with precomputed LB basis
%   S2: The target shape S2 with precomputed LB basis
%   T12: the initial point-wise map from S1 to S2
%   para: a structure stores the following parameters
%       k_init: the initial dimension of functional map
%       k_step: the step size of chaning the functional map
%       k_final: the final dimension (zoomOut) of the functional map
%       num_samples: the sampling density to speed-up the zoomOut refinement
%   seed: the seed for the random sampling
% Output:
%   T12_refined: the refined point-wise map
```

Comments
-------------------------
- The script ```eg1_cat_wolf.m``` shows how to use the above function (to reproduce the Fig.2 in our paper). 
- The script ```eg2_self_symmetry.m``` shows how to compute self-symmetric maps for human shapes.
- You can also find the code for initial map computation (both for self-symmetric maps or maps between two shapes) [here](https://github.com/llorz/SGA18_orientation_BCICP_code)
- Please let us know (jing.ren@kaust.edu.sa, melzismn@gmail.com) if you have any question regarding the algorithms/paper ୧((〃•̀ꇴ•〃))૭⁺✧ or you find any bugs in the code ʃ͠ʘɷʘ͠ƪ
- **(update 2019-10-14)** As pointed out by **@jayis520**, you might get different results from ours for the **cat-dog** example. The reason might be that your computed LB basis of the cat and the wolf does not align in the same way as ours. In this case, please change line 27 in the script ```eg1_cat_wolf.m``` to ```C12_ini = diag([1, -1, -1, -1])``` for example to align the first four basis of the cat with the wolf that are computed in your machine. Please check the [issue](https://github.com/llorz/SGA19_zoomOut/issues/1) for more explanation. Sorry for the confusion :(




[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)

This work is licensed under a [Creative Commons Attribution-NonCommercial 4.0 International License](http://creativecommons.org/licenses/by-nc/4.0/). For any commercial uses or derivatives, please contact us (jing.ren@kaust.edu.sa, peter.wonka@kaust.edu.sa, maks@lix.polytechnique.fr).
