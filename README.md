# ShockFrontedFourthOrder
MATLAB and Mathematica notebooks to compute shock-fronted travelling waves of the RND PDE with fourth-order spatial regularisation $-\varepsilon^2 \partial x^4$



README

…/ShockFrontedFourthOrder




File types:

*.m  -  MATLAB notebooks
*.mat - MATLAB data arrays
*.nb - Mathematica notebooks
*.cpp - C++ code used with MATLAB implementation of DOP853 ODE solver


SOFTWARE VERSIONS USED BY AUTHOR:
MATLAB 2021a 
Mathematica 12.0.0.0

YMMV.



WARNING: Running MATLAB files with .cpp dependencies requires DOP853 code to be available and executable in the working directory.
See e.g. http://www.unige.ch/~hairer/software.html. But these parts are easily modified to use MATLAB ODE solvers instead.



  
File descriptions and dependencies:

-- interpolredyifei.m 

Requires delwhcparams2.mat and delparams3.mat (parameter values at which singular heteroclinics at eps = 0 exist)

Requires interpred2.cpp to integrate two-dimensional reduced problem

Output: reducedwave.mat, coordinates for reduced singular heteroclinic for equal area rule shock-fronted travelling waves


-- interpolreg_bvp_shock.m

Requires reducedwave.mat (data for reduced singular heteroclinic orbit)

Output: 4-dimensional shock-fronted travelling wave for eps > 0



