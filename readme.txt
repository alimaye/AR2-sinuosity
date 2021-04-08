April 7, 2021
Ajay B. Limaye (ajay@virginia.edu)
University of Virginia

This repository contains files for modeling river centerlines and calculating fit parameters for natural river centerlines using the second-order autogressive model described by Ferguson (1976) and reproduced by Limaye et al. (submitted). All codes are written in MATLAB and have been tested for version R2020a. Ajay B. Limaye authored all code, except where indicated. This software is provided as-is, with no performance guarantees, and is subject to the included license file. Please submit any feedback, including issues running the code, to Ajay B. Limaye (ajay@virginia.edu). If any part of this code is used in preparing a publication, please cite Limaye et al. (submitted) using the reference below.

The folder "AR2 model" contains the following codes that generate modeled centerlines:

-wrapper_modelCenterlineAR2.m: A wrapper script to define parameters and execute centerline generation using modelCenterlineAR2.m.

-modelCenterlineAR2.m: generates channel planforms using the second-order
autoregressive (AR2) model. Calls pathBanks.m, loopLength.m, and my_intersections.m.

-pathBanks.m: Creates river channel bank coordinates from a centerline ('path') and an associated channel width ('w'). Bank geometry is output for each bank as well as the full channel boundary (a polygon).

-loopLength.m: Checks for loops as neck cutoffs in the input centerline ('path') with a fixed channel width ('w') and returns the centerline with loops removed ('path') and the fraction of nodes in the original centerline that belonged to loops ('fractionNodesCutoff').

-my_intersections.m: Identifies self-intersection of the channel banks. This code is slightly modified from the original third-party software intersections.m by Douglas M. Schwartz, accessed via the MATLAB File Exchange (https://www.mathworks.com/matlabcentral/fileexchange/11837-fast-and-robust-curve-intersections). 


The folder "AR2 fit" contains the following files that calculate fit parameters for the autoregressive statistical model, applied to natural centerlines: 

- BeattonRiverCenterline.mat: an example data file for a natural river channel, the Beatton River, Alberta, Canada (57.1°N, 121.1°W).  

- wrapper_ARestimator_movingWindow.m: wrapper script that defines parameters and executes statistical fit using ARestimator.m.

- AR2stimator.m: function to estimate parameters for an autoregressive model with Gaussian disturbance, for an input angle series theta.

- getAzimuthSeries.m: Translates in input coordinate series (x,y) to a direction series (azimuth). 


References: 

Ferguson (1976), Disturbed periodic model for river meanders, Earth Surface Processes 1(4), 337-347. 

Limaye, A. B., Lazarus, E. D., Li, Y., and Schwenk, J., submitted, River sinuosity describes a continuum between randomness and ordered growth.