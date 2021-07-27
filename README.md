# AR2-sinuosityJuly 26, 2021
Ajay B. Limaye (ajay@virginia.edu)
University of Virginia
ajay@virginia.edu

This repository contains files for modeling and calculating fit parameters for river centerlines using the second-order autogressive model described by Ferguson (1976) and reproduced by Limaye et al. (2021). All codes are written in MATLAB and have been tested for version R2020a. Ajay B. Limaye authored all code, except where indicated. This software is provided as-is, with no performance guarantees, and is subject to the included license file. Please submit any feedback, including issues running the code, via email to ajay@virginia.edu. If any part of this code is used in preparing a publication, please cite Limaye et al. (2021) using the reference below.

The folder "AR2 model" contains the following codes that generate modeled centerlines:

-wrapper_modelCenterlineAR2.m: A wrapper script to define parameters and execute centerline generation using modelCenterlineAR2.m.

-modelCenterlineAR2.m: generates channel planforms using the second-order autoregressive (AR2) model. Calls pathBanks.m, removeLoops.m, and my_intersections.m.

-pathBanks.m: Creates river channel bank coordinates from a centerline and an associated channel width. 

-removeLoops.m: Checks for loops as neck cutoffs in the input centerline with a fixed channel width and returns the centerline with loops removed  and the fraction of nodes in the original centerline that belonged to loops.

-my_intersections.m: Identifies self-intersection of the channel banks. This code is slightly modified from the original third-party software intersections.m by Douglas M. Schwartz, accessed via the MATLAB File Exchange (https://www.mathworks.com/matlabcentral/fileexchange/11837-fast-and-robust-curve-intersections). 

-windowedSinuosity.m: Calculates channel sinuosity using a moving average. 

The folder "AR2 fit" contains the following files that calculate fit parameters for the autoregressive statistical model: 

-demoCenterline.mat: an example data file a channel centerline coordinates (x,y), with the centerline sampled at unit spacing.

-wrapper_ARestimator.m: a wrapper script that defines parameters and executes statistical fitting using ARestimator.m.

-AR2estimator.m: function to estimate parameters for an autoregressive model for an input angle series theta.

-getAzimuthSeries.m: Translates an input coordinate series (x,y) to a direction series.


References: 

Ferguson, R. I., 1976, Disturbed periodic model for river meanders, Earth Surface Processes 1(4), 337-347, https://doi.org/10.1002/esp.3290010403

Limaye, A. B., Lazarus, E. D., Li, Y., and Schwenk, J., in press, River sinuosity describes a continuum between randomness and ordered growth, Geology, https://doi.org/10.1130/G49153.1
