% wrapper_modelCenterlineAR2.m: wrapper script to generate centerlines 
% using a second-order autoregressive model adapted from Ferguson (1976), 
% Disturbed periodic model for river meanders, Earth Surface Processes
% 1(4), 337-347, doi:10.1002/esp.3290010403.
% Centerlines are constructed using a fixed node spacing and fixed centerline
% length in the x-direction. Three parameters are systematically varied to
% create the centerlines: b1 (the lag-1 coefficient), b2 (the lag-2 
% coefficient), and sigma (the standard deviation of "innovation", i.e.,
% direction noise). For each parameter set, several realizations are made
% using different random noise series. The output includes the original
% centerline ("centerlineRaw"), a centerline version with any self-intersections
% removed ("centerlineNoIntersections"), a flag to identify unstable model
% runs ("unstableThetaFlag") and the fraction of the original centerline length that
% is composed by cutoffs ("fractionNodesCutoff").
% Created June 29, 2020 by Ajay B. Limaye, University of Virginia
% (ajay@virginia.edu). 
% Last edited June 26, 2021 by Ajay B. Limaye.

clear,clc % clear variables and command window
dbstop if error % pause execution if an error is thrown

% Set file output directory and name; create the directory if it does not
% exist
outDir = [pwd,'\data\'];
if ~exist(outDir)
    mkdir(outDir)
end

% Set fixed parameters for generated modeled centerlines
w = 1; % channel width, in widths
deltaS = 1; % step length, in channel widths
maxSteps = round(1.5*100/deltaS); % set maximum number of steps (This case equates 
% to 100 channel widths of straight-line distance for a natural channel with meander wavelenghth = 10 channel widths and sinuosity = 1.5)

% Set range of values to use for each parameter in the AR2 model.
b1 = 0.1:0.1:2.0; % Values of the first coefficient, b1. Skipping 0.
b2 = -0.1:-0.1:-1; % Values of the second coefficient, b2. Skipping 0.
sigma = 0.1:0.1:0.5; % Values of the standard deviation of "innovation," sigma.
nReplicateSets = 10; % number of replicate parameter sets to use, each with a different random seed

% Set name for output file
outName = ['AR2_modelData_allParameterSets_',num2str(nReplicateSets),'replicates.mat']; 
outNameNoRaw = ['AR2_modelData_allParameterSets_',num2str(nReplicateSets),'replicates_noRaw.mat']; 
outNameStatsOnly = ['AR2_modelData_allParameterSets_',num2str(nReplicateSets),'replicates_statsOnly.mat']; 

centerlines = struct; % Initialize structure array to store centerlines
nCenterlines = nReplicateSets*numel(b1)*numel(b2)*numel(sigma); % Calculate number of centerlines to be generated, for tracking progress of iteration loops
count = 0; % counter for number of centerlines generated

% Generate centerlines by loop through the specified model parameter
% combinations. The iteration index 'l' creates replicate centerlines with
% the same parameters (b1, b2, sigma) but a different random seed
% (epsilonUnscaled). 

fprintf('Running AR-2 model and calculating sinuosity...\n');
for l=1:nReplicateSets
    epsilonUnscaled = normrnd(0,1,[maxSteps,1]); % use same random number series for each replicate set
    for k=1:numel(sigma)
        epsilon = epsilonUnscaled*sigma(k); % scale the epsilon values by sigma
        for j=1:numel(b2)
            for i=1:numel(b1)
                count = count+1; % iterate counter

                % Store the parameter values in the 'centerlines' structure
                % array
                centerlines(i,j,k,l).parameters.AR1_coeff = b1(i);
                centerlines(i,j,k,l).parameters.AR2_coeff = b2(j);
                centerlines(i,j,k,l).parameters.sigma = sigma(k);
                % Create and store the "raw" centerline
                [centerlineRaw,unstableFlag] = modelCenterlineAR2(b1(i),b2(j),epsilon,deltaS);
                centerlines(i,j,k,l).centerlineRaw = centerlineRaw;
                centerlines(i,j,k,l).unstableFlag = unstableFlag;
                % Check for cutoffs in the "raw" centerline; save a version
                % of the centerline with no cutoffs as well as the fraction
                % of nodes in the "raw" centerline that belong to cutoff
                % loops.
                if ~centerlines(i,j,k,l).unstableFlag
                    [centerlineNoIntersections,fractionNodesCutoff]=removeLoops(centerlineRaw,w,deltaS); % removeLoops.m creates a version of the centerline with no intersections (cutoffs), and records the fraction of the raw centerline nodes bound in cutoff loops
                    centerlines(i,j,k,l).centerlineNoIntersections = centerlineNoIntersections; % save the variables from the previous line to the output structure array
                    centerlines(i,j,k,l).fractionNodesCutoff = fractionNodesCutoff;
                    
                    % Calculate sinuosity using a moving window.
                    % In model runs, the centerline node spacing is 1 and represents the
                    % channel width. This is the same node spacing used for all
                    % analyses, so the spacing is proper to calculate the
                    % windowed sinuosity using windowedSinuosity.m. That code
                    % uses a fixed window of 50 elements along-stream 
                    % (equivalent to 50 channel widths) to calculate the
                    % moving average sinuosity.
                    path.xWidths= centerlines(i,j,k,l).centerlineNoIntersections.x; % x-coordinates, which are already in units of channel width
                    path.yWidths = centerlines(i,j,k,l).centerlineNoIntersections.y; % y-coordinates, which are already in units of channel width
                    centerlines(i,j,k,l).windowedSinuosity = windowedSinuosity(path); % average sinuosity using moving window
                else
                    centerlines(i,j,k,l).centerlineNoIntersections = [];
                    centerlines(i,j,k,l).fractionNodesCutoff = [];
                    centerlines(i,j,k,l).windowedSinuosity = [];
                end
            end
        end
    end
    % Print progress of model execution to screen
    fprintf('Progress:%2.0f percent\n',100*count/nCenterlines)
end

% Write the full set of modeled parameters to structure array
parameters.b1 = b1;
parameters.b2 = b2;
parameters.sigma = sigma;
parameters.nReplicateSets = nReplicateSets;

% Export the modeled centerlines and parameters to a .mat file
save([outDir,outName],'centerlines','parameters','-v7.3') % version 7.3 enables saving larger file sizes