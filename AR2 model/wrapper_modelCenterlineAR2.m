% wrapper_modelCenterlineAR2.m: wrapper script to generate centerlines 
% a second-order autoregressive model (the disturbed periodic model of 
% Ferguson (1976). For a fixed centerline node spacing (step length), three 
% parameters are systematically varied: the AR-1 coefficient, the AR-2 
% coefficient, and sigma (the standard deviation of "innovation", i.e.,
% direction noise). Output centerlines are statistically characterized by 
% measuring their sinuosity and the  proportion of centerline length that
% is composed by cutoffs.
% Created June 29, 2020 by Ajay B. Limaye, University of Virginia
% (ajay@virginia.edu). 
% Last edited April 7, 2021 by Ajay B. Limaye.

clear,clc
dbstop if error % pause execution if an error is thrown

% Set file output directory and name; create the directory if it does not
% exist
outDir = [pwd,'\data\'];
if ~exist(outDir)
    mkdir(outDir)
end

% Set name for output file
outName = 'AR2_modelDataFinal.mat'; 

% Set fixed parameters for generated modeled centerlines
xMax = 100; % maximum domain length x-direction, in channel widths
w = 1; % channel width, in widths
deltaS = 1; % step length, in channel widths

% Set range of values to use for each parameter in the AR2 model.
b1 = 1:0.1:2; % Values of the first coefficient, b1.
b2 = -0.5:-0.1:-1; % Values of the second coefficient, b2.
sigma = [0.05 0.1 0.2 0.3]; % Values of the standard deviation of "innovation", sigma.
nReplicateSets = 10; % number of replicate parameter sets to use, each with a different random seed

centerlines = struct; % Initialize structure array to store centerlines
nCenterlines = nReplicateSets*numel(b1)*numel(b2)*numel(sigma); % Calculate number of centerlines to be generated, for tracking progress of iteration loops
count = 0; % counter for number of centerlines generated

% Generate centerlines by loop through the specified model parameter
% combinations. The iteration index 'l' creates replicate centerlines with
% the same parameters (b1, b2, sigma) but a different random seed
% (epsilonUnscaled). 
for l=1:nReplicateSets
    epsilonUnscaled = normrnd(0,1,[1e4,1]); % use same random number series for this parameter set
    for k=1:numel(sigma)
        epsilon = epsilonUnscaled*sigma(k);
        for j=1:numel(b2)
            for i=1:numel(b1)
                % Store the parameter values in the 'centerlines' structure
                % array
                centerlines(i,j,k,l).parameters.AR1_coeff = b1(i);
                centerlines(i,j,k,l).parameters.AR2_coeff = b2(j);
                centerlines(i,j,k,l).parameters.sigma = sigma(k);
                % Create and store the "raw" centerline
                [centerlineRaw,unstableThetaFlag] = modelCenterlineAR2(b1(i),b2(j),epsilon,deltaS,xMax);
                centerlines(i,j,k,l).centerlineRaw = centerlineRaw;
                centerlines(i,j,k,l).unstableThetaFlag = true;
                % Check for cutoffs in the "raw" centerline; save a version
                % of the centerline with no cutoffs as well as the fraction
                % of nodes in the "raw" centerline that belong to cutoff
                % loops.
                [centerlineNoIntersections,fractionNodesCutoff]=loopLength(centerlineRaw,w);
                centerlines(i,j,k,l).centerlineNoIntersections = centerlineNoIntersections;
                centerlines(i,j,k,l).fractionNodesCutoff = fractionNodesCutoff;
                count = count+1; % iterate counter
            end
            % Print progress of model execution to screen
            progressPct = 100*count/nCenterlines;
            fprintf('Progress:%2.0f percent\n',progressPct)
        end
    end
end

% Write the full set of modeled parameters to structure array
parameters.b1 = b1;
parameters.b2 = b2;
parameters.sigma = sigma;
parameters.nReplicateSets = nReplicateSets;

% Export the modeled centerlines and parameters to a .mat file
save([outDir,outName],'centerlines','parameters')