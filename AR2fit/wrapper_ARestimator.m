% wrapper_ARestimator.m: Estimates coefficients (b1,b2,sigma) for
% second-order autoregressive model for natural rivers. Input coordinates for channel
% centerlines should be interpolated to a consistent node spacing equal to
% the mean channel width.
% Created July 1, 2020 by Ajay B. Limaye, University of Virginia
% (ajay@virginia.edu). 
% Last edited July 26, 2021 by Ajay B. Limaye.

% Define input and output file names
inputFilename = 'demoCenterline.mat'; % a demo centerline 
outputFilename = 'AR2fits.mat';

% Load the input centerline coordinates
load(inputFilename)

fit = struct; % initialize structure array to store fit parameters 
fh = waitbar(0,'Fitting AR parameters for centerlines...'); % create status bar
for i=1:numel(centerline)
    % convert coordinates to a direction series 
    try
        theta = getAzimuthSeries(centerline(i).x,centerline(i).y);
    catch
        theta = NaN; % can happen if only 1 coordinate
    end
    
    ARparameters=ARestimator(theta); % fit AR parameters to full direction series
    fit(i).ARparameters.b1 = ARparameters.b1;
    fit(i).ARparameters.b2 = ARparameters.b2;
    fit(i).ARparameters.sigma = ARparameters.sigma;
end
close(fh)

% Export fit parameters to file
save(outputFilename,'fit')