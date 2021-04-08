% wrapperARestimator_movingWindow.m: Estimates autoregression coefficients (b1,b2,sigma)
% for fixed step length (deltaS) for natural rivers using a moving window.
% Created July 1, 2020 by Ajay B. Limaye, University of Virginia
% (ajay@virginia.edu). 
% Last edited April 8, 2021 by Ajay B. Limaye.

deltaS = 1; % centerline point spacing, in channel widths; set for 
% consistency with spacing used for AR2 model centerlines 

fitWindow = 10; % number of points in moving window for AR2 model fit

% Define output directory and create it if it doesn't exist
outputDir = [pwd,'\fitParameters\'];
if ~exist(outputDir)
    mkdir(outputDir);
end

% Define output file name
outputFilename = 'naturalAR2fits.mat';

% Define input file(s) with river (x,y) coordinates
sourceDir = pwd;
locationName = {'BeattonRiver'};
inputFilename{1} = [sourceDir,'\BeattonRiverCenterline.mat'];


fit = struct; % initialize structure array to store fit parameters 

% For each location (and corresponding coordinate file)...
for i=1:numel(locationName)

    % Load the centerline coordinates and metadata
    load(inputFilename{i},'centerline','metadata')
    
    % Subtract offsets from coordinates (x,y) coordinates 
    % and divide both by channel width
    path.original.xWidths = (centerline.X-centerline.X(1))/metadata(1).widthMean;
    path.original.yWidths = (centerline.Y-centerline.Y(1))/metadata(1).widthMean;
    
    % Calculate node-to-node distances, in channel widths
    path.original.distanceWidths = [0;cumsum(sqrt(sum(diff([path.original.xWidths,path.original.yWidths],1,1).^2,2)))];  
    
    % calculate direction series, with angle wrapping to remove step changes in direction series greater than pi/2
    path.original.theta = getAzimuthSeries(path.original.xWidths,path.original.yWidths);
    
    % Interpolate for the direction (theta) at the defined centerline node spacing. Note that
    % spline interpolation yields slight variation in the sampling
    % distance. 
    distanceSample = (0:deltaS:max(path.original.distanceWidths))'; % vector of distances at which to interpolate for direction
    path.interpolated.xWidths = interp1(path.original.distanceWidths,path.original.xWidths,...
                                    distanceSample,'spline'); % interpolate for x-coordinates, in units of channel width
    path.interpolated.yWidths = interp1(path.original.distanceWidths,path.original.yWidths,...
                                    distanceSample,'spline'); % interpolate for y-coordinates, in units of channel width
    path.interpolated.distanceWidths = [0;cumsum(sqrt(sum(diff([path.interpolated.xWidths,path.interpolated.yWidths],1,1).^2,2)))];  % calculate the actual distance between centerlrine nodes (in channel widths) after interpolation
    straightLineLength =  sqrt(diff(path.interpolated.xWidths([1 end]))^2+diff(path.interpolated.yWidths([1 end]))^2); % Calculate the straight-line length of the channel
    path.interpolated.sinuosity =path.interpolated.distanceWidths(end)/straightLineLength; % sinuosity of the full centerline is along-stream length divided by straight-line length

    % Convert the (x,y) coordinates to a direction series
    path.interpolated.theta = getAzimuthSeries(path.interpolated.xWidths,path.interpolated.yWidths);
    
    % Fit AR parameters to the centerline using a moving window   
    
    % initialize arrays to store fit parameters as NaN, with one value for
    % each node of the centerline
    N = numel(path.interpolated.theta); % number of nodes in centerline
    b1 = nan(N,1); % fit coefficient b1
    b2 = nan(N,1); % fit coefficient b2
    sigma = nan(N,1); % fit coefficient sigma
    fh = waitbar(0,['Fitting AR parameters for ',locationName{i}]); % create status bar
    for j=1:N % for each centerline node...
        waitbar(j/N) % update states bar
        ind = (j-fitWindow/2):(j+fitWindow/2); % set node indices to use in AR fit
        if all(ind>0) && all(ind<=N) % only fit the AR model if window around the node fits within the length of the centerline
            theta = path.interpolated.theta(ind);
            ARparameters=ARestimator(theta); % fits AR parameters to direction series
            % Extract b1, b2, sigma from structure array
            b1(j) = ARparameters.AR1_coeff;
            b2(j) = ARparameters.AR2_coeff;
            sigma(j) = ARparameters.sigma;
        end
    end 
    close(fh)
    
    % AR parameters are now fit to every node in the centerline, with
    % non-fitted nodes remaining as NaN. Calculate the median of each fit
    % parameter, excluding the NaN values.
    b1  = nanmedian(b1);
    b2  = nanmedian(b2);
    sigma  = nanmedian(sigma);
   
    % Save autoregressive fit parameters to structure array
    fit(i).deltaS = deltaS;
    fit(i).ARparameters.b1 = b1;
    fit(i).ARparameters.b2 = b2;
    fit(i).ARparameters.sigma = sigma;
    fit(i).locationName = locationName{i};
    fit(i).path = path;
end

% Export fit parameters to file
save([outputDir,outputFilename],'fit')