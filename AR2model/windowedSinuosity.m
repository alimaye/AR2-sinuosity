function sinuosityMean = windowedSinuosity(in)
% windowedSinuosity.m: Calculates average sinuosity using a moving average.
% This calculation assumes a fixed spacing between centerline
% nodes. For centerlines that have fewer nodes than the window length, the
% sinuosity is calculated using the full centerline with no moving average.
% Created July 17, 2020 by Ajay B. Limaye, University of Virginia
% (ajay@virginia.edu). 
% Last edited July 26, 2021 by Ajay B. Limaye.

windowLength = 50; % window length (number of array elements)
x = in.xWidths; % x-coordinates of input centerline
y = in.yWidths; % y-coordinates of input centerline
N = numel(x); % number of nodes in the centerline
  
if N<=windowLength
      length = sum(sqrt(sum(diff([x,y],1,1).^2,2))); % along-stream channel length
      straightLineDistance  = sqrt(diff(x([1 end]))^2+diff(y([1 end]))^2); % straight-line channel length
      sinuosityMean = length/straightLineDistance; % sinuosity 
else    
    halfWindow = round(windowLength/2); % half of window length
    sinuosity = nan(N,1); % preallocate an array to store sinuosity locally along the centerline
    for i=1:numel(x) % for each element of the centerline...
        if and(i> halfWindow, i+halfWindow<=N) % if the window fits fully around the current node
            xWindow = x((i-halfWindow):(i+halfWindow)); % get the x-coordinates within the window..
            yWindow = y((i-halfWindow):(i+halfWindow)); % .. and the y-coordiantes too
            length = sum(sqrt(sum(diff([xWindow,yWindow],1,1).^2,2))); % along-stream channel length
            straightLineDistance = sqrt(diff(xWindow([1 end]))^2+diff(yWindow([1 end]))^2); % straight-line channel length
            sinuosity(i) =  length/straightLineDistance; % sinuosity for this node
        end
    end
    sinuosityMean = nanmean(sinuosity); % calculate the mean of all the windowed sinuosity values, excluding NaN
end