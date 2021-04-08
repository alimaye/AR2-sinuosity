function [path,unstableThetaFlag] = modelCenterlineAR2(b1,b2,epsilon,deltaS,xMax)
% modelCenterlineAR2.m: generates channel planforms using the second-order
% autoregressive (AR2) model introduced in Ferguson (1976), Disturbed periodic
% model for meanders. 
% Inputs: AR2 model parameters (b1, b2, epsilon), centerline node spacing
% (deltaS), and length of centerline (xMax).
% Outputs: path (the centerline) and unstableThetaFlag (a logical flag to
% indicate whether centerline construction was halted due to crossing a
% threshold number of iterations). 
% Created June 29, 2020 by Ajay B, Limaye, University of Virginia (ajay@virginia.edu). 
% Last edited April 7, 2021 by Ajay Limaye. 

maxIt = 10*xMax/deltaS; % maximum number of iterations; this effectively corresponds to sinuosities greater than 10, which correspond to highly tangled loops.

unstableThetaFlag = false; % Initialize flag to note whether the direction series (theta) yields an unstable centerline. 
% initialize direction of first two steps as 0 (i.e., along x-axis in
% positive direction). 
theta(1) = 0;
theta(2) = 0;

% Calculate direction series (theta) using second-order autogressive formula and based
% on the input parameters.
for i=3:numel(epsilon)
    theta(i) = b1*theta(i-1)+b2*theta(i-2)+epsilon(i);
end

% Some combinations of parameters cause theta to spiral off, so note those 
% cases using unstableThetaFlag and truncate direction series if necessary
if any(isinf(theta))
    ind = find(isinf(theta),1,'first');
    theta = theta(1:ind-1);
    unstableThetaFlag = true;
end

% Initialize centerline as structure array, 'path'. The first step of the
% path is along the x-axis in the positive direction.
path.x = [0;deltaS];
path.y = [0;0];

% Calculate each step of coordinate series using the direction series, theta.
for i=3:numel(theta)
    path.x(i) = path.x(i-1)+deltaS*cos(theta(i-1));
    path.y(i) = path.y(i-1)+deltaS*sin(theta(i-1));
    if path.x(i) >= xMax || i>maxIt % stop building path once x-coordinate exceed x_max or maximum number of iterations reached
        break
    end
end

% As a check, throw an error if any NaNs were created in the centerline coordinates
if any(isnan(path.x))
    error('NaN created in centerline coordinates');
end

end