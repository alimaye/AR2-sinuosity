function [centerline,unstableFlag] = modelCenterlineAR2(b1,b2,epsilon,deltaS)
% modelCenterlineAR2.m: generates channel planforms using the second-order
% autoregressive (AR2) model introduced in Ferguson (1976), Disturbed periodic
% model for meanders, Earth Surface Processes 1(4), 337-347,
% doi:10.1002/esp.3290010403.
% Inputs: AR2 model parameters (b1, b2, epsilon) and centerline node spacing
% (deltaS).
% Outputs: centerline (the channel centerline) and unstableFlag (flag to
% identify cases with unstable direction series). 
% Created June 29, 2020 by Ajay B, Limaye, University of Virginia (ajay@virginia.edu). 
% Last edited June 22, 2021 by Ajay B. Limaye. 

% initialize direction of first two steps as 0 (i.e., along x-axis in
% positive direction). 
theta(1) = 0;
theta(2) = 0;

% Calculate direction series (theta) using second-order autogressive formula and based
% on the input parameters.
for i=3:numel(epsilon)
    theta(i) = b1*theta(i-1)+b2*theta(i-2)+epsilon(i);
end

unstableFlag=any(isinf(theta)); % Some combinations of parameters cause theta to increase without bound (i.e., become unstable). Identify those cases with this flag.
%%% 0 --> stable direction series
%%% 1 --> unstable direction series

% Initialize centerline as structure array, 'path'. The first step of the
% path is along the x-axis in the positive direction.
centerline.x = [0;deltaS];
centerline.y = [0;0];

% Calculate each step of coordinate series using the direction series, theta.
for i=3:numel(theta)
    centerline.x(i) = centerline.x(i-1)+deltaS*cos(theta(i-1));
    centerline.y(i) = centerline.y(i-1)+deltaS*sin(theta(i-1));
end

% As a check, throw an error if any NaNs were created in the centerline coordinates
if any(isnan(centerline.x))
    error('NaN created in centerline coordinates');
end

end