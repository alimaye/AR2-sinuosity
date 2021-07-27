function azimuth = getAzimuthSeries(x,y)
% getAzimuthSeries.m: Translates in input coordinate series (x,y) to a
% direction series (azimuth). 
% Created July 1, 2020 by Ajay B. Limaye, University of Virginia
% (ajay@virginia.edu). 
% Last edited April 8, 2021 by Ajay B. Limaye.

% Construct an Nx2 array of vectors from one centerline node to the next
v=[[x(2:end)-x(1:end-1);NaN],[y(2:end)-y(1:end-1);NaN]];

% set vector for endpoint equal to previous 
v(end,:) = v(end-1,:);

% Normalize v to get unit vectors
ind0=v==0;
v = v./(repmat(sqrt(sum(v.^2,2)),1,2));
v(ind0)=0; % avoids getting NaN from dividing by zero

% Convert vector to azimuth value
azimuth=atan2(v(:,2),v(:,1));

% Remove differences greater than pi by wrapping azimuth values to +/- 2*pi
for i=1:(numel(azimuth)-1)
    if (azimuth(i+1)-azimuth(i)) > pi
        azimuth((i+1):end) = azimuth((i+1):end) - 2*pi;
    elseif azimuth(i+1)-azimuth(i) < -pi
        azimuth((i+1):end) = azimuth((i+1):end) + 2*pi;
    end
end

end