function [leftBank,rightBank,bankAll] = pathBanks(path,w)
% pathBanks.m: Creates river channel bank coordinates from a centerline ('path') and an associated channel width ('w').
% Bank geometry is output for the each bank (leftBank,rightBank) as well as the full channel boundary (a polygon).
% Created April 9, 2020 by Ajay B. Limaye, University of Virginia
% (ajay@virginia.edu).
% Last edited April 7, 2021 by Ajay B. Limaye.

% Create vectors to represent channel azimuth
v1=[[NaN;path.x(2:end)-path.x(1:end-1)],[NaN;path.y(2:end)-path.y(1:end-1)]];
v2=[[path.x(2:end)-path.x(1:end-1);NaN],[path.y(2:end)-path.y(1:end-1);NaN]];

% Normalize v1 and v2 to get unit vectors
ind0=v1==0;
v1 = v1./(repmat(sqrt(sum(v1.^2,2)),1,2));
v1(ind0)=0; % avoids getting NaN from divide by zero

ind0=v2==0;
v2 = v2./(repmat(sqrt(sum(v2.^2,2)),1,2));
v2(ind0)=0; % avoids getting NaN from divide by zero

% Calculate azimuth tangent to each point on the cennel centerline
v3=v1+v2;
tangent_az=atan2(v3(:,2),v3(:,1));

% Set azimuth for endpoints of the centerline
tangent_az(1) = 0;
tangent_az(end) = tangent_az(end-1);

% Set coordinates of left bank and right bank
leftBank.x = path.x+(w/2)*cos(tangent_az+pi/2);
leftBank.y = path.y+(w/2)*sin(tangent_az+pi/2); 
rightBank.x= path.x+(w/2)*cos(tangent_az-pi/2);
rightBank.y= path.y+(w/2)*sin(tangent_az-pi/2);

% Save ordered sequence of all bank coordinates to an additional structure
% array. These coordinates define a polygon that corresponds to the
% "footprint" of the channel.
bankAll.x=[leftBank.x(end:-1:1);rightBank.x]; 
bankAll.y=[leftBank.y(end:-1:1);rightBank.y];
end