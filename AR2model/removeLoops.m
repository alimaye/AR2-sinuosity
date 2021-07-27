function [centerline,fractionNodesCutoff]=removeLoops(centerline,w,deltaS)
% removeLoops.m: Checks for loops as neck cutoffs in the input centerline
% ('centerline') with a fixed channel width ('w') and node spacing (deltaS) and returns the centerline 
% with loops removed ('centerline') and the fraction of nodes in the
% original centerline that belonged to loops ('fractionNodesCutoff'). Any 
% gaps in the centerline with loops removed are filled in using linear
% interpolation with the original node spacing.
% Created June 29, 2020 by Ajay B. Limaye, University of Virginia (ajay@virginia.edu).
% Last edited July 26, 2021 by Ajay B. Limaye.

centerline_in=centerline;% save in case needed for debugging 
fractionNodesCutoff = 0; % initialize as zero, will be overwritten if cutoffs occur        
[leftBank,rightBank,~] = pathBanks(centerline,w); % find banks

% Locate points of self-intersection of channel banks
[xo1,yo1,u1,d1]=my_intersections(leftBank.x,leftBank.y); % indices are upstream, downstream indices of self-intersection
[xo2,yo2,u2,d2]=my_intersections(rightBank.x,rightBank.y);

    if ~isempty(u1) || ~isempty(u2) % If intersections are located...
        % Order the node index and (x,y) coordinates for each intersection
        u=[u1;u2];
        d=[d1;d2];
        x_int=[xo1;xo2]; 
        y_int=[yo1;yo2];

        % Check for and remove and NaN values
        if any(isnan(u))
            indRemove = find(isnan(u));
            u(indRemove)=[];
            d(indRemove)=[];
            x_int(indRemove)=[];
            y_int(indRemove)=[];
        end

        % Enact the cutoff by removing nodes contained within the cutoff
        % loop
        centerline_nodes_remove=zeros(100,1);
        count_temp=1;
        for m=1:numel(u)
            nodes_remove = ceil(u(m)):floor(d(m));
            nn_remove = numel(nodes_remove);
            centerline_nodes_remove(count_temp:(count_temp+nn_remove-1))=nodes_remove;
            count_temp=count_temp+nn_remove;
        end
        centerline_nodes_remove(count_temp:end)=[];
        centerline_nodes_remove=sort(centerline_nodes_remove);
        centerline_nodes_remove = unique(centerline_nodes_remove); % in this model, some nodes can be part of multiple cutoffs
        
        % nodes within a cutoff have consecutive indices, so number of
        % cutoffs equals diff(centerline_nodes_remove)>1
       
        % Calculate the fraction of nodes in the centerline that belong to
        % a cutoff loop. Because the nodes spacing is constant, this
        % quantity is equal to the fraction length of the channel in cutoff
        % loops.
        fractionNodesCutoff = numel(centerline_nodes_remove)/numel(centerline.x);
        
        % remove cutoff nodes from the centerline
        centerline.x(centerline_nodes_remove)=[]; 
        centerline.y(centerline_nodes_remove)=[]; 
        
        if any(centerline_nodes_remove)
            % check for gaps in the trimmed centerline and fill in using linear
            % interpolation
            i = 1;
            tol = 2*deltaS; % set gap tolerance to twice the original node spacing
            while i<numel(centerline.x)
                % check distance to next node
                distToNextNode = sqrt((centerline.x(i+1)-centerline.x(i))^2 + (centerline.y(i+1)-centerline.y(i))^2);
                if distToNextNode >= tol % if distance to next node is greater than the tolerance...
                    % insert node(s) at original node spacing using linear interpolation
                    xInsert = interp1([0;distToNextNode],[centerline.x(i);centerline.x(i+1)],(deltaS:deltaS:distToNextNode-deltaS)','linear');
                    yInsert = interp1([0;distToNextNode],[centerline.y(i);centerline.y(i+1)],(deltaS:deltaS:distToNextNode-deltaS)','linear');
                    centerline.x = [centerline.x(1:i);xInsert;centerline.x(i+1:end)];
                    centerline.y = [centerline.y(1:i);yInsert;centerline.y(i+1:end)];
                    i = i+numel(xInsert); % advance the iteration beyond the inserted nodes
                end
                i=i+1; % advance to next node
            end
        end
    end
end