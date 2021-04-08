function [path,fractionNodesCutoff]=loopLength(path,w)
% loopLength.m: Checks for loops as neck cutoffs in the input centerline
% ('path') with a fixed channel width ('w') and returns the centerline 
% with loops removed ('path') and the fraction of nodes in the
% original centerline that belonged to loops ('fractionNodesCutoff').
% Created June 29, 2020 by Ajay B. Limaye, University of Virginia (ajay@virginia.edu).
% Last edited April 8, 2021 by Ajay B. Limaye.

path_in=path;% save in case needed for debugging 
fractionNodesCutoff = 0; % initialize as zero, will be overwritten if cutoffs occur        
[leftBank,rightBank,~] = pathBanks(path,w); % find banks

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
       
        % Calcualate the fraction of nodes in the centerline that belong to
        % a cutoff loop
        fractionNodesCutoff = numel(centerline_nodes_remove)/numel(path.x);
        
        % remove cutoff nodes from the centerline
        path.x(centerline_nodes_remove)=[]; 
        path.y(centerline_nodes_remove)=[]; 
    end
end