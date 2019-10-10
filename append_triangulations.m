function [ TR ] = append_triangulations( varargin )
%append_triangulations appends a list of triangulation objects into one
%   This function appends a list of triangulation objects and merges them
%   in a single triangulation
    fprintf('Total number of inputs = %d\n',nargin);

    nVarargs = length(varargin);
    fprintf('Inputs in varargin(%d):\n', nVarargs)
    for k = 1:nVarargs
        varargin{k};
    end
    
    % extract point and connectivity from the first triangulation object
    pts = varargin{1}.Points;
    con = varargin{1}.ConnectivityList;
    
    
    % loop over all remaining, supplied triangulations and add them
    for i=2:nVarargs      
        % deal with points first
        start_pt = size(pts,1) + 1;
        end_pt   = start_pt + size(varargin{i}.Points,1) - 1;    
        pts(start_pt:end_pt,:) = varargin{i}.Points;
                      
        % deal with connectivity next
        start_con  = size(con, 1) + 1;
        end_con    = start_con + size(varargin{i}.ConnectivityList,1) - 1;        
        con(start_con:end_con,:) = varargin{i}.ConnectivityList + start_pt-1;        
    end
      
    TR  = triangulation(con, pts);

   
    fprintf('END of Function: append_triangulations\n');
    
end


