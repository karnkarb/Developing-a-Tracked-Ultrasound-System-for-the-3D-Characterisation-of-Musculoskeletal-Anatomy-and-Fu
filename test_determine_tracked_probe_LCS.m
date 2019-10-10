% test function to establish probe LCS ------------------------------------

% test data
marker_coords = [ 1.0  10.0 1.0; ...
                 -2.0 -2.0 1.5; ...
                  5.0  1.0 0.75; ...
                  0.5 -0.5 0.25 ];

% call the function
[ T, R, ori, success_status ] = determine_tracked_probe_LCS( marker_coords )

% check where marker locations in LCS are 
% --> do they end up where we expect them?
local_marker_coords = apply_4by4_tmatrix( T, marker_coords)