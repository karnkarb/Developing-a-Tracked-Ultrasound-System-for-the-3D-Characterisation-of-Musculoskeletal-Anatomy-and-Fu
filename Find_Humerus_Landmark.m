clear all;
close all;
% loading raw data
load ('H:\ResearchProject\calibration data\September 12\Humerus_Landmark\postlat_shaft_t1.mat')

%Filtering data
[~,~,Marker_Stylus,Marker_Stylus_ref,~,~,Marker_Humerus,Marker_Humerus_ref,~,...
    Quality_Infor] = clean_frame_OPA (data,3,12);
filename = 'Stylus_postlat_shaft_t1_12.vtk';

  %Transform to Humerus LCS
for i = 1:size(Marker_Stylus,3)
    %Determine LCS of the Humerus
    Marker_Stylus_tip = [-199.2522447	-28.81881573	-0.0];
    [d,Z(:,:,i),tr{i}] = procrustes(Marker_Humerus_ref,Marker_Humerus(:,:,i)','reflection',false,'scaling',false);   
    [d_st,Z_st(:,:,i),tr_st{i}] = procrustes(Marker_Stylus(:,:,i)',Marker_Stylus_ref,'reflection',false,'scaling',false);
    Stylus_Tip_S(:,:,i) = tr_st{i}.b * Marker_Stylus_tip * tr_st{i}.T + tr_st{i}.c(1,:);
    Stylus_Tip_H(:,:,i) = tr{i}.b * Stylus_Tip_S(:,:,i) * tr{i}.T + tr{i}.c(1,:);
    Tip_polyline(i,:) = Stylus_Tip_H(:,:,i);

end
% write vtk file for visualizing 
polylines(1).vertices = Tip_polyline;
polylines(1).is_closed = 0;
[ success_code ] = write_vtk_multipolyline( filename, polylines )
