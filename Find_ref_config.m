
%Filtering data

[Marker_Probe_T,Marker_Stylus_T,Marker_Phantom_T,Marker_Humerus_T,~,...
        ~] = clean_frame(data,5,14);
    
for i = 1:size(Marker_Humerus_T,3)
        [ Thum(:,:,i), Rhum(:,:,i), orihum(:,:,i), success_status ] = determine_Humerus_LCS(Marker_Humerus_T(:,:,i)');
            Marker_Humerus_tmp(:,:,i) = apply_4by4_tmatrix( Thum(:,:,i),Marker_Humerus_T(:,:,i)');
        [ Tpro(:,:,i), Rpro(:,:,i), oripro(:,:,i), success_status ] = determine_tracked_probe_LCS(Marker_Probe_T(:,:,i)');
            Marker_Probe_tmp(:,:,i) = apply_4by4_tmatrix( Tpro(:,:,i),Marker_Probe_T(:,:,i)');
        [ Tpha(:,:,i), Rpha(:,:,i), oripha(:,:,i), success_status ] = determine_tracked_probe_LCS(Marker_Phantom_T(:,:,i)');
            Marker_Phantom_tmp(:,:,i) = apply_4by4_tmatrix( Tpha(:,:,i),Marker_Phantom_T(:,:,i)');
        [ Tsty(:,:,i), Rsty(:,:,i), oristy(:,:,i), success_status ] = determine_Humerus_LCS(Marker_Stylus_T(:,:,i)');
            Marker_Stylus_tmp(:,:,i) = apply_4by4_tmatrix( Tsty(:,:,i),Marker_Stylus_T(:,:,i)');
end