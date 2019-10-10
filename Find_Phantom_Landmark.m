clear all;
close all;

%set directory
directory = 'H:\ResearchProject\calibration data\September 12\';
%Landmark point index
for p = 1:16
    %Load and Filter data with OPA
    input_data = sprintf('P%d_t2.mat', p);
    load ([directory input_data])
    output_data = sprintf('P%d_t2_analysed.mat', p);
    
    [~,~,Marker_Stylus_T,Marker_Stylus_ref,...
        Marker_Phantom_T,Marker_Phantom_ref,~,~,...
        ~,Quality_Infor] = clean_frame_OPA (data,1,12);
    
    % Create transformation of ref. config. into tracker frame
    for i = 1:size(Marker_Phantom_T,3)
        [d,Z(:,:,i),tr{i}] = procrustes(Marker_Phantom_T(:,:,i)',Marker_Phantom_ref(:,:),'reflection',false,'scaling',false);
        [d_st,Z_st(:,:,i),tr_st] = procrustes(Marker_Stylus_T(:,:,i)',Marker_Stylus_ref,'reflection',false,'scaling',false);
        local_Marker_Stylus(:,:,i) = tr{i}.b * Marker_Stylus_T(:,:,i)' * tr{i}.T + tr{i}.c(1:size(Marker_Stylus_T,2),:);
        
    end
    %Calculate CoR with SARA
    [ ~,Residual, Scaled_Residual, ...
        ~, Range_Of_Motion,Local_Centres_SARA, ~,~, ...
        Global_Centres_1, Global_Centres_2, Mean_Global_Centres, ...
        ~, ~, ~, ~, ~, ~ ] ...
        = SARA ( permute(Z,[2 1 3]), permute(Z_st,[2 1 3]));
    
    %Transform centre to Phantom LCS
    for j = 1:size(Marker_Stylus_T,3)
        [d_st,Z_st(:,:,j),tr_st] = procrustes(Marker_Stylus_ref,Marker_Stylus_T(:,:,j)','reflection',false,'scaling',false);
        local_tip(:,:,j) = tr_st.b * Global_Centres_1(:,1)' * tr_st.T + tr_st.c(1,:);
        [d,Z_ph(:,:,j),tr_ph{j}] = procrustes(Marker_Phantom_ref(:,:),Marker_Phantom_T(:,:,i)','reflection',false,'scaling',false);
        Center_SCoRE(:,:,j) = tr_ph{j}.b * Global_Centres_1(:,j)' * tr_ph{j}.T + tr_ph{j}.c(1:size(Global_Centres_1(:,j),2),:);
    end
    
    %Calculate CoR with Sphere-Fit
    for n = 1:size(Marker_Phantom_T,3)
        local_marker_sphere1(n,:) = local_Marker_Stylus(1,:,n);
    end
    [Center_LSE(:,:),Radius_LSE(:,:)] = sphereFit(local_marker_sphere1);
    clear data
    %save data for further analysis 
    save(output_data)
    
    clear local_Marker_Phantom local_Marker_Stylus Center_SCoRE Z Z_st Z_ph
    
    
end
