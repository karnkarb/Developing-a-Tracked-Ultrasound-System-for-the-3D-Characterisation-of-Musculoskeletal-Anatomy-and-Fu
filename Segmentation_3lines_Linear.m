
%--------------------------------------------------------------------------
% automatic_fiducial_finding.m
%--------------------------------------------------------------------------
% loads US images from probe calibration tests and picks the 6 N-fiducials
% marking the location of the wires in the image.

%--------------------------------------------------------------------------

% get rid of everything currently on the workspace, close all figures ...
clear all;
close all;
warning('off','all');

%Set file directory
%Linear Probe
directory = 'H:\ResearchProject\calibration data\September 17\USscan\';
%Curvilinear Probe
% directory = '\\filestore.soton.ac.uk\users\ww1u18\mydocuments\ResearchProject\calibration data\September 12\Calibration\';

% prepare to load the images;


Set = {'arbitary_t1.mat';'ftb_t1.mat';'ftf_t1.mat';'btb_t1.mat';'arbitary_t2.mat'};

n_tests  = size(Set,2);

XI_total = [];
XI_1 = [];XI_2 = [];XI_3 = [];
XH_total = [];
XH_1 = [];XH_2 = [];XH_3 = [];
XH_X = [];
XH_P = [];
XH_P_1 = [];XH_P_2 = [];XH_P_3 = [];
XH_T = [];
XH_ac = [];
XH_ac_1 = [];XH_ac_2 = [];XH_ac_3 = [];
tr_reg_phantom_col = {};
tr_reg_probe_col = {};
Construct = [];
tmp = 1;

%Markers (Long-Short) in LCS of the Probe
Marker_Probe_P = [48.874,107.0689,109.6146; 48.86750,7.06890,109.6146;...
    90.3675,78.5633,109.6146;15.3675,78.5746,109.6146];

%Markers (Long-Short) in LCS of SolidWorks
Marker_Phantom_Sol = [ -7.8 117.19 -77.58; -7.8 167.17 9.03; -7.8 167.36 -73.64; -7.8 102.41 -36.15];

%Landmarks (#ascending 1-15) in LCS of SolidWorks
% ([5 6 13] 138.21 30 -5;138.21 30 -55;120 34 -5;)
Landmark_Sol = [1.79 30 -55; 1.79 30 -15; 1.79 5 -5; 1.79 5 -55;...
    138.21 10 -5;138.21 5 -55;...
    30 34 -55; 75 34 -55; 120 34 -55; 130 34 -40;...
    55 34 -5; 30 34 -5;20 34 -20];

%Landmark (#ascending 1-16) in LCS of Phantom by SCoRE
Local_Landmark_SC_tmp = load('H:\ResearchProject\calibration data\September 12\Phantom_Landmark\Analysed points\t2\Allpoints_analysed.mat');
Local_Landmark_SC = Local_Landmark_SC_tmp.Center_SCoRE;
clear Local_Landmark_SC_tmp
Local_Landmark_SC(:,[5 6 13]) = [];
% %Find the transformation between actual Phantom and CAD model Phantom
[ d_SS, Z_SS, tr_SS ] ...
    = procrustes ( Local_Landmark_SC',Landmark_Sol,'reflection',false);
[ ~, ~, tr_SS_reg ] ...
    = procrustes ( Landmark_Sol,Local_Landmark_SC','reflection',false);

%
for j = 2
    
    load([directory Set{j,:}])
    
    % load US images data
    [Marker_Probe_T,Marker_Probe_ref,~,~,Marker_Phantom_T,Marker_Phantom_ref,~,~,...
        tmp_USimage,Quality_Infor{j}] = clean_frame_OPA (data,2,17);
    
    clear image_no data
    image_no = 1:size(Marker_Phantom_T,3);
    n_images(j,:) = size(image_no,2)
    USim = tmp_USimage(:,:,image_no(1));
    
    % show the original image
    figure(100);
    imshow(USim, [16 234]);
    h = gca;
    h.Visible = 'On';
    crop(:,:) = [394,312,190,132;377,278,186,126;383,294,179,128;395,331,182,128;403,289,186,132];
    %         crop(j,:) = getrect;
    XH_size = size(XH_total,2);
    
    for i = 1:size(image_no,2)
        
        % load US images
        USim = tmp_USimage(:,:,image_no(i));
        
        %         figure(1)
        %         imshow(USim, [16 234]);
        % show the original image
        figure(101)
        subplot(2,2,1),imshow(USim, [16 234]);
        h = gca;
        h.Visible = 'On';
        
        if size(XH_total,2)-XH_size >= 2
            x = min(coordinates(:,1)) - 15;
            y = min(coordinates(:,2)) - 15;
            width_2 = width + 10;
            height_2 = height + 10;
            img = imcrop(USim, [x y width_2 height_2]);
        else
            % crop image to show only the actual US image
            x = crop(j,1);  % the x coordinate of the top left most corner of the desired image
            y = crop(j,2);  % the y coordinate of the top left most corner of the desired image
            width = crop(j,3);  % the width of the desired image
            height = crop(j,4);  % the height of the desired image
            img = imcrop(USim, [x y width height]);
        end
        
        subplot(2,2,2),imshow(img, [16 234]);
        
        %         % remove background noise and smooth image
        img = noise_cleaning_linear(img, [3 1], 1);
        subplot(2,2,3),imshow(img, [16 234]);
        
        %find optimal thresholding
        
        thre = multithresh(img)/255.0;
        
        %normalise the image
        I = img/255.0;
        
        %binarize the image
        img = im2bw(I, thre);
        
        subplot(2,2,4),imshow(img, [0 1]);
        
        % find the location of the centroids
        
        centroids = get_fiducial_centroids(img);
        
        % find the geometries of the lines
        [slopes, intercepts, fiducial_indices] = get_lines(centroids);
        
        % display the lines on the image
        if(size(slopes, 2) < 3)
            fprintf('Image no.%.0f NOT All 6 fiducials were detected\n',i);
        elseif (size(slopes, 2) > 30)
            fprintf('Image no.%.0f TOO MUCH noise were detected\n',i);
        else
            best_lines = get_parallel_three_lines(slopes, fiducial_indices);
            
            if best_lines == [0, 0, 0]
                fprintf('Image no.%.0f NOT All 6 fiducials were detected\n',i);
            else
                % line 1
                x0 = 0;
                y0 = slopes(best_lines(1)) * x0 + intercepts(best_lines(1));
                x1 = size(img, 2);
                y1 = slopes(best_lines(1)) * x1 + intercepts(best_lines(1));
                
                % line 2
                x2 = 0;
                y2 = slopes(best_lines(2)) * x2 + intercepts(best_lines(2));
                x3 = size(img, 2);
                y3 = slopes(best_lines(2)) * x3 + intercepts(best_lines(2));
                
                % line 3
                x4 = 0;
                y4 = slopes(best_lines(3)) * x4 + intercepts(best_lines(3));
                x5 = size(img, 2);
                y5 = slopes(best_lines(3)) * x5 + intercepts(best_lines(3));
                
                hold on;
                imshow(img);
                plot([x0, x1], [y0, y1],'b', 'LineWidth', 0.5);
                plot([x2, x3], [y2, y3],'r','LineWidth', 0.5);
                plot([x4, x5], [y4, y5],'g','LineWidth', 0.5);
                hold off;
                
                %Define the reference slope for each Set
                if size(XH_total,2)-XH_size == 2  && size(slopes,1)~=0 && best_lines(1)~=0
                    ref_slope(j) = slopes(best_lines(1));
                end
                
                % output the final 6 coordinates found
                [coordinates] = get_three_coordinates(best_lines, fiducial_indices, centroids);
                
                coordinates(:,1) = coordinates(:,1) + x;
                coordinates(:,2) = coordinates(:,2) + y;
                
                % determine which coordinates line is the top row wire
                sort_line = sort([coordinates(1,2), coordinates(4,2), coordinates(7,2)]);
                [sorted_cordinates] = sorting_lines(sort_line,coordinates);
                X1 = sorted_cordinates(1, :);
                X2 = sorted_cordinates(2, :);
                X3 = sorted_cordinates(3, :);
                X4 = sorted_cordinates(4, :);
                X5 = sorted_cordinates(5, :);
                X6 = sorted_cordinates(6, :);
                X7 = sorted_cordinates(7, :);
                X8 = sorted_cordinates(8, :);
                X9 = sorted_cordinates(9, :);
                
                for k = 1:9
                    sorted_cordinates(k, 3) = USim(round(coordinates(k, 2)), round(coordinates(k, 1)));
                end
                
                %                 coordinates
                
                % calculate alpha
                X12 = [X1;X2]; d12 = pdist(X12,'euclidean');
                X13 = [X1;X3]; d13 = pdist(X13,'euclidean');
                
                alpha = d12/d13;
                
                % calculate 2nd alpha
                X45 = [X4;X5]; d45 = pdist(X45,'euclidean');
                X46 = [X4;X6]; d46 = pdist(X46,'euclidean');
                
                alpha2 = d45/d46;
                
                % calculate 3rd alpha
                X78 = [X7;X8]; d78 = pdist(X78,'euclidean');
                X79 = [X7;X9]; d79 = pdist(X79,'euclidean');
                
                alpha3 = d78/d79;
                
                % coordinates from the phantom geometry for top row (1)
                H1 = [68.53, 24.75, -10];
                h1 = [68.53, 24.75, -53.67];
                O1 = [103.53, 24.75, -1.88];
                o1 = [103.53, 24.75, -53.67];
                % coordinates from the phantom geometry for middle row (3)
                H3 = [68.53, 14.75, -10];
                h3 = [68.53, 14.75, -53.67];
                O3 = [103.53, 14.75, -1.88];
                o3 = [103.53, 14.75, -53.67];
                % coordinates from the phantom geometry for bottom row (5)
                H5 = [68.53, 4.75, -6.33];
                h5 = [68.53, 4.75, -50];
                O5 = [103.53, 4.75, -6.33];
                o5 = [103.53, 4.75, -50];
                % Image plane matrix
                XI = [];
                XI(1, 1) = X1(1); XI(2, 1) = X2(1); XI(3, 1) = X3(1);
                XI(4, 1) = X4(1); XI(5, 1) = X5(1); XI(6, 1) = X6(1);
                XI(7, 1) = X7(1); XI(8, 1) = X8(1); XI(9, 1) = X9(1);
                XI(1, 2) = X1(2); XI(2, 2) = X2(2); XI(3, 2) = X3(2);
                XI(4, 2) = X4(2); XI(5, 2) = X5(2); XI(6, 2) = X6(2);
                XI(7, 2) = X7(2); XI(8, 2) = X8(2); XI(9, 2) = X9(2);
                
                XI(:, 3) = 0;
                
                XI = XI';
                
                
                % XH(ri,c1,c2) = A(ri,c1)+alphai(D(ri,c2)-A(ri,c1))
                % Phantom plane matrix
                %% Transform Phantom to Vicon capture
                %Calculate middlde point position using N-wire geometry in
                %SolidWorks frame
                XH1 = h1 +(1-alpha)*(O1-h1);
                XH1 = XH1';
                %Transform to actual Phantom frame
                Sol_XH1 = tr_SS.b * XH1' * tr_SS.T + tr_SS.c(1:size(XH1',1),:);
                %Determine the LCS of Markers in each time frame
                
                %Determine transformation of the markers to Vicon frame
                [ d_Phantom_tmp, Z_Phantom_tmp, tr_Phantom_tmp ] ...
                    = procrustes ( Marker_Phantom_T(:,2:4,image_no(i))', Marker_Phantom_ref(2:4,:),'reflection',false,'scaling',false);
                [ ~, ~, tr_Phantom_reg ] ...
                    = procrustes (  Marker_Phantom_ref(2:4,:),Marker_Phantom_T(:,2:4,image_no(i))','reflection',false,'scaling',false);
                
                XH1_T = tr_Phantom_tmp.b * Sol_XH1 * tr_Phantom_tmp.T + tr_Phantom_tmp.c(1:size(Sol_XH1,1),:);
                
                %%
                %Calculate middlde point position using N-wire geometry in
                %SolidWorks frame
                XH3 = h3 +(1-alpha2)*(O3-h3);
                XH3 = XH3';
                %Transform to actual Phantom frame
                Sol_XH3 = tr_SS.b * XH3' * tr_SS.T+tr_SS.c(1:size(XH3',1),:);
                
                XH3_T = tr_Phantom_tmp.b * Sol_XH3 * tr_Phantom_tmp.T + tr_Phantom_tmp.c(1:size(Sol_XH3,1),:);
                
                %%
                %Calculate middlde point position using N-wire geometry in
                %SolidWorks frame
                XH5 = H5 +(1-alpha3)*(o5-H5);
                XH5 = XH5';
                %Transform to actual Phantom frame
                Sol_XH5 = tr_SS.b * XH5' * tr_SS.T+tr_SS.c(1:size(XH5',1),:);
                
                XH5_T = tr_Phantom_tmp.b * Sol_XH5 * tr_Phantom_tmp.T + tr_Phantom_tmp.c(1:size(Sol_XH5,1),:);
                
                %% Transform Vicon capture to Probe
                
                [ d_Probe_tmp, Z_Probe_tmp, tr_Probe_tmp ] ...
                    = procrustes ( Marker_Probe_ref(2:4,:),Marker_Probe_T(:,2:4,image_no(i))','reflection',false,'scaling',false);
                [ ~, ~, tr_Probe_reg ]...
                    = procrustes ( Marker_Probe_T(:,2:4,image_no(i))',Marker_Probe_ref(2:4,:),'reflection',false,'scaling',false);
                XH1_P = tr_Probe_tmp.b * XH1_T * tr_Probe_tmp.T + tr_Probe_tmp.c(1:size(XH1_T,1),:);
                XH3_P = tr_Probe_tmp.b * XH3_T * tr_Probe_tmp.T + tr_Probe_tmp.c(1:size(XH3_T,1),:);
                XH5_P = tr_Probe_tmp.b * XH5_T * tr_Probe_tmp.T + tr_Probe_tmp.c(1:size(XH5_T,1),:);
                
                %Thresholding locations of N-fiducials
                if size(XH_total,2)-XH_size >= 2
                    thresh1 = norm (XH1(:) - XH_total(:,size(XH_total,2)-2)); %location of XH1
                    
                    thresh2 = norm (XH3(:) - XH_total(:,size(XH_total,2)-1)); %location of XH3
                    
                    thresh3 = norm (XH5(:) - XH_total(:,size(XH_total,2))); %location of XH3
                    
                    thresh4 = norm (XI(:,2) - XI(:,5)); %Distance between parallel lines
                    thresh5 = norm (XI(:,5) - XI(:,8)); %Distance between parallel lines
                    
                    if thresh1>10
                        cprintf ('Magenta','Image no.%.0f miss-locating XH1, %.2f\n',i,thresh1);
                        XI_total = XI_total;
                        XH_total = XH_total;
                    elseif thresh2>10
                        cprintf ('Magenta','Image no.%.0f miss-locating XH3, %.2f\n',i,thresh2);
                        XI_total = XI_total;
                        XH_total = XH_total;
                    elseif thresh3>10
                        cprintf ('Magenta','Image no.%.0f miss-locating XH5, %.2f\n',i,thresh3);
                        XI_total = XI_total;
                        XH_total = XH_total;
                    elseif thresh4>45
                        cprintf ('Magenta','Image no.%.0f miss-locating parallel lines, %.2f\n',i,thresh3);
                        XI_total = XI_total;
                        XH_total = XH_total;
                    else
                        XI_total = [XI_total XI(:,2) XI(:,5) XI(:,8)];
                        XI_1 = [XI_1 XI(:,2) XI(:,5)];
                        XI_2 = [XI_2 XI(:,2) XI(:,8)];
                        XI_3 = [XI_3 XI(:,5) XI(:,8)];
                        XH_total = [XH_total XH1 XH3 XH5];
                        XH_1 = [XH_1 XH1 XH3];
                        XH_2 = [XH_2 XH1 XH5];
                        XH_3 = [XH_3 XH3 XH5];
                        XH_ac = [XH_ac; Sol_XH1; Sol_XH3; Sol_XH5];
                        XH_ac_1 = [XH_ac_1; Sol_XH1; Sol_XH3];
                        XH_ac_2 = [XH_ac_2; Sol_XH1; Sol_XH5];
                        XH_ac_3 = [XH_ac_3; Sol_XH3; Sol_XH5];
                        XH_X = [XH_X; X2; X5; X8];
                        XH_P = [XH_P; XH1_P; XH3_P; XH5_P ];
                        XH_P_1 = [XH_P_1; XH1_P; XH3_P ];
                        XH_P_2 = [XH_P_2; XH1_P; XH5_P ];
                        XH_P_3 = [XH_P_3; XH3_P; XH5_P ];
                        XH_T = [XH_T; XH1_T; XH3_T; XH5_T];
                        Construct = [Construct XI(:,1) XI(:,3) XI(:,4) XI(:,6) XI(:,7) XI(:,9)];
                        tr_reg_phantom_col{tmp} = tr_Phantom_reg;
                        tr_reg_probe_col{tmp} = tr_Probe_reg;
                        %%
                        [d_cross_wire(j,i),Z_cross_wire_plot,tr_cross_wire_plot{tmp}] = ...
                            procrustes(XH_P,XI_total');
                        diff_plot = Z_cross_wire_plot - XH_P;
                        error_norm_plot{j,i} = vecnorm(diff_plot,2,2);
                        Mean_error_plot(j,i) = mean(error_norm_plot{j,i});
                        display(Mean_error_plot(j,i))
                        Max_error_plot(j,i) = max(error_norm_plot{j,i});
                        SD_error_plot(j,i) = std(error_norm_plot{j,i});
                        
                        [d_cross_wire_1(j,i),Z_cross_wire_plot_1,tr_cross_wire_plot_1] = ...
                            procrustes(XH_P_1,XI_1');
                        diff_plot_1 = Z_cross_wire_plot_1 - XH_P_1;
                        error_norm_plot_1{j,i} = vecnorm(diff_plot_1,2,2);
                        Mean_error_plot_1(j,i) = mean(error_norm_plot_1{j,i});
                        
                        Max_error_plot_1(j,i) = max(error_norm_plot_1{j,i});
                        SD_error_plot_1(j,i) = std(error_norm_plot_1{j,i});
                        
                        [d_cross_wire_2(j,i),Z_cross_wire_plot_2,tr_cross_wire_plot_2] = ...
                            procrustes(XH_P_2,XI_2');
                        diff_plot_2 = Z_cross_wire_plot_2 - XH_P_2;
                        error_norm_plot_2{j,i} = vecnorm(diff_plot_2,2,2);
                        Mean_error_plot_2(j,i) = mean(error_norm_plot_2{j,i});
                        
                        Max_error_plot_2(j,i) = max(error_norm_plot_2{j,i});
                        SD_error_plot_2(j,i) = std(error_norm_plot_2{j,i});
                        
                        [d_cross_wire_3(j,i),Z_cross_wire_plot_3,tr_cross_wire_plot_3] = ...
                            procrustes(XH_P_3,XI_3');
                        diff_plot_3 = Z_cross_wire_plot_3 - XH_P_3;
                        error_norm_plot_3{j,i} = vecnorm(diff_plot_3,2,2);
                        Mean_error_plot_3(j,i) = mean(error_norm_plot_3{j,i});
                        
                        Max_error_plot_3(j,i) = max(error_norm_plot_3{j,i});
                        SD_error_plot_3(j,i) = std(error_norm_plot_3{j,i});
                        
                        %% Calculate FRE
                        [FRE_mean(tmp),FRE_max(tmp),FRE_SD(tmp)] = ...
                            FRE_Feedback(Z_cross_wire_plot,XH_ac,tr_reg_probe_col,tr_reg_phantom_col);
                        [FRE_mean_1(tmp),FRE_max_1(tmp),FRE_SD_1(tmp)] = ...
                            FRE_Feedback_2lines(Z_cross_wire_plot_1, XH_ac_1, tr_reg_probe_col, tr_reg_phantom_col);
                        [FRE_mean_2(tmp),FRE_max_2(tmp),FRE_SD_2(tmp)] = ...
                            FRE_Feedback_2lines(Z_cross_wire_plot_2, XH_ac_2, tr_reg_probe_col, tr_reg_phantom_col);
                        [FRE_mean_3(tmp),FRE_max_3(tmp),FRE_SD_3(tmp)] = ...
                            FRE_Feedback_2lines(Z_cross_wire_plot_3, XH_ac_3, tr_reg_probe_col, tr_reg_phantom_col);
                        
                        tmp = tmp+1;
                        
                    end
                else
                    XI_total = [XI_total XI(:,2) XI(:,5) XI(:,8)];
                    XI_1 = [XI_1 XI(:,2) XI(:,5)];
                    XI_2 = [XI_2 XI(:,2) XI(:,8)];
                    XI_3 = [XI_3 XI(:,5) XI(:,8)];
                    XH_total = [XH_total XH1 XH3 XH5];
                    XH_1 = [XH_1 XH1 XH3];
                    XH_2 = [XH_2 XH1 XH5];
                    XH_3 = [XH_3 XH3 XH5];
                    XH_ac = [XH_ac; Sol_XH1; Sol_XH3; Sol_XH5];
                    XH_ac_1 = [XH_ac_1; Sol_XH1; Sol_XH3];
                    XH_ac_2 = [XH_ac_2; Sol_XH1; Sol_XH5];
                    XH_ac_3 = [XH_ac_3; Sol_XH3; Sol_XH5];
                    XH_X = [XH_X; X2; X5; X8];
                    XH_P = [XH_P; XH1_P; XH3_P; XH5_P ];
                    XH_P_1 = [XH_P_1; XH1_P; XH3_P ];
                    XH_P_2 = [XH_P_2; XH1_P; XH5_P ];
                    XH_P_3 = [XH_P_3; XH3_P; XH5_P ];
                    XH_T = [XH_T; XH1_T; XH3_T; XH5_T];
                    Construct = [Construct XI(:,1) XI(:,3) XI(:,4) XI(:,6) XI(:,7) XI(:,9)];
                    tr_reg_phantom_col{tmp} = tr_Phantom_reg;
                    tr_reg_probe_col{tmp} = tr_Probe_reg;
                    %%
                    [d_cross_wire(j,i),Z_cross_wire_plot,tr_cross_wire_plot{tmp}] = ...
                        procrustes(XH_P,XI_total');
                    diff_plot = Z_cross_wire_plot - XH_P;
                    error_norm_plot{j,i} = vecnorm(diff_plot,2,2);
                    Mean_error_plot(j,i) = mean(error_norm_plot{j,i});
                    display(Mean_error_plot(j,i));
                    Max_error_plot(j,i) = max(error_norm_plot{j,i});
                    SD_error_plot(j,i) = std(error_norm_plot{j,i});
                    
                    [d_cross_wire_1(j,i),Z_cross_wire_plot_1,tr_cross_wire_plot_1] = ...
                        procrustes(XH_P_1,XI_1');
                    diff_plot_1 = Z_cross_wire_plot_1 - XH_P_1;
                    error_norm_plot_1{j,i} = vecnorm(diff_plot_1,2,2);
                    Mean_error_plot_1(j,i) = mean(error_norm_plot_1{j,i});
                    
                    Max_error_plot_1(j,i) = max(error_norm_plot_1{j,i});
                    SD_error_plot_1(j,i) = std(error_norm_plot_1{j,i});
                    
                    [d_cross_wire_2(j,i),Z_cross_wire_plot_2,tr_cross_wire_plot_2] = ...
                        procrustes(XH_P_2,XI_2');
                    diff_plot_2 = Z_cross_wire_plot_2 - XH_P_2;
                    error_norm_plot_2{j,i} = vecnorm(diff_plot_2,2,2);
                    Mean_error_plot_2(j,i) = mean(error_norm_plot_2{j,i});
                    
                    Max_error_plot_2(j,i) = max(error_norm_plot_2{j,i});
                    SD_error_plot_2(j,i) = std(error_norm_plot_2{j,i});
                    
                    [d_cross_wire_3(j,i),Z_cross_wire_plot_3,tr_cross_wire_plot_3] = ...
                        procrustes(XH_P_3,XI_3');
                    diff_plot_3 = Z_cross_wire_plot_3 - XH_P_3;
                    error_norm_plot_3{j,i} = vecnorm(diff_plot_3,2,2);
                    Mean_error_plot_3(j,i) = mean(error_norm_plot_3{j,i});                    
                    Max_error_plot_3(j,i) = max(error_norm_plot_3{j,i});
                    SD_error_plot_3(j,i) = std(error_norm_plot_3{j,i});
                    
                    
                    %% Calculate FRE
                    [FRE_mean(tmp),FRE_max(tmp),FRE_SD(tmp)] = ...
                        FRE_Feedback(Z_cross_wire_plot,XH_ac,tr_reg_probe_col,tr_reg_phantom_col);
                    [FRE_mean_1(tmp),FRE_max_1(tmp),FRE_SD_1(tmp)] = ...
                        FRE_Feedback_2lines(Z_cross_wire_plot_1, XH_ac_1, tr_reg_probe_col, tr_reg_phantom_col);
                    [FRE_mean_2(tmp),FRE_max_2(tmp),FRE_SD_2(tmp)] = ...
                        FRE_Feedback_2lines(Z_cross_wire_plot_2, XH_ac_2, tr_reg_probe_col, tr_reg_phantom_col);
                    [FRE_mean_3(tmp),FRE_max_3(tmp),FRE_SD_3(tmp)] = ...
                        FRE_Feedback_2lines(Z_cross_wire_plot_3, XH_ac_3, tr_reg_probe_col, tr_reg_phantom_col);
                                      
                    tmp = tmp+1;
                end
                
            end
            
        end
    end
    
    n_images(j,2) = size(XH_total,2)
    cprintf('Red','Select new Set %d\n',j+1);
    %     clear Marker_Probe_T Marker_Stylus_T Marker_Phantom_T Marker_Humerus_T tmp_USimage;
end

%Create Combinations of Double-N from Triple-N

XI_line1 = XI_total(:,1:3:end);
XI_line2 = XI_total(:,2:3:end);
XI_line3 = XI_total(:,3:3:end);
XI_Set1 = [XI_line1 XI_line2];
XI_Set2 = [XI_line1 XI_line3];
XI_Set3 = [XI_line2 XI_line3];
XH_line1 = XH_P(1:3:end,:);
XH_line2 = XH_P(2:3:end,:);
XH_line3 = XH_P(3:3:end,:);
XH_Set1 = [XH_line1; XH_line2];
XH_Set2 = [XH_line1; XH_line3];
XH_Set3 = [XH_line2; XH_line3];

%Procruste to transform points in Y to X
[d_cross_wire_total,Z_cross_wire_total,tr_cross_wire_total] = ...
    procrustes(XH_P,XI_total');
diff_total = Z_cross_wire_total - XH_P;
error_norm_total = vecnorm(diff_total,2,2);
Mean_error_total = mean(error_norm_total);
Max_error_total = max(error_norm_total);
SD_error_total = std(error_norm_total);
% Calculate Double-N 1st Comb.
[d_cross_wire_comb{1},Z_cross_wire_comb{1},tr_cross_wire_comb{1}] = ...
    procrustes(XH_Set1,XI_Set1');
diff_comb{1} = Z_cross_wire_comb{1} - XH_Set1;
error_norm_comb{1} = vecnorm(diff_comb{1},2,2);
Mean_error_comb{1} = mean(error_norm_comb{1});
Max_error_comb{1} = max(error_norm_comb{1});
SD_error_comb{1} = std(error_norm_comb{1});

% Calculate Double-N 2nd Comb.
[d_cross_wire_comb{2},Z_cross_wire_comb{2},tr_cross_wire_comb{2}] = ...
    procrustes(XH_Set2,XI_Set2');
diff_comb{2} = Z_cross_wire_comb{2} - XH_Set2;
error_norm_comb{2} = vecnorm(diff_comb{2},2,2);
Mean_error_comb{2} = mean(error_norm_comb{2});
Max_error_comb{2} = max(error_norm_comb{2});
SD_error_comb{2} = std(error_norm_comb{2});
% Calculate Double-N 3rd Comb.
[d_cross_wire_comb{3},Z_cross_wire_comb{3},tr_cross_wire_comb{3}] = ...
    procrustes(XH_Set3,XI_Set3');
diff_comb{3} = Z_cross_wire_comb{3} - XH_Set3;
error_norm_comb{3} = vecnorm(diff_comb{3},2,2);
Mean_error_comb{3} = mean(error_norm_comb{3});
Max_error_comb{3} = max(error_norm_comb{3});
SD_error_comb{3} = std(error_norm_comb{3});


% [d_cross_wire,Z_cross_wire,tr_cross_wire] = procrustes(XH_P,XI_total','reflection',false);
% B_ang = asin(tr_cross_wire.T(1,3));
% A_ang = acos(tr_cross_wire.T(3,3)/cos(B_ang));
% C_ang = acos(tr_cross_wire.T(1,1)/cos(B_ang));
% diff = Z_cross_wire - XH_P;

%Preserve empty cell array
n = size(XI_total,2);
XI_random = {}; XH_random = {};
d_cross_wire = {}; Z_cross_wire = {}; tr_cross_wire = {};
diff = {}; error_norm = {};

%%
% Eliminate points with error_norm higher than 10
eliminate = find(error_norm_total> 10);
XH_eli(:,:) = XH_P;
XH_eli(eliminate,:) = [];
XI_eli(:,:) = XI_total;
XI_eli(:,eliminate) = [];
[d_cross_wire_eli,Z_cross_wire_eli,tr_cross_wire_eli] = ...
    procrustes(XH_eli,XI_eli');
diff_eli = Z_cross_wire_eli - XH_eli;
error_norm_eli = vecnorm(diff_eli,2,2);
Mean_error_eli = mean(error_norm_eli);
Max_error_eli = max(error_norm_eli);
SD_error_eli = std(error_norm_eli);

figure();
subplot(1,3,1);scatter3 (XI_total(1,1:3:end),XI_total(2,1:3:end),XI_total(3,1:3:end))
hold on; scatter3 (XI_total(1,2:3:end),XI_total(2,2:3:end),XI_total(3,2:3:end))
hold on; scatter3 (XI_total(1,3:3:end),XI_total(2,3:3:end),XI_total(3,3:3:end))
title ('X_I');xlabel('x'); ylabel('y'); zlabel('z')
axis equal
subplot(1,3,2); scatter3 (XH_total(1,1:3:end),XH_total(2,1:3:end),XH_total(3,1:3:end))
hold on; scatter3 (XH_total(1,2:3:end),XH_total(2,2:3:end),XH_total(3,2:3:end))
hold on; scatter3 (XH_total(1,3:3:end),XH_total(2,3:3:end),XH_total(3,3:3:end))
title ('X_H');xlabel('x'); ylabel('y'); zlabel('z')
axis equal
XH_P = XH_P';
subplot(1,3,3);scatter3 (XH_P(1,1:3:end),XH_P(2,1:3:end),XH_P(3,1:3:end))
hold on; scatter3 (XH_P(1,2:3:end),XH_P(2,2:3:end),XH_P(3,2:3:end))
hold on; scatter3 (XH_P(1,3:3:end),XH_P(2,3:3:end),XH_P(3,3:3:end))
title ('X_P');xlabel('x'); ylabel('y'); zlabel('z')
scatter3(Marker_Probe_ref(:,1),Marker_Probe_ref(:,2),Marker_Probe_ref(:,3),'filled')
axis equal
clear Local_Landmark_SC_tmp USim error_norm_plot error_norm_plot_1 error_norm_plot_2 error_norm_plot_3 tmp_USimage data

tmp = 1;
for m =1:size(Set,1)
    for n = 1:size(Mean_error_plot,2)
        if Mean_error_plot(m,n)~=0
            conv(tmp) = Mean_error_plot(m,n);
            conv_1(tmp) = Mean_error_plot_1(m,n);
            conv_2(tmp) = Mean_error_plot_2(m,n);
            conv_3(tmp) = Mean_error_plot_3(m,n);
            
            tmp = tmp+1;
        end
    end
end
n(1) = n_images(1,2);
n(2) = n_images(2,2)-n_images(1,2);
n(3) = n_images(3,2)-n_images(2,2);
n(4) = n_images(4,2)-n_images(3,2);
n(5) = n_images(5,2)-n_images(4,2);

figure();hold on;plot(conv);plot(conv_1);plot(conv_2);plot(conv_3);
title ('');xlabel('number of input image'); ylabel('|X_P-^P_UT*X_U| [mm]');
yL = ylim;
plot([n(1) n(1)], yL,'k--');
plot([n(2) n(2)], yL,'k--');
plot([n(3) n(3)], yL,'k--');
plot([n(4) n(4)], yL,'k--');
plot([n(5) n(5)], yL,'k--');
legend ('Triple-N','Comb1','Comb2','Comb3')

figure();hold on;
plot(FRE_mean);plot(FRE_mean_1);plot(FRE_mean_2);plot(FRE_mean_3)
plot([n(1) n(1)], yL,'k--');
plot([n(2) n(2)], yL,'k--');
plot([n(3) n(3)], yL,'k--');
plot([n(4) n(4)], yL,'k--');
plot([n(5) n(5)], yL,'k--');
title ('');xlabel('number of input image'); ylabel('FRE [mm]');
legend ('Triple-N','Comb1','Comb2','Comb3')
mean(abs(diff_plot),1)

dim1 = size(XI_total,2);
dim2 = dim1*2/3;

try_box = [error_norm_total;error_norm_comb{1,1};error_norm_comb{1,2};error_norm_comb{1,3}];
grp1 = repmat(['Triple-N'],dim1,1);
grp2 = repmat(['Comb1   '],dim2,1);
grp3 = repmat(['Comb2   '],dim2,1);
grp4 = repmat(['Comb3   '],dim2,1);
grp = [grp1;grp2;grp3;grp4];
figure()
boxplot(try_box,grp)
hold on;
h = findobj(gca,'Tag','Median');
set(h,'Visible','off');
scatter([1;2;3;4],[Mean_error_total;Mean_error_comb{1};Mean_error_comb{2};Mean_error_comb{3}],'filled')
ylabel('FRE [mm]');