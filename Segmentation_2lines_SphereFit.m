
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

%set working directory for Double-N with Curvilinear probe
directory = 'H:\ResearchProject\calibration data\August 15\Calibration\';

% prepare to load the images; 'Tilting_part1.mat'
Set = { 'Front_to_Back_1.mat'; 'Back_to_Front_1.mat';'m_static_1.mat'; 'f_static_1.mat'; 'b_static_1.mat'; ...
    'Upper_layer_1.mat'; 'Lower_layer_1.mat';'Tilting_part2.mat'  };

n_tests  = size(set,2);

XI_total = [];

XH_total = [];
XH_X = [];
XH_P = [];
XH_T = [];
XH_ac = [];
Construct = [];
tr_reg_phantom_col = {};
tr_reg_probe_col = {};
tmp = 1;

%Markers (Long-Short) in LCS of the Probe
Marker_Probe_P = [48.874,107.0689,109.6146; 48.86750,7.06890,109.6146;...
    90.3675,78.5633,109.6146;15.3675,78.5746,109.6146];

%Markers (Long-Short) in LCS of SolidWorks
Marker_Phantom_Sol = [ -7.8 117.19 -77.58; -7.8 167.17 9.03; -7.8 167.36 -73.64; -7.8 102.41 -36.15];

%Landmarks (#ascending 1-16) in LCS of SolidWorks
%eliminate point 5 7 16 [138.21 30 -5;138.21 10 -5;; 20 34 -20];
Landmark_Sol = [1.79 30 -55; 1.79 30 -15; 1.79 5 -5; 1.79 5 -55;...
    138.21 30 -55; 138.21 5 -55;...
    30 34 -55; 75 34 -55; 120 34 -55; 130 34 -40;...
    120 34 -5; 55 34 -5; 30 34 -5];

%Landmark (#ascending 1-16) in LCS of Phantom by SphereFit
Local_Landmark_SC_tmp = load('H:\ResearchProject\calibration data\August 15\Phantom_Landmark\Analyzed\allpoints.mat');
Local_Landmark_SC = Local_Landmark_SC_tmp.Center_SphereFit;
clear Local_Landmark_SC_tmp
Local_Landmark_SC(:,[5 7 16]) = [];
% %Find the transformation between actual Phantom and CAD model Phantom
[ d_SS, Z_SS, tr_SS ] ...
    = procrustes ( Local_Landmark_SC',Landmark_Sol,'reflection',false);
%1:

for j = 1:size(set,1)
    
    load([directory set{j,:}])
    
    % load US images data
    [Marker_Probe_T,Marker_Probe_ref,~,~,Marker_Phantom_T,Marker_Phantom_ref,~,~,...
        tmp_USimage,Quality_Infor{j}] = clean_frame_OPA (data,2,17);
    
    clear image_no
    if j == 1
        image_no = 1:140;
    else
        image_no = 1:size(Marker_Phantom_T,3);
    end
    n_images(j,:) = size(image_no,2)
    USim = tmp_USimage(:,:,image_no(1));
    
    % show the original image
    figure(100);
    imshow(USim, [16 234]);
    h = gca;
    h.Visible = 'On';
    crop = [392,374,224,111;392,374,224,111;324,432,224,110;415,363,217,122;459,432,237,131;371,374,251,120;337,320,234,114;491,365,219,130;290,347,229,146];
    
    XH_size = size(XH_total,2);
    
    for i = 1:size(image_no,2)
        
        % load US images
        USim = tmp_USimage(:,:,image_no(i));
        
        
        % show the original image
        figure(101)
        subplot(2,2,1),imshow(USim, [16 234]);
        
        if size(XH_total,2)-XH_size >= 2
            x = min(coordinates(:,1)) - 30;
            y = min(coordinates(:,2)) - 20;
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
        img = noise_cleaning(img, [3 1], 1);
        subplot(2,2,3),imshow(img, [16 234]);
        
        %find optimal thresholding
        
        thre = multithresh(img)/255.0;
        
        %normalise the image
        I = img/255.0;
        
        %binarize the image
        img = im2bw(I, thre-0.05);
        
        subplot(2,2,4),imshow(img, [0 1]);
        
        % find the location of the centroids
        centroids = get_fiducial_centroids(img);
        
        % find the geometries of the lines
        [slopes, intercepts, fiducial_indices] = get_lines(centroids);
        
        % display the lines on the image
        if(size(slopes, 2) < 2)
            fprintf('Image no.%.0f NOT All 6 fiducials were detected\n',i);
        elseif (size(slopes, 2) > 20)
            fprintf('Image no.%.0f TOO MUCH noise were detected\n',i);
        else
            best_lines = get_parallel_lines(slopes, fiducial_indices);
            
            if best_lines == [0, 0]
                fprintf('Image no.%.0f All 6 fiducials were not detected\n',i);
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
                
                hold on;
                imshow(img);
                plot([x0, x1], [y0, y1],'b', 'LineWidth', 0.5);
                plot([x2, x3], [y2, y3],'r','LineWidth', 0.5);
                hold off;
                
                %Define the reference slope for each set
                if size(XH_total,2)-XH_size == 2  && size(slopes,1)~=0 && best_lines(1)~=0
                    ref_slope(j) = slopes(best_lines(1));
                end
                
                % output the final 6 coordinates found
                [coordinates] = get_coordinates(best_lines, fiducial_indices, centroids);
                
                coordinates(:,1) = coordinates(:,1) + x;
                coordinates(:,2) = coordinates(:,2) + y;
                
                % determine which coordinates line is the top row wire
                topline = min(coordinates(1,2), coordinates(4,2));
                
                % label N-fiducials according to chen et als ordering
                if topline == coordinates(1,2)
                    X3 = coordinates(1, :);
                    X2 = coordinates(2, :);
                    X1 = coordinates(3, :);
                    X6 = coordinates(4, :);
                    X5 = coordinates(5, :);
                    X4 = coordinates(6, :);
                else
                    X3 = coordinates(4, :);
                    X2 = coordinates(5, :);
                    X1 = coordinates(6, :);
                    X6 = coordinates(1, :);
                    X5 = coordinates(2, :);
                    X4 = coordinates(3, :);
                end
                
                coordinates(1, :) = X3;
                coordinates(2, :) = X2;
                coordinates(3, :) = X1;
                coordinates(4, :) = X6;
                coordinates(5, :) = X5;
                coordinates(6, :) = X4;
                
                for k = 1:6
                    coordinates(k, 3) = USim(round(coordinates(k, 2)), round(coordinates(k, 1)));
                end
                
                %                 coordinates
                
                % calculate alpha
                X12 = [X1;X2];d12 = pdist(X12,'euclidean');
                
                X13 = [X1;X3];d13 = pdist(X13,'euclidean');
                
                alpha = d12/d13;
                
                % calculate 2nd alpha
                X45 = [X4;X5];
                d45 = pdist(X45,'euclidean');
                
                X46 = [X4;X6];
                d46 = pdist(X46,'euclidean');
                
                alpha2 = d45/d46;
                
                % coordinates from the phantom geometry for top row (2)
                H1 = [68.53, 24.75, -1.88];
                h1 = [68.53, 24.75, -58.74];
                O1 = [103.53, 24.75, -1.88];
                o1 = [103.53, 24.75, -58.74];
                % coordinates from the phantom geometry for bottom row (5)
                H4 = [68.53, 9.75, -1.88];
                h4 = [68.53, 9.75, -58.74];
                O4 = [103.53, 9.75, -1.88];
                o4 = [103.53, 9.75, -58.74];
                
                % Image plane matrix
                XI = [];
                XI(1, 1) = X1(1); XI(2, 1) = X2(1); XI(3, 1) = X3(1);
                XI(4, 1) = X4(1); XI(5, 1) = X5(1); XI(6, 1) = X6(1);
                XI(1, 2) = X1(2); XI(2, 2) = X2(2); XI(3, 2) = X3(2);
                XI(4, 2) = X4(2); XI(5, 2) = X5(2); XI(6, 2) = X6(2);
                
                XI(:, 3) = 0;
                XI = XI';
%                 subplot(3,2,6)
%                 hold on;scatter(X2(1),X2(2),'filled')
%                 hold on;scatter(X5(1),X5(2),'filled')
%                 hold off;
%                 
                % XH(ri,c1,c2) = A(ri,c1)+alphai(D(ri,c2)-A(ri,c1))
                % Phantom plane matrix
                %% Transform Phantom to Vicon capture
                %Calculate middlde point position using N-wire geometry in
                %SolidWorks frame
                XH1 = h1 +alpha*(O1-h1);
                XH1 =XH1';
                %Transform to actual Phantom frame
                Act_XH1 = tr_SS.b * XH1' * tr_SS.T + tr_SS.c(1:size(XH1,2),:);
                
                %Determine transformation of the markers to Vicon frame
                [ d_Phantom_tmp, Z_Phantom_tmp, tr_Phantom_tmp ] ...
                    = procrustes ( Marker_Phantom_T(:,:,image_no(i))', Marker_Phantom_ref(:,:),'reflection',false,'scaling',false);
                [ ~, ~, tr_Phantom_reg ] ...
                    = procrustes (  Marker_Phantom_ref(:,:),Marker_Phantom_T(:,:,image_no(i))','reflection',false,'scaling',false);
                
                XH1_T = tr_Phantom_tmp.b * Act_XH1 * tr_Phantom_tmp.T + tr_Phantom_tmp.c(1:size(Act_XH1,1),:);
                
                %%
                %Calculate middlde point position using N-wire geometry in
                %SolidWorks frame
                XH4 = H4 +alpha2*(o4-H4);
                XH4 = XH4';
                %Transform to actual Phantom frame
                Act_XH4 = tr_SS.b * XH4' * tr_SS.T+tr_SS.c(1:size(XH4,2),:);
                
                XH4_T = tr_Phantom_tmp.b * Act_XH4 * tr_Phantom_tmp.T + tr_Phantom_tmp.c(1:size(Act_XH4,1),:);
                %% Transform Vicon capture to Probe
                [ d_Probe_tmp, Z_Probe_tmp, tr_Probe_tmp ] ...
                    = procrustes ( Marker_Probe_ref(:,:),Marker_Probe_T(:,:,image_no(i))','reflection',false,'scaling',false);
                [ ~, ~, tr_Probe_reg ]...
                    = procrustes ( Marker_Probe_T(:,:,image_no(i))',Marker_Probe_ref(:,:),'reflection',false,'scaling',false);
                
                XH1_P = tr_Probe_tmp.b * XH1_T * tr_Probe_tmp.T + tr_Probe_tmp.c(1:size(XH1_T,1),:);
                XH4_P = tr_Probe_tmp.b * XH4_T * tr_Probe_tmp.T + tr_Probe_tmp.c(1:size(XH4_T,1),:);
                
                %Thresholding locations of N-fiducials
                if size(XH_total,2)-XH_size >= 2
                    thresh1 = norm (XH1(:) - XH_total(:,size(XH_total,2)-1)); %location of XH1
                    thresh2 = norm (XH4(:) - XH_total(:,size(XH_total,2))); %location of XH4
                    thresh3 = norm (XI(:,2) - XI(:,5)); %Distance between parallel lines
                    
                    if thresh1>15
                        cprintf ('Magenta','Image no.%.0f miss-locating XH1, %.2f\n',i,thresh1);
                        XI_total = XI_total;
                        XH_total = XH_total;
                    elseif thresh2>15
                        cprintf ('Magenta','Image no.%.0f miss-locating XH4, %.2f\n',i,thresh2);
                        XI_total = XI_total;
                        XH_total = XH_total;
                    elseif ((thresh3<65) || (abs(slopes(best_lines(1))-ref_slope(j))>0.5))
                        cprintf ('Magenta','Image no.%.0f miss-locating lines, %.2f\n',i,thresh3);
                        XI_total = XI_total;
                        XH_total = XH_total;
                    else
                        XI_total = [XI_total XI(:,2) XI(:,5)];
                        XH_total = [XH_total XH1 XH4];
                        XH_X = [XH_X; X2; X5;];
                        XH_ac = [XH_ac; Act_XH1; Act_XH4];                        
                        XH_P = [XH_P; XH1_P; XH4_P];
                        XH_T = [XH_T; XH1_T; XH4_T];
                        Construct = [Construct XI(:,1) XI(:,3) XI(:,4) XI(:,6)];
                        tr_reg_phantom_col{tmp} = tr_Phantom_reg;
                        tr_reg_probe_col{tmp} = tr_Probe_reg;
                        [d_cross_wire(j,i),Z_cross_wire_plot,tr_cross_wire_plot{tmp}] = ...
                            procrustes(XH_P,XI_total');
                        diff_plot = Z_cross_wire_plot - XH_P;
                        error_norm_plot{j,i} = vecnorm(diff_plot,2,2);
                        Mean_error_plot(j,i) = mean(error_norm_plot{j,i});
                        display(Mean_error_plot(j,i))
                        Max_error_plot(j,i) = max(error_norm_plot{j,i});
                        SD_error_plot(j,i) = std(error_norm_plot{j,i});
                        
                        [FRE_mean(tmp),FRE_max(tmp),FRE_SD(tmp)] = ...
                            FRE_Feedback_2lines(Z_cross_wire_plot,XH_ac,tr_reg_probe_col,tr_reg_phantom_col);
                        
                        tmp = tmp+1
                    end
                else
                    XI_total = [XI_total XI(:,2) XI(:,5)];
                    XH_total = [XH_total XH1 XH4];
                    XH_X = [XH_X; X2; X5;];
                    XH_ac = [XH_ac; Act_XH1; Act_XH4];
                    XH_P = [XH_P; XH1_P; XH4_P];
                    XH_T = [XH_T; XH1_T; XH4_T];
                    Construct = [Construct XI(:,1) XI(:,3) XI(:,4) XI(:,6)];
                    
                    tr_reg_phantom_col{tmp} = tr_Phantom_reg;
                    tr_reg_probe_col{tmp} = tr_Probe_reg;
                    [d_cross_wire(j,i),Z_cross_wire_plot,tr_cross_wire_plot{tmp}] = ...
                        procrustes(XH_P,XI_total');
                    diff_plot = Z_cross_wire_plot - XH_P;
                    error_norm_plot{j,i} = vecnorm(diff_plot,2,2);
                    Mean_error_plot(j,i) = mean(error_norm_plot{j,i});
                    display(Mean_error_plot(j,i))
                    Max_error_plot(j,i) = max(error_norm_plot{j,i});
                    SD_error_plot(j,i) = std(error_norm_plot{j,i});
                    
                    [FRE_mean(tmp),FRE_max(tmp),FRE_SD(tmp)] = ...
                        FRE_Feedback_2lines(Z_cross_wire_plot,XH_ac,tr_reg_probe_col,tr_reg_phantom_col);
                    
                    tmp = tmp+1;
                end
                
            end
        end
    end
    clear image_no;
    n_images(j,2) = size(XH_total,2)
    cprintf('Red','Select new set %d\n',j+1);
    %     clear Marker_Probe_T Marker_Stylus_T Marker_Phantom_T Marker_Humerus_T tmp_USimage;
end

% %Procruste to transform points in Y to X
[ OPA_cross_wire, Rotation_cross_wire, Translation_cross_wire ] ...
    = OPA ( XI_total,XH_P' );

[d_cross_wire_total,Z_cross_wire_total,tr_cross_wire_total] = ...
    procrustes(XH_P,XI_total');
diff_total = Z_cross_wire_total - XH_P;
error_norm_total = vecnorm(diff_total,2,2);
Mean_error_total = mean(error_norm_total);
Max_error_total = max(error_norm_total);
SD_error_total = std(error_norm_total);

% [d_cross_wire,Z_cross_wire,tr_cross_wire] = procrustes(XH_P,XI_total','reflection',false);
B_ang = asin(tr_cross_wire_total.T(1,3));
A_ang = acos(tr_cross_wire_total.T(3,3)/cos(B_ang));
C_ang = acos(tr_cross_wire_total.T(1,1)/cos(B_ang));
% diff = Z_cross_wire - XH_P;

%Plot intersection points in different frames
figure();
subplot(1,3,1);scatter3 (XI_total(1,1:2:end),XI_total(2,1:2:end),XI_total(3,1:2:end))
hold on; scatter3 (XI_total(1,2:2:end),XI_total(2,2:2:end),XI_total(3,2:2:end))
title ('X_I','FontSize',17);xlabel('x','FontSize',17); ylabel('y','FontSize',17); zlabel('z','FontSize',17);
axis equal
subplot(1,3,2); scatter3 (XH_total(1,1:2:end),XH_total(2,1:2:end),XH_total(3,1:2:end))
hold on; scatter3 (XH_total(1,2:2:end),XH_total(2,2:2:end),XH_total(3,2:2:end))
title ('X_H','FontSize',17);xlabel('x','FontSize',17); ylabel('y','FontSize',17); zlabel('z','FontSize',17);
axis equal
XH_P = XH_P';
subplot(1,3,3); scatter3 (XH_P(1,1:2:end),XH_P(2,1:2:end),XH_P(3,1:2:end))
hold on; scatter3 (XH_P(1,2:2:end),XH_P(2,2:2:end),XH_P(3,2:2:end))
title ('X_P','FontSize',17);xlabel('x','FontSize',17); ylabel('y','FontSize',17); zlabel('z','FontSize',17);
scatter3(Marker_Probe_ref(:,1),Marker_Probe_ref(:,2),Marker_Probe_ref(:,3),'filled')
axis equal

clear Local_Landmark_SC_tmp USim error_norm_plot tmp_USimage data I

tmp = 1;
for m =1:size(Set,1)
    for n = 1:size(Mean_error_plot,2)
        if Mean_error_plot(m,n)~=0
%             conv(tmp) = Mean_error_plot(m,n);
%             conv_1(tmp) = Mean_error_plot_1(m,n);
%             conv_2(tmp) = Mean_error_plot_2(m,n);
%             conv_3(tmp) = Mean_error_plot_3(m,n);
            
            tmp = tmp+1;
        end
    end
end
n(1) = n_images(1,2)/2;
n(2) = n(1)+(n_images(2,2)-n_images(1,2))/2;
n(3) = n(2)+(n_images(3,2)-n_images(2,2))/2;
n(4) = n(3)+(n_images(4,2)-n_images(3,2))/2;
n(5) = n(4)+(n_images(5,2)-n_images(4,2))/2;
n(6) = n(5)+(n_images(6,2)-n_images(5,2))/2;
n(7) = n(6)+(n_images(7,2)-n_images(6,2))/2;
n(8) = n(7)+(n_images(8,2)-n_images(7,2))/2;


figure();hold on;plot(conv);
title ('');xlabel('number of input image','FontSize',17); ylabel('|X_P-^P_UT*X_U| [mm]','FontSize',17);
yL = ylim;
plot([n(1) n(1)], yL,'k--');
plot([n(2) n(2)], yL,'k--');
plot([n(3) n(3)], yL,'k--');
plot([n(4) n(4)], yL,'k--');
plot([n(5) n(5)], yL,'k--');
plot([n(6) n(6)], yL,'k--');
plot([n(7) n(7)], yL,'k--');
plot([n(8) n(8)], yL,'k--');
legend ('Triple-N','FontSize',17);

figure();hold on;
plot(FRE_mean);
yL = ylim;
plot([n(1) n(1)], yL,'k--');
plot([n(2) n(2)], yL,'k--');
plot([n(3) n(3)], yL,'k--');
plot([n(4) n(4)], yL,'k--');
plot([n(5) n(5)], yL,'k--');
plot([n(6) n(6)], yL,'k--');
plot([n(7) n(7)], yL,'k--');
plot([n(8) n(8)], yL,'k--');
title ('');xlabel('number of input image','FontSize',17); ylabel('FRE [mm]','FontSize',17);
legend ('Triple-N','Comb1','Comb2','Comb3','FontSize',10)
mean(abs(diff_plot),1)