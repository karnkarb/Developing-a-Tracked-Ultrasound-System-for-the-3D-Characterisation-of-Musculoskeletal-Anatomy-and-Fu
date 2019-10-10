
%--------------------------------------------------------------------------
% automatic_fiducial_finding.m
%--------------------------------------------------------------------------
% loads US images from probe calibration tests and picks the 6 N-fiducials
% marking the location of the wires in the image.
% point 1: top left corner; point 3: top right, point 6: bottom right
%--------------------------------------------------------------------------

% get rid of everything currently on the workspace, close all figures ...
clear all;
close all;
warning('off','all');

% prepare to load the images
set = [ "Calibration\MIddle_static_1" "Calibration\Front_static_1" ...
    "Calibration\Upper_layer_1" "Calibration\Lower_layer_1" "Calibration\Tilting_part1"...
    "Calibration\Tilting_part2" "Calibration\Front_to_Back_1" "Calibration\Back_to_Front_1" ];

n_tests  = size(set,2);

XI_total = [];
XH_total = [];
XH_X = [];
XH_P = [];
XH_T = [];
Construct = [];

%Markers (Long-Short) in LCS of the Probe 
Marker_Probe_P = [48.874,107.0689,109.6146; 48.86750,7.06890,109.6146;...
    90.3675,78.5633,109.6146;15.3675,78.5746,109.6146];

%Markers (Long-Short) in LCS of SolidWorks
Marker_Phantom_Sol = [ -7.8 117.19 -77.58; -7.8 167.17 9.03; -7.8 167.36 -73.64; -7.8 102.41 -36.15];

%Landmarks (#ascending 1-16) in LCS of SolidWorks
Landmark_Sol = [1.79 30 -55; 1.79 30 -15; 1.79 5 -5; 1.79 5 -55;...
    138.21 30 -5; 138.21 30 -55; 138.21 10 -5; 138.21 5 -55;...
    30 34 -55; 75 34 -55; 120 34 -55; 130 34 -40;...
    120 34 -5; 55 34 -5; 30 34 -5; 20 34 -20];

[ Tsol(:,:), Rsol(:,:), orisol(:,:), success_status ] = determine_tracked_probe_LCS(Marker_Phantom_Sol);
%apply Phantom LCS Transform to matrix (nx3) 
Local_Sol(:,:) = apply_4by4_tmatrix(Tsol(:,:),Marker_Phantom_Sol);
Local_Landmark_Sol(:,:) = apply_4by4_tmatrix( Tsol(:,:),Landmark_Sol);

%Landmark (#ascending 1-16) in LCS of Phantom by SCoRE
Local_Landmark_SC_tmp = load('SCoRE_Landmarks.mat');
Local_Landmark_SC = Local_Landmark_SC_tmp.Center_SCoRE;

%Find the transformation between actual Phantom and CAD model Phantom
[ d_SS, Z_SS, tr_SS ] ...
   = procrustes ( Local_Landmark_SC',Local_Landmark_Sol);

% 1:size(set,2)
for j = 1
    
    % load US images data
    [Marker_Probe_T,~,Marker_Phantom_T,~,tmp_USimage] = clean_frame (set(j),1,11);
    
%     if j == 2 
%         image_no(:) = 1:40;
%     elseif j == 4
%         image_no(:) = 1:75;
%     elseif j == 7 || j == 8 
%         image_no(:) = 1:140;
%     else
%         image_no(:) = 1:size(Marker_Probe_T,3);
%     end

        image_no(:) = 1:2;
        n_images(j,:) = size(image_no,2)
        
        USim = tmp_USimage(:,:,image_no(1));
        
        % show the original image
        figure(1);
        imshow(USim, [16 234]);
        h = gca;
        h.Visible = 'On';
%         crop(j,:) = getrect;
        crop = [319,428,230,122;414,358,225,129;362,367,252,127;331,321,236,116;...
            490,362,225,136;285,360,228,145;390,373,231,113;433,368,230,120];
        XH_size = size(XH_total,2);
           
    for i = 1:n_images(j)        
       
        % load US images
        USim = tmp_USimage(:,:,image_no(i));
        
%         figure(1)
%         imshow(USim, [16 234]);
        % show the original image
        figure(2)
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
        
        % remove background noise and smooth image
        img = noise_cleaning(img, [2 1], 1);
        subplot(2,2,3),imshow(img, [16 234]);
        
        % binarize the original image
%          img = binary_thresholding(img, 0);
        %find optimal thresholding

        thre = (multithresh(img) - 0)/255.0;

        %normalise the image
        I = img/255.0;

        %binarize the image
        img = im2bw(I, thre);

        subplot(2,2,4),imshow(img, [0 1]);
        [xi,yi] = getpts;
        % remove the rectangular waterbath noise
%         [K_img, K_img2, ctr] = kmeans_clustering(img, 6);
%         figure(3)
%         subplot(2,1,1),imshow(K_img);
%         subplot(2,1,2),imshow(K_img2);

        % find the location of the centroids
        centroids = get_fiducial_centroids(img);

        % find the geometries of the lines 
        [slopes, intercepts, fiducial_indices] = get_lines(centroids);

        % display the lines on the image
        if(size(slopes, 2) < 2)
            fprintf('Image no.%.0f All 6 fiducials were not detected\n',i);
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
                X1 = [xi(1)+x,yi(1)+y];
                X2 = [xi(2)+x,yi(2)+y];
                X3 = [xi(3)+x,yi(3)+y];
                X4 = [xi(4)+x,yi(4)+y];
                X5 = [xi(5)+x,yi(5)+y];
                X6 = [xi(6)+x,yi(6)+y];
                % calculate alpha
                X12 = [X1;X2];
                d12 = pdist(X12,'euclidean');
                
                X13 = [X1;X3];
                d13 = pdist(X13,'euclidean');
                
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
                XI(:, 4) = 1;
                XI = XI';
                
                % XH(ri,c1,c2) = A(ri,c1)+alphai(D(ri,c2)-A(ri,c1)) 
                % Phantom plane matrix
%% Transform Phantom to Vicon capture
                %Calculate middlde point position using N-wire geometry in
                %SolidWorks frame
                XH1 = h1 +alpha*(O1-h1);
                XH1 = XH1';
                %Transform to LCS of Phantom model frame
                XH1_H = apply_4by4_tmatrix( Tsol(:,:),XH1);
                %Transform to actual Phantom frame               
                Sol_XH1 = tr_SS.b * XH1_H' * tr_SS.T + tr_SS.c(1:size(XH1_H',1),:);
                %Determine the LCS of Markers in each time frame
                [ T(:,:,i), R(:,:,i), ori(:,:,i), success_status ] = determine_tracked_probe_LCS(Marker_Phantom_T(:,:,image_no(i))');
                
                %apply Phantom LCS Transform to matrix (nx3) 
                Marker_Phantom_tmp(:,:,i) = apply_4by4_tmatrix( T(:,:,i),Marker_Phantom_T(:,:,image_no(i))');
                
                %Determine transformation of the markers to Vicon frame
                [ d_Phantom_tmp, Z_Phantom_tmp, tr_Phantom_tmp ] ...
                    = procrustes ( Marker_Phantom_T(:,:,image_no(i))', Marker_Phantom_tmp(:,:,i));
                    
                                           
                XH1_T = tr_Phantom_tmp.b * Sol_XH1 * tr_Phantom_tmp.T + tr_Phantom_tmp.c(1:size(Sol_XH1,1),:);
                
%%               
                %Calculate middlde point position using N-wire geometry in
                %SolidWorks frame
                XH4 = H4 +alpha2*(o4-H4);
                XH4 = XH4';
                %Transform to LCS of Phantom model frame
                XH4_H = apply_4by4_tmatrix( Tsol(:,:),XH4);
                %Transform to actual Phantom frame                
                Sol_XH4 = tr_SS.b * XH4_H' * tr_SS.T+tr_SS.c(1:size(XH4_H',1),:);
                
                XH4_T = tr_Phantom_tmp.b * Sol_XH4 * tr_Phantom_tmp.T + tr_Phantom_tmp.c(1:size(Sol_XH4,1),:);
%% Transform Vicon capture to Probe
    [ Tpro(:,:,i), Rpro(:,:,i), oripro(:,:,i), success_status ] = determine_tracked_probe_LCS(Marker_Probe_T(:,:,image_no(i))');
    Marker_Probe_tmp(:,:,i) = apply_4by4_tmatrix( Tpro(:,:,i),Marker_Probe_T(:,:,image_no(i))');
    [ d_Probe_tmp, Z_Probe_tmp, tr_Probe_tmp ] ...
                    = procrustes ( Marker_Probe_tmp(:,:,i),Marker_Probe_T(:,:,image_no(i))');
                
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
                        XH_P = [XH_P; XH1_P; XH4_P];
                        XH_T = [XH_T; XH1_T; XH4_T];
                        Construct = [Construct XI(:,1) XI(:,3) XI(:,4) XI(:,6)];
                        diff_transHT = Z_Phantom_tmp - Marker_Phantom_T(:,:,image_no(i))';
                        error_norm_transHT(i,:,j) = vecnorm(diff_transHT,2,2);
                        Mean_transHT(i,:,j) = mean(error_norm_transHT(i,:,j));
                        diff_transTP = Z_Probe_tmp - Marker_Probe_T(:,:,image_no(i))';
                        error_norm_transTP(i,:,j) = vecnorm(diff_transTP,2,2);
                        Mean_transTP(i,:,j) = mean(error_norm_transTP(i,:,j));
                        
                    end
                else
                XI_total = [XI_total XI(:,2) XI(:,5)];
                XH_total = [XH_total XH1 XH4];
                XH_X = [XH_X; X2; X5;];
                XH_P = [XH_P; XH1_P; XH4_P];
                XH_T = [XH_T; XH1_T; XH4_T];
                
                Construct = [Construct XI(:,1) XI(:,3) XI(:,4) XI(:,6)];
                diff_transHT = Z_Phantom_tmp - Marker_Phantom_T(:,:,image_no(i))';
                error_norm_transHT(i,:,j) = vecnorm(diff_transHT,2,2);
                Mean_transHT(i,:,j) = mean(error_norm_transHT(i,:,j));
                diff_transTP = Z_Probe_tmp - Marker_Probe_T(:,:,image_no(i))';
                error_norm_transTP(i,:,j) = vecnorm(diff_transTP,2,2);
                Mean_transTP(i,:,j) = mean(error_norm_transTP(i,:,j));
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
% [d,Z,tr] = procrustes(XI_proc,XH_proc);
XI_total (4,:) = [];
[ OPA_cross_wire, Rotation_cross_wire, Translation_cross_wire ] ...
   = OPA ( XI_total,XH_P' );

[d_cross_wire,Z_cross_wire,tr_cross_wire] = procrustes(XH_P,XI_total','reflection',false);
% B_ang = asin(tr_cross_wire.T(1,3));
% A_ang = acos(tr_cross_wire.T(3,3)/cos(B_ang));
% C_ang = acos(tr_cross_wire.T(1,1)/cos(B_ang));
% diff = Z_cross_wire - XH_P;

%Preserve empty cell array
n = size(XI_total,2);
XI_random = {}; XH_random = {};
d_cross_wire = {}; Z_cross_wire = {}; tr_cross_wire = {};
diff = {}; error_norm = {}; 

%% Sampling data from 10-100% and apply Procrustes analysis
% error is calculated from the differnce between points in Probe's LCS and
% points transform from procrustes (Z)

index1 = 1;
for perc = 0.1:0.05:1
    
for samp = 1:10000
    random{index1,samp} = randsample(n,round(size(XI_total,2)*perc));
    XI_random{index1,samp} = [(XI_total(:,random{index1,samp}))'];
    XH_random{index1,samp} = [XH_P(random{index1,samp},:)];
    
    [d_cross_wire{index1,samp},Z_cross_wire{index1,samp},tr_cross_wire{index1,samp}] = ...
        procrustes(XH_random{index1,samp},XI_random{index1,samp});  
    
    diff{index1,samp} = Z_cross_wire{index1,samp} - XH_random{index1,samp};
    error_norm{index1,samp} = vecnorm(diff{index1,samp},2,2);
    Mean_error(index1,samp) = mean(error_norm{index1,samp});
    min_max_mean_sd(index1,1) = min(Mean_error(index1,:));
    min_max_mean_sd(index1,2) = max(Mean_error(index1,:));
    min_max_mean_sd(index1,3) = mean(Mean_error(index1,:));
    min_max_mean_sd(index1,4) = std(Mean_error(index1,:));
    
        %Calculate Volume of the random samples
    
        width_volume{index1,samp} = max(XH_random{index1,samp}(:,1)) - min(XH_random{index1,samp}(:,1));
        length_volume{index1,samp} = max(XH_random{index1,samp}(:,2)) - min(XH_random{index1,samp}(:,2));
        depth_volume{index1,samp} = max(XH_random{index1,samp}(:,3)) - min(XH_random{index1,samp}(:,3));
        Volume(index1,samp)  =  width_volume{index1,samp}*length_volume{index1,samp}*depth_volume{index1,samp};
end
index1 = index1+1;
end



sorted = sort(Mean_error(1,:));
cap = sorted(50);
points = find(Mean_error(1,:)<=cap);
for m = 1:50
    no_Points(:,m) = random{1,points(m)};
    Volume_err(:,m) = Volume(1,points(m));
end

sorted_no_Points = sort(no_Points,1);
for o = 1:size(sorted_no_Points,2)
    
count(1,o) = size(find(sorted_no_Points(:,o) <= n_images(1,2)),1);
count(2,o) = size(find(sorted_no_Points(:,o) <= n_images(2,2)),1);
count(3,o) = size(find(sorted_no_Points(:,o) <= n_images(3,2)),1);
count(4,o) = size(find(sorted_no_Points(:,o) <= n_images(4,2)),1);
count(5,o) = size(find(sorted_no_Points(:,o) <= n_images(5,2)),1);
count(6,o) = size(find(sorted_no_Points(:,o) <= n_images(6,2)),1);
count(7,o) = size(find(sorted_no_Points(:,o) <= n_images(7,2)),1);

end

figure()
xx = [1:size(Mean_error,1)];
yy = min_max_mean_sd(:,3);
err = min_max_mean_sd(:,4)/sqrt(1000);
errorbar(xx,yy,err)

% Eliminate points with error_norm higher than 10
eliminate = find(error_norm{19,1}> 10);
XH_eli(:,:) = XH_random{19, 1};
XH_eli(eliminate,:) = [];
XI_eli(:,:) = XI_random{19, 1};
XI_eli(eliminate,:) = [];
[d_cross_wire_eli,Z_cross_wire_eli,tr_cross_wire_eli] = ...
        procrustes(XH_eli,XI_eli); 
    diff_eli = Z_cross_wire_eli - XH_eli;
    error_norm_eli = vecnorm(diff_eli,2,2);
    Mean_error_eli = mean(error_norm_eli);


figure(); scatter3 (XH_total(1,1:2:end),XH_total(2,1:2:end),XH_total(3,1:2:end))
hold on; scatter3 (XH_total(1,2:2:end),XH_total(2,2:2:end),XH_total(3,2:2:end))
xlabel('x')
ylabel('y')
zlabel('z')

XH_P = XH_P';
figure(); scatter3 (XH_P(1,1:2:end),XH_P(2,1:2:end),XH_P(3,1:2:end))
hold on; scatter3 (XH_P(1,2:2:end),XH_P(2,2:2:end),XH_P(3,2:2:end))
xlabel('x')
ylabel('y')
zlabel('z')

figure(); scatter3 (XI_total(1,1:2:end),XI_total(2,1:2:end),XI_total(3,1:2:end))
hold on; scatter3 (XI_total(1,2:2:end),XI_total(2,2:2:end),XI_total(3,2:2:end))
xlabel('x')
ylabel('y')
zlabel('z')

0.0177297298392205	0.989808075634794	0.141299787996921
0.967085759633955	0.0189009618831251	-0.253747684034964
-0.253832208738457	0.141147890696387	-0.956894185769209

0.195378150873029

58.9268631036768	-89.3226785040787	-82.6307962602305
% 
% figure(); scatter3 (Z_cross_wire(1:2:end,1),Z_cross_wire(1:2:end,2),Z_cross_wire(1:2:end,2))
% hold on; scatter3 (Z_cross_wire(2:2:end,1),Z_cross_wire(2:2:end,2),Z_cross_wire(2:2:end,2))
% 
% figure(); scatter3 (OPA_cross_wire(1,1:2:end),OPA_cross_wire(2,1:2:end),OPA_cross_wire(3,1:2:end))
% hold on; scatter3 (OPA_cross_wire(1,2:2:end),OPA_cross_wire(2,2:2:end),OPA_cross_wire(3,2:2:end))