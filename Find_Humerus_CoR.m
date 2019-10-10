clear all;
close all;

% set working directory
directory = 'H:\ResearchProject\calibration data\September 12\USscan\Humerus\';
%;'pta_t1.mat';'saggital_t2.mat';'atp_t1.mat';;'coronal_t4.mat';'pta_t1.mat'
set = {'humerus_slice_t3.mat'};

n_tests  = size(set,1); % we consider 5 tests (size of set)

%reserve ecpry array
curve_total = [];
curve_P_total = [];
curve_T_total = [];
curve_H_total = [];
curve_H_c1_total = [];
curve_H_c2_total = [];
curve_H_c3_total = [];

%load calibrated transformation matrix
Transformation_all = load('H:\ResearchProject\calibration data\September 12\12Workspace\all_3m_FRE.mat');
Transformation = Transformation_all.tr_cross_wire_total;
Transformation_c1 = Transformation_all.tr_cross_wire_comb{1, 1};
Transformation_c2 = Transformation_all.tr_cross_wire_comb{1, 2};
Transformation_c3 = Transformation_all.tr_cross_wire_comb{1, 3};
clear Transformation_all

% start locating the humerus 
for j = 1:n_tests
        load([directory set{j,:}])
    %Filtering data
    [Marker_Probe_T,Marker_Probe_ref,~,~,~,~,Marker_Humerus_T,Marker_Humerus_ref,...
        USimage,Quality_Infor{j}] = clean_frame_OPA (data,4,12);       
    divide = round(size(Marker_Humerus_T,3)/20);   
    image_no{j} = [1:divide:size(Marker_Humerus_T,3)];
    if divide ~= 0
    for i = 1:size(image_no{j},2)
        % show the original image
        USim = USimage(:,:,image_no{j}(i));
        
        figure(10);
        imshow(USim, [16 234]);
        %mouse-click for points selection
        [xi,yi] = getpts;
        curve = [xi yi];
        if xi ~= 0
        spline = cscvn(curve');
        hold on;
        fnplt(cscvn(curve'),'r',2)
        end
        curve(:,3) = 0;
            
        %Transform curve to Probe's LCS    
        curve_P = Transformation.b * curve * Transformation.T+ Transformation.c(1:size(curve,1),:);
        curve_P_c1 = Transformation_c1.b * curve * Transformation_c1.T+ Transformation_c1.c(1:size(curve,1),:);
        curve_P_c2 = Transformation_c2.b * curve * Transformation_c2.T+ Transformation_c2.c(1:size(curve,1),:);
        curve_P_c3 = Transformation_c3.b * curve * Transformation_c3.T+ Transformation_c3.c(1:size(curve,1),:);
        
        number(j,i) = size(curve,1)
        %Transform curve to Vicon
        [ d_Probe_tmp, Z_Probe_tmp, tr_Probe_tmp ] ...
                = procrustes ( Marker_Probe_T(:,2:4,image_no{j}(i))',Marker_Probe_ref(2:4,:),'scaling',false,'reflection',false);
        curve_T = tr_Probe_tmp.b * curve_P * tr_Probe_tmp.T + repmat(tr_Probe_tmp.c(1,:),size(curve_P,1),1); 
        curve_T_c1 = tr_Probe_tmp.b * curve_P_c1 * tr_Probe_tmp.T + repmat(tr_Probe_tmp.c(1,:),size(curve_P,1),1); 
        curve_T_c2 = tr_Probe_tmp.b * curve_P_c2 * tr_Probe_tmp.T + repmat(tr_Probe_tmp.c(1,:),size(curve_P,1),1); 
        curve_T_c3 = tr_Probe_tmp.b * curve_P_c3 * tr_Probe_tmp.T + repmat(tr_Probe_tmp.c(1,:),size(curve_P,1),1);
        
        scale(:,:,i) = tr_Probe_tmp.b;
        rotation(:,:,i) = tr_Probe_tmp.T;
        translation(:,:,i) = tr_Probe_tmp.c(1,:);
        

        %Determine LCS of the Humerus model
        [ d_Humerus_tmp, Z_Humerus_tmp, tr_Humerus_tmp ] ...
                = procrustes ( Marker_Humerus_ref,Marker_Humerus_T(:,:,image_no{j}(i))','scaling',false,'reflection',false);
        curve_H = tr_Humerus_tmp.b * curve_T * tr_Humerus_tmp.T + repmat(tr_Humerus_tmp.c(1,:),size(curve_T,1),1);  
        curve_H_c1 = tr_Humerus_tmp.b * curve_T_c1 * tr_Humerus_tmp.T + repmat(tr_Humerus_tmp.c(1,:),size(curve_T,1),1);  
        curve_H_c2 = tr_Humerus_tmp.b * curve_T_c2 * tr_Humerus_tmp.T + repmat(tr_Humerus_tmp.c(1,:),size(curve_T,1),1);  
        curve_H_c3 = tr_Humerus_tmp.b * curve_T_c3 * tr_Humerus_tmp.T + repmat(tr_Humerus_tmp.c(1,:),size(curve_T,1),1);  
% plot probe trajectories and their corresponding curves (same frame)
        figure(29);scatter3(curve_T(:,1),curve_T(:,2),curve_T(:,3))
        hold on;
        if xi ~= 0
        fnplt(cscvn(curve_T'),'r',2)
        end
        hold on;
        xlabel('X') 
        ylabel('Y') 
        zlabel('Z')
        
        figure(30);scatter3(curve_T(:,1),curve_T(:,2),curve_T(:,3))
        xlabel('X') 
        ylabel('Y') 
        zlabel('Z')
        if xi ~= 0
        hold on; fnplt(cscvn(curve_T'),'r',2)
        end
        hold on;
        probe_mark1 = Marker_Probe_T(:,1,image_no{j}(i));
        probe_mark2 = Marker_Probe_T(:,2,image_no{j}(i));
        probe_mark3 = Marker_Probe_T(:,3,image_no{j}(i));
        probe_mark4 = Marker_Probe_T(:,4,image_no{j}(i));
        
        scatter3 (probe_mark1(1),probe_mark1(2),probe_mark1(3))  
        scatter3 (probe_mark2(1),probe_mark2(2),probe_mark2(3))
        scatter3 (probe_mark3(1),probe_mark3(2),probe_mark3(3))
        scatter3 (probe_mark4(1),probe_mark4(2),probe_mark4(3))
        line1 = [probe_mark1';probe_mark2'];
        line2 = [probe_mark4';probe_mark3'];
        plot3 (line1(:,1),line1(:,2),line1(:,3))
        plot3 (line2(:,1),line2(:,2),line2(:,3))
        

        curve_total = [curve_total; curve];
        curve_P_total = [curve_P_total; curve_P];
        curve_T_total = [curve_T_total; curve_T];
        curve_H_total = [curve_H_total; curve_H];
        curve_H_c1_total = [curve_H_c1_total; curve_H_c1];
        curve_H_c2_total = [curve_H_c2_total; curve_H_c2];
        curve_H_c3_total = [curve_H_c3_total; curve_H_c3];
        
    end
    end
        
    end
%plot curves in different spaces
figure();scatter3(curve_T_total(:,1),curve_T_total(:,2),curve_T_total(:,3),'filled')
title ('curve_T');xlabel('x'); ylabel('y'); zlabel('z')
axis equal
figure();scatter3(curve_H_total(:,1),curve_H_total(:,2),curve_H_total(:,3),'filled')
title ('curve_H');xlabel('x'); ylabel('y'); zlabel('z')
axis equal

%% calculate CoR and radius using sphere-fitting method
%total transformation
figure();clf;
plot3(curve_H_total(:,1),curve_H_total(:,2),curve_H_total(:,3),'r.')
hold on;
[Center_LSE,Radius_LSE] = sphereFit(curve_H_total)
[Base_X,Base_Y,Base_Z] = sphere(20);
surf(Radius_LSE*Base_X+Center_LSE(1),...
    Radius_LSE*Base_Y+Center_LSE(2),...
    Radius_LSE*Base_Z+Center_LSE(3),'faceAlpha',0.3,'Facecolor','b')
scatter3(Marker_Humerus_ref(:,1),Marker_Humerus_ref(:,2),Marker_Humerus_ref(:,3),'filled')
axis equal

figure();clf;
plot3(curve_H_c1_total(:,1),curve_H_c1_total(:,2),curve_H_c1_total(:,3),'r.')
hold on;
[Center_LSE_c1,Radius_LSE_c1] = sphereFit(curve_H_c1_total)
[Base_X,Base_Y,Base_Z] = sphere(20);
surf(Radius_LSE_c1*Base_X+Center_LSE_c1(1),...
    Radius_LSE_c1*Base_Y+Center_LSE_c1(2),...
    Radius_LSE_c1*Base_Z+Center_LSE_c1(3),'faceAlpha',0.3,'Facecolor','g')
scatter3(Marker_Humerus_ref(:,1),Marker_Humerus_ref(:,2),Marker_Humerus_ref(:,3),'filled')
axis equal

figure();clf;
plot3(curve_H_c2_total(:,1),curve_H_c2_total(:,2),curve_H_c2_total(:,3),'r.')
hold on;
[Center_LSE_c2,Radius_LSE_c2] = sphereFit(curve_H_c2_total)
[Base_X,Base_Y,Base_Z] = sphere(20);
surf(Radius_LSE_c2*Base_X+Center_LSE_c2(1),...
    Radius_LSE_c2*Base_Y+Center_LSE_c2(2),...
    Radius_LSE_c2*Base_Z+Center_LSE_c2(3),'faceAlpha',0.3,'Facecolor','c')
scatter3(Marker_Humerus_ref(:,1),Marker_Humerus_ref(:,2),Marker_Humerus_ref(:,3),'filled')
axis equal
figure();clf;
plot3(curve_H_c3_total(:,1),curve_H_c2_total(:,2),curve_H_c2_total(:,3),'r.')
hold on;
[Center_LSE_c3,Radius_LSE_c3] = sphereFit(curve_H_c2_total)
[Base_X,Base_Y,Base_Z] = sphere(20);
surf(Radius_LSE_c3*Base_X+Center_LSE_c3(1),...
    Radius_LSE_c3*Base_Y+Center_LSE_c3(2),...
    Radius_LSE_c3*Base_Z+Center_LSE_c3(3),'faceAlpha',0.3,'Facecolor','y')
scatter3(Marker_Humerus_ref(:,1),Marker_Humerus_ref(:,2),Marker_Humerus_ref(:,3),'filled')
axis equal
clear USim USimage data

load('\\filestore.soton.ac.uk\users\ww1u18\mydocuments\ResearchProject\calibration data\September 17\Humerus_Landmark\Fitted.mat')