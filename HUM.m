clear all;
close all;

directory = ['H:\ResearchProject\calibration data\August 15\Humerus_USimage\'];
set = [ "Stylus_HumeralHead_Slide_2" ];

%     "P1_ROM.mat" "P2_ROM.mat" "P3_ROM.mat" "P4_ROM.mat"
%     "P5_ROM.mat" "P6_ROM.mat" "P7_ROM.mat" "P8_ROM.mat" 
%     "P9_ROM.mat" "P11_ROM.mat" "P13_ROM.mat" "P15_ROM.mat"
%     "P5_ROM_2" "P6_ROM_2" "P7_ROM_2" "P8_ROM_2" "P9_ROM_2" "P13_ROM_2" "P15_ROM_2" "P10_ROM_2" 
%     "P9_3.mat" "P10_ROM_2" "P12_ROM_2" "P16_ROM_2" "P15_1" "P16_1"

Marker_Probe = [];
Marker_Stylus = [];
Marker_Phantom = [];

for j = 1:size(set,2)
    load( directory+set(j))
    number_of_frame = size(data.FrameNumber,2);
    
    temp_frame = data.FrameNumber(1):(number_of_frame+data.FrameNumber(1)-1);
    subtract = data.FrameNumber - temp_frame;
    
    tmp_end_frame = find(subtract<0);
    
    if ( isempty(tmp_end_frame))
        end_frame = number_of_frame;
    else
        end_frame = tmp_end_frame(1)-1;
    end
        
for i = 1:end_frame   
    
    
    Marker_Probe (:,1,i) = data.Marker1(:,i);
    Marker_Probe (:,2,i) = data.Marker3(:,i);
    Marker_Probe (:,3,i) = data.Marker4(:,i);
    Marker_Probe (:,4,i) = data.Marker2(:,i);
    
    Marker_Stylus (:,1,i) = data.Marker6(:,i);
    Marker_Stylus (:,2,i) = data.Marker5(:,i);
    Marker_Stylus (:,3,i) = data.Marker7(:,i);
     
    Marker_Phantom (:,1,i) = data.Marker8(:,i);
    Marker_Phantom (:,2,i) = data.Marker10(:,i);
    Marker_Phantom (:,3,i) = data.Marker9(:,i);
    Marker_Phantom (:,4,i) = data.Marker11(:,i);
    
    Marker_Humerus (:,1,i) = data.Marker12(:,i);
    Marker_Humerus (:,2,i) = data.Marker13(:,i);
    Marker_Humerus (:,3,i) = data.Marker14(:,i);
    
end

[ GPA_Marker_Probe, conv ] = GPA ( Marker_Probe );
[ GPA_Marker_Phantom, conv ] = GPA ( Marker_Phantom );
[ GPA_Marker_Stylus, conv ] = GPA ( Marker_Stylus );
   
%Calculate CoR
[ Singular_Values, ...% Reduced_Singular_Values, ...
           Residual, Scaled_Residual(j), ...
           Angular_Variance, Range_Of_Motion(j,:), ...
           Local_Centres, Local_Axis_Points, Local_Axis_Directions, ...
           Global_Centres_1, Global_Centres_2, Mean_Global_Centres, ...
           Global_Axis_Points_1, Global_Axis_Directions_1, ... 
           Global_Axis_Points_2, Global_Axis_Directions_2, ...
           Global_Mean_Axis_Points, Global_Mean_Axis_Directions ] ...
    = SARA ( GPA_Marker_Phantom, GPA_Marker_Stylus); 

Tip = Mean_Global_Centres;



for k = 1:end_frame

figure(1)
probe_mark1 = GPA_Marker_Probe(:,1,k);
probe_mark2 = GPA_Marker_Probe(:,2,k);
probe_mark3 = GPA_Marker_Probe(:,3,k);
probe_mark4 = GPA_Marker_Probe(:,4,k);
scatter3 (probe_mark1(1),probe_mark1(2),probe_mark1(3))
hold on;
scatter3 (probe_mark2(1),probe_mark2(2),probe_mark2(3))
hold on;
scatter3 (probe_mark3(1),probe_mark3(2),probe_mark3(3))
hold on;
scatter3 (probe_mark4(1),probe_mark4(2),probe_mark4(3))

line1 = [probe_mark1';probe_mark4'];
line2 = [probe_mark2';probe_mark3'];
hold on;
plot3 (line1(:,1),line1(:,2),line1(:,3))
hold on;
plot3 (line2(:,1),line2(:,2),line2(:,3))


figure(2)
scatter3 (Tip(1,k),Tip(2,k),Tip(3,k))
hold on;
phantom_mark1 = GPA_Marker_Phantom(:,1,k);
phantom_mark2 = GPA_Marker_Phantom(:,2,k);
phantom_mark3 = GPA_Marker_Phantom(:,3,k);
phantom_mark4 = GPA_Marker_Phantom(:,4,k);

scatter3 (phantom_mark1(1),phantom_mark1(2),phantom_mark1(3))
hold on;
scatter3 (phantom_mark2(1),phantom_mark2(2),phantom_mark2(3))
hold on;
scatter3 (phantom_mark3(1),phantom_mark3(2),phantom_mark3(3))
hold on;
scatter3 (phantom_mark4(1),phantom_mark4(2),phantom_mark4(3))

line4 = [phantom_mark1';phantom_mark4'];
line5 = [phantom_mark2';phantom_mark3'];
hold on;
plot3 (line4(:,1),line4(:,2),line4(:,3))
hold on;
plot3 (line5(:,1),line5(:,2),line5(:,3))

    if k>1
        figure(3)
        tail  = [GPA_Marker_Stylus(:,1,k-1)';GPA_Marker_Stylus(:,1,k)'];
        mid = [GPA_Marker_Stylus(:,2,k-1)';GPA_Marker_Stylus(:,2,k)'];
        upper = [GPA_Marker_Stylus(:,3,k-1)';GPA_Marker_Stylus(:,3,k)']; 
        vector = mid - tail;
        unit_vec(1,:) = vector(1,:)/norm(vector(1,:));
        unit_vec(2,:) = vector(2,:)/norm(vector(2,:));
        tip = unit_vec*153 + tail;

        plot3 (tip(:,1),tip(:,2),tip(:,3))
        hold on;
        % plot3 (tail(:,1),tail(:,2),tail(:,3))
        % hold on;
        % plot3 (mid(:,1),mid(:,2),mid(:,3))
        % hold on;
        % plot3 (upper(:,1),upper(:,2),upper(:,3))

    end
end
end
