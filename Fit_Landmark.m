clear all;
close all;

%set file directory 
Direct = 'H:\ResearchProject\calibration data\September 12\Analysed points\02-10\';
%Load all 16 points data
for p=1:16
    input_data = sprintf('P%d_t2_analysed.mat', p);
    pnt{p} = load ([Direct input_data]);
    Center_SCoRE(:,p) = pnt{p}.Center_SCoRE (:,:,1);
    Center_SphereFit(:,p) = pnt{p}.Center_LSE (:,:)';
    Range_of_Motion (:,p) = pnt{p}.Range_Of_Motion(:,:);
    Scaled_residual (:,p) = pnt{p}.Scaled_Residual(:);
end

%Landmarks (#ascending 1-16) in LCS of SolidWorks 
Center_Sol = [1.79 30 -55; 1.79 30 -15; 1.79 5 -5; 1.79 5 -55;...
     138.21 30 -5; 138.21 30 -55; 138.21 10 -5;138.21 5 -55;...
    30 34 -55; 75 34 -55; 120 34 -55; 130 34 -40;...
    120 34 -5; 55 34 -5; 30 34 -5; 20 34 -20];

%Fit Centers to SolidWorks model with Procrustes analysis
[ d_SC, Z_SC, tr_SC ] ...
   = procrustes ( Center_Sol,Center_SCoRE','reflection',false);
[ d_SF, Z_SF, tr_SF ] ...
   = procrustes ( Center_Sol,Center_SphereFit','reflection',false);

%Calculate error 
Difference_SCoRE = Center_Sol - Z_SC;
Difference_SphereFit = Center_Sol - Z_SF;
Norm_SCoRE = vecnorm(Difference_SCoRE,2,2);
Norm_SphereFit = vecnorm(Difference_SphereFit,2,2);
Mean_SCoRE = mean(Norm_SCoRE);
Mean_SphereFit = mean(Norm_SphereFit);
SD_SCoRE = std(Norm_SCoRE);
SD_SphereFit = std(Norm_SphereFit);
%Plot error from all fitted points
figure();plot(Norm_SCoRE)
hold on; plot(Norm_SphereFit)
xlabel ('points')
ylabel ('error from procrustes transformed points[mm]')
legend('SCoRE','Spherefit')
savefig('Error_allpoints.fig')

for r = 1:16
motion_vol(r) = (Range_of_Motion(1,r)*Range_of_Motion(2,r)*Range_of_Motion(3,r));
end
figure();scatter(motion_vol,Norm_SCoRE)
%eliminate high error point (#16) and repeat the Procrustes analysis
Center_Sol([5 6 13],:) = [];
Center_SCoRE(:,[5 6 13]) = [];
Center_SphereFit(:,[5 6 13]) = [];
[ d_SC_2, Z_SC_2, tr_SC_2 ] ...
   = procrustes ( Center_Sol,Center_SCoRE','reflection',false);
[ d_SF_2, Z_SF_2, tr_SF_2 ] ...
   = procrustes ( Center_Sol,Center_SphereFit','reflection',false);
Difference_SCoRE_2 = Center_Sol - Z_SC_2;
Difference_SphereFit_2 = Center_Sol - Z_SF_2;
Norm_SCoRE_2 = vecnorm(Difference_SCoRE_2,2,2);
Norm_SphereFit_2 = vecnorm(Difference_SphereFit_2,2,2);
Mean_SCoRE_2 = mean(Norm_SCoRE_2);
Mean_SphereFit_2 = mean(Norm_SphereFit_2);
SD_SCoRE_2 = std(Norm_SCoRE_2);
SD_SphereFit_2 = std(Norm_SphereFit_2);
%Plot error from remain fitted point
figure();plot(Norm_SCoRE_2)
hold on; plot(Norm_SphereFit_2)
xlabel ('points')
ylabel ('error from procrustes transformed points[mm]')
legend('SCoRE','Spherefit')
%save data
savefig('Error_allpoints.fig')
savefig('Error_eliminatedpoints.fig')
