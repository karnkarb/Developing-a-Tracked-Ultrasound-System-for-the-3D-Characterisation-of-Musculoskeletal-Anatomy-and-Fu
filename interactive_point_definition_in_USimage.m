%--------------------------------------------------------------------------
% interactive_line_definition_in_USimage.m
%--------------------------------------------------------------------------
% load US images from probe calibration tests and interactively pick sets
% of 6 points marking the location of the wires in the image
% point 1: top left corner; point 3: top right, point 6: bottom right
%--------------------------------------------------------------------------
% v1.0 16.04.2017 (c) Markus O. Heller m.o.heller@soton.ac.uk
%--------------------------------------------------------------------------

% get rid of everything currently on the workspace, close all figures ...
clear all; clc; close all;

% prepare for loading the image data for the calibration trials 
% -- define the input_dir, and 
% -- the list of file_prefixes 
% NOTE: these might have to be adjusted!

input_dir    = '/Users/rebecca/Documents/Uni Work/Year 3/Individual Project/MATLAB code/ImageBank/calibration data';
file_prefix  = { 'StraightUpAndDown', 'AxisLeft', 'AxisRight_2', 'RotatedLeft', 'RotatedRight', 'FigEight'};
file_postfix = '.mat';

%-------------------------------------------------------------------
% define the image_no(no_tests, no_images) to consider
% Here the 1st index refers to the test, 
% and the 2nd index refers to the respective test image number 
%-------------------------------------------------------------------
image_no = [  1 21 117; ...
             50 72 100; ...
              7 54  90; ...
              5 19  49; ...
             29  9  92 ];

n_tests  = 5; % we consider 5 tests - i.e. we currently skip the 6th test
              %(figure of eight, "FigEight")
n_images = 3; % 3 test images for each test


% create a container for the figure
fig1 = figure(1);

for i = 1:n_tests
    
    % load the stack of images for the current test
    filename     = [ input_dir file_prefix{i} file_postfix ];
    load(filename);
    
    % loop over all images we want to assess for the current test
    for j = 1:n_images
        
        % set the image to the number we want to work with - this might ---
        % need to be adjusted! --------------------------------------------
        
        % now display the image
        image(data.USimages(:,:,image_no(i,j)))
        hold on
 
        % use a grey-scale color map
        colormap('gray')

        % now the magic: interactively select 6 points - do this 3 times --
        % for each line: left click to add points, double click to end 
        [x1,y1] = getpts(fig1)
        [x2,y2] = getpts(fig1)
        [x3,y3] = getpts(fig1)
        %------------------------------------------------------------------

        % convert the line data (doubles) to actual pixel values (integers) 
        % and get the intensities associated to each of the point sets 
        % (row/column swap !) ---------------------------------------------
        [ pts1 ] = round([ y1(1:6), x1(1:6) ]);
        [ pts2 ] = round([ y2(1:6), x2(1:6) ]);
        [ pts3 ] = round([ y3(1:6), x3(1:6) ]);
        
        for k=1:6 % size(pts1,1)
            intensity1(k) = data.USimages(pts1(k,1), pts1(k,2),image_no(i,j));   
        end
        for k=1:6 % size(pts2,1)
            intensity2(k) = data.USimages(pts2(k,1), pts2(k,2),image_no(i,j));   
        end
        for k=1:6 % size(pts3,1)
            intensity3(k) = data.USimages(pts3(k,1), pts3(k,2),image_no(i,j));   
        end

        % plot the point sets as plus es ----------------------------------
        plot(pts1(:,2), pts1(:,1), 'r+')
        plot(pts2(:,2), pts2(:,1), 'y+')
        plot(pts3(:,2), pts3(:,1), 'c+')
        %------------------------------------------------------------------

        % save the current figure to a file using a variety of formats ----
        
        outfile_prefix  = [ input_dir file_prefix{i} '_ImageNo_' num2str(image_no(i,j), '%03d') '.Points' ];
      
        outfile_postfix = '.emf';
        outfile = [ outfile_prefix outfile_postfix];   
        saveas(fig1, outfile);

        outfile_postfix = '.fig';
        outfile = [ outfile_prefix outfile_postfix];   
        saveas(fig1, outfile);

        outfile_postfix = '.eps';
        outfile = [ outfile_prefix outfile_postfix];   
        saveas(fig1, outfile);

        % outfile_postfix = '.pdf';
        % outfile = [ outfile_prefix outfile_postfix];   
        % saveas(fig1, outfile);

        outfile_postfix = '.png';
        outfile = [ outfile_prefix outfile_postfix];   
        saveas(fig1, outfile);

        outfile_postfix = '.svg';
        outfile = [ outfile_prefix outfile_postfix];   
        saveas(fig1, outfile);
        %------------------------------------------------------------------
        
        % collect the data in a variable where we can hold all info -------
        pts{i,j,1} = [ pts1 intensity1' ];
        pts{i,j,2} = [ pts2 intensity2' ];
        pts{i,j,3} = [ pts3 intensity3' ];        
        %------------------------------------------------------------------
                       
    end % for j
    
end % for i

% remove the USimage data object from the work space ----------------------
clear data;

% remove the current figure from the work space ---------------------------
clear fig1;

% calculate mean and SD data for the 3 repetitions of each image ----------
for i=1:n_tests
    for j=1:n_images
        mean_res{i,j} = (pts{i,j,1} + pts{i,j,2} + pts{i,j,3})/3.0;
        tmp_res = cat(3, pts{i,j,1}, pts{i,j,2}, pts{i,j,3});
        sd_res{i,j} = std(tmp_res, 0, 3); % 0 - compute SD using N-1
    end
end

% now save all data on the workspace to a .mat file -----------------------
outfile_data_postfix = '.mat';
outfile = [ input_dir 'US_CalibrationImage_UserDefined_WirePoints' outfile_data_postfix ];
save(outfile)
%--------------------------------------------------------------------------

% save the data to an EXCEL file ------------------------------------------
outfile_data_postfix = '.xlsx';
outfile = [ input_dir 'US_CalibrationImage_UserDefined_WirePoints' outfile_data_postfix ];
for i=1:n_tests
    for j=1:n_images
        sheet_name = [ 'RawData_Test' num2str(i) '_Image' num2str(j) ];  
        tmp_array  = [ pts{i,j,1}, pts{i,j,2}, pts{i,j,3} ];
        T = table( tmp_array(:,1), tmp_array(:,2), tmp_array(:,3), ...
                   tmp_array(:,4), tmp_array(:,5), tmp_array(:,6), ...
                   tmp_array(:,7), tmp_array(:,8), tmp_array(:,9) );
        T.Properties.VariableNames = {'x1' 'y1' 'i1' 'x2' 'y2' 'i2' 'x3' 'y3' 'i3'};
        writetable(T, outfile, 'Sheet', sheet_name);
        
        sheet_name = [ 'MeanData_Test' num2str(i) '_Image' num2str(j) ];
        tmp_array  = mean_res{i,j};
        T = table( tmp_array(:,1), tmp_array(:,2), tmp_array(:,3) );
        T.Properties.VariableNames = {'Mean_x' 'Mean_y' 'Mean_i'};
        writetable(T, outfile, 'Sheet', sheet_name);
        
        sheet_name = [ 'SDData_Test' num2str(i) '_Image' num2str(j) ];
        tmp_array  = sd_res{i,j};
        T = table( tmp_array(:,1), tmp_array(:,2), tmp_array(:,3) );
        T.Properties.VariableNames = {'SD_x' 'SD_y' 'SD_i'};
        writetable(T, outfile, 'Sheet', sheet_name);
    end
end

%--------------------------------------------------------------------------