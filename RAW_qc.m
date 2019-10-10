marker_number_start = 1;
marker_number_end = 14;
    tmp_column = [];
    
    %collect sufficient frame
    number_of_frame = size(data.FrameNumber,2);
    
    temp_frame = data.FrameNumber(1):(number_of_frame+data.FrameNumber(1)-1);
    subtract = data.FrameNumber - temp_frame;
    
    tmp_end_frame = find(subtract<0);
    
    if ( isempty(tmp_end_frame))
        end_frame = number_of_frame;
    else
        end_frame = tmp_end_frame(1)-1;
    end
    
for i = 1:end_frame-8   
    
    %Arrage markers in 3D array (x y z:number:timeframe)
    Marker_Probe_raw (:,1,i) = data.Marker1(:,i);
    Marker_Probe_raw (:,2,i) = data.Marker2(:,i);
    Marker_Probe_raw (:,3,i) = data.Marker3(:,i);
    Marker_Probe_raw (:,4,i) = data.Marker4(:,i);
    
    Marker_Stylus_raw (:,1,i) = data.Marker6(:,i);
    Marker_Stylus_raw (:,2,i) = data.Marker5(:,i);
    Marker_Stylus_raw (:,3,i) = data.Marker7(:,i);
     
    Marker_Phantom_raw (:,1,i) = data.Marker8(:,i);
    Marker_Phantom_raw (:,3,i) = data.Marker9(:,i);
    Marker_Phantom_raw (:,4,i) = data.Marker10(:,i);    
    Marker_Phantom_raw (:,2,i) = data.Marker11(:,i);
    
    Marker_Humerus_raw (:,1,i) = data.Marker12(:,i);
    Marker_Humerus_raw (:,2,i) = data.Marker13(:,i);
    Marker_Humerus_raw (:,3,i) = data.Marker14(:,i);
    
    tmp_USimage(:,:,i) = data.USimages(:,:,i);
end

%Filtering noise from RAW data with GPA   
    [ Marker_Probe, ~ ] = GPA ( Marker_Probe_raw );
    [ Marker_Phantom, ~ ] = GPA ( Marker_Phantom_raw );
    [ Marker_Stylus, ~ ] = GPA ( Marker_Stylus_raw );
    [ Marker_Humerus, ~ ] = GPA ( Marker_Humerus_raw );
%Calculate mean different between RAW data and GPA data    
    Difference = {'Probe','Phantom','Stylus','Stylus';...
        Marker_Probe_raw-Marker_Probe,Marker_Phantom_raw-Marker_Phantom,...
        Marker_Stylus_raw-Marker_Stylus,Marker_Humerus_raw-Marker_Humerus};
    Norm = {'Probe','Phantom','Stylus','Humerus';...
        vecnorm(Difference{2, 1}),vecnorm(Difference{2, 2}),...
        vecnorm(Difference{2, 3}),vecnorm(Difference{2, 4})};
    Mean = {'Probe','Phantom','Stylus','Humerus';...
        mean(Norm{2, 1}),mean(Norm{2, 2}),...
        mean(Norm{2, 3}),mean(Norm{2, 4})};
    Mean_different = {'Probe','Phantom','Stylus','Humerus';...
        mean(Mean{2, 1}),mean(Mean{2, 2}),...
        mean(Mean{2, 3}),mean(Mean{2, 4})};
    SD_different = {'Probe','Phantom','Stylus','Humerus';...
        std(Mean{2, 1}),std(Mean{2, 2}),...
        std(Mean{2, 3}),std(Mean{2, 4})};
    Max_different = {'Probe','Phantom','Stylus','Humerus';...
        max(max(Norm{2, 1})),max(max(Norm{2, 2})),...
        max(max(Norm{2, 3})),max(max(Norm{2, 4}))};
    Quality_Infor = {' ','Probe','Phantom','Stylus','Humerus';...
        'Mean', mean(Mean{2, 1}),mean(Mean{2, 2}),mean(Mean{2, 3}),mean(Mean{2, 4});...
        'Max', max(max(Norm{2, 1})),max(max(Norm{2, 2})),max(max(Norm{2, 3})),max(max(Norm{2, 4}));...
        'SD', std(Mean{2, 1}),std(Mean{2, 2}),std(Mean{2, 3}),std(Mean{2, 4})};
    
    figure()
    subplot(2,2,1);
    plot(squeeze(Norm{2,1}(:,1,:)))
    hold on;
    plot(squeeze(Norm{2,1}(:,2,:)))
    plot(squeeze(Norm{2,1}(:,3,:)))
    plot(squeeze(Norm{2,1}(:,4,:)))
    title('Probe Markers')
    subplot(2,2,2);
    plot(squeeze(Norm{2,2}(:,1,:)))
    hold on;
    plot(squeeze(Norm{2,2}(:,2,:)))
    plot(squeeze(Norm{2,2}(:,3,:)))
    plot(squeeze(Norm{2,2}(:,4,:)))  
    title('Phantom Markers')

    subplot(2,2,3);
    plot(squeeze(Norm{2,3}(:,1,:)))
    hold on;
    plot(squeeze(Norm{2,3}(:,2,:)))
    plot(squeeze(Norm{2,3}(:,3,:)))
    title('Stylus Markers')

    subplot(2,2,4);
    plot(squeeze(Norm{2,4}(:,1,:)))
    hold on;
    plot(squeeze(Norm{2,4}(:,2,:)))
    plot(squeeze(Norm{2,4}(:,3,:)))
    title('Humerus Markers')
