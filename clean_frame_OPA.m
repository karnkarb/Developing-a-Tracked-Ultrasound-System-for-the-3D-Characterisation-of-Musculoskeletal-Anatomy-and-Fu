function [Marker_Probe_tmp,Probe_ref_conf,Marker_Stylus_tmp,Stylus_ref_conf,...
    Marker_Phantom_tmp,Phantom_ref_conf,Marker_Humerus_tmp,Humerus_ref_conf,...
    tmp_USimage,Quality_Infor] = clean_frame_OPA (data,mode,date)
    
%%input 
    %MODE 1: For data of Phantom registration
    %MODE 2: For data of Probe Calibration
    %MODE 3: For data of Humerus Registration
    %MODE 4: For data of Humerus Scan
    
    %date refers to operations
    %15: Curvilinear with Double-N
    %12: Curvilinear with Triple-N
    %17: Linear with Triple-N    
    
    p = mode;
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
    
    if date == 12
        Probe_ref_conf = [27.97711802	0	0;-71.3224611	0	0;...
        0.037077693	33.47267899	0.54430272;-0.037077693	-41.48390416	0.54430272];
        Phantom_ref_conf = [29.67995873	0	0;-71.96701007	0	0;...
        0.173924245	41.23693575	-0.730449091;-0.173924245	-34.42010967	-0.730449091];
        Humerus_ref_conf = [0	71.67691992	0;42.17769813	0	0;-33.02098741	0	0];
        Stylus_ref_conf = [0	8.34973872	0;-49.3835141	0	0;-110.1456908	0	0];
    elseif date == 17
        Probe_ref_conf = [28.75766431	0	0;-71.39933267	0	0;...
        0.12682366	32.80644912	-3.595924997;-0.12682366	-42.44972975	-3.595924997];
        Phantom_ref_conf = [28.93722547	0	0;-71.7407264	0	0;...
        0.079355238	41.89126961	-0.686342974;-0.079355238	-33.09824517	-0.686342974];
        Humerus_ref_conf = [0	28.77849425	0;33.07472256	0	0;-40.90671253	0	0];
        Stylus_ref_conf = [0	8.34973872	0;-49.3835141	0	0;-110.1456908	0	0];
    else 
        Probe_ref_conf = [28.42968646	0	0;-71.56360411	0	0;...
        0.360758739	33.32798025	-0.282002703;-0.360758739	-41.91083722	-0.282002703];     
        Phantom_ref_conf = [29.16818472	0	0;-71.98496373	0	0;...
        -0.449394774	42.56890067	-1.187737256;0.449394774	-33.43624802	-1.187737256];
        Humerus_ref_conf = [0	21.22298192	0;34.26591485	0	0;-52.75612594	0	0];
        Stylus_ref_conf = [0	7.654770404	0;-49.86169581	0	0;-110.1688029	0	0];
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
for t = 1:size(Marker_Humerus_raw,3)
    [ OPA_Marker_Probe(:,:,t), Rotation_Probe, Translation_Probe ] = OPA ( Marker_Probe_raw(:,:,t), Probe_ref_conf' );
    [ OPA_Marker_Probe_2(:,:,t), Rotation_Probe_2, Translation_Probe_2 ] = OPA ( Marker_Probe_raw(:,1:3,t), Probe_ref_conf(1:3,:)' );
    [ OPA_Marker_Probe_3(:,:,t), Rotation_Probe_3, Translation_Probe_3 ] = OPA ( Marker_Probe_raw(:,[1 2 4],t), Probe_ref_conf([1 2 4],:)' );
    [ OPA_Marker_Probe_4(:,:,t), Rotation_Probe_4, Translation_Probe_4 ] = OPA ( Marker_Probe_raw(:,[1 3 4],t), Probe_ref_conf([1 3 4],:)' );
    [ OPA_Marker_Probe_5(:,:,t), Rotation_Probe_5, Translation_Probe_3 ] = OPA ( Marker_Probe_raw(:,2:4,t), Probe_ref_conf(2:4,:)' );    
    [ OPA_Marker_Phantom(:,:,t), Rotation_Phantom(:,:,t), Translation_Phantom(:,:,t) ] = OPA ( Marker_Phantom_raw(:,:,t), Phantom_ref_conf'  );
    [ OPA_Marker_Phantom_2(:,:,t), Rotation_Phantom_2, Translation_Phantom_2 ] = OPA ( Marker_Phantom_raw(:,1:3,t), Phantom_ref_conf(1:3,:)' );
    [ OPA_Marker_Phantom_3(:,:,t), Rotation_Phantom_3, Translation_Phantom_3 ] = OPA ( Marker_Phantom_raw(:,[1 2 4],t), Phantom_ref_conf([1 2 4],:)' );
    [ OPA_Marker_Phantom_4(:,:,t), Rotation_Phantom_4, Translation_Phantom_4 ] = OPA ( Marker_Phantom_raw(:,[1 3 4],t), Phantom_ref_conf([1 3 4],:)' );
    [ OPA_Marker_Phantom_5(:,:,t), Rotation_Phantom_5, Translation_Phantom_5 ] = OPA ( Marker_Phantom_raw(:,2:4,t), Phantom_ref_conf(2:4,:)' );    
    [ OPA_Marker_Stylus(:,:,t), Rotation_Stylus(:,:,t), Translation_Stylus(:,:,t) ] = OPA ( Marker_Stylus_raw(:,:,t), Stylus_ref_conf'  );
    [ OPA_Marker_Humerus(:,:,t), Rotation_Humerus(:,:,t), Translation_Humerus(:,:,t) ] = OPA ( Marker_Humerus_raw(:,:,t), Humerus_ref_conf'  );
end
    
%Calculate mean different between RAW data and Ref. conf. data    
    Difference = {'Probe','Phantom','Stylus','Humerus';...
        repmat(Probe_ref_conf',[1 1 size(Marker_Humerus_raw,3)])-OPA_Marker_Probe,...
        repmat(Phantom_ref_conf',[1 1 size(Marker_Humerus_raw,3)])-OPA_Marker_Phantom,...
        repmat(Stylus_ref_conf',[1 1 size(Marker_Humerus_raw,3)])-OPA_Marker_Stylus,...
        repmat(Humerus_ref_conf',[1 1 size(Marker_Humerus_raw,3)])-OPA_Marker_Humerus;
        repmat(Probe_ref_conf(1:3,:)',[1 1 size(Marker_Humerus_raw,3)])-OPA_Marker_Probe_2,...
        repmat(Probe_ref_conf([1 2 4],:)',[1 1 size(Marker_Humerus_raw,3)])-OPA_Marker_Probe_3,...
        repmat(Probe_ref_conf([1 3 4],:)',[1 1 size(Marker_Humerus_raw,3)])-OPA_Marker_Probe_4,...
        repmat(Probe_ref_conf(2:4,:)',[1 1 size(Marker_Humerus_raw,3)])-OPA_Marker_Probe_5;...
        repmat(Phantom_ref_conf(1:3,:)',[1 1 size(Marker_Humerus_raw,3)])-OPA_Marker_Phantom_2,...
        repmat(Phantom_ref_conf([1 2 4],:)',[1 1 size(Marker_Humerus_raw,3)])-OPA_Marker_Phantom_3,...
        repmat(Phantom_ref_conf([1 3 4],:)',[1 1 size(Marker_Humerus_raw,3)])-OPA_Marker_Phantom_4,...
        repmat(Phantom_ref_conf(2:4,:)',[1 1 size(Marker_Humerus_raw,3)])-OPA_Marker_Phantom_5};
    
    Norm = {'Probe','Phantom','Stylus','Humerus';...
        vecnorm(Difference{2, 1}),vecnorm(Difference{2, 2}),...
        vecnorm(Difference{2, 3}),vecnorm(Difference{2, 4});...
        vecnorm(Difference{3, 1}),vecnorm(Difference{3, 2}),...
        vecnorm(Difference{3, 3}),vecnorm(Difference{3, 4});...
        vecnorm(Difference{4, 1}),vecnorm(Difference{4, 2}),...
        vecnorm(Difference{4, 3}),vecnorm(Difference{4, 4});};
    Mean = {'Probe','Phantom','Stylus','Humerus';...
        mean(Norm{2, 1}),mean(Norm{2, 2}),...
        mean(Norm{2, 3}),mean(Norm{2, 4})};
    Quality_Infor = {' ','Probe','Phantom','Stylus','Humerus';...
        'Mean', mean(Mean{2, 1}),mean(Mean{2, 2}),mean(Mean{2, 3}),mean(Mean{2, 4});...
        'Max', max(max(Norm{2, 1})),max(max(Norm{2, 2})),max(max(Norm{2, 3})),max(max(Norm{2, 4}));...
        'SD', std(Mean{2, 1}),std(Mean{2, 2}),std(Mean{2, 3}),std(Mean{2, 4})};
    %Eliminate untrackable markers frame
    if p == 1 %Phantom_Landmark
        Norm_ph_tmp = squeeze(Norm{4,4})';
        Norm_st_tmp = squeeze(Norm{2,3})';
        [row1,~] = find(Norm_ph_tmp>1.5);
        [row2,~] = find(Norm_st_tmp>1.5);
        row = unique([row1;row2]);
    elseif p == 2 %Probe_Calibration
        Norm_ph_tmp = squeeze(Norm{4,4})';
        Norm_pr_tmp = squeeze(Norm{3,4})';
        [row1,~] = find(Norm_ph_tmp>1.5);
        [row2,~] = find(Norm_pr_tmp>1.5);
        row = unique([row1;row2]);
    elseif p == 3 %Humerus_Landmark
        Norm_hm_tmp = squeeze(Norm{2,4})';
        Norm_st_tmp = squeeze(Norm{2,3})';
        [row1,~] = find(Norm_hm_tmp>1.5);
        [row2,~] = find(Norm_st_tmp>1.5);
        row = unique([row1;row2]);
    else %Humerus_CoR
        Norm_pr_tmp = squeeze(Norm{3,4})';
        Norm_hm_tmp = squeeze(Norm{2,4})';
        [row1,~] = find(Norm_pr_tmp>1.5);
        [row2,~] = find(Norm_hm_tmp>3);
        row = unique([row1;row2]);
    end
    Marker_Probe_tmp = Marker_Probe_raw;
    Marker_Stylus_tmp = Marker_Stylus_raw;
    Marker_Phantom_tmp = Marker_Phantom_raw;
    Marker_Humerus_tmp = Marker_Humerus_raw;
    Marker_Probe_tmp(:,:,row) = [];
    Marker_Stylus_tmp(:,:,row) = [];
    Marker_Phantom_tmp(:,:,row) = [];
    Marker_Humerus_tmp(:,:,row) = [];
    tmp_USimage = data.USimages;
    tmp_USimage(:,:,row) = [];
   
    
        figure()
        title('Probe with 3 Markers')
    subplot(2,2,1);
    plot(squeeze(Norm{3,1}(:,1,:)))
    hold on;
    plot(squeeze(Norm{3,1}(:,2,:)))
    plot(squeeze(Norm{3,1}(:,3,:)))
    subplot(2,2,2);
    plot(squeeze(Norm{3,2}(:,1,:)))
    hold on;
    plot(squeeze(Norm{3,2}(:,2,:)))
    plot(squeeze(Norm{3,2}(:,3,:)))
    subplot(2,2,3);
    plot(squeeze(Norm{3,3}(:,1,:)))
    hold on;
    plot(squeeze(Norm{3,3}(:,2,:)))
    plot(squeeze(Norm{3,3}(:,3,:)))
    subplot(2,2,4);
    plot(squeeze(Norm{3,4}(:,1,:)))
    hold on;
    plot(squeeze(Norm{3,4}(:,2,:)))
    plot(squeeze(Norm{3,4}(:,3,:)))
     figure()
        title('Phantom with 3 Markers')
    subplot(2,2,1);
    plot(squeeze(Norm{4,1}(:,1,:)))
    hold on;
    plot(squeeze(Norm{4,1}(:,2,:)))
    plot(squeeze(Norm{4,1}(:,3,:)))
    subplot(2,2,2);
    plot(squeeze(Norm{4,2}(:,1,:)))
    hold on;
    plot(squeeze(Norm{4,2}(:,2,:)))
    plot(squeeze(Norm{4,2}(:,3,:)))
    subplot(2,2,3);
    plot(squeeze(Norm{4,3}(:,1,:)))
    hold on;
    plot(squeeze(Norm{4,3}(:,2,:)))
    plot(squeeze(Norm{4,3}(:,3,:)))
    subplot(2,2,4);
    plot(squeeze(Norm{4,4}(:,1,:)))
    hold on;
    plot(squeeze(Norm{4,4}(:,2,:)))
    plot(squeeze(Norm{4,4}(:,3,:)))
    
    %plot error OPA VS RAW
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
    
    figure()
    scatter3(OPA_Marker_Probe(1,:), OPA_Marker_Probe(2,:), OPA_Marker_Probe(3,:))
    hold on;
    scatter3(Probe_ref_conf(:,1),Probe_ref_conf(:,2),Probe_ref_conf(:,3))
    axis equal

end