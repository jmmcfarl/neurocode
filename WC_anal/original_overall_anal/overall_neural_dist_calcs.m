clear all

load('C:\WC_Germany\overall_calcs\overall_dir.mat')

dsf = 8;
Fsd = 2016/dsf;
niqf = 2016/2;
lcf = 0.05/niqf;
hcf= 40/niqf;
[b,a] = butter(2,[lcf hcf]);
vxgrid = [-4; 4; 400];

for d = 1:length(over_dir)
    
    disp(sprintf('session %d',d))
    cd(over_dir{d});
    d
    
    load used_data wcv_minus_spike lf8 lf3 lf2
    
    %bandlimit signals
    down_w = filtfilt(b,a,wcv_minus_spike);
    down_8 = filtfilt(b,a,lf8);
    down_3 = filtfilt(b,a,lf3);
    down_2 = filtfilt(b,a,lf2);
    
    down_w = downsample(down_w,dsf);
    down_8 = downsample(down_8,dsf);
    down_3 = downsample(down_3,dsf);
    down_2 = downsample(down_2,dsf);
    
    %zscore
    down_w = zscore(down_w);
    down_8 = zscore(down_8);
    down_3 = zscore(down_3);
    down_2 = zscore(down_2);
    
    temp_cor = corrcoef(down_3,down_8);
    hc_cor = temp_cor(2,1)
    down3_sub = down_3-hc_cor*down_8;
    
    [wcv_dist(d,:),xgrid] = gpkde(down_w,-3,vxgrid);
    [lf8_dist(d,:),xgrid] = gpkde(down_8,-3,vxgrid);
    [lf3_dist(d,:),xgrid] = gpkde(down_3,-3,vxgrid);
    [lf2_dist(d,:),xgrid] = gpkde(down_2,-3,vxgrid);
    [lf3s_dist(d,:),xgrid] = gpkde(down3_sub,-3,vxgrid);


    plot(xgrid,wcv_dist(d,:),'linewidth',2)
    hold on
    plot(xgrid,lf8_dist(d,:),'r','linewidth',2)
    plot(xgrid,lf3_dist(d,:),'g','linewidth',2)
    plot(xgrid,lf2_dist(d,:),'k','linewidth',2)
    plot(xgrid,lf3s_dist(d,:),'m','linewidth',2)
    legend('WCV','LF8','LF3','LF2','LF3orth')
    xlim([-3 3])
    tname = ['C:\WC_Germany\overall_calcs\neural_dist\' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    close
  
    clear down* wcv* lf8 lf3 lf2
    
end


save C:\WC_Germany\overall_calcs\neural_dist\neural_dist_data xgrid *dist