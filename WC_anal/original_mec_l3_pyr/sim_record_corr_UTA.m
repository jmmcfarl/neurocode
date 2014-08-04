%%calculate overall correlograms and up triggered MP averages

dsf = 8;
Fsd = 2016/dsf;
niqf = 2016/2;
hcf = 40/niqf;
lcf = 0.05/niqf;
[b,a] = butter(2,[lcf hcf]);
maxlag = round(2*Fsd);
tlags = [-maxlag:maxlag]/Fsd;

    load C:\WC_Germany\JMM_Analysis_pyr\sim_record_UDS_dur_2
    load used_data lf15 lf8
 
    wcv = filtfilt(b,a,lf15);
    wcv_d = downsample(wcv,dsf);
    wcv_d = zscore(wcv_d);
    
    lf8 = filtfilt(b,a,lf8);
    lf8_d = downsample(lf8,dsf);
    lf8_d = zscore(lf8_d);
    
    up_trans8(up_trans8<maxlag) = [];
    up_trans8(up_trans8>length(wcv_d)-maxlag) = [];
    
    num_8ups = length(up_trans8);
    
    up_trans_mat = zeros(num_8ups,2*maxlag+1);
    up_trans_mat8 = up_trans_mat;
    
    for t = 1:num_8ups
       
        up_trans_mat(t,:) = wcv_d(up_trans8(t)-maxlag:up_trans8(t)+maxlag);
        up_trans_mat8(t,:) = lf8_d(up_trans8(t)-maxlag:up_trans8(t)+maxlag);
    end
    
    up_trig_wcv = mean(up_trans_mat);
    up_trig_lfp = mean(up_trans_mat8);
    
    down_trans8(down_trans8<maxlag) = [];
    down_trans8(down_trans8>length(wcv_d)-maxlag) = [];
    
    num_8downs = length(down_trans8);
    
    down_trans_mat = zeros(num_8downs,2*maxlag+1);
    down_trans_mat8 = down_trans_mat;
    
    for t = 1:num_8downs
       
        down_trans_mat(t,:) = wcv_d(down_trans8(t)-maxlag:down_trans8(t)+maxlag);
        down_trans_mat8(t,:) = lf8_d(down_trans8(t)-maxlag:down_trans8(t)+maxlag);
        
    end
    
    down_trig_wcv = mean(down_trans_mat);
    down_trig_lfp = mean(down_trans_mat8);
    
    
    
    %calculate overall correlograms
    [l_acorr,lags] = xcov(lf8_d,maxlag,'coeff');
    [w_acorr,lags] = xcov(wcv_d,maxlag,'coeff');
    [x_corr,lags] = xcov(wcv_d,lf8_d,maxlag,'coeff');
    
    
    plot(lags/Fsd,x_corr,'linewidth',2)
    xlim([-0.2 0.2])
    grid 
    
plot(tlags,up_trig_wcv,'linewidth',2)
hold on
plot(tlags,up_trig_lfp,'r','linewidth',2)

plot(tlags,down_trig_wcv,'linewidth',2)
hold on
plot(tlags,down_trig_lfp,'r','linewidth',2)


