clear all
close all
cd C:\WC_Germany\Cortical_analysis

load C:\WC_Germany\Cortical_analysis\cortical_dir

Fs = 2016;
niqf = Fs/2;

lcf = 5/niqf;
% hcf = 250/niqf;

% [b,a] = butter(2,[lcf hcf]);
[b,a] = butter(2,lcf,'high');
dsf = 1;
Fsd = Fs/dsf;

maxlag = round(0.5*Fsd);
lags = -maxlag:maxlag;


for d = 1:length(sess_data)
    
    cd(sess_data(d).directory)
    pwd   
    load spike_time_jmm  
    load used_data lf8 lf5 lf7 lf6 
    spkid = round(spkid/dsf);
    
    lf8_f = filtfilt(b,a,lf8);
    lf7_f = filtfilt(b,a,lf7);
    lf6_f = filtfilt(b,a,lf6);
    lf5_f = filtfilt(b,a,lf5);
    
    lf8_d = downsample(lf8_f,dsf);
    lf7_d = downsample(lf8_f,dsf);
    lf6_d = downsample(lf7_f,dsf);
    lf5_d = downsample(lf5_f,dsf);
 
        lf8_d = zscore(lf8_f);
    lf7_d = zscore(lf8_f);
    lf6_d = zscore(lf7_f);
    lf5_d = zscore(lf5_f);

    
    lf8_sta_mat = zeros(length(spkid),length(lags));
    lf7_sta_mat = zeros(length(spkid),length(lags));
    lf6_sta_mat = zeros(length(spkid),length(lags));
    lf5_sta_mat = zeros(length(spkid),length(lags));
 
    datalen = length(lf8_d);
    for i = 1:length(spkid)
       if spkid(i) > maxlag & spkid(i)+maxlag < datalen
           lf8_sta_mat(i,:) = lf8_d(spkid(i)-maxlag:spkid(i)+maxlag);
            lf7_sta_mat(i,:) = lf7_d(spkid(i)-maxlag:spkid(i)+maxlag);
           lf6_sta_mat(i,:) = lf6_d(spkid(i)-maxlag:spkid(i)+maxlag);
           lf5_sta_mat(i,:) = lf5_d(spkid(i)-maxlag:spkid(i)+maxlag);
       else
           lf8_sta_mat(i,:) = nan;
           lf7_sta_mat(i,:) = nan;
           lf6_sta_mat(i,:) = nan;
           lf5_sta_mat(i,:) = nan;
       end     
    end
    
    lf8_sta = nanmean(lf8_sta_mat);
    lf7_sta = nanmean(lf7_sta_mat);
    lf6_sta = nanmean(lf6_sta_mat);
    lf5_sta = nanmean(lf5_sta_mat);
    
    plot(lags/Fsd,lf8_sta)
    hold on
    plot(lags/Fsd,lf7_sta,'r')
    plot(lags/Fsd,lf6_sta,'g')
    plot(lags/Fsd,lf5_sta,'k')
    legend('LF8','LF7','LF6','LF5')
    title(['Number of Spikes = ' num2str(length(spkid))])
               cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
t_names = ['C:\WC_Germany\Cortical_analysis\sta\' cell_name];
   print('-dpng',t_names);
   close
   
    
end

cd C:\WC_Germany\Cortical_analysis\sta
save cortical_sta_data *sta lags Fsd