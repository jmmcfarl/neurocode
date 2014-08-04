clear all
close all

load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\overall_calcs\UDS_synch_state_dur\UDS_synch_state_dur_data
load C:\WC_Germany\overall_calcs\overall_info_file_data

dsf = 1;
Fs = 2016;
Fsd = Fs/dsf;
niqf = Fs/2;

[b,a] = butter(2,[1/niqf 15/niqf]);
[b2,a2] = butter(2,[20/niqf 80/niqf]);
[b2,a2] = butter(2,[0.5/niqf 15/niqf]);
maxLag = round(2*Fsd);
lags = -maxLag:maxLag;
freq = 1:0.25:15;
for d = 1:length(over_dir)
    d=19
    cd(over_dir{d})
    pwd
    
    load used_data wcv_minus_spike
    
    wcv_f = filtfilt(b,a,wcv_minus_spike/mp_gain(d));
    wcv_f2 = filtfilt(b2,a2,wcv_minus_spike/mp_gain(d));
    wcv_down = downsample(wcv_minus_spike,dsf);
    wcv_r = filtfilt(b2,a2,wcv_minus_spike/mp_gain(d));
    wcv_r = downsample(wcv_r,dsf);
    wcv_d = downsample(wcv_f,dsf);
    wcv_fd = downsample(wcv_f2,dsf);
    wcv_fd = jmm_smooth_1d(wcv_fd.^2,10);
%     wcv_d = zscore(wcv_d);
%     wcv_down = zscore(wcv_down);
%     wcv_r = zscore(wcv_r);
    
%% old up and down transitions
    up_trans_pts = round(up_trans{d}(synch_ups{d})*2);
    down_trans_pts = round(down_trans{d}(synch_ups{d})*2);
 
%% calculate new state transitions based on high frequency

thresh = 0.005;

temp = wcv_fd;
temp(wcv_fd > thresh) = 1;
temp(wcv_fd < thresh) = 0;
dtemp = [0;diff(temp)];
ucross = find(dtemp == 1);
dcross = find(dtemp == -1);

%makes sure first crossing is for an up transition
if ucross(1) > dcross(1)
    dcross(1) = [];
end
%make sure crossings end on a down transition
if ucross(end) > dcross(end)
    ucross(end) = [];
end

%get rid of supposed up states that dont have a point exceeding 0.05
badcross = [];
for i = 1:length(ucross)
    curmax = max(wcv_f(ucross(i):dcross(i)));
    if curmax < 0.05
        badcross = [badcross i];
    end
end
ucross(badcross) = [];
dcross(badcross) = [];

%% get rid of down states that are less than 0.5 seconds
shortdowns = [];
for i = 1:length(dcross)-1
    downdur = (ucross(i+1)-dcross(i))/Fsd;
    if downdur < 0.5
        shortdowns = [shortdowns i];
    end
end
dcross(shortdowns) = [];
ucross(shortdowns+1) = [];

shortups = [];
for i = 1:length(ucross)
    updur = (dcross(i)-ucross(i))/Fsd;
    if updur < 0.3
        shortups = [shortups i];
    end
end
dcross(shortups) = [];
ucross(shortups) = [];


up_state_dur = (dcross-ucross)/Fsd;
down_state_dur = (ucross(2:end)-dcross(1:end-1))/Fsd;

    bad_ups = find(up_state_dur < 2.0);
    bad_downs = find(down_state_dur < 2.0);

    up_state_dur(bad_ups) = [];
    down_state_dur(bad_downs) = [];
    
    good_up_st = ucross;
    good_up_st(bad_ups) = [];
    good_up_end = dcross;
    good_up_end(bad_ups) = [];
    good_down_st = dcross(1:end-1);
    good_down_end = ucross(2:end);
    good_down_st(bad_downs) = [];
    good_down_end(bad_downs) = [];
    
%%
    
    
%     bad_ups = find(up_state_dur{d}(synch_ups{d}) < 2.0);
%     bad_downs = find(down_state_dur{d}(synch_downs{d}) < 2.0);
%     
%     good_up_st = up_trans_pts;
%     good_up_st(bad_ups) = [];
%     good_up_end = down_trans_pts;
%     good_up_end(bad_ups) = [];
%     
    up_state_cor_mat = zeros(length(good_up_st),length(lags));
    up_state_pow = zeros(length(good_up_st),length(freq));
%     
%     good_down_st = down_trans_pts(1:end-1);
%     good_down_end = up_trans_pts(2:end);
%     good_down_st(bad_downs) = [];
%     good_down_end(bad_downs) = [];
%   

    down_state_cor_mat = zeros(length(good_down_st),length(lags));
     down_state_pow = zeros(length(good_down_st),length(freq));
  
    for i = 1:length(good_up_st)
        up_state_cor_mat(i,:) = xcov(wcv_d(good_up_st(i)+round(0.1*Fsd):good_up_end(i)-round(0.1*Fsd)),maxLag,'coeff');
        up_state_pow(i,:) = pwelch(wcv_r(good_up_st(i)+round(0.1*Fsd):good_up_end(i)-round(0.1*Fsd)),[],[],freq,Fsd);
    i
    end
    
    for i = 1:length(good_down_st)
        down_state_cor_mat(i,:) = xcov(wcv_d(good_down_st(i)+round(.1*Fsd):good_down_end(i)-round(0.1*Fsd)),maxLag,'coeff');
            down_state_pow(i,:) = pwelch(wcv_r(good_down_st(i)+round(0.1*Fsd):good_down_end(i)-round(0.1*Fsd)),[],[],freq,Fsd);
   i
    end
    
%     good_ups = synch_ups{d};
%     good_ups(bad_ups) = [];
%     good_downs = synch_downs{d};
%     good_downs(bad_downs) = [];
%    
%     [dummy,up_sort] = sort(up_state_dur{d}(good_ups));
%         [dummy,down_sort] = sort(down_state_dur{d}(good_downs));

%     figure
%     subplot(2,1,1)
%     pcolor(1:length(good_up_st),lags/Fsd,up_state_cor_mat(up_sort,:)');shading flat;colorbar
%     ylim([0 2])
%     caxis([-0.4 0.4])
%     title('Up state')
%     subplot(2,1,2)
%     pcolor(1:length(good_down_st),lags/Fsd,down_state_cor_mat(down_sort,:)');shading flat;colorbar
%     ylim([0 2])
%         caxis([-0.4 0.4])
%     title('Down state')
%     t_names = ['C:\WC_Germany\overall_calcs\up_down_corpow\wb_corr_sort_' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',t_names);
%     close

    
%     figure
%     subplot(2,1,1)
%     pcolor(1:length(good_up_st),freq,10*log10(up_state_pow(up_sort,:)'));shading flat;colorbar
%         title('Up state')
%         subplot(2,1,2)
%     pcolor(1:length(good_down_st),freq,10*log10(down_state_pow(down_sort,:)'));shading flat;colorbar
%         title('Down state')
%         caxis([-40 -5])
%     t_names = ['C:\WC_Germany\overall_calcs\up_down_corpow\short_filt_pow_sort_' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',t_names);
%     close

        figure
    subplot(2,1,1)
    pcolor(1:length(good_up_st),lags/Fsd,up_state_cor_mat');shading flat;colorbar
    ylim([0 2])
    caxis([-0.4 0.4])
    title('Up state')
    subplot(2,1,2)
    pcolor(1:length(good_down_st),lags/Fsd,down_state_cor_mat');shading flat;colorbar
    ylim([0 2])
        caxis([-0.4 0.4])
    title('Down state')
    t_names = ['C:\WC_Germany\overall_calcs\up_down_corpow\wb_short_corr_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',t_names);
    close

    
%     figure
%     subplot(2,1,1)
%     pcolor(1:length(good_up_st),freq,10*log10(up_state_pow'));shading flat;colorbar
%         title('Up state')
%         subplot(2,1,2)
%     pcolor(1:length(good_down_st),freq,10*log10(down_state_pow'));shading flat;colorbar
%         title('Down state')
%         caxis([-40 -5])
%     t_names = ['C:\WC_Germany\overall_calcs\up_down_corpow\short_filt_pow_' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',t_names);
%     close

    
end
