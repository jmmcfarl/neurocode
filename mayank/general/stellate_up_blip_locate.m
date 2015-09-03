clear all
close all

load('C:\WC_Germany\JMM_analysis_ste\stellate_heka_dir.mat')
% load('C:\WC_Germany\JMM_analysis_pyr\pyr_heka_dir.mat')

for cdir =1;
    load(f_loc{cdir})
    dat_name = f_loc{cdir}(25:end);

    eval(['data = ' dat_name '_MP;'])
    eval(['time = ' dat_name '_sampletimes;'])
    Fs = 5000;
    dsf = 10;
    Fsd = Fs/dsf;
    niqf = Fs/2;
    [b,a] = butter(2,3/niqf,'low');
    [b2,a2] = butter(2,[1.5/niqf 6/niqf]);

    sweep_t = [0;find(diff(time) < 0)];
    low = min(data);
    high = max(data)-0.3;
    dataseg = data(1:sweep_t(1));
    filtseg = filtfilt(b,a,dataseg);

    maxnoise = 0.04;
    mindowndur = 0.75; %min dur in sec
    minupdur = 1.0;
    minpeakheight = 0.2;

    amp_diff_thresh = 0.01;
    
    down_cnt = 0;
    up_cnt = 0;


    for i = 1:length(sweep_t)-1

        dataseg = data(sweep_t(i)+1:sweep_t(i+1));
        filtseg = filtfilt(b,a,dataseg);
        filtseg2 = filtfilt(b2,a2,dataseg);
        down_filt = downsample(filtseg,dsf);
        diff_filt = [0;diff(down_filt)]*Fsd;

        if max(diff_filt > minpeakheight)
            [dummy,cur_peaks] = findpeaks(diff_filt,'minpeakheight',minpeakheight);
        else
            cur_peaks = [];
        end
        cur_down_start = [];
        cur_down_stop = [];
        cur_up_start = [];
        cur_up_stop = [];


        %locate peaks in filt2 signal
        [temp_peaks,temp_peak_locs] = findpeaks(filtseg2);
        [temp_npeaks,temp_npeak_locs] = findpeaks(-filtseg2);
        temp_npeaks = -temp_npeaks;
        
        %make sure you start on up blip and end in down blip
        if temp_peak_locs(1) > temp_npeak_locs(1)
            temp_npeak_locs(1) = [];
            temp_npeaks(1) = [];
        end
        if temp_peak_locs(end) > temp_npeak_locs(end)
            temp_peak_locs(end) = [];
            temp_peaks(end) = [];
        end
        
        num_peaks = length(temp_npeak_locs);
        pk_locs = zeros(2*num_peaks,1);
        pk_amps = zeros(2*num_peaks,1);
        for nm = 1:num_peaks
            pk_locs(nm*2-1)=temp_peak_locs(nm);
            pk_locs(nm*2)=temp_npeak_locs(nm);
            pk_amps(nm*2-1) = temp_peaks(nm);
            pk_amps(nm*2) = temp_npeaks(nm);
        end

        dpeak = [amp_diff_thresh+1; diff(pk_amps)];
        
        bad_peaks = find(abs(dpeak) < amp_diff_thresh);
               
        bad_peaks = [bad_peaks bad_peaks-1];
        pk_id = zeros(1,length(pk_locs));
        pk_id(mod((1:length(pk_locs)),2)==1) = 1;
        pk_locs(bad_peaks) = [];
        pk_amps(bad_peaks) = [];
        pk_id(bad_peaks) = [];
        pk_id = logical(pk_id);
        up_locs = pk_locs(pk_id);
        down_locs = pk_locs(~pk_id);
        
        time_vec = time(sweep_t(i)+1:sweep_t(i+1));




               plot(time_vec,data(sweep_t(i)+1:sweep_t(i+1)))
%                hold on
%                   plot(time_vec,filtseg2,'r','linewidth',2)
                  hold on
%                   plot(time_vec(cur_down_start),dataseg(cur_down_start),'go','markersize',12,'linewidth',2)
%                   plot(time_vec(cur_down_stop),dataseg(cur_down_stop),'ro','markersize',12,'linewidth',2)
                  plot(time_vec(up_locs),dataseg(up_locs),'ko','markersize',8,'linewidth',2)
                  plot(time_vec(down_locs),dataseg(down_locs),'mo','markersize',8,'linewidth',2)

                  %                   plot(time_vec(temp_npeak_locs),filtseg2(temp_npeak_locs),'mo','markersize',8,'linewidth',2)

%                   plot(time_vec(cur_up_start),dataseg(cur_up_start),'g*')
%                   plot(time_vec(cur_up_stop),dataseg(cur_up_stop),'r*')
        
%                ylim([low high])
        
           grid
           pause
           clf

    end


%     maxlag = 0.75*Fs;
%     lags = -maxlag:maxlag;
%     
%     %% now go through down states and create xcov matrix
%     niqf = Fs/2;
%     [b,a] = butter(2,2/niqf,'high');
%     xcov_mat = zeros(length(down_state_array),length(lags));
%     freqd = linspace(0,niqf,mindowndur*Fs);
%     pow_mat = zeros(length(down_state_array),length(freqd));
%     for i = 1:length(down_state_array)
%         mean_mp(i) = mean(down_state_array{i});
%         xcov_mat(i,:) = xcov(filtfilt(b,a,down_state_array{i}),maxlag,'coeff');
%         pow_mat(i,:) = periodogram(down_state_array{i}-mean(down_state_array{i}),[],freqd,Fs);
%         i
%     end
%     [dummy,mp_order] = sort(mean_mp);
% %     figure
% %     pcolor(lags/Fs,1:length(down_state_array),xcov_mat);
% %     shading flat
% %     colorbar
% %     xlim([0 maxlag/Fs])
%     figure
%     pcolor(freqd,1:length(down_state_array),10*log10(pow_mat));
%     shading flat
%     xlim([0 50])
%     caxis([-100 -60])
%     colorbar
%     t_names = ['C:\WC_Germany\down_spec_mat_' f_names{cdir}];
%     print('-dpng',t_names);
%     close 
%     
%     %% now go through up states and create xcov matrix
%     niqf = Fs/2;
%     [b,a] = butter(2,2/niqf,'high');
%     xcov_mat_u = zeros(length(up_state_array),length(lags));
%     frequ = linspace(0,niqf,minupdur*Fs);
%     pow_mat_u = zeros(length(up_state_array),length(frequ));
%     for i = 1:length(up_state_array)
%         mean_mp(i) = mean(up_state_array{i});
%         xcov_mat_u(i,:) = xcov(filtfilt(b,a,up_state_array{i}),maxlag,'coeff');
%         pow_mat_u(i,:) = periodogram(up_state_array{i}-mean(up_state_array{i}),[],frequ,Fs);
%     end
%     [dummy,mp_order] = sort(mean_mp);
% %     figure
% %     pcolor(lags/Fs,1:length(up_state_array),xcov_mat_u);
% %     shading flat
% %     colorbar
% %     xlim([0 maxlag/Fs])
%     figure
%     pcolor(frequ,1:length(up_state_array),10*log10(pow_mat_u));
%     shading flat
%     xlim([0 50])
%     caxis([-70 -10])
%     colorbar
%         t_names = ['C:\WC_Germany\up_spec_mat_' f_names{cdir}];
%     print('-dpng',t_names);
% close
    
end