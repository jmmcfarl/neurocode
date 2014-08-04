clear all

load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\overall_info_file_data
load C:\WC_Germany\overall_calcs\eyeball_theta_times

Fs = 2016;
niqf = Fs/2;

lcf = 10/niqf;
[b,a] = butter(2,lcf,'high');
maxLag = round(0.1*Fs);
lags = -maxLag:maxLag;
for d = 1:33

    cd(over_dir{d});
    disp(['session ' num2str(d)]);

    load used_data lf8 lf3 lf2 wcv_minus_spike

    wcv_f = filtfilt(b,a,wcv_minus_spike/mp_gain(d));
    lf8_f = filtfilt(b,a,lf8/lf8_gain(d));
    lf3_f = filtfilt(b,a,lf3/lf3_gain(d));
    lf2_f = filtfilt(b,a,lf2/lf2_gain(d));
    
    if exist('spike_time_br.mat')
        load spike_time_br
    else
        load spike_time
    end


    if ~isempty(theta_times{d}) && ~isempty(nt_times{d})

        tot_theta = sum(theta_times{d}(:,2) - theta_times{d}(:,1));
        tot_nt = sum(nt_times{d}(:,2) - nt_times{d}(:,1));
        theta_perc(d) = tot_theta/(tot_theta+tot_nt);

        theta_spikes = [];
        nt_spikes = [];

        for i = 1:size(theta_times{d},1);
            theta_spikes = [theta_spikes; spkid(spkid > theta_times{d}(i,1)*Fs & spkid < theta_times{d}(i,2)*Fs)];
        end
        for i = 1:size(nt_times{d},1);
            nt_spikes = [nt_spikes; spkid(spkid > nt_times{d}(i,1)*Fs & spkid < nt_times{d}(i,2)*Fs)];
        end

        theta_sta_mat_8 = zeros(length(theta_spikes),length(lags));
        nt_sta_mat_8 = zeros(length(nt_spikes),length(lags));
        theta_sta_mat_3 = zeros(length(theta_spikes),length(lags));
        nt_sta_mat_3 = zeros(length(nt_spikes),length(lags));
        theta_sta_mat_2 = zeros(length(theta_spikes),length(lags));
        nt_sta_mat_2 = zeros(length(nt_spikes),length(lags));

        total_length = length(wcv_f);

        for i = 1:length(theta_spikes)
            if theta_spikes(i) > maxLag && total_length - theta_spikes(i) > maxLag
                theta_sta_mat_8(i,:) = lf8_f(theta_spikes(i)-maxLag:theta_spikes(i)+maxLag);
                theta_sta_mat_3(i,:) = lf3_f(theta_spikes(i)-maxLag:theta_spikes(i)+maxLag);
                                theta_sta_mat_2(i,:) = lf2_f(theta_spikes(i)-maxLag:theta_spikes(i)+maxLag);

            else
                theta_sta_mat_8(i,:) = nan;
                theta_sta_mat_3(i,:) = nan;
                                theta_sta_mat_2(i,:) = nan;

            end
        end
        for i = 1:length(nt_spikes)
            if nt_spikes(i) > maxLag && total_length - nt_spikes(i) > maxLag

                nt_sta_mat_8(i,:) = lf8_f(nt_spikes(i)-maxLag:nt_spikes(i)+maxLag);
                nt_sta_mat_3(i,:) = lf3_f(nt_spikes(i)-maxLag:nt_spikes(i)+maxLag);
                                nt_sta_mat_2(i,:) = lf2_f(nt_spikes(i)-maxLag:nt_spikes(i)+maxLag);

            else
                nt_sta_mat_8(i,:) = nan;
                nt_sta_mat_3(i,:) = nan;
                                nt_sta_mat_2(i,:) = nan;

            end
        end

        if ~isempty(theta_spikes)
            theta_sta_8(d,:) = nanmean(theta_sta_mat_8);
            theta_sta_3(d,:) = nanmean(theta_sta_mat_3);
                        theta_sta_2(d,:) = nanmean(theta_sta_mat_2);

        end
        if ~isempty(nt_spikes)
            nt_sta_8(d,:) = nanmean(nt_sta_mat_8);
            nt_sta_3(d,:) = nanmean(nt_sta_mat_3);
                        nt_sta_2(d,:) = nanmean(nt_sta_mat_2);

        end

        theta_rate(d) = length(theta_spikes)/tot_theta;
        nt_rate(d) = length(nt_spikes)/tot_nt;
        
        if ~isempty(theta_spikes) && ~isempty(nt_spikes)
        plot(lags/Fs,theta_sta_8(d,:),'r')
        hold on
        plot(lags/Fs,nt_sta_8(d,:),'r--')
        plot(lags/Fs,theta_sta_3(d,:),'k')
        plot(lags/Fs,nt_sta_3(d,:),'k--')
        plot(lags/Fs,theta_sta_2(d,:),'g')
        plot(lags/Fs,nt_sta_2(d,:),'g--')
        xlim([-0.1 0.1])
        legend('8theta','8','3theta','3','2theta','2')
        title([num2str(theta_perc(d)) 'Percent theta.  Theta Rate: ' num2str(theta_rate(d)) ' NT Rate: ' num2str(nt_rate(d))]);
        tname = ['C:\WC_Germany\overall_calcs\state_dep_sta\' num2str(cell_type(d)) '_' over_names{d}];
        print('-dpng',tname);
        close all
        end
        
        
    end


end