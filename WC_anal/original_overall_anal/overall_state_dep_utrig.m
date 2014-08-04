clear all

load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\overall_info_file_data
load C:\WC_Germany\overall_calcs\eyeball_theta_times
load C:\WC_Germany\overall_calcs\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\overall_calcs\UDS_synch_state_dur\UDS_synch_state_dur_data

Fs = 2016;
niqf = Fs/2;
dsf = 8;
Fsd = Fs/dsf;
lcf = 2/niqf;
hcf = 6/niqf;
[b,a] = butter(1,[lcf hcf]);
maxLag = round(2*Fsd);
lags = -maxLag:maxLag;

for d = 1:33

    cd(over_dir{d});
    disp(['session ' num2str(d)]);

    load used_data lf8 lf3 lf2 wcv_minus_spike

    if ~isempty(theta_times{d})



        wcv_f = filtfilt(b,a,wcv_minus_spike);
        lf8_f = filtfilt(b,a,lf8);
        lf3_f = filtfilt(b,a,lf3);
        lf2_f = filtfilt(b,a,lf2);

        down_w = downsample(wcv_f,dsf);
        down_8 = downsample(lf8_f,dsf);
        down_3 = downsample(lf3_f,dsf);
        down_2 = downsample(lf2_f,dsf);

        down_w = zscore(down_w);
        down_8 = zscore(down_8);
        down_3 = zscore(down_3);
        down_2 = zscore(down_2);

        tot_theta = sum(theta_times{d}(:,2) - theta_times{d}(:,1));
        if ~isempty(nt_times{d})
            tot_nt = sum(nt_times{d}(:,2) - nt_times{d}(:,1));
            theta_perc(d) = tot_theta/(tot_theta+tot_nt);
        end
        total_length = length(wcv_f);

        theta_ups = [];
        nt_ups = [];

        for i = 1:size(theta_times{d},1)
            theta_ups = [theta_ups synch_ups{d}(up_trans{d}(synch_ups{d}) > theta_times{d}(i,1)*Fsd ...
                & up_trans{d}(synch_ups{d}) < theta_times{d}(i,2)*Fsd)];
        end
        for i = 1:size(nt_times{d},1)
            nt_ups = [nt_ups synch_ups{d}(up_trans{d}(synch_ups{d}) > nt_times{d}(i,1)*Fsd & ...
                up_trans{d}(synch_ups{d}) < nt_times{d}(i,2)*Fsd)];
        end


        mp_utrig_mp_mat_theta = zeros(length(theta_ups),length(lags));
        mp_utrig_lf8_mat_theta = zeros(length(theta_ups),length(lags));
        mp_utrig_lf3_mat_theta = zeros(length(theta_ups),length(lags));
        mp_utrig_lf2_mat_theta = zeros(length(theta_ups),length(lags));

        mp_utrig_mp_mat_nt = zeros(length(nt_ups),length(lags));
        mp_utrig_lf8_mat_nt = zeros(length(nt_ups),length(lags));
        mp_utrig_lf3_mat_nt = zeros(length(nt_ups),length(lags));
        mp_utrig_lf2_mat_nt = zeros(length(nt_ups),length(lags));


        for i = 1:length(theta_ups)

            if up_trans{d}(theta_ups(i)) > maxLag && ...
                    length(down_w) - up_trans{d}(theta_ups(i)) > maxLag

                mp_utrig_mp_mat_theta(i,:) = down_w(up_trans{d}(theta_ups(i))-maxLag:...
                    up_trans{d}(theta_ups(i))+maxLag);
                mp_utrig_lf8_mat_theta(i,:) = down_8(up_trans{d}(theta_ups(i))-maxLag:...
                    up_trans{d}(theta_ups(i))+maxLag);
                mp_utrig_lf3_mat_theta(i,:) = down_3(up_trans{d}(theta_ups(i))-maxLag:...
                    up_trans{d}(theta_ups(i))+maxLag);
                mp_utrig_lf2_mat_theta(i,:) = down_2(up_trans{d}(theta_ups(i))-maxLag:...
                    up_trans{d}(theta_ups(i))+maxLag);

            else

                mp_utrig_mp_mat_theta(i,:) = nan;
                mp_utrig_lf8_mat_theta(i,:) = nan;
                mp_utrig_lf3_mat_theta(i,:) = nan;
                mp_utrig_lf2_mat_theta(i,:) = nan;

            end

        end

        for i = 1:length(nt_ups)

            if up_trans{d}(nt_ups(i)) > maxLag && ...
                    length(down_w) - up_trans{d}(nt_ups(i)) > maxLag

                mp_utrig_mp_mat_nt(i,:) = down_w(up_trans{d}(nt_ups(i))-maxLag:...
                    up_trans{d}(nt_ups(i))+maxLag);
                mp_utrig_lf8_mat_nt(i,:) = down_8(up_trans{d}(nt_ups(i))-maxLag:...
                    up_trans{d}(nt_ups(i))+maxLag);
                mp_utrig_lf3_mat_nt(i,:) = down_3(up_trans{d}(nt_ups(i))-maxLag:...
                    up_trans{d}(nt_ups(i))+maxLag);
                mp_utrig_lf2_mat_nt(i,:) = down_2(up_trans{d}(nt_ups(i))-maxLag:...
                    up_trans{d}(nt_ups(i))+maxLag);

            else

                mp_utrig_mp_mat_nt(i,:) = nan;
                mp_utrig_lf8_mat_nt(i,:) = nan;
                mp_utrig_lf3_mat_nt(i,:) = nan;
                mp_utrig_lf2_mat_nt(i,:) = nan;

            end

        end

        mp_utrig_mp_theta(d,:) = nanmean(mp_utrig_mp_mat_theta);
        mp_utrig_lf8_theta(d,:) = nanmean(mp_utrig_lf8_mat_theta);
        mp_utrig_lf3_theta(d,:) = nanmean(mp_utrig_lf3_mat_theta);
        mp_utrig_lf2_theta(d,:) = nanmean(mp_utrig_lf2_mat_theta);
        if ~isempty(nt_times{d})
            mp_utrig_mp_nt(d,:) = nanmean(mp_utrig_mp_mat_nt);
            mp_utrig_lf8_nt(d,:) = nanmean(mp_utrig_lf8_mat_nt);
            mp_utrig_lf3_nt(d,:) = nanmean(mp_utrig_lf3_mat_nt);
            mp_utrig_lf2_nt(d,:) = nanmean(mp_utrig_lf2_mat_nt);
        else
            mp_utrig_mp_nt(d,:) = nan;
            mp_utrig_lf8_nt(d,:) = nan;
            mp_utrig_lf3_nt(d,:) = nan;
            mp_utrig_lf2_nt(d,:) = nan;
        end


        %         plot(lags/Fsd,mp_utrig_mp_theta)
        %         hold on
        %         plot(lags/Fsd,mp_utrig_mp_nt,'r')
        %         legend('theta','no theta')
        %         if ~isempty(nt_times{d})
        %         title([num2str(theta_perc(d)) 'Percent theta']);
        %         else
        %             title('all theta')
        %         end
        %         tname = ['C:\WC_Germany\overall_calcs\state_dep_trig\theta_MP_' num2str(cell_type(d)) '_' over_names{d}];
        %         print('-dpng',tname);
        %         close all
        %
        %           plot(lags/Fsd,mp_utrig_lf8_theta)
        %         hold on
        %         plot(lags/Fsd,mp_utrig_lf8_nt,'r')
        %         legend('theta','no theta')
        %         if ~isempty(nt_times{d})
        %         title([num2str(theta_perc(d)) 'Percent theta']);
        %         else
        %             title('all theta')
        %         end
        %         tname = ['C:\WC_Germany\overall_calcs\state_dep_trig\theta_LF8' num2str(cell_type(d)) '_' over_names{d}];
        %         print('-dpng',tname);
        %         close all
        %
        %                   plot(lags/Fsd,mp_utrig_lf3_theta)
        %         hold on
        %         plot(lags/Fsd,mp_utrig_lf3_nt,'r')
        %         legend('theta','no theta')
        %         if ~isempty(nt_times{d})
        %         title([num2str(theta_perc(d)) 'Percent theta']);
        %         else
        %             title('all theta')
        %         end
        %         tname = ['C:\WC_Germany\overall_calcs\state_dep_trig\theta_LF3' num2str(cell_type(d)) '_' over_names{d}];
        %         print('-dpng',tname);
        %         close all
        %
        %                   plot(lags/Fsd,mp_utrig_lf2_theta)
        %         hold on
        %         plot(lags/Fsd,mp_utrig_lf2_nt,'r')
        %         legend('theta','no theta')
        %         if ~isempty(nt_times{d})
        %         title([num2str(theta_perc(d)) 'Percent theta']);
        %         else
        %             title('all theta')
        %         end
        %         tname = ['C:\WC_Germany\overall_calcs\state_dep_trig\theta_LF2' num2str(cell_type(d)) '_' over_names{d}];
        %         print('-dpng',tname);
        %         close all

    end

end