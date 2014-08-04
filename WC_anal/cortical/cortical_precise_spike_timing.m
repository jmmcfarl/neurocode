clear all
close all
cd C:\WC_Germany\Cortical_analysis

load C:\WC_Germany\Cortical_analysis\cortical_dir
load C:\WC_Germany\Cortical_analysis\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\Cortical_analysis\UDS_synch_state_dur\UDS_synch_state_dur_data
load C:\WC_Germany\Cortical_analysis\sigmoid_fit\sig_fit_alldata

Fs = 2016;
niqf = Fs/2;
lcf = 0.05/niqf;
hcf = 40/niqf;
[b,a] = butter(2,[lcf hcf]);

maxlag = 0.5*Fs;
lags = -maxlag:maxlag;

for d = 1:length(sess_data)


    cd(sess_data(d).directory)
    disp(num2str(d))

    load used_data lf8 wcv_minus_spike
    load spike_time_jmm

    lf8_f = filtfilt(b,a,lf8);
    lf8_f = zscore(lf8_f);

    wcv_f = filtfilt(b,a,wcv_minus_spike);
    wcv_f = zscore(wcv_f);

    t_50_wup{d} = rlid_wup{d};

    %initialize
    mp_utrig_lf8_mat = zeros(length(synch_ups{d}),length(lags));
    mp_utrig_mp_mat = zeros(length(synch_ups{d}),length(lags));


    datalen = length(wcv_f);

    %calculate mp utrigs
    for i = 1:length(synch_ups{d})

        if t_10_wup{d}(i) > maxlag && datalen - t_90_wup{d}(i) > maxlag && ~isnan(t_50_wup{d}(i))

            mp_utrig_mp_mat(i,:) = wcv_f(t_50_wup{d}(i)-maxlag:t_50_wup{d}(i)+maxlag);
            mp_utrig_lf8_mat(i,:) = lf8_f(t_50_wup{d}(i)-maxlag:t_50_wup{d}(i)+maxlag);
        else
            mp_utrig_mp_mat(i,:) = nan;
            mp_utrig_lf8_mat(i,:) = nan;
        end

    end


    %% plot mp up trig matrices
    Fig = figure(1)
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    imagesc(lags/Fs,(1:length(synch_ups{d}))/length(synch_ups{d}),mp_utrig_mp_mat);
    shading flat; colorbar; hold on
    caxis([-3 3]);colorbar
    line([0 0],[0 1],'Color','k')
    num_trans = size(mp_utrig_mp_mat,1);
    for i = 1:num_trans
        cur_spikes = find((spkid > t_50_wup{d}(i)-maxlag) & (spkid < t_50_wup{d}(i)+maxlag));
        cur_spikes = spkid(cur_spikes) - t_50_wup{d}(i);
        plot(cur_spikes/Fs,ones(size(cur_spikes))*i/num_trans,'k.')
    end
    xlim([-0.2 0.5])
    cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
    t_names = ['C:\WC_Germany\Cortical_analysis\spike_timing\mp_utrig_mp_' cell_name];
    print('-dpng',t_names);
    close

        Fig = figure(1)
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    imagesc(lags/Fs,(1:length(synch_ups{d}))/length(synch_ups{d}),mp_utrig_lf8_mat);
    shading flat; colorbar; hold on
    caxis([-3 3]);colorbar
    line([0 0],[0 1],'Color','k')
    num_trans = size(mp_utrig_mp_mat,1);
    for i = 1:num_trans
        cur_spikes = find((spkid > t_50_wup{d}(i)-maxlag) & (spkid < t_50_wup{d}(i)+maxlag));
        cur_spikes = spkid(cur_spikes) - t_50_wup{d}(i);
        plot(cur_spikes/Fs,ones(size(cur_spikes))*i/num_trans,'k.')
    end
    xlim([-0.2 0.5])
    cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
    t_names = ['C:\WC_Germany\Cortical_analysis\spike_timing\mp_utrig_lf8_' cell_name];
    print('-dpng',t_names);
    close

end

