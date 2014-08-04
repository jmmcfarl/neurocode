clear all
close all
cd C:\WC_Germany\Cortical_analysis

load C:\WC_Germany\Cortical_analysis\cortical_dir
load C:\WC_Germany\Cortical_analysis\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\Cortical_analysis\UDS_synch_state_dur\UDS_synch_state_dur_data
load C:\WC_Germany\Cortical_analysis\sigmoid_fit\sig_fit_alldata

Fs = 2016;

maxlag = 0.5;
maxlagid = maxlag*Fs;
xrange = linspace(-maxlag,maxlag,100);

max_up_timing_diff = 1*Fs;

for d = 1:length(sess_data)


    cd(sess_data(d).directory)
    disp(num2str(d))

    load spike_time_jmm
    t_50_wup{d} = rlid_wup{d};
    t_50_lup{d} = rlid_lup{d};


    %initialize
    tL10_tM10{d} = zeros(1,length(t_50_lup{d}));
    tL10_tFS{d} = zeros(1,length(t_50_lup{d}));
    L10_rel_spike_prob(d,:) = zeros(length(xrange),1);

    num_used_ups = 0;

    for i = 1:length(t_50_lup{d})

        if ~isnan(t_50_lup{d}(i))
            [dummy,nearest_mp_up] = min(abs(t_50_lup{d}(i) - t_50_wup{d}));
            if abs(t_50_wup{d}(nearest_mp_up) - t_50_lup{d}(i)) < max_up_timing_diff
                next_spike = find(spkid > t_50_wup{d}(nearest_mp_up),1,'first');
                tL10_tM10{d}(i) = t_10_lup{d}(i) - t_10_wup{d}(nearest_mp_up);
                used_spikes = spkid(spkid > t_10_wup{d}(nearest_mp_up));
                cur_spikes = find(used_spikes > t_10_lup{d}(i)-maxlagid & used_spikes < t_10_lup{d}(i) + maxlagid);
                cur_spikes = (t_10_lup{d}(i) - used_spikes(cur_spikes))/Fs;
                cur_spikes = cur_spikes';
                cur_hist = histc(cur_spikes,xrange);
                if ~isempty(cur_hist)
                    L10_rel_spike_prob(d,:) = L10_rel_spike_prob(d,:)+histc(cur_spikes,xrange);
                end
                num_used_ups = num_used_ups+1;
                if ~isempty(next_spike)
                    tL10_tFS{d}(i) = t_10_lup{d}(i) - spkid(next_spike);
                else
                    tL10_tFS{d}(i) = nan;
                end
            else
                tL10_tFS{d}(i) = nan;
                tL10_tM10{d}(i) = nan;
                tL10_tFS{d}(i) = nan;
            end
        else
            tL10_tM10{d}(i) = nan;
            tL10_tFS{d}(i) = nan;
        end

    end

    L10_rel_spike_prob(d,:) = L10_rel_spike_prob(d,:)/num_used_ups;

    tL10_tM10{d} = tL10_tM10{d}/Fs;
    tL10_tFS{d} = tL10_tFS{d}/Fs;


    n_tL10_tM10(d,:) = histc(tL10_tM10{d},xrange);
    n_tL10_tM10(d,:) = n_tL10_tM10(d,:)/sum(n_tL10_tM10(d,:));

end

