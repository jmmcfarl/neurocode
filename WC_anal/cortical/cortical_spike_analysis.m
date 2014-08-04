clear all
close all
cd C:\WC_Germany\Cortical_analysis
load cortical_dir

Fs = 2016;

for d = 1:length(sess_data)

    cd(sess_data(d).directory)
    pwd


    load spike_time_jmm

    num_spikes = length(spkid);
    if num_spikes > 0
        if length(spkdur) > num_spikes
            spkdur(num_spikes+1:end) = [];
            spkamp(num_spikes+1:end) = [];
        end

        spkdur = spkdur*1000;
        isis = [0;diff(spkid)]/Fs;
        isis2 = [0;isis(1:end-1)];
        [n,x] = log_hist_non_norm(isis,[0.001 2],100);
        n(end) = 0;
        bar(x,n)
        set(gca,'xscale','log')
        xlim([0.001 2])
        cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
        t_names = ['C:\WC_Germany\Cortical_analysis\spike_analysis\isi_hist_' cell_name];
        print('-dpng',t_names);
        close

        [n,x] = hist(spkamp,100);
        figure
        bar(x,n)
        xlabel('Spike Amplitude (AU)')
        cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
        t_names = ['C:\WC_Germany\Cortical_analysis\spike_analysis\amp_hist_' cell_name];
        print('-dpng',t_names);
        close

        [n,x] = hist(spkdur,100);
        bar(x,n)
        xlabel('Spike duration (ms)')
        cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
        t_names = ['C:\WC_Germany\Cortical_analysis\spike_analysis\duration_hist_' cell_name];
        print('-dpng',t_names);
        close

        figure
        plot(isis,spkamp,'.')
        set(gca,'xscale','log')
        cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
        t_names = ['C:\WC_Germany\Cortical_analysis\spike_analysis\isi_v_amp_' cell_name];
        print('-dpng',t_names);
        close

        figure
        plot(spkamp,spkdur,'.')
        cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
        t_names = ['C:\WC_Germany\Cortical_analysis\spike_analysis\amp_v_dur_' cell_name];
        print('-dpng',t_names);
        close

        figure
        plot(isis,spkdur,'.')
        set(gca,'xscale','log')
        cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
        t_names = ['C:\WC_Germany\Cortical_analysis\spike_analysis\isi_v_dur_' cell_name];
        print('-dpng',t_names);
        close
    end
    clear isis spkdur spkamp n x
end