clear all

load G:\WC_Germany\overall_EC\overall_EC_dir
addpath('G:\WC_Germany\parietal_cortical_2010\')
addpath('G:\WC_Germany\hsmm_state_detection\')
drive_letter = 'G';

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);
Fsd = 2016/8;

% Calculate distributions of up and down states
up_range = [0.1 10];
down_range = [0.1 10];
% up_range = [0.1 3];
% down_range = [0.1 3];
numBins = 60;

lin_grid = linspace(up_range(1),up_range(2),numBins+1);
log_grid = logspace(log10(up_range(1)),log10(up_range(2)),numBins+1);

for d = 1:length(sess_data)
    
    fprintf('session %d\n',d)
    cdir = sess_data(d).directory;
    cdir(1) = drive_letter;
    cd(cdir)
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);
    
    load ./pa_hsmm_state_seq_new2.mat
    load ./pa_hsmm_state_seq8_new2.mat
    
    mp_state_seq =  hsmm_bbstate_seq;
    lf8_state_seq = hsmm_bbstate_seq8;
    
    [mp_state_durations{d}] = compute_state_durations_seg(mp_state_seq,Fsd);
    [lfp_state_durations{d}] = compute_state_durations_seg(lf8_state_seq,Fsd);
     
    up_hist(d,:) = histc(mp_state_durations{d}{2},lin_grid);
    up_hist8(d,:) = histc(lfp_state_durations{d}{2},lin_grid);
    
    down_hist(d,:) = histc(mp_state_durations{d}{1},lin_grid);
    down_hist8(d,:) = histc(lfp_state_durations{d}{1},lin_grid);
    
    up_loghist(d,:) = histc(mp_state_durations{d}{2},log_grid);
    up_loghist8(d,:) = histc(lfp_state_durations{d}{2},log_grid);
    
    down_loghist(d,:) = histc(mp_state_durations{d}{1},log_grid);
    down_loghist8(d,:) = histc(lfp_state_durations{d}{1},log_grid);
    
    %%
    net_up_time(d) = nansum(mp_state_durations{d}{2});
    net_down_time(d) = nansum(mp_state_durations{d}{1});
    net_up_time8(d) = nansum(lfp_state_durations{d}{2});
    net_down_time8(d) = nansum(lfp_state_durations{d}{1});
    
    mean_up_dur(d) = nanmean(mp_state_durations{d}{2});
    mean_down_dur(d) = nanmean(mp_state_durations{d}{1});
    mean_up_dur8(d) = nanmean(lfp_state_durations{d}{2});
    mean_down_dur8(d) = nanmean(lfp_state_durations{d}{1});
 
    max_up_dur(d) = nanmax(mp_state_durations{d}{2});
    max_down_dur(d) = nanmax(mp_state_durations{d}{1});
    max_up_dur8(d) = nanmax(lfp_state_durations{d}{2});
    max_down_dur8(d) = nanmax(lfp_state_durations{d}{1});

    total_number_up_states(d) = length(mp_state_durations{d}{2});
    total_number_up_states8(d) = length(lfp_state_durations{d}{2});
    total_number_down_states(d) = length(mp_state_durations{d}{1});
    total_number_down_states8(d) = length(lfp_state_durations{d}{1});
    
end

cd G:\WC_Germany\persistent_9_27_2010
save pa_state_dur_stats_new2 mp_* lfp_* *hist* *time* *dur* total_* lin_grid log_grid