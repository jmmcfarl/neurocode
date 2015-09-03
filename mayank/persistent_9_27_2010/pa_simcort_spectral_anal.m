clear all

load('F:\WC_Germany\persistent_9_27_2010\pa_simcort_dir.mat')
addpath('F:\WC_Germany\parietal_cortical_2010\')
addpath('F:\WC_Germany\hsmm_state_detection\')
addpath('F:\WC_Germany\persistent_2010\')
addpath('F:\Code\smoothing\software\')
addpath('F:\Code\general\')

dsf = 8;
Fsd = 2016/dsf;
params.Fs = Fsd;
params.err = [2 .05];
params.fpass = [0 45];
window = [60 60];
niqf = 2016/2;
hcf = 100/niqf;
lcf = 0.05/niqf;
[b1,a1] = butter(2,[lcf hcf]);

for d = 1:length(dir_array)
    
    cdir = dir_array{d};
    cdir(1) = 'F';
    disp(sprintf('session %d',d))
    cd(cdir);
    s_name = f_names{d};
    
    load ./used_data lf8 lf4
    load ./desynch_times_lf8
    
    %bandlimit signals
    down_w = filtfilt(b1,a1,lf4);
    down_8 = filtfilt(b1,a1,lf8);
    
    down_w = downsample(down_w,dsf);
    down_8 = downsample(down_8,dsf);
        
    %zscore
    zdown_w = zscore(down_w);
    zdown_8 = zscore(down_8);
    
    %compute markers indicating segments of data to be used
    t_axis = (1:length(down_w))/Fsd;
    if ~isempty(desynch_times_lf8)
        desynch_start = round(interp1(t_axis,1:length(t_axis),desynch_times_lf8(:,1)));
        desynch_stop = round(interp1(t_axis,1:length(t_axis),desynch_times_lf8(:,2)));
    else
        desynch_start = [];
        desynch_stop = [];
    end
    desynch_ind = zeros(size(down_w));
    for i = 1:length(desynch_start)
        desynch_ind(desynch_start(i):desynch_stop(i)) = 1;
    end
    synch_starts = find(desynch_ind(1:end-1)==1 & desynch_ind(2:end)==0)+1;
    if desynch_ind(1) == 0
        synch_starts = [1; synch_starts];
    end
    synch_stops = find(desynch_ind(1:end-1)==0 & desynch_ind(2:end)==1)+1;
    if desynch_ind(end) == 0
        synch_stops = [synch_stops; length(down_w)];
    end
    sMarkers = [synch_starts(:) synch_stops(:)];
    
    %
    params.tapers = [5 9];
    [Cmn(d,:),Phimn(d,:),Smn,Smm,f,ConfC(d),PhiStd,Cerr(d,:,:)] = ...
        coherencyc_unequal_length_trials([down_w down_8],window, params, sMarkers);
    
    params.tapers = [3 5];
    [Pww(d,:), f, ~]= mtspectrumc_unequal_length_trials(down_w, window, params, sMarkers );
    [P88(d,:), f, ~]= mtspectrumc_unequal_length_trials(down_8, window, params, sMarkers );
    [zPww(d,:), f, ~]= mtspectrumc_unequal_length_trials(zdown_w, window, params, sMarkers );
    [zP88(d,:), f, ~]= mtspectrumc_unequal_length_trials(zdown_8, window, params, sMarkers );
    
    clear down_w down_8 lf4 lf8
    
end

cd F:\WC_Germany\persistent_9_27_2010\
save simcort_spectral_data Cmn Pww P88 zPww zP88 f Phimn Cerr ConfC