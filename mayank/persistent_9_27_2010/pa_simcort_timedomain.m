clear all

load('C:\WC_Germany\persistent_9_27_2010\pa_simcort_dir.mat')
addpath('C:\WC_Germany\parietal_cortical_2010\')
addpath('C:\WC_Germany\hsmm_state_detection\')
addpath('C:\WC_Germany\persistent_2010\')
addpath('C:\Code\smoothing\software\')
addpath('C:\Code\general\')

dsf = 8;
Fsd = 2016/dsf;
minSegLength = 60;
maxLag = 10*Fsd;
niqf = 2016/2;
lcf = .05/niqf;
hcf = 40/niqf;
[b,a] = butter(2,[lcf hcf]);

for d = 1:length(dir_array)
    
    cdir = dir_array{d};
    cdir(1) = 'C';
    disp(sprintf('session %d',d))
    cd(cdir);
    s_name = f_names{d};
    
    load ./used_data lf8 lf4
    load ./desynch_times_lf8
    
    %bandlimit signals
    down_w = filtfilt(b,a,lf4);
    down_8 = filtfilt(b,a,lf8);
    
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
    
    seg_durs = diff(sMarkers')'/Fsd;
    too_short = find(seg_durs < minSegLength);
    sMarkers(too_short,:) = [];
    seg_durs(too_short) = [];
    
    cnt = 0;
    for i = 1:size(sMarkers,1)
        cnt = cnt+1;
        
        [wcv_acorr(cnt,:),lags] = xcov(down_w(sMarkers(i,1):sMarkers(i,2)),maxLag,'coeff');
        [lf8_acorr(cnt,:),lags] = xcov(down_8(sMarkers(i,1):sMarkers(i,2)),maxLag,'coeff');
        
        [w8_x(cnt,:),lags] = xcov(down_w(sMarkers(i,1):sMarkers(i,2)),down_8(sMarkers(i,1):sMarkers(i,2)),maxLag,'coeff');
    end
    
    %compute weighted averages (weighted by relative duration
    weighting = seg_durs/sum(seg_durs);
    tot_wcv_acorr(d,:) = sum(wcv_acorr.*repmat(weighting(:),1,length(lags)),1);
    tot_lf8_acorr(d,:) = sum(lf8_acorr.*repmat(weighting(:),1,length(lags)),1);
    tot_w8_x(d,:) = sum(w8_x.*repmat(weighting(:),1,length(lags)),1);
    
    clear down_w down_8 lf4 lf8 wcv_acorr lf8_acorr w8_x
    
end

cd C:\WC_Germany\persistent_9_27_2010\
save simcort_timedomain_data lags tot*