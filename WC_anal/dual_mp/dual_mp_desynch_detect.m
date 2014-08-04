clear all
close all
cd C:\WC_Germany\DUAL-MP_recordings
load dual_mp_dir

dsf = 8;
Fsd = 2016/dsf;

niqf = 2016/2;
hcf1 = 100/niqf;
[b,a] = butter(2,hcf1,'low');
[b2,a2] = butter(2,[0.1/niqf 10/niqf]);

params.Fs = Fsd;
params.err = 0;
params.tapers = [3 5];
params.fpass = [0 10];

winlength = 20;
winslide = 2.5;
movingwin = [winlength winslide];

desynch_thresh = -10;
removal_window = 6; %number of window slides

for d = 1:length(dir_array)

    disp(sprintf('session %d',d))
    cd(dir_array{d});
    load used_data lf8 

    %bandlimit signals
    down_8 = filtfilt(b,a,lf8);
    down_8f = filtfilt(b2,a2,lf8);

    down_8 = downsample(down_8,dsf);
    down_8f = downsample(down_8f,dsf);
    
    %zscore
    down_8 = zscore(down_8);
    down_8f = zscore(down_8f);
    
    [P8{d},t{d},f{d}]=mtspecgramc(down_8,movingwin,params);
   
    f_beg_so = find(f{d} > 0.2,1,'first');
    f_end_so = find(f{d} > 1.5,1,'first');
    f_theta_beg = find(f{d} > 3.5,1,'first');
    f_theta_end = find(f{d} > 6,1,'first');

    max_theta_8 = max(10*log10(P8{d}(:,f_theta_beg:f_theta_end)),[],2);
    max_so_8 = max(10*log10(P8{d}(:,f_beg_so:f_end_so)),[],2);
     pow_diff_8 = max_theta_8 - max_so_8;
    
     desynch_ind = zeros(size(pow_diff_8));
    desynch_ind(pow_diff_8 > desynch_thresh) = 1;
    desynch_diff = [0;diff(desynch_ind)];
    desynch_start_ids = find(desynch_diff == 1);
    desynch_stop_ids = find(desynch_diff == -1);
    
    %make sure you start and stop in desynchronized epochs in correct order
    if ~isempty(desynch_start_ids)
        if isempty(desynch_stop_ids)
            desynch_stop_ids = length(t{d});
        else
            if desynch_start_ids(1) > desynch_stop_ids(1)
                desynch_start_ids = [1 desynch_start_ids];
            end
            if desynch_start_ids(end) > desynch_stop_ids(end)
                desynch_stop_ids = [desynch_stop_ids length(t{d})];
            end
        end
    end
     
    if length(desynch_start_ids) ~= length(desynch_stop_ids)
        disp('error start and stop not equal!')
    end
    
    %now make a window around desynchronized times
    for w = 1:length(desynch_start_ids)
        if desynch_start_ids(w) <= removal_window
           desynch_start_ids(w) = 1; 
        else
            desynch_start_ids(w) = desynch_start_ids(w)-removal_window;
        end
        
        if length(t{d})-desynch_stop_ids(w) <= removal_window
            desynch_stop_ids(w) = length(t{d});
        else
            desynch_stop_ids(w) = desynch_stop_ids(w)+removal_window;
        end
    end
     
    %now make sure there are no overlapping windows
    bad_desynch_start = [];
    for w = 2:length(desynch_start_ids)
        if desynch_start_ids(w) < desynch_stop_ids(w-1)
            bad_desynch_start = [bad_desynch_start w];           
        end
    end
    
    desynch_start_ids(bad_desynch_start) = [];
    desynch_stop_ids(bad_desynch_start-1) = [];
    
    desynch_start_times{d} = t{d}(desynch_start_ids);
    desynch_stop_times{d} = t{d}(desynch_stop_ids);
    
    
    Fig = figure(1);
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    subplot(2,1,1)
    pcolor(t{d},f{d},10*log10(P8{d}'));shading flat;
    caxis([-25 2]);
    ylim([0 7])
    subplot(2,1,2)
%     pcolor(At{d},lags/Fsd,A8{d}');shading flat; colorbar
%     ylim([0 maxLag/Fsd])
%     caxis([-0.4 0.4])
plot(t{d},max_theta_8,'r')
hold on
plot(t{d},max_so_8,'k')
plot(t{d},max_theta_8-max_so_8)
plot(t{d}(desynch_start_ids),ones(size(desynch_start_ids))*-10,'go')
plot(t{d}(desynch_stop_ids),ones(size(desynch_stop_ids))*-10,'ro')
xlim([t{d}(1) t{d}(end)])
    tname = ['C:\WC_Germany\DUAL-MP_recordings\desynch_detect\' f_names{d}];
    print('-dpng',tname);
    close

end

cd C:\WC_Germany\DUAL-MP_recordings\desynch_detect
save desynch_times desynch*times
