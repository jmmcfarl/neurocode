clear all
close all

load('C:\WC_Germany\persistent_revised\pers_revised_dir.mat')

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

desynch_thresh = -6;
mp_theta_thresh = -8;
removal_window = 6; %number of window slides

for d = 1:28

    disp(sprintf('session %d',d))
    cd(dir_array{d});
    load used_data wcv_minus_spike lf8 


    %bandlimit signals
    down_w = filtfilt(b,a,wcv_minus_spike);
    down_8 = filtfilt(b,a,lf8);
    down_wf = filtfilt(b2,a2,wcv_minus_spike);
    down_8f = filtfilt(b2,a2,lf8);

    down_w = downsample(down_w,dsf);
    down_8 = downsample(down_8,dsf);
    down_wf = downsample(down_wf,dsf);
    down_8f = downsample(down_8f,dsf);
    
    %zscore
    down_w = zscore(down_w);
    down_8 = zscore(down_8);
    down_wf = zscore(down_wf);
    down_8f = zscore(down_8f);
    
     [Pw{d},t{d},f{d}]=mtspecgramc(down_w,movingwin,params);
    [P8{d},t{d},f{d}]=mtspecgramc(down_8,movingwin,params);
   
    f_beg_so = find(f{d} > 0.2,1,'first');
    f_end_so = find(f{d} > 1.0,1,'first');
    f_theta_beg = find(f{d} > 1.5,1,'first');
    f_theta_end = find(f{d} > 4.5,1,'first');

    max_theta_w = max(10*log10(Pw{d}(:,f_theta_beg:f_theta_end)),[],2);
    max_so_w = max(10*log10(Pw{d}(:,f_beg_so:f_end_so)),[],2);
     pow_diff_w = max_theta_w - max_so_w;

         max_theta_8 = max(10*log10(P8{d}(:,f_theta_beg:f_theta_end)),[],2);
    max_so_8 = max(10*log10(P8{d}(:,f_beg_so:f_end_so)),[],2);
     pow_diff_8 = max_theta_8 - max_so_8;

     desynch_ind = zeros(size(pow_diff_w));
    desynch_ind(max_theta_w > mp_theta_thresh) = 1;
    desynch_diff = [0 diff(desynch_ind)'];
    desynch_start_ids = find(desynch_diff == 1);
    desynch_stop_ids = find(desynch_diff == -1);
 
         desynch_ind_8 = zeros(size(pow_diff_8));
    desynch_ind_8(pow_diff_8 > desynch_thresh) = 1;
    desynch_diff_8 = [0 diff(desynch_ind_8)'];
    desynch_start_ids_8 = find(desynch_diff_8 == 1);
    desynch_stop_ids_8 = find(desynch_diff_8 == -1);

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
    
  
    
        %make sure you start and stop in desynchronized epochs in correct order
    if ~isempty(desynch_start_ids_8)
        if isempty(desynch_stop_ids_8)
            desynch_stop_ids_8 = length(t{d});
        else
            if desynch_start_ids_8(1) > desynch_stop_ids_8(1)
                desynch_start_ids_8 = [1 desynch_start_ids_8];
            end
            if desynch_start_ids_8(end) > desynch_stop_ids_8(end)
                desynch_stop_ids_8 = [desynch_stop_ids_8 length(t{d})];
            end
        end
    end
     
    if length(desynch_start_ids_8) ~= length(desynch_stop_ids_8)
        disp('error start and stop not equal!')
    end
    
    %now make a window around desynchronized times
    for w = 1:length(desynch_start_ids_8)
        if desynch_start_ids_8(w) <= removal_window
           desynch_start_ids_8(w) = 1; 
        else
            desynch_start_ids_8(w) = desynch_start_ids_8(w)-removal_window;
        end
        
        if length(t{d})-desynch_stop_ids_8(w) <= removal_window
            desynch_stop_ids_8(w) = length(t{d});
        else
            desynch_stop_ids_8(w) = desynch_stop_ids_8(w)+removal_window;
        end
    end
     
    %now make sure there are no overlapping windows
    bad_desynch_start = [];
    for w = 2:length(desynch_start_ids_8)
        if desynch_start_ids_8(w) < desynch_stop_ids_8(w-1)
            bad_desynch_start = [bad_desynch_start w];           
        end
    end
    
    desynch_start_ids_8(bad_desynch_start) = [];
    desynch_stop_ids_8(bad_desynch_start-1) = [];
    
    desynch_start_times_8{d} = t{d}(desynch_start_ids_8);
    desynch_stop_times_8{d} = t{d}(desynch_stop_ids_8);

    
    
    
    Fig = figure(1);
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    subplot(2,1,1)
    pcolor(t{d},f{d},10*log10(Pw{d}'));shading flat;
    caxis([-25 2]);
    ylim([0 7])
    subplot(2,1,2)
%     pcolor(At{d},lags/Fsd,A8{d}');shading flat; colorbar
%     ylim([0 maxLag/Fsd])
%     caxis([-0.4 0.4])
plot(t{d},max_theta_w,'r')
hold on
plot(t{d},max_so_w,'k')
plot(t{d},max_theta_w-max_so_w)
plot(t{d}(desynch_start_ids),ones(size(desynch_start_ids))*-10,'go')
plot(t{d}(desynch_stop_ids),ones(size(desynch_stop_ids))*-10,'ro')
xlim([t{d}(1) t{d}(end)])
    tname = ['C:\WC_Germany\persistent_revised\desynch_detect_2\mp_' f_names{d}];
    print('-dpng',tname);
    close

 
    
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
plot(t{d}(desynch_start_ids_8),ones(size(desynch_start_ids_8))*-10,'go')
plot(t{d}(desynch_stop_ids_8),ones(size(desynch_stop_ids_8))*-10,'ro')
xlim([t{d}(1) t{d}(end)])
    tname = ['C:\WC_Germany\persistent_revised\desynch_detect_2\lf8_' f_names{d}];
    print('-dpng',tname);
    close


end

cd C:\WC_Germany\persistent_revised\desynch_detect_2
save desynch_times_2 desynch*times*
