close all
clear all

drive_letter = 'G';
addpath(strcat(drive_letter,':\Code\wavelet_tools\'))
save_dir = strcat(drive_letter,':\WC_Germany\overall_EC\ripple_analysis\');
cd(strcat(drive_letter,':\WC_Germany\overall_EC\'))
load overall_EC_dir

%%
dsf = 2;
Fs = 2016;
Fsd = Fs/dsf;
niqf = Fs/2;

spw_lcf = 1;
spw_hcf = 40;
[b_spw,a_spw] = butter(2,[spw_lcf spw_hcf]/niqf);

rip_lcf = 40;
rip_hcf = 250;
rip_sm = round(Fsd*.040);

[b_rip,a_rip] = butter(2,[rip_lcf rip_hcf]/niqf);

bb_lcf = 0.5;
bb_hcf = 250;
[b_bb,a_bb] = butter(2,[bb_lcf bb_hcf]/niqf);

spw_win = round(Fsd*0.075);
maxlag = spw_win;
lags = -maxlag:maxlag;

%%
for d = 1:length(sess_data)
    cdir = sess_data(d).directory;
    cdir(1) = 'G';
    disp(sprintf('session %d',d))
    cd(cdir);
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);
    
    
    load used_data lf3 lf2
    lf3_spw = zscore(downsample(filtfilt(b_spw,a_spw,lf3),dsf));
    lf3_rip = zscore(downsample(filtfilt(b_rip,a_rip,lf3),dsf));
    lf3_bb = zscore(downsample(filtfilt(b_bb,a_bb,lf3),dsf));
    lf3_rip_pow = zscore(abs(hilbert(lf3_rip)));
    [lf3_rip_amps,lf3_rip_locs] = findpeaks(lf3_rip_pow,'minpeakheight',6,'minpeakdistance',round(0.1*Fsd));
    lf2_spw = zscore(downsample(filtfilt(b_spw,a_spw,lf2),dsf));
    lf2_rip = zscore(downsample(filtfilt(b_rip,a_rip,lf2),dsf));
    lf2_bb = zscore(downsample(filtfilt(b_bb,a_bb,lf2),dsf));
    lf2_rip_pow = zscore(abs(hilbert(lf2_rip)));
    [lf2_rip_amps,lf2_rip_locs] = findpeaks(lf2_rip_pow,'minpeakheight',6,'minpeakdistance',round(0.1*Fsd));
    
    bad_rips = find(lf3_rip_locs < maxlag | lf3_rip_locs > length(lf3_rip)-maxlag);
    lf3_rip_locs(bad_rips) = [];
    lf3_rip_amps(bad_rips) = [];
    bad_rips = find(lf2_rip_locs < maxlag | lf2_rip_locs > length(lf2_rip)-maxlag);
    lf2_rip_locs(bad_rips) = [];
    lf2_rip_amps(bad_rips) = [];
    
    n_lf3_rips = length(lf3_rip_locs);
    n_lf2_rips = length(lf2_rip_locs);
    
    t_axis = (1:length(lf3_bb))/Fsd;
    
    lf3_spw_amp = nan(n_lf3_rips,1);
    lf3_spw_mat = nan(n_lf3_rips,length(lags));
    for i = 1:n_lf3_rips
        cur_win = (lf3_rip_locs(i)-spw_win):(lf3_rip_locs(i)+spw_win);
        [~,cur_spw_maxloc] = max(abs(lf3_spw(cur_win)));
        lf3_spw_amp(i) = lf3_spw(cur_win(cur_spw_maxloc));
        lf3_spw_mat(i,:) = lf3_spw(cur_win);
    end
    
    lf2_spw_amp = nan(n_lf2_rips,1);
    lf2_spw_mat = nan(n_lf2_rips,length(lags));
    for i = 1:n_lf2_rips
        cur_win = (lf2_rip_locs(i)-spw_win):(lf2_rip_locs(i)+spw_win);
        [~,cur_spw_maxloc] = max(abs(lf2_spw(cur_win)));
        lf2_spw_amp(i) = lf2_spw(cur_win(cur_spw_maxloc));
        lf2_spw_mat(i,:) = lf2_spw(cur_win);
    end
    
    lf3_med_rip_amp(d) = median(lf3_rip_amps);
    lf3_med_spw_amp(d) = median(lf3_spw_amp);
    lf2_med_rip_amp(d) = median(lf2_rip_amps);
    lf2_med_spw_amp(d) = median(lf2_spw_amp);
    
    %%
    figure('visible','off')
    [lf3_rip_dens,cur_axis] = ksdensity(lf3_rip);
    plot(cur_axis,lf3_rip_dens), hold on
    pred_dens = 1/sqrt(2*pi)*exp(-cur_axis.^2/2);
    plot(cur_axis,pred_dens,'r')
    [lf2_rip_dens,cur_axis] = ksdensity(lf2_rip);
    plot(cur_axis,lf2_rip_dens,'k')
    set(gca,'yscale','log')
    yl = ylim();
    ylim([1e-6 yl(2)])
    legend('LF3','Gaussian','LF2')
    line([0 0],yl,'color','k')
    savename = strcat(save_dir,'rip_amp_',s_name);
    print('-dpng',savename), close
    
    figure('visible','off')
    [lf3_spw_dens,cur_axis] = ksdensity(lf3_spw);
    plot(cur_axis,lf3_spw_dens), hold on
    pred_dens = 1/sqrt(2*pi)*exp(-cur_axis.^2/2);
    plot(cur_axis,pred_dens,'r')
    [lf2_spw_dens,cur_axis] = ksdensity(lf2_spw);
    plot(cur_axis,lf2_spw_dens,'k')
    set(gca,'yscale','log')
    yl = ylim();
    ylim([1e-6 yl(2)])
    legend('LF3','Gaussian','LF2')
    line([0 0],yl,'color','k')
    savename = strcat(save_dir,'spw_amp_',s_name);
    print('-dpng',savename), close
    
    figure('visible','off')
    plot(lf3_rip_amps,lf3_spw_amp,'o'), hold on
    plot(lf2_rip_amps,lf2_spw_amp,'ko')
    legend('LF3','LF2')
    xlabel('Ripple amplitude (z)')
    ylabel('SPW amplitude (z)')
    savename = strcat(save_dir,'scatter_',s_name);
    print('-dpng',savename), close
    
    
    figure('visible','off')
    if n_lf3_rips > 1
    errorbar(lags/Fsd,nanmean(lf3_spw_mat),nanstd(lf3_spw_mat)/sqrt(length(lf3_rip_locs))), hold on
    end
    if n_lf2_rips > 1
    errorbar(lags/Fsd,nanmean(lf2_spw_mat),nanstd(lf2_spw_mat)/sqrt(length(lf2_rip_locs)),'r'), hold on
    end
    legend('LF3','LF2')
    xlabel('Lag (s)')
    xlim([-spw_win spw_win]/Fsd)
    ylabel('SPW amplitude (z)')
    yl = ylim();
    line([0 0],yl,'color','k')
    title(sprintf('%d LF3 Rips  %d LF2 Rips',n_lf3_rips,n_lf2_rips))
    savename = strcat(save_dir,'avg_spw_',s_name);
    print('-dpng',savename), close
    
end
drive_letter = 'G';
cd(strcat(drive_letter,':\WC_Germany\overall_EC\'))
save spw_rip_data lf3_med_rip_amp lf3_med_spw_amp lf2_med_rip_amp lf2_med_spw_amp

