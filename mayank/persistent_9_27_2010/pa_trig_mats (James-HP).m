clear all
close all

load F:\WC_Germany\overall_EC\overall_EC_dir
addpath('F:\WC_Germany\parietal_cortical_2010\')
addpath('F:\WC_Germany\hsmm_state_detection\')

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

drive_letter = 'F';

Fs = 2016;
niqf = Fs/2;
lcf = 0.05/niqf;
hcf = 10/niqf;
[b,a] = butter(2,[lcf hcf]);
lcf2 = 1/niqf;
hcf2 = 10/niqf;
[b2,a2] = butter(2,[lcf2 hcf2]);
lcf3 = 15/niqf;
hcf3 = 80/niqf;
[b3,a3] = butter(2,[lcf3 hcf3]);
dsf = 8;
Fsd = Fs/dsf;
pow_smooth = round(Fsd*0.05);

backlag = 4*Fsd;
forwardlag = 10*Fsd;
lags = (-backlag:forwardlag)/Fsd;


for d = 7:10
    
    cdir = sess_data(d).directory;
    cdir(1) = 'F';
    disp(sprintf('session %d',d))
    cd(cdir);
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);
    
    load ./used_data lf8 wcv_minus_spike lf3 lf2 lf5
    
    wcv_f = filtfilt(b,a,wcv_minus_spike);
    lf8_f = filtfilt(b,a,lf8);
    wcv_f = downsample(wcv_f,dsf)/sess_data(d).gains(1);
    lf8_f = downsample(lf8_f,dsf)/sess_data(d).gains(8);
    
    wcv_f = zscore(wcv_f);
    lf8_f = zscore(lf8_f);
    
    t_axis = (1:length(lf8_f))/Fsd;
    
    %% extract up and down transition times for MP and LF8
    load pa_hsmm_state_seq
    load pa_hsmm_state_seq8
    
    mp_state_seq_c =  hsmm_bbstate_seq;
    lf8_state_seq_c = hsmm_bbstate_seq8;
    
    [new_mp_seg_inds] = round(resample_uds_seg_inds(hmm.UDS_segs,50.4,Fsd,length(t_axis)));
    [new_lf8_seg_inds] = round(resample_uds_seg_inds(hmm8.UDS_segs,50.4,Fsd,length(t_axis)));
    
    mp_state_seq = nan(size(t_axis));
    lf8_state_seq = nan(size(t_axis));
    
    mp_utrans = [];
    mp_dtrans = [];
    for n = 1:hmm.Nsegs
        mp_state_seq(new_mp_seg_inds(n,1):new_mp_seg_inds(n,2)) = mp_state_seq_c{n};
        cur_mp_utrans = new_mp_seg_inds(n,1) + find(mp_state_seq_c{n}(1:end-1) == 1 & mp_state_seq_c{n}(2:end) == 2);
        cur_mp_dtrans = new_mp_seg_inds(n,1) + find(mp_state_seq_c{n}(1:end-1) == 2 & mp_state_seq_c{n}(2:end) == 1);
        cur_mp_dtrans(cur_mp_dtrans < cur_mp_utrans(1)) = [];
        cur_mp_utrans(cur_mp_utrans > cur_mp_dtrans(end)) = [];
        mp_utrans = [mp_utrans; cur_mp_utrans];
        mp_dtrans = [mp_dtrans; cur_mp_dtrans];
    end
    
    lf8_utrans = [];
    lf8_dtrans = [];
    for n = 1:hmm8.Nsegs
        lf8_state_seq(new_lf8_seg_inds(n,1):new_lf8_seg_inds(n,2)) = lf8_state_seq_c{n};
        cur_lf8_utrans = new_lf8_seg_inds(n,1) + find(lf8_state_seq_c{n}(1:end-1) == 1 & lf8_state_seq_c{n}(2:end) == 2);
        cur_lf8_dtrans = new_lf8_seg_inds(n,1) + find(lf8_state_seq_c{n}(1:end-1) == 2 & lf8_state_seq_c{n}(2:end) == 1);
        cur_lf8_dtrans(cur_lf8_dtrans < cur_lf8_utrans(1)) = [];
        cur_lf8_utrans(cur_lf8_utrans > cur_lf8_dtrans(end)) = [];
        lf8_utrans = [lf8_utrans; cur_lf8_utrans];
        lf8_dtrans = [lf8_dtrans; cur_lf8_dtrans];
    end
    
    n_mp_ups = length(mp_utrans);
    n_lf8_ups = length(lf8_utrans);
    
    %% initialize
    mp_utrig_mp_mat = nan(n_mp_ups,length(lags));
    mp_utrig_lf8_mat = nan(n_mp_ups,length(lags));
    mp_dtrig_mp_mat = nan(n_mp_ups,length(lags));
    mp_dtrig_lf8_mat = nan(n_mp_ups,length(lags));
    
    %     calculate mp utrigs
    for i = 1:n_mp_ups
        if mp_utrans(i) > backlag && length(wcv_f) - mp_utrans(i) > forwardlag
            mp_utrig_mp_mat(i,:) = wcv_f(mp_utrans(i)-backlag:mp_utrans(i)+forwardlag);
            mp_utrig_lf8_mat(i,:) = lf8_f(mp_utrans(i)-backlag:mp_utrans(i)+forwardlag);
        end
    end
    
%     calculate mp dtrigs
    for i = 1:n_mp_ups
        if mp_dtrans(i) > backlag && length(wcv_f) - mp_dtrans(i) > forwardlag
            mp_dtrig_mp_mat(i,:) = wcv_f(mp_dtrans(i)-backlag:mp_dtrans(i)+forwardlag);
            mp_dtrig_lf8_mat(i,:) = lf8_f(mp_dtrans(i)-backlag:mp_dtrans(i)+forwardlag);
        end
    end
    
    mp_utrig_mp(d,:) = nanmean(mp_utrig_mp_mat);
    mp_utrig_lf8(d,:) = nanmean(mp_utrig_lf8_mat);
    mp_dtrig_mp(d,:) = nanmean(mp_dtrig_mp_mat);
    mp_dtrig_lf8(d,:) = nanmean(mp_dtrig_lf8_mat);
    
    mp_updur = (mp_dtrans-mp_utrans)/Fsd;
    mp_downdur = (mp_utrans(2:end)-mp_dtrans(1:end-1))/Fsd;
    
    [dummy,up_order] = sort(mp_updur);
    
    % plot mp up trig matrices
    Fig= figure('visible','off');
    set(Fig,'PaperUnits','centimeters');
    set(Fig, 'PaperSize', [15 15]);
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))]);
    subplot(2,1,1)
    pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_utrig_mp_mat(up_order,:));
    shading flat; colorbar; hold on
    plot(mp_updur(up_order),(1:n_mp_ups)/n_mp_ups,'w','linewidth',2)
    caxis([-2.3 2.]);colorbar
    xlim([-2 8])
    line([0 0],[0 1],'Color','k')
    xlabel('Time since MP up-transition (s)','fontsize',14)
    ylabel('Up state number','fontsize',14)
    title('MP')
    subplot(2,1,2)
    pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_utrig_lf8_mat(up_order,:));
    shading flat; colorbar; hold on
    plot(mp_updur(up_order),(1:n_mp_ups)/n_mp_ups,'w','linewidth',2)
    caxis([-1.5 2.5]);colorbar
    xlim([-2 8])
    line([0 0],[0 1],'Color','k')
    xlabel('Time since MP up-transition (s)','fontsize',14)
    ylabel('Up state number','fontsize',14)
    title('Cortical LFP')
    t_names = ['F:\WC_Germany\persistent_9_27_2010\trig_mats2\mp_utrig_' s_name];
    print('-dpng',t_names);
    close

        % plot mp up trig matrices
    Fig= figure('visible','off');
    set(Fig,'PaperUnits','centimeters');
    set(Fig, 'PaperSize', [15 15]);
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))]);
    subplot(2,1,1)
    pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_utrig_mp_mat);
    shading flat; colorbar; hold on
    caxis([-2.3 2.]);colorbar
    xlim([-2 8])
    line([0 0],[0 1],'Color','k')
    xlabel('Time since MP up-transition (s)','fontsize',14)
    ylabel('Up state number','fontsize',14)
    title('MP')
    subplot(2,1,2)
    pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_utrig_lf8_mat);
    shading flat; colorbar; 
    caxis([-1.5 2.5]);colorbar
    xlim([-2 8])
    line([0 0],[0 1],'Color','k')
    xlabel('Time since MP up-transition (s)','fontsize',14)
    ylabel('Up state number','fontsize',14)
    title('Cortical LFP')
    t_names = ['F:\WC_Germany\persistent_9_27_2010\trig_mats2\mp_utrig_un_' s_name];
    print('-dpng',t_names);
    close

    [dummy,down_order] = sort(mp_downdur);
    
    Fig= figure('visible','off');
    set(Fig,'PaperUnits','centimeters');
    set(Fig, 'PaperSize', [15 15]);
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))]);
    subplot(2,1,1)
    pcolor(lags,(1:length(down_order))/length(down_order),mp_dtrig_mp_mat(down_order,:));
    shading flat; colorbar; hold on
    plot(mp_downdur(down_order),(1:length(down_order))/length(down_order),'w','linewidth',2)
    caxis([-2.3 2.]);
    xlim([-4 4])
    line([0 0],[0 1],'Color','k')
    xlabel('Time since MP down-transition (s)','fontsize',14)
    ylabel('Down state number','fontsize',14)
    title('MP')
    subplot(2,1,2)
    pcolor(lags,(1:length(down_order))/length(down_order),mp_dtrig_lf8_mat(down_order,:));
    shading flat; colorbar; hold on
    plot(mp_downdur(down_order),(1:length(down_order))/length(down_order),'w','linewidth',2)
    caxis([-1.5 2.5]);
    xlim([-4 4])
    line([0 0],[0 1],'Color','k')
    xlabel('Time since MP down-transition (s)','fontsize',14)
    ylabel('Down state number','fontsize',14)
    title('Cortical LFP')
    t_names = ['F:\WC_Germany\persistent_9_27_2010\trig_mats2\mp_dtrig_' s_name];
    print('-dpng',t_names);
    close

    Fig= figure('visible','off');
    set(Fig,'PaperUnits','centimeters');
    set(Fig, 'PaperSize', [15 15]);
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))]);
    subplot(2,1,1)
    pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_dtrig_mp_mat);
    shading flat; colorbar;
    caxis([-2.3 2.]);colorbar
    xlim([-4 4])
    line([0 0],[0 1],'Color','k')
    xlabel('Time since MP down-transition (s)','fontsize',14)
    ylabel('Down state number','fontsize',14)
    title('MP')
    subplot(2,1,2)
    pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_dtrig_lf8_mat);
    shading flat; colorbar;
    caxis([-1.5 2.5]);colorbar
    xlim([-4 4])
    line([0 0],[0 1],'Color','k')
    xlabel('Time since MP down-transition (s)','fontsize',14)
    ylabel('Down state number','fontsize',14)
    title('Cortical LFP')
    t_names = ['F:\WC_Germany\persistent_9_27_2010\trig_mats2\mp_dtrig_un_' s_name];
    print('-dpng',t_names);
    close
    
    
    %% For LF8 triggered analysis initialize
    lf8_utrig_mp_mat = nan(n_lf8_ups,length(lags));
    lf8_utrig_lf8_mat = nan(n_lf8_ups,length(lags));
    lf8_dtrig_mp_mat = nan(n_lf8_ups,length(lags));
    lf8_dtrig_lf8_mat = nan(n_lf8_ups,length(lags));
    
    %calculate mp utrigs
    for i = 1:n_lf8_ups
        if lf8_utrans(i) > backlag && length(wcv_f) - lf8_utrans(i) > forwardlag
            lf8_utrig_mp_mat(i,:) = wcv_f(lf8_utrans(i)-backlag:lf8_utrans(i)+forwardlag);
            lf8_utrig_lf8_mat(i,:) = lf8_f(lf8_utrans(i)-backlag:lf8_utrans(i)+forwardlag);
        end
    end
    
    %calculate mp dtrigs
    for i = 1:n_lf8_ups
        if lf8_dtrans(i) > backlag && length(wcv_f) - lf8_dtrans(i) > forwardlag
            lf8_dtrig_mp_mat(i,:) = wcv_f(lf8_dtrans(i)-backlag:lf8_dtrans(i)+forwardlag);
            lf8_dtrig_lf8_mat(i,:) = lf8_f(lf8_dtrans(i)-backlag:lf8_dtrans(i)+forwardlag);
        end
    end
    
    lf8_utrig_mp(d,:) = nanmean(lf8_utrig_mp_mat);
    lf8_utrig_lf8(d,:) = nanmean(lf8_utrig_lf8_mat);
    lf8_dtrig_mp(d,:) = nanmean(lf8_dtrig_mp_mat);
    lf8_dtrig_lf8(d,:) = nanmean(lf8_dtrig_lf8_mat);
    
    lf8_updur = (lf8_dtrans-lf8_utrans)/Fsd;
    lf8_downdur = (lf8_utrans(2:end)-lf8_dtrans(1:end-1))/Fsd;
    
    
    %% plot lf8 up trig matrices
    [dummy,up_order] = sort(lf8_updur);

    figure('visible','off')
    subplot(2,1,1)
    pcolor(lags,(1:n_lf8_ups)/n_lf8_ups,lf8_utrig_lf8_mat(up_order,:));
    shading flat; colorbar; hold on
    plot(lf8_updur(up_order),(1:n_lf8_ups)/n_lf8_ups,'w','linewidth',2)
    caxis([-1.5 2.5]);colorbar
    xlim([-2 3])
    line([0 0],[0 1],'Color','k')
    subplot(2,1,2)
    pcolor(lags,(1:n_lf8_ups)/n_lf8_ups,lf8_utrig_mp_mat(up_order,:));
    shading flat; colorbar; hold on
    plot(lf8_updur(up_order),(1:n_lf8_ups)/n_lf8_ups,'w','linewidth',2)
    caxis([-2.3 2]);colorbar
    xlim([-2 3])
    line([0 0],[0 1],'Color','k')
    t_names = ['F:\WC_Germany\persistent_9_27_2010\trig_mats2\lf8_utrig_' s_name];
    print('-dpng',t_names);
    close

    figure('visible','off')
    subplot(2,1,1)
    pcolor(lags,(1:n_lf8_ups)/n_lf8_ups,lf8_utrig_lf8_mat);
    shading flat; colorbar; hold on
    caxis([-1.5 2.5]);colorbar
    xlim([-2 3])
    line([0 0],[0 1],'Color','k')
    subplot(2,1,2)
    pcolor(lags,(1:n_lf8_ups)/n_lf8_ups,lf8_utrig_mp_mat);
    shading flat; colorbar; hold on
    caxis([-2.3 2]);colorbar
    xlim([-2 3])
    line([0 0],[0 1],'Color','k')
    t_names = ['F:\WC_Germany\persistent_9_27_2010\trig_mats2\lf8_utrig_un_' s_name];
    print('-dpng',t_names);
    close

    [dummy,down_order] = sort(lf8_downdur);
    figure('visible','off')
    subplot(2,1,1)
    pcolor(lags,(1:length(down_order))/length(down_order),lf8_dtrig_lf8_mat(down_order,:));
    shading flat; colorbar; hold on
    plot(lf8_downdur(down_order),(1:length(down_order))/length(down_order),'w','linewidth',2)
    caxis([-1.5 2.5]);colorbar
    xlim([-3 3])
    line([0 0],[0 1],'Color','k')
    subplot(2,1,2)
    pcolor(lags,(1:length(down_order))/length(down_order),lf8_dtrig_mp_mat(down_order,:));
    shading flat; colorbar; hold on
    plot(lf8_downdur(down_order),(1:length(down_order))/length(down_order),'w','linewidth',2)
    caxis([-2.3 2]);colorbar
    xlim([-3 3])
    line([0 0],[0 1],'Color','k')
    t_names = ['F:\WC_Germany\persistent_9_27_2010\trig_mats2\lf8_dtrig_' s_name];
    print('-dpng',t_names);
    close
    
    figure('visible','off')
    subplot(2,1,1)
    pcolor(lags,(1:n_lf8_ups)/n_lf8_ups,lf8_dtrig_lf8_mat);
    shading flat; colorbar;
    caxis([-1.5 2.5]);colorbar
    xlim([-3 3])
    line([0 0],[0 1],'Color','k')
    subplot(2,1,2)
    pcolor(lags,(1:n_lf8_ups)/n_lf8_ups,lf8_dtrig_mp_mat);
    shading flat; colorbar;
    caxis([-2.3 2]);colorbar
    xlim([-3 3])
    line([0 0],[0 1],'Color','k')
    t_names = ['F:\WC_Germany\persistent_9_27_2010\trig_mats2\lf8_dtrig_un_' s_name];
    print('-dpng',t_names);
    close
    
end

