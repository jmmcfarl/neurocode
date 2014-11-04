clear all
close all

cur_dir_letter = 'C';
addpath('C:\WC_Germany\persistent_9_27_2010\')
cd C:\WC_Germany\sven_thomas_combined\
load ./combined_dir.mat
addpath('C:\WC_Germany\new_mec\')
uset = sort([l3mec l3lec]);
all_cells = 1:61;
l3mec = find(ismember(all_cells(uset),l3mec));
l3lec = find(ismember(all_cells(uset),l3lec));
combined_dir = combined_dir(uset);
hpc_mua = hpc_mua(uset);
hpc_lfp = hpc_lfp(uset);
ctx_lfp = ctx_lfp(uset);

%%
raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
niqf = raw_Fs/2;
hcf = 10;
lcf = 0.05;
lcf_hf = 15;
hcf_hf = 80;
hcf_sm = 0.05;
rate_sm = round(Fsd*0.1);

backlag = 4*Fsd;
forwardlag = 10*Fsd;
lags = (-backlag:forwardlag)/Fsd;

for d = 59
% d = 61;
cur_dir = combined_dir{d};
    cd(cur_dir)
    pwd
    last_slash = find(combined_dir{d} == '\',1,'last');
    s_name = combined_dir{d}(last_slash+1:end);
    
    load ./used_data lf8 wcv_minus_spike lf7
    if ctx_lfp(d) == 7
        lf8 = lf7;
    end
    [lf8_lf,t_axis] = get_lf_features(lf8,raw_Fs,Fsd,[lcf hcf]);
    wcv_lf = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[lcf hcf]);
%     if hpc_lfp(d) == 3
%         load ./used_data lf3 lf5
%         if ismember(d,old_data_inds)
%             lf3 = lf3 + lf5; %redefine LF3 wrt gnd
%         end
%         hpc_lf = get_lf_features(lf3,raw_Fs,Fsd,[lcf hcf]);
%         hpc_hf = get_hf_features(lf3,raw_Fs,Fsd,[lcf_hf hcf_hf],hcf_sm);
%     else
%         load ./used_data lf2
%         hpc_lf = get_lf_features(lf2,raw_Fs,Fsd,[lcf hcf]);
%         hpc_hf = get_hf_features(lf2,raw_Fs,Fsd,[lcf_hf hcf_hf],hcf_sm);
%     end
    if ~isnan(hpc_mua(d))
        load ./mua_data2
        load ./sync_times.mat
        synct_d = downsample(synct,dsf);
        hpc_mua_times = mua_times{hpc_mua(d)};
          hpc_mua_times(hpc_mua_times > synct_d(end) | hpc_mua_times < synct_d(1)) = [];
        hpc_mua_rate =hist(hpc_mua_times,synct_d)*Fsd;
        hpc_mua_rate = jmm_smooth_1d_cor(hpc_mua_rate,rate_sm);
        mua_rate = zscore(hpc_mua_rate');
        
%         [log_mua,offset(d)] = log_transform_sig(mua_rate);
%         mua_rate = zscore(log_mua);
        
        if length(mua_rate) > length(t_axis)
            mua_rate = mua_rate(1:length(t_axis));
        end
    else
        mua_rate = nan(size(lf8_lf));
    end
    

    %% extract up and down transition times for MP and LF8
    load ./pa_hsmm_state_seq_combined_fin.mat
    load ./pa_hsmm_state_seq7_combined_fin.mat
    hsmm_bbstate_seq8 = hsmm_bbstate_seq7;
    
    [new_seg_inds] = resample_uds_seg_inds_v2(hmm.UDS_segs,hmm.Fs,Fsd,hsmm_bbstate_seq);
    
    dur_uds = sum(diff(new_seg_inds,[],2))/Fsd;
    [mp_utrans,mp_dtrans] = compute_state_transitions_seg(new_seg_inds,hsmm_bbstate_seq);
    [lf8_utrans,lf8_dtrans] = compute_state_transitions_seg(new_seg_inds,hsmm_bbstate_seq8);
    
    n_mp_ups = length(mp_utrans);
    n_lf8_ups = length(lf8_utrans);
    
    %% initialize
    mp_utrig_mp_mat = nan(n_mp_ups,length(lags));
    mp_utrig_lf8_mat = nan(n_mp_ups,length(lags));
%     mp_utrig_hpch_mat = nan(n_mp_ups,length(lags));
    mp_utrig_mua_mat = nan(n_mp_ups,length(lags));
    mp_dtrig_mp_mat = nan(n_mp_ups,length(lags));
    mp_dtrig_lf8_mat = nan(n_mp_ups,length(lags));
%     mp_dtrig_hpch_mat = nan(n_mp_ups,length(lags));
    mp_dtrig_mua_mat = nan(n_mp_ups,length(lags));
    
    %     calculate mp utrigs
    for i = 1:n_mp_ups
        if mp_utrans(i) > backlag && length(wcv_lf) - mp_utrans(i) > forwardlag
            mp_utrig_mp_mat(i,:) = wcv_lf(mp_utrans(i)-backlag:mp_utrans(i)+forwardlag);
            mp_utrig_lf8_mat(i,:) = lf8_lf(mp_utrans(i)-backlag:mp_utrans(i)+forwardlag);
%             if ~isnan(hpc_lfp(d))
%                 mp_utrig_hpch_mat(i,:) = hpc_hf(mp_utrans(i)-backlag:mp_utrans(i)+forwardlag);
%             end
            if ~isnan(hpc_mua(d))
                 mp_utrig_mua_mat(i,:) = mua_rate(mp_utrans(i)-backlag:mp_utrans(i)+forwardlag);               
            end
        end
    end
    
    %     calculate mp dtrigs
    for i = 1:n_mp_ups
        if mp_dtrans(i) > backlag && length(wcv_lf) - mp_dtrans(i) > forwardlag
            mp_dtrig_mp_mat(i,:) = wcv_lf(mp_dtrans(i)-backlag:mp_dtrans(i)+forwardlag);
            mp_dtrig_lf8_mat(i,:) = lf8_lf(mp_dtrans(i)-backlag:mp_dtrans(i)+forwardlag);
%             if ~isnan(hpc_lfp(d))
%                 mp_dtrig_hpch_mat(i,:) = hpc_hf(mp_dtrans(i)-backlag:mp_dtrans(i)+forwardlag);
%             end
            if ~isnan(hpc_mua(d))
                mp_dtrig_mua_mat(i,:) = mua_rate(mp_dtrans(i)-backlag:mp_dtrans(i)+forwardlag);
            end
        end
    end
        
    mp_updur = (mp_dtrans-mp_utrans)/Fsd;
    mp_downdur = (mp_utrans(2:end)-mp_dtrans(1:end-1))/Fsd;  
    [dummy,up_order] = sort(mp_updur);
    
    %%
    % plot mp up trig matrices
%     Fig= figure('visible','off');
%     set(Fig,'PaperUnits','centimeters');
%     set(Fig, 'PaperSize', [15 15]);
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))]);
%     subplot(2,1,1)
%     pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_utrig_mp_mat(up_order,:));
%     shading flat; colorbar; hold on
%     plot(mp_updur(up_order),(1:n_mp_ups)/n_mp_ups,'w','linewidth',2)
%     caxis([-2.3 2.]);colorbar
%     xlim([-2 8])
%     line([0 0],[0 1],'Color','k')
%     xlabel('Time since MP up-transition (s)','fontsize',14)
%     ylabel('Up state number','fontsize',14)
%     title('MP')
%     subplot(2,1,2)
%     pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_utrig_lf8_mat(up_order,:));
%     shading flat; colorbar; hold on
%     plot(mp_updur(up_order),(1:n_mp_ups)/n_mp_ups,'w','linewidth',2)
%     caxis([-1.5 2.5]);colorbar
%     xlim([-2 8])
%     line([0 0],[0 1],'Color','k')
%     xlabel('Time since MP up-transition (s)','fontsize',14)
%     ylabel('Up state number','fontsize',14)
%     title('Cortical LFP')
%     t_names = ['C:\WC_Germany\sven_thomas_combined\trig_mats\mp_utrig_lf8_' s_name];
%     print('-dpng',t_names);
%     close

%     if ~isnan(hpc_lfp(d))
%     Fig= figure('visible','off');
%     set(Fig,'PaperUnits','centimeters');
%     set(Fig, 'PaperSize', [15 15]);
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))]);
%     subplot(2,1,1)
%     pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_utrig_mp_mat(up_order,:));
%     shading flat; colorbar; hold on
%     plot(mp_updur(up_order),(1:n_mp_ups)/n_mp_ups,'w','linewidth',2)
%     caxis([-2.3 2.]);colorbar
%     xlim([-2 8])
%     line([0 0],[0 1],'Color','k')
%     xlabel('Time since MP up-transition (s)','fontsize',14)
%     ylabel('Up state number','fontsize',14)
%     title('MP')
%     subplot(2,1,2)
%     pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_utrig_hpch_mat(up_order,:));
%     shading flat; colorbar; hold on
%     plot(mp_updur(up_order),(1:n_mp_ups)/n_mp_ups,'w','linewidth',2)
%     caxis([-1.5 1.5]);colorbar
%     xlim([-2 8])
%     line([0 0],[0 1],'Color','k')
%     xlabel('Time since MP up-transition (s)','fontsize',14)
%     ylabel('Up state number','fontsize',14)
%     title('Hpc LFP')
%     t_names = ['C:\WC_Germany\sven_thomas_combined\trig_mats\mp_utrig_hpch_' s_name];
%     print('-dpng',t_names);
%     close
%     end
    if ~isnan(hpc_mua(d))
    Fig= figure('visible','off');
    set(Fig,'PaperUnits','centimeters');
    set(Fig, 'PaperSize', [15 25]);
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))]);
    subplot(3,1,1)
    pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_utrig_mp_mat(up_order,:));
    shading flat; colorbar; hold on
    plot(mp_updur(up_order),(1:n_mp_ups)/n_mp_ups,'w','linewidth',2)
    caxis([-2.3 2.]);colorbar
    xlim([-1 7])
    line([0 0],[0 1],'Color','k')
    xlabel('Time since MP up-transition (s)','fontsize',14)
    ylabel('Up state number','fontsize',14)
    title('MP')
    subplot(3,1,2)
%     for ii = 1:n_mp_ups
%        cur_spks = hpc_mua_times(hpc_mua_times > synct_d(mp_utrans(ii)) - 8e6 & hpc_mua_times < synct_d(mp_utrans(ii)) + 8e6);
%        cur_spks = (cur_spks - synct_d(mp_utrans(ii)))/1e6;
%        plot(cur_spks,ones(size(cur_spks))*ii/n_mp_ups,'k.','markersize',0.1); hold on
%     end
    pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_utrig_mua_mat(up_order,:));
    shading flat; colorbar; hold on
    plot(mp_updur(up_order),(1:n_mp_ups)/n_mp_ups,'w','linewidth',2)
    caxis([-1. 1.]);colorbar
    set(gca,'yscale','log')
    xlim([-1 7])
    line([0 0],[0 1],'Color','k')
    xlabel('Time since MP up-transition (s)','fontsize',14)
    ylabel('Up state number','fontsize',14)
    subplot(3,1,3)
    pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_utrig_lf8_mat(up_order,:));
    shading flat; colorbar; hold on
    plot(mp_updur(up_order),(1:n_mp_ups)/n_mp_ups,'w','linewidth',2)
    caxis([-1.5 1.5]);colorbar
    xlim([-1 7])
    line([0 0],[0 1],'Color','k')
    xlabel('Time since MP up-transition (s)','fontsize',14)
    ylabel('Up state number','fontsize',14)
    title('Neocortical LFP')
    t_names = ['C:\WC_Germany\sven_thomas_combined\trig_mats\fin_mp_utrig_mua_' s_name];
    print('-dpng',t_names);
    close
    end
    
    %%
    [dummy,down_order] = sort(mp_downdur);
    
%     Fig= figure('visible','off');
%     set(Fig,'PaperUnits','centimeters');
%     set(Fig, 'PaperSize', [15 15]);
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))]);
%     subplot(2,1,1)
%     pcolor(lags,(1:length(down_order))/length(down_order),mp_dtrig_mp_mat(down_order,:));
%     shading flat; colorbar; hold on
%     plot(mp_downdur(down_order),(1:length(down_order))/length(down_order),'w','linewidth',2)
%     caxis([-2.3 2.]);
%     xlim([-4 4])
%     line([0 0],[0 1],'Color','k')
%     xlabel('Time since MP down-transition (s)','fontsize',14)
%     ylabel('Down state number','fontsize',14)
%     title('MP')
%     subplot(2,1,2)
%     pcolor(lags,(1:length(down_order))/length(down_order),mp_dtrig_lf8_mat(down_order,:));
%     shading flat; colorbar; hold on
%     plot(mp_downdur(down_order),(1:length(down_order))/length(down_order),'w','linewidth',2)
%     caxis([-1.5 2.5]);
%     xlim([-4 4])
%     line([0 0],[0 1],'Color','k')
%     xlabel('Time since MP down-transition (s)','fontsize',14)
%     ylabel('Down state number','fontsize',14)
%     title('Cortical LFP')
%     t_names = ['C:\WC_Germany\sven_thomas_combined\trig_mats\mp_dtrig_lf8_' s_name];
%     print('-dpng',t_names);
%     close

%     if ~isnan(hpc_lfp(d))
%     Fig= figure('visible','off');
%     set(Fig,'PaperUnits','centimeters');
%     set(Fig, 'PaperSize', [15 15]);
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))]);
%     subplot(2,1,1)
%     pcolor(lags,(1:length(down_order))/length(down_order),mp_dtrig_mp_mat(down_order,:));
%     shading flat; colorbar; hold on
%     plot(mp_downdur(down_order),(1:length(down_order))/length(down_order),'w','linewidth',2)
%     caxis([-2.3 2.]);
%     xlim([-4 4])
%     line([0 0],[0 1],'Color','k')
%     xlabel('Time since MP down-transition (s)','fontsize',14)
%     ylabel('Down state number','fontsize',14)
%     title('MP')
%     subplot(2,1,2)
%     pcolor(lags,(1:length(down_order))/length(down_order),mp_dtrig_hpch_mat(down_order,:));
%     shading flat; colorbar; hold on
%     plot(mp_downdur(down_order),(1:length(down_order))/length(down_order),'w','linewidth',2)
%     caxis([-1.5 2.5]);
%     xlim([-4 4])
%     line([0 0],[0 1],'Color','k')
%     xlabel('Time since MP down-transition (s)','fontsize',14)
%     ylabel('Down state number','fontsize',14)
%     title('Hpc LFP')
%     t_names = ['C:\WC_Germany\sven_thomas_combined\trig_mats\mp_dtrig_hpch_' s_name];
%     print('-dpng',t_names);
%     close
%     end
    if ~isnan(hpc_mua(d))
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
    pcolor(lags,(1:length(down_order))/length(down_order),mp_dtrig_mua_mat(down_order,:));
    shading flat; colorbar; hold on
    plot(mp_downdur(down_order),(1:length(down_order))/length(down_order),'w','linewidth',2)
    caxis([-1.5 2.5]);
    xlim([-4 4])
    line([0 0],[0 1],'Color','k')
    xlabel('Time since MP down-transition (s)','fontsize',14)
    ylabel('Down state number','fontsize',14)
    title('Hpc MUA')
    t_names = ['C:\WC_Germany\sven_thomas_combined\trig_mats\fin_mp_dtrig_mua_' s_name];
    print('-dpng',t_names);
    close
    end
end

