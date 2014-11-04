clear all
close all

%%
cd C:\WC_Germany\sven_thomas_combined\
load ./combined_dir.mat

raw_Fs = 2016;
win = 25;

lcf = 0.05;
hcf = 100;
dsf = 8;
Fsd = raw_Fs/dsf;
params.Fs = Fsd;
params.fpass = [0 100];
params.tapers = [3 5];

rate_sm = round(raw_Fs*0.05);

freqs = linspace(0.05,80,5000);
%%
for d = 1:length(combined_dir)
    cd(combined_dir{d})
    pwd
    load ./used_data 
    
     if ctx_lfp(d) == 7
        lf8 = lf7;
    end
    [lf8_lf,t_axis] = get_lf_features(lf8,raw_Fs,Fsd,[lcf hcf]);
    wcv_lf = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[lcf hcf]);

    [desynch_times,desynch_inds,P_lf8,f,t] = locate_desynch_times_individual(lf8);
    %compute markers indicating segments of data to be used
    if ~isempty(desynch_times)
        desynch_start = round(interp1(t_axis,1:length(t_axis),desynch_times(:,1)));
        desynch_stop = round(interp1(t_axis,1:length(t_axis),desynch_times(:,2)));
    else
        desynch_start = [];
        desynch_stop = [];
    end
    desynch_ind = zeros(size(lf8_lf));
    for i = 1:length(desynch_start)
        desynch_ind(desynch_start(i):desynch_stop(i)) = 1;
    end
    desynch_fract(d) = sum(desynch_ind)/length(desynch_ind);
    
    synch_starts = find(desynch_ind(1:end-1)==1 & desynch_ind(2:end)==0)+1;
    if desynch_ind(1) == 0
        synch_starts = [1; synch_starts];
    end
    synch_stops = find(desynch_ind(1:end-1)==0 & desynch_ind(2:end)==1)+1;
    if desynch_ind(end) == 0
        synch_stops = [synch_stops; length(lf8_lf)];
    end
    sMarkers = [synch_starts(:) synch_stops(:)];
        
    if ~isnan(hpc_mua(d))
        load ./mua_data
        load ./sync_times.mat
        hpc_mua_times = mua_times{hpc_mua(d)};
        temp_mua_rate =hist(hpc_mua_times,synct)*raw_Fs;
        mua_rate = jmm_smooth_1d_cor(temp_mua_rate,rate_sm);
        mua_rate = downsample(mua_rate,dsf);
        mua_rate = zscore(mua_rate);
        if length(mua_rate) > length(t_axis)
            mua_rate = mua_rate(1:length(t_axis));
        end
    end
        
    [Sw(d,:),f]= mtspectrumc_unequal_length_trials(wcv_lf,[win win],params,sMarkers);
    [S8(d,:),f]=mtspectrumc_unequal_length_trials(lf8_lf,[win win],params,sMarkers);
    if ~isnan(hpc_mua(d))
    [Smua(d,:),f]=mtspectrumc_unequal_length_trials(mua_rate(:),[win win],params,sMarkers);
    else
        Smua(d,:) = nan(1,length(f));
    end
end

%%
n_recs = size(Sw,1);
n_mec = length(l3mec);
n_lec = length(l3lec);

uds_freqs = find(f > 0.2 & f < 1);
cd C:\WC_Germany\sven_thomas_combined\
save combined_spectra_short_z S* f

%%
l3mec_m = l3mec(~isnan(hpc_mua(l3mec)));

figure; set(gca,'fontname','arial','fontsize',14)
plot(f,mean(Sw(l3mec,:)),'r','linewidth',2)
hold on
plot(f,mean(Sw(l3lec,:)),'b','linewidth',2)
plot(f,mean(S8),'color',[0.2 0.2 0.2],'linewidth',2)
plot(f,mean(Smua(l3mec_m,:)),'color',[0.2 0.8 0.2],'linewidth',2)
legend('MEC MP','LEC MP','Ctx LFP','Hpc MUA')
shadedErrorBar(f,mean(S8),std(S8)/sqrt(length(combined_dir)),{'color',[0.2 0.2 0.2]});
shadedErrorBar(f,mean(Sw(l3lec,:)),std(Sw(l3lec,:))/sqrt(length(l3lec)),{'b'});
shadedErrorBar(f,mean(Sw(l3mec,:)),std(Sw(l3mec,:))/sqrt(length(l3mec)),{'r'});
shadedErrorBar(f,mean(Smua(l3mec_m,:)),std(Smua(l3mec_m,:))/sqrt(length(l3mec_m)),{'color',[0.2 0.8 0.2]});
xlim([0 1.])
xlabel('Frequency (Hz)','fontsize',16)
ylabel('Relative Power','fontsize',16)

uds_freqs = find(f > 0.1 & f < 0.8);
[mp_p,mp_f] = max(Sw(:,uds_freqs),[],2);
[lfp_p,lfp_f] = max(S8(:,uds_freqs),[],2);
[mua_p,mua_f] = max(Smua(:,uds_freqs),[],2);
mp_m = mean(Sw(:,uds_freqs),2);
lfp_m = mean(S8(:,uds_freqs),2);
mua_m = mean(Smua(:,uds_freqs),2);
%% PLOT AVERAGE SPECTRA
n_recs = size(S,1);
n_mec = length(l3mec);
n_lec = length(l3lec);
green = [0.1 0.9 0.2];
xl = [0 2];

m_ctx_lS = mean(ctx_lS);
u_ctx_lS = m_ctx_lS + std(ctx_lS)/sqrt(n_recs);
l_ctx_lS = m_ctx_lS - std(ctx_lS)/sqrt(n_recs);
m_hpc_lS = mean(hpc_lS(~isnan(hpc_lfp),:));
u_hpc_lS = m_hpc_lS + std(hpc_lS(~isnan(hpc_lfp),:))/sqrt(sum(~isnan(hpc_lfp)));
l_hpc_lS = m_hpc_lS - std(hpc_lS(~isnan(hpc_lfp),:))/sqrt(sum(~isnan(hpc_lfp)));

shadedErrorBar(f,squeeze(mean(lS(l3mec,1,:))),squeeze(std(lS(l3mec,1,:)))/sqrt(n_mec),{'r'});
hold on
shadedErrorBar(f,squeeze(mean(lS(l3lec,1,:))),squeeze(std(lS(l3lec,1,:)))/sqrt(n_lec),{'b'});
xlim(xl)
ylabel('Log MP Power','fontsize',16)
ax1 = gca; set(gca,'box','off')
ax2 = axes('Position',get(ax1,'Position'),...
    'XAxisLocation','bottom',...
    'YAxisLocation','right',...
    'Color','none',...
    'XColor','k','YColor','k');
line(f,m_ctx_lS,'color',[0.2 0.2 0.2],'linewidth',2,'parent',ax2);
line(f,m_hpc_lS,'color',[0.2 0.8 0.2],'linewidth',2,'parent',ax2);
legend('Ctx LFP','Hpc LFP')
line(f,u_ctx_lS,'color',[0.2 0.2 0.2],'linestyle','--','parent',ax2);
line(f,l_ctx_lS,'color',[0.2 0.2 0.2],'linestyle','--','parent',ax2);
line(f,u_hpc_lS,'color',[0.2 0.8 0.2],'linestyle','--','parent',ax2);
line(f,l_hpc_lS,'color',[0.2 0.8 0.2],'linestyle','--','parent',ax2);
xlim(xl)

xlabel('Frequency (Hz)','fontsize',16)
ylabel('Log LFP Power','fontsize',16)

%% PLOT AVG RELATIVE POWER VS DEPTH
rel_peak_S = bsxfun(@minus,peak_S,peak_S(:,8));
figure
h = errorbar(2:8,nanmean(peak_S(:,2:8)),nanstd(peak_S(:,2:8))/sqrt(n_recs));
xlabel('Electrode Number','fontsize',16)
ylabel('Log UDS Power','fontsize',16)

figure
h = errorbar(2:8,nanmean(peak_uds_freq(:,2:8)),nanstd(peak_uds_freq(:,2:8))/sqrt(n_recs));
xlabel('Electrode Number','fontsize',16)
ylabel('Peak UDS Frequency (Hz)','fontsize',16)

%% VISUALIZE ALL LFP SPECTRA
cmap = colormap(jet(7));
for i = 1:n_recs
    fprintf('Rec %d\n',i);
    temp = real(squeeze(lS(i,2:8,:)));
    imagesc(f,2:8,temp); caxis([-90 -70]); colorbar
% for i = 1:7
%     plot(f,temp(i,:),'color',cmap(i,:),'linewidth',2); hold on
% end
xlim([0 2])
    pause
    clf
end
