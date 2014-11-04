clear all
close all

%%
cd C:\WC_Germany\sven_thomas_combined\
load ./distal_dir.mat

raw_Fs = 2016;
params.Fs = raw_Fs;
params.fpass = [0 100];
params.tapers = [3 5];
win = 50;

lcf_hf = 15;
hcf_hf = 80;
hcf_sm = 0.025;

rate_sm = round(raw_Fs*0.05);

freqs = linspace(0.05,80,5000);
%%
for d = 1:length(distal_dir)
    cd(distal_dir{d})
    pwd
        load ./used_data
    
        load ./used_data lf2
        hpc_hf = get_hf_features(lf2,raw_Fs,raw_Fs,[lcf_hf hcf_hf],hcf_sm);

    [desynch_times,desynch_inds,P_lf8,f,t] = locate_desynch_times_individual_v2(lf8);
    %compute markers indicating segments of data to be used
    t_axis = (1:length(lf8))/raw_Fs;
    if ~isempty(desynch_times)
        desynch_start = round(interp1(t_axis,1:length(t_axis),desynch_times(:,1)));
        desynch_stop = round(interp1(t_axis,1:length(t_axis),desynch_times(:,2)));
    else
        desynch_start = [];
        desynch_stop = [];
    end
    desynch_ind = zeros(size(lf8));
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
        synch_stops = [synch_stops; length(lf8)];
    end
    sMarkers = [synch_starts(:) synch_stops(:)];
    
    [S(d,1,:),f]=mtspectrumc_unequal_length_trials(wcv_minus_spike,[win win],params,sMarkers);
    [S(d,8,:),f]=mtspectrumc_unequal_length_trials(lf8,[win win],params,sMarkers);
    [S(d,7,:),f]=mtspectrumc_unequal_length_trials(lf7,[win win],params,sMarkers);
    [S(d,6,:),f]=mtspectrumc_unequal_length_trials(lf6,[win win],params,sMarkers);
    [S(d,5,:),f]=mtspectrumc_unequal_length_trials(lf5,[win win],params,sMarkers);
    [S(d,4,:),f]=mtspectrumc_unequal_length_trials(lf4,[win win],params,sMarkers);
    [S(d,3,:),f]=mtspectrumc_unequal_length_trials(lf3,[win win],params,sMarkers);
    [S(d,2,:),f]=mtspectrumc_unequal_length_trials(lf2,[win win],params,sMarkers);
    [Sh(d,:),f]=mtspectrumc_unequal_length_trials(hpc_hf,[win win],params,sMarkers);
end

%%
cd C:\WC_Germany\sven_thomas_combined
save distal_spectra_fin_nd S* f
%%
n_recs = size(S,1);

uds_freqs = find(f > 0.2 & f < 1);
lS = real(10*log10(S));
[peak_S,peakf] = max(lS(:,:,uds_freqs),[],3);
peak_uds_freq = f(uds_freqs(peakf));

high_freqs = find(f >= 15);
hg_pow = trapz(lS(:,:,high_freqs),3);

ctx_lS = squeeze(lS(:,8,:));
hpc_lS = nan(size(ctx_lS));
for i = 1:n_recs
        hpc_lS(i,:) = squeeze(lS(i,2,:));
end

%%
mp_S = squeeze(S(:,1,:));
for i = 1:size(S,1)
ctx_S(i,:) = squeeze(S(i,ctx_lfp(i),:));
end
l3mec_m = l3mec(~isnan(hpc_mua(l3mec)));
l3mec_mlfp = l3mec(~isnan(hpc_lfp(l3mec)));
lmp_S = log10(mp_S);
lhpchf_S = log10(Sh);
lctx_S = log10(ctx_S);
lSmua = log10(Smua);
figure; set(gca,'fontname','arial','fontsize',14)
plot(f,mean(mp_S(l3mec,:)),'r','linewidth',2)
hold on
plot(f,mean(mp_S(l3lec,:)),'b','linewidth',2)
plot(f,mean(ctx_S),'color',[0.2 0.2 0.2],'linewidth',2)
% plot(f,mean(Smua(l3mec_m,:)),'color',[0.2 0.8 0.2],'linewidth',2)
legend('MEC MP','LEC MP','Ctx LFP','Hpc MUA')
shadedErrorBar(f,mean(ctx_S),std(ctx_S)/sqrt(length(combined_dir)),{'color',[0.2 0.2 0.2]});
shadedErrorBar(f,mean(mp_S(l3lec,:)),std(mp_S(l3lec,:))/sqrt(length(l3lec)),{'b'});
shadedErrorBar(f,mean(mp_S(l3mec,:)),std(mp_S(l3mec,:))/sqrt(length(l3mec)),{'r'});
% shadedErrorBar(f,mean(Smua(l3mec_m,:)),std(Smua(l3mec_m,:))/sqrt(length(l3mec_m)),{'color',[0.2 0.8 0.2]});
xlim([0 1.])
xlabel('Frequency (Hz)','fontsize',16)
ylabel('Relative Power','fontsize',16)

uds_freqs = find(f > 0.1 & f < 0.8);
[mp_p,mp_f] = max(mp_S(:,uds_freqs),[],2);
[lfp_p,lfp_f] = max(ctx_S(:,uds_freqs),[],2);
[mua_p,mua_f] = max(Smua(:,uds_freqs),[],2);
mp_m = mean(mp_S(:,uds_freqs),2);
lfp_m = mean(ctx_S(:,uds_freqs),2);
mua_m = mean(Smua(:,uds_freqs),2);
ctx_uds_freq = f(uds_freqs(lfp_f));
mp_uds_freq = f(uds_freqs(mp_f));
mua_uds_freq = f(uds_freqs(mua_f));

noise_freqs = find(f < 0.01);
ctx_offset = sum(ctx_S(:,noise_freqs),2);


%%
df = f(2)-f(1);
norm_mp_S = bsxfun(@rdivide,mp_S,trapz(f,mp_S,2));
norm_ctx_S = bsxfun(@rdivide,ctx_S,trapz(f,ctx_S,2));
norm_hpc_mua = bsxfun(@rdivide,Smua,trapz(f,Smua,2));
norm_hpc_lfp = bsxfun(@rdivide,Sh,trapz(f,Sh,2));
figure; set(gca,'fontname','arial','fontsize',14)
hold on
shadedErrorBar(f,mean(norm_ctx_S),std(norm_ctx_S)/sqrt(length([l3lec l3mec])),{'color',[0.2 0.2 0.2]});
shadedErrorBar(f,mean(norm_mp_S(l3lec,:)),std(norm_mp_S(l3lec,:))/sqrt(length(l3lec)),{'b'});
shadedErrorBar(f,mean(norm_mp_S(l3mec,:)),std(norm_mp_S(l3mec,:))/sqrt(length(l3mec)),{'r'});
shadedErrorBar(f,mean(norm_hpc_mua(l3mec_m,:)),std(norm_hpc_mua(l3mec_m,:))/sqrt(length(l3mec_m)),{'g'});
% shadedErrorBar(f,mean(norm_hpc_lfp(l3mec_mlfp,:)),std(norm_hpc_lfp(l3mec_mlfp,:))/sqrt(length(l3mec_mlfp)),{'c'});
xlim([0 1.])
xlabel('Frequency (Hz)','fontsize',16)
ylabel('Relative Power','fontsize',16)


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

new_set = l3mec(23:end);
old_set = [l3mec(1:22) l3lec];
figure
h = errorbar(2:8,nanmean(peak_S(new_set,2:8)),nanstd(peak_S(new_set,2:8))/sqrt(length(new_set)),'k');
hold on
h = errorbar(2:8,nanmean(peak_S(old_set,2:8)),nanstd(peak_S(old_set,2:8))/sqrt(length(old_set)),'r');
legend('Sven','Thomas')
xlabel('Electrode Number','fontsize',14)
ylabel('Relative power','fontsize',14)

figure
h = errorbar(2:8,nanmean(rel_peak_S(new_set,2:8)),nanstd(rel_peak_S(new_set,2:8))/sqrt(length(new_set)),'k');
hold on
h = errorbar(2:8,nanmean(rel_peak_S(old_set,2:8)),nanstd(rel_peak_S(old_set,2:8))/sqrt(length(old_set)),'r');
legend('Sven','Thomas')
xlabel('Electrode Number','fontsize',14)
ylabel('Relative power','fontsize',14)

good_hpc_recs = new_set(~isnan(hpc_mua(poss_set)));
bad_hpc_recs = new_set(isnan(hpc_mua(poss_set)));
figure
h = errorbar(2:8,nanmean(peak_S(new_set,2:8)),nanstd(peak_S(new_set,2:8))/sqrt(length(new_set)),'k');
hold on
h = errorbar(2:8,nanmean(peak_S(good_hpc_recs,2:8)),nanstd(peak_S(good_hpc_recs,2:8))/sqrt(length(good_hpc_recs)));
h = errorbar(2:8,nanmean(peak_S(bad_hpc_recs,2:8)),nanstd(peak_S(bad_hpc_recs,2:8))/sqrt(length(bad_hpc_recs)),'r');
xlabel('Electrode Number','fontsize',16)
ylabel('Log UDS Power','fontsize',16)

figure
h = errorbar(2:8,nanmean(peak_uds_freq(:,2:8)),nanstd(peak_uds_freq(:,2:8))/sqrt(n_recs));
xlabel('Electrode Number','fontsize',16)
ylabel('Peak UDS Frequency (Hz)','fontsize',16)

%% VISUALIZE ALL LFP SPECTRA
% cmap = colormap(jet(7));
for i = 1:n_recs
    fprintf('Rec %d\n',i);
    temp = squeeze(lS(i,2:8,:));
    imagesc(f,2:8,temp);caxis([-90 -70]); colorbar
% for i = 1:7
%     plot(f,temp(i,:),'color',cmap(i,:),'linewidth',2); hold on
% end
xlim([0 2])
    pause
    clf
end
