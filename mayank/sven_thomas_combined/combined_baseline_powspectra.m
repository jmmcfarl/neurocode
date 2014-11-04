clear all
close all

%%
cd C:\WC_Germany\sven_thomas_combined\
load ./combined_dir.mat

raw_Fs = 2016;
params.Fs = raw_Fs;
params.fpass = [0 100];
params.tapers = [3 5];
win = 25;

rate_sm = round(raw_Fs*0.05);

freqs = linspace(0.05,80,5000);
%%
for d = 1:length(combined_dir)
    cd(combined_dir{d})
    pwd
    if exist('./used_data2.mat','file')
        load ./used_data2
    else
        load ./used_data
    end
        
    %for old data, create approximate LF3 wrt ground
    if ismember(d,old_data_inds)
        lf3 = lf3 + lf5;
    end
    
    if ~isnan(hpc_mua(d))
        load ./mua_data
        load ./sync_times.mat
        hpc_mua_times = mua_times{hpc_mua(d)};
        temp_mua_rate =hist(hpc_mua_times,synct)*raw_Fs;
        mua_rate = jmm_smooth_1d_cor(temp_mua_rate,rate_sm);
        mua_rate = zscore(mua_rate);
        if length(mua_rate) > length(lf8)
            mua_rate = mua_rate(1:length(t_axis));
        end
    end
        
    [S(d,1,:),f]=mtspectrumsegc(wcv_minus_spike,win,params);
    [S(d,8,:),f]=mtspectrumsegc(lf8,win,params);
    [S(d,7,:),f]=mtspectrumsegc(lf7,win,params);
    [S(d,6,:),f]=mtspectrumsegc(lf6,win,params);
    [S(d,5,:),f]=mtspectrumsegc(lf5,win,params);
    [S(d,4,:),f]=mtspectrumsegc(lf4,win,params);
    [S(d,3,:),f]=mtspectrumsegc(lf3,win,params);
    [S(d,2,:),f]=mtspectrumsegc(lf2,win,params);
    if ~isnan(hpc_mua(d))
    [Smua(d,:),f]=mtspectrumsegc(mua_rate(:),win,params);
    else
        Smua(d,:) = nan(1,length(f));
    end
end

%%
n_recs = size(S,1);
n_mec = length(l3mec);
n_lec = length(l3lec);

uds_freqs = find(f > 0.2 & f < 1);
lS = real(10*log10(S));
[peak_S,peakf] = max(lS(:,:,uds_freqs),[],3);
peak_uds_freq = f(uds_freqs(peakf));

high_freqs = find(f >= 15);
hg_pow = trapz(lS(:,:,high_freqs),3);
cd C:\WC_Germany\sven_thomas_combined\
save combined_baseline_spectra S* f

%%
mp_S = squeeze(S(:,1,:));
ctx_S = squeeze(S(:,7,:));
l3mec_m = l3mec(~isnan(hpc_mua(l3mec)));
% lmp_S = log10(mp_S);
% lctx_S = log10(ctx_S);
% lSmua = log10(Smua);
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
