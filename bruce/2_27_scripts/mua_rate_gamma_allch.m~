clear all
% close all

cd ~/Data/bruce/2_27_12/
load Blocks.mat
accept_window = [-5 5;-5 5];
sac_eyespeed = 10;
thresh_eyespeed = 2.5;
sac_buffer = 0;
min_fix_dur = 0.2;
blink_thresh = 5;
min_blink_dur = 0.05;
% nlags = 4;


Pix2Deg = 0.018837;

% down-sampling fraction for image
dsfrac = 2;
Fsd = 1/Pix2Deg/dsfrac;

Nyp = 1024;
Nxp = 1280;
Ny = round(Nyp/dsfrac);
Nx = round(Nxp/dsfrac);
xax = linspace(-Nx/2,Nx/2,Nx)/Fsd; yax = linspace(-Ny/2,Ny/2,Ny)/Fsd;
% [XX,YY] = meshgrid(xax,yax);

% RF_patch = [3.5 6; -3.5 -1]; %location of RFs in degrees [x1 x2;y1 y2]
% RF_patch = [2.5 6.5; -4 0]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch = [2.75 6.25; -3.5 0]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch_pix = Fsd*RF_patch;
RF_patch_cent = mean(RF_patch,2);
RF_patch_width = diff(RF_patch,[],2);

% patch_inds = find(XX >= RF_patch(1,1) & XX <= RF_patch(1,2) & YY >= RF_patch(2,1) & YY <= RF_patch(2,2));
xpatch_inds = find(xax >= RF_patch(1,1) & xax <= RF_patch(1,2));
ypatch_inds = find(yax >= RF_patch(2,1) & yax <= RF_patch(2,2));

cd /Users/James/Data/bruce/2_27_12/stimrecon
load fixation_data_v3

cd /Users/James/Data/bruce/2_27_12/stimrecon
load sing_eye_pgabortrack_d2
load sing_eye_pgabortrack_fin_est_d2
cd ~/Data/bruce/2_27_12/ExptA

SDIM = length(xpatch_inds);
cd /Users/James/Data/bruce/2_27_12/stimrecon
load fixation_image_patches_corrected
cd /Users/James/James_scripts/bruce/modelfits/
% load fixation_gabor_models
load gabor_tempmodfits_en_lin
%%
Fs = 1000;
dsf = 2;
Fsd = Fs/dsf;
niqf = Fs/2;
[b,a] = butter(2,[1 250]/niqf);

scales = logspace(log10(3),log10(70),50);
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);

%%
cd ~/Data/bruce/2_27_12

avg_lfp_amps = nan(24,length(all_fix_start_times),length(scales));
avgb_lfp_amps = nan(24,length(all_fix_start_times),length(scales));
avge_lfp_amps = nan(24,length(all_fix_start_times),length(scales));
net_spk_cnts = nan(length(all_fix_start_times),1);
netb_spk_cnts = nan(length(all_fix_start_times),1);
nete_spk_cnts = nan(length(all_fix_start_times),1);
spk_cnts = nan(length(all_fix_start_times),24);
beg_spk_cnts = nan(length(all_fix_start_times),24);
fin_spk_cnts = nan(length(all_fix_start_times),24);
used_pts = [];

ind_bins_edges = round(Fsd*(0:0.05:0.5));
ind_bins_cents = 0.5*ind_bins_edges(1:end-1) + 0.5*ind_bins_edges(2:end);
all_spk_phases = cell(24,1);
init_spk_phases = cell(24,1);
fin_spk_phases = cell(24,1);
bin_spk_phases = cell(24,length(ind_bins_edges)-1);
for blockid = 1:3;
    
    fprintf('Processing block %d of %d...\n',blockid,3);
    
    block_times = Blocks{blockid}.blocktimes;
    stim_times = Blocks{blockid}.stimtime;
    stimID = Blocks{blockid}.stimids;
    
    %%
    cd /Users/James/Data/bruce/2_27_12
    load(sprintf('lemM232A.5%d.lfp.mat',blockid));
    %get start times of each LFP trial
    n_lfp_trials = length(LFP.Trials);
    lfp_trial_start = nan(n_lfp_trials,1);
    lfp_trial_stop = nan(n_lfp_trials,1);
    for i = 1:n_lfp_trials
        lfp_trial_start(i) = LFP.Trials(i).Start/1e4; %originally in tenths of ms
        lfp_trial_stop(i) = LFP.Trials(i).End/1e4;
        lfp_dur(i) = size(LFP.Trials(i).LFP,1)/1000;
    end
    
    lfp_time = [];
    lfp_samps = [];
    for i = 1:n_lfp_trials
        lfp_time = [lfp_time linspace(lfp_trial_start(i),lfp_trial_start(i)+lfp_dur(i),size(LFP.Trials(i).LFP,1))];
        lfp_samps = [lfp_samps; LFP.Trials(i).LFP];
    end
    % lfp_samps = zscore(lfp_samps);
    
    lfp_samps = filtfilt(b,a,lfp_samps);
    
    lfp_sampsd = downsample(lfp_samps,dsf);
    lfp_sampsd = zscore(lfp_sampsd);
    lfp_timed = downsample(lfp_time,dsf);
    
    %%
    ampgram = zeros(24,length(lfp_timed),length(wfreqs));
    phasegram = zeros(24,length(lfp_timed),length(wfreqs));
    for cc = 1:24
        fprintf('Channel %d of %d\n',cc,24);
        temp = cwt(lfp_sampsd(:,cc),scales,'cmor1-1');
        tampgram = abs(temp);
        tampgram = bsxfun(@minus,tampgram,mean(tampgram,2));
        tampgram = bsxfun(@rdivide,tampgram,std(tampgram,[],2));
        ampgram(cc,:,:) = tampgram';
        
        phasegram(cc,:,:) = angle(temp');
    end
    
    %%
    spikes_binned = zeros(24,length(lfp_timed));
    for c = 1:10
        temp = round(interp1(lfp_timed,1:length(lfp_timed),Blocks{blockid}.spktimes{c}));
        temp(isnan(temp)) = [];
        all_spike_lfp_inds{c} = temp;
        
        cur_hist = hist(Blocks{blockid}.spktimes{c},lfp_timed);
        cur_hist(end) = 0;
        spikes_binned(c,:) = smooth(cur_hist,2);
    end
    for c = 1:14
        temp = round(interp1(lfp_timed,1:length(lfp_timed),Blocks{blockid}.mutimes{c}));
        temp(isnan(temp)) = [];
        all_mua_lfp_inds{c} = temp;
        
        cur_hist = hist(Blocks{blockid}.mutimes{c},lfp_timed);
        cur_hist(end) = 0;
        spikes_binned(c+10,:) = smooth(cur_hist,2);
    end
    
    spike_rates = spikes_binned'*Fsd;
    spike_rates = zscore(spike_rates)';
    net_rate = mean(spike_rates);
    spike_rates = spike_rates';
    %%
    
    buff_win = round(Fsd*0.15);
    temp_fix_set = find(blockids == blockid & all_fix_start_times > lfp_timed(1));
    cur_fix_start_times = all_fix_start_times(temp_fix_set);
    cur_fix_stop_times = all_fix_stop_times(temp_fix_set);
    cur_fix_start_inds = round(interp1(lfp_timed,1:length(lfp_timed),cur_fix_start_times));
    cur_fix_stop_inds = round(interp1(lfp_timed,1:length(lfp_timed),cur_fix_stop_times));
    bad_inds = find(isnan(cur_fix_start_inds) | isnan(cur_fix_stop_inds));
    
    nlags = round(0.2*Fsd);
    
    %     fix_trig_avg = zeros(24,nlags+1);
    temp_fix_set(bad_inds) = [];
    all_used_inds = [];
    all_init_inds = [];
    all_fin_inds = [];
    all_inds_bin = cell(length(ind_bins_edges)-1,1);
    for i = 1:length(temp_fix_set)
        start_ind = cur_fix_start_inds(i);
        end_ind = cur_fix_stop_inds(i);
        %         fix_trig_avg = fix_trig_avg + lfp_sampsd(start_ind:(start_ind+nlags),:)';
        cur_range1 = start_ind:end_ind;
        cur_range2 = (start_ind + buff_win):end_ind;
        cur_range3 = start_ind:(start_ind + buff_win);
        net_spk_cnts(temp_fix_set(i)) = mean(net_rate(cur_range1));
        netb_spk_cnts(temp_fix_set(i)) = mean(net_rate(cur_range3));
        nete_spk_cnts(temp_fix_set(i)) = mean(net_rate(cur_range2));
        
        avg_lfp_amps(:,temp_fix_set(i),:) = mean(ampgram(:,cur_range1,:),2);
        avgb_lfp_amps(:,temp_fix_set(i),:) = mean(ampgram(:,cur_range3,:),2);
        avge_lfp_amps(:,temp_fix_set(i),:) = mean(ampgram(:,cur_range2,:),2);
        
        spk_cnts(temp_fix_set(i),:) = mean(spike_rates(cur_range1,:));
        beg_spk_cnts(temp_fix_set(i),:) = mean(spike_rates(cur_range3,:));
        fin_spk_cnts(temp_fix_set(i),:) = mean(spike_rates(cur_range2,:));
        
        all_used_inds = [all_used_inds cur_range1];
        all_init_inds = [all_init_inds cur_range3];
        all_fin_inds = [all_fin_inds cur_range2];
        for ii = 1:length(ind_bins_edges)-1
            cur_range = start_ind + (ind_bins_edges(ii):ind_bins_edges(ii+1));
            all_inds_bin{ii} = [all_inds_bin{ii} cur_range];
        end
    end
    %     fix_trig_avg = fix_trig_avg/length(temp_fix_set);
    
%     for cellid = 1:24;
    for cellid = 1:10;
        if cellid > 10
            cur_spktimes = Blocks{blockid}.mutimes{cellid-10};
        else
            cur_spktimes = Blocks{blockid}.spktimes{cellid};
        end
        cur_spk_inds = round(interp1(lfp_timed,1:length(lfp_timed),cur_spktimes));
        cur_all_spk = cur_spk_inds(ismember(cur_spk_inds,all_used_inds));
        cur_init_spk = cur_spk_inds(ismember(cur_spk_inds,all_init_inds));
        cur_fin_spk = cur_spk_inds(ismember(cur_spk_inds,all_fin_inds));
        all_spk_phases{cellid} = cat(2,all_spk_phases{cellid},phasegram(:,cur_all_spk,:));
        init_spk_phases{cellid} = cat(2,init_spk_phases{cellid},phasegram(:,cur_init_spk,:));
        fin_spk_phases{cellid} = cat(2,fin_spk_phases{cellid},phasegram(:,cur_fin_spk,:));
        for ii = 1:length(ind_bins_edges)-1
            cur_spk = cur_spk_inds(ismember(cur_spk_inds,all_inds_bin{ii}));
            bin_spk_phases{cellid,ii} = cat(2,bin_spk_phases{cellid,ii},phasegram(:,cur_spk,:));
        end
    end
    
end
%%
cd /Users/James/James_scripts/bruce/modelfits
save lfp_spike_calcs_allch3 all_* init_* fin_* *spk_cnts avg*_lfp* wfreqs Fsd *spk_phases
%%
used_set = find(~isnan(net_spk_cnts) & ~isnan(squeeze(avge_lfp_amps(1,:,1))'));

%% Analyze relationship between net spiking and LFP power
for cc = 1:24
    anet_alfp(cc,:) = corr(net_spk_cnts(used_set),squeeze(avg_lfp_amps(cc,used_set,:)));
    anet_blfp(cc,:) = corr(net_spk_cnts(used_set),squeeze(avgb_lfp_amps(cc,used_set,:)));
    anet_elfp(cc,:) = corr(net_spk_cnts(used_set),squeeze(avge_lfp_amps(cc,used_set,:)));
    bnet_alfp(cc,:) = corr(netb_spk_cnts(used_set),squeeze(avg_lfp_amps(cc,used_set,:)));
    bnet_blfp(cc,:) = corr(netb_spk_cnts(used_set),squeeze(avgb_lfp_amps(cc,used_set,:)));
    bnet_elfp(cc,:) = corr(netb_spk_cnts(used_set),squeeze(avge_lfp_amps(cc,used_set,:)));
    enet_alfp(cc,:) = corr(nete_spk_cnts(used_set),squeeze(avg_lfp_amps(cc,used_set,:)));
    enet_blfp(cc,:) = corr(nete_spk_cnts(used_set),squeeze(avgb_lfp_amps(cc,used_set,:)));
    enet_elfp(cc,:) = corr(nete_spk_cnts(used_set),squeeze(avge_lfp_amps(cc,used_set,:)));
end

subplot(3,1,1)
% plot(wfreqs,anet_alfp)
hold on
plot(wfreqs,anet_blfp,'r')
plot(wfreqs,anet_elfp,'k')
xlim(wfreqs([end 1]))
% set(gca,'xscale','log')
title('Overall spiking','fontsize',14)

subplot(3,1,2)
% plot(wfreqs,bnet_alfp)
hold on
plot(wfreqs,bnet_blfp,'r')
plot(wfreqs,bnet_elfp,'k')
xlim(wfreqs([end 1]))
% set(gca,'xscale','log')
title('Initial spiking','fontsize',14)

subplot(3,1,3)
% plot(wfreqs,enet_alfp)
hold on
plot(wfreqs,enet_blfp,'r')
plot(wfreqs,enet_elfp,'k')
xlim(wfreqs([end 1]))
% set(gca,'xscale','log')
title('Late spiking','fontsize',14)

%% Compute correlations between SU spk cnts and LFP power
for i = 1:24
    fprintf('Cell %d of %d\n',i,10);
    for cc = 1:24
        asu_alfp(i,cc,:) = corr(spk_cnts(used_set,i),squeeze(avg_lfp_amps(cc,used_set,:)));
        asu_blfp(i,cc,:) = corr(spk_cnts(used_set,i),squeeze(avgb_lfp_amps(cc,used_set,:)));
        asu_elfp(i,cc,:) = corr(spk_cnts(used_set,i),squeeze(avge_lfp_amps(cc,used_set,:)));
        bsu_alfp(i,cc,:) = corr(beg_spk_cnts(used_set,i),squeeze(avg_lfp_amps(cc,used_set,:)));
        bsu_blfp(i,cc,:) = corr(beg_spk_cnts(used_set,i),squeeze(avgb_lfp_amps(cc,used_set,:)));
        bsu_elfp(i,cc,:) = corr(beg_spk_cnts(used_set,i),squeeze(avge_lfp_amps(cc,used_set,:)));
        esu_alfp(i,cc,:) = corr(fin_spk_cnts(used_set,i),squeeze(avg_lfp_amps(cc,used_set,:)));
        esu_blfp(i,cc,:) = corr(fin_spk_cnts(used_set,i),squeeze(avgb_lfp_amps(cc,used_set,:)));
        esu_elfp(i,cc,:) = corr(fin_spk_cnts(used_set,i),squeeze(avge_lfp_amps(cc,used_set,:)));
    end
end
%%
close all
cd /Users/James/James_scripts/bruce/modelfits
for c = 1:24
    if c <= 10
        cur_lfp = Blocks{1}.suprobes(c);
    else
        cur_lfp = Blocks{1}.muprobes(c-10);
    end
    cur_max = max([max(max(bsu_blfp(c,:,:))) max(max(bsu_elfp(c,:,:))) ...
        max(max(esu_blfp(c,:,:))) max(max(esu_elfp(c,:,:)))]);
    cur_min = min([min(min(bsu_blfp(c,:,:))) min(min(bsu_elfp(c,:,:))) ...
        min(min(esu_blfp(c,:,:))) min(min(esu_elfp(c,:,:)))]);
    subplot(2,2,1)
    pcolor(wfreqs,1:24,squeeze(bsu_blfp(c,:,:)));shading flat
    xl = xlim();
    line(xl,[cur_lfp cur_lfp],'color','w','linewidth',2)
    caxis([cur_min cur_max]);
    title('Early spikes, Early LFP')
    xlabel('Frequency (Hz)','fontsize',14)
    ylabel('LFP Channel','fontsize',14)
    colorbar
    subplot(2,2,2)
    pcolor(wfreqs,1:24,squeeze(bsu_elfp(c,:,:)));shading flat
    xl = xlim();
    line(xl,[cur_lfp cur_lfp],'color','w','linewidth',2)
    caxis([cur_min cur_max]);
    title('Early spikes, Late LFP')
    xlabel('Frequency (Hz)','fontsize',14)
    ylabel('LFP Channel','fontsize',14)
    colorbar
    subplot(2,2,3)
    pcolor(wfreqs,1:24,squeeze(esu_blfp(c,:,:)));shading flat
    xl = xlim();
    line(xl,[cur_lfp cur_lfp],'color','w','linewidth',2)
    caxis([cur_min cur_max]);
    title('Late spikes, Early LFP')
    xlabel('Frequency (Hz)','fontsize',14)
    ylabel('LFP Channel','fontsize',14)
    colorbar
    subplot(2,2,4)
    pcolor(wfreqs,1:24,squeeze(esu_elfp(c,:,:)));shading flat
    xl = xlim();
    line(xl,[cur_lfp cur_lfp],'color','w','linewidth',2)
    caxis([cur_min cur_max]);
    title('Late spikes, Late LFP')
    xlabel('Frequency (Hz)','fontsize',14)
    ylabel('LFP Channel','fontsize',14)
     colorbar
   
    for i = 1:4
        subplot(2,2,i)
        set(gca,'xscale','log')
    end
    if c <= 10
        fname = sprintf('SU_%d_fullLFPpow_corr',c);
    else
        fname = sprintf('MUA_%d_fullLFPpow_corr',c-10);
    end
    fillPage(gcf,'margins',[0 0 0 0],'papersize',[25 15]);
    print('-dpng',fname);
    close
end
%% Compute relationships to Saccade amplitude
[sacamp_asus,sacamp_asus_p] = corr(spk_cnts(used_set,:),all_sac_amps(used_set));
[sacamp_bsus,sacamp_bsus_p] = corr(beg_spk_cnts(used_set,:),all_sac_amps(used_set));
[sacamp_esus,sacamp_esus_p] = corr(fin_spk_cnts(used_set,:),all_sac_amps(used_set));
[sacamp_anet,sacamp_anet_p] = corr(net_spk_cnts(used_set,:),all_sac_amps(used_set));
[sacamp_bnet,sacamp_bnet_p] = corr(netb_spk_cnts(used_set,:),all_sac_amps(used_set));
[sacamp_enet,sacamp_enet_p] = corr(nete_spk_cnts(used_set,:),all_sac_amps(used_set));

for cc = 1:24
    [sacamp_alfp(cc,:),sacamp_lfp_p] = corr(all_sac_amps(used_set),squeeze(avg_lfp_amps(cc,used_set,:)));
    [sacamp_blfp(cc,:),sacamp_blfp_p] = corr(all_sac_amps(used_set),squeeze(avgb_lfp_amps(cc,used_set,:)));
    [sacamp_elfp(cc,:),sacamp_elfp_p] = corr(all_sac_amps(used_set),squeeze(avge_lfp_amps(cc,used_set,:)));
end

% figure
% % plot(wfreqs,sacamp_alfp,'b')
% hold on
% plot(wfreqs,sacamp_blfp,'r')
% plot(wfreqs,sacamp_elfp,'k')

close all
figure
pcolor(wfreqs,1:24,sacamp_blfp)
shading flat
set(gca,'xscale','log')
colorbar
xlabel('Frequency (Hz)','fontsize',14)
ylabel('Correlation','fontsize',14)

figure
pcolor(wfreqs,1:24,sacamp_elfp)
set(gca,'xscale','log')
shading flat
colorbar
xlabel('Frequency (Hz)','fontsize',14)
ylabel('Correlation','fontsize',14)


%% ANALYZE SU PHASE_LOCKING

cmean_all = zeros(10,24,length(wfreqs));
ckapp_all = zeros(10,24,length(wfreqs));
cmean_init = zeros(10,24,length(wfreqs));
ckapp_init = zeros(10,24,length(wfreqs));
cmean_fin = zeros(10,24,length(wfreqs));
ckapp_fin = zeros(10,24,length(wfreqs));
% cmean_bin = zeros(10,24,length(wfreqs),length(ind_bins_edges)-1);
% ckapp_bin = zeros(10,24,length(wfreqs),length(ind_bins_edges)-1);
for cellid = 1:24
    fprintf('Cell %d of %d\n',cellid,24);
    for cc = 1:24
        for wbin = 1:length(wfreqs)
            %     fprintf('Freq %d of %d\n',wbin,length(wfreqs));
            cmean_all(cellid,cc,wbin) = circ_mean(squeeze(all_spk_phases{cellid}(cc,:,wbin))');
            ckapp_all(cellid,cc,wbin) = circ_kappa(squeeze(all_spk_phases{cellid}(cc,:,wbin))');
            cmean_init(cellid,cc,wbin) = circ_mean(squeeze(init_spk_phases{cellid}(cc,:,wbin))');
            ckapp_init(cellid,cc,wbin) = circ_kappa(squeeze(init_spk_phases{cellid}(cc,:,wbin))');
            cmean_fin(cellid,cc,wbin) = circ_mean(squeeze(fin_spk_phases{cellid}(cc,:,wbin))');
            ckapp_fin(cellid,cc,wbin) = circ_kappa(squeeze(fin_spk_phases{cellid}(cc,:,wbin))');
            %         cpval(cellid,wbin) = circ_otest(all_spk_phases{cellid}(:,wbin));
            
%             for ii = 1:length(ind_bins_edges)-1
%                 cmean_bin(cellid,cc,wbin,ii) = circ_mean(bin_spk_phases{cellid,ii}(cc,:,wbin)');
%                 ckapp_bin(cellid,cc,wbin,ii) = circ_kappa(bin_spk_phases{cellid,ii}(cc,:,wbin)');
%             end
        end
    end
end

%% PLOT OF ALL PHASE-LOCKING
for i = 1:5
    cur_probe = Blocks{1}.suprobes(i)+0.5;
    subplot(5,4,(i-1)*4+1)
    pcolor(wfreqs,1:24,squeeze(ckapp_init(i,:,:)));shading flat
    xlim([10 100]);
    xl = xlim();
    line(xl,[cur_probe cur_probe],'color','w','linewidth',2)
    caxis([0 0.5])
    xlabel('Frequency (Hz)','fontsize',16)
    ylabel('Channel Number','fontsize',16)
    title(sprintf('Cell%d_early',i),'fontsize',20)
    subplot(5,4,(i-1)*4+2)
    pcolor(wfreqs,1:24,squeeze(ckapp_fin(i,:,:)));shading flat
    xlim([10 100])
     caxis([0 0.5])
    xl = xlim();
    line(xl,[cur_probe cur_probe],'color','w','linewidth',2)
    xlabel('Frequency (Hz)','fontsize',16)
    ylabel('Channel Number','fontsize',16)
    title(sprintf('Cell%d Late',i),'fontsize',20)
end   
for i = 1:5
    cur_probe = Blocks{1}.suprobes(i+5)+0.5;
    subplot(5,4,(i-1)*4+3)
    pcolor(wfreqs,1:24,squeeze(ckapp_init(i+5,:,:)));shading flat
    xlim([10 100])
    xl = xlim();
    line(xl,[cur_probe cur_probe],'color','w','linewidth',2)
    caxis([0 0.5])
    xlabel('Frequency (Hz)','fontsize',16)
    ylabel('Channel Number','fontsize',16)
    title(sprintf('Cell%d Early',i+5),'fontsize',20)
    subplot(5,4,(i-1)*4+4)
    pcolor(wfreqs,1:24,squeeze(ckapp_fin(i+5,:,:)));shading flat
    xlim([10 100])
    xl = xlim();
    line(xl,[cur_probe cur_probe],'color','w','linewidth',2)
     caxis([0 0.5])
    xlabel('Frequency (Hz)','fontsize',16)
    ylabel('Channel Number','fontsize',16)
    title(sprintf('Cell%d Late',i+5),'fontsize',20)
end   
fillPage(gcf,'margins',[0 0 0 0],'papersize',[30 30])

%% PLOT OF ALL POWER CORRS
for i = 1:5
    cur_probe = Blocks{1}.suprobes(i)+0.5;
    subplot(5,4,(i-1)*4+1)
    pcolor(wfreqs,1:24,squeeze(bsu_blfp(i,:,:)));shading flat
    xlim([10 100]);
    xl = xlim();
    line(xl,[cur_probe cur_probe],'color','w','linewidth',2)
    caxis([0 0.3])
    xlabel('Frequency (Hz)','fontsize',16)
    ylabel('Channel Number','fontsize',16)
    title(sprintf('Cell%d_early',i),'fontsize',20)
    subplot(5,4,(i-1)*4+2)
    pcolor(wfreqs,1:24,squeeze(esu_elfp(i,:,:)));shading flat
    xlim([10 100])
    caxis([0 0.3])
    xl = xlim();
    line(xl,[cur_probe cur_probe],'color','w','linewidth',2)
    xlabel('Frequency (Hz)','fontsize',16)
    ylabel('Channel Number','fontsize',16)
    title(sprintf('Cell%d Late',i),'fontsize',20)
end   
for i = 1:5
    cur_probe = Blocks{1}.suprobes(i+5)+0.5;
    subplot(5,4,(i-1)*4+3)
    pcolor(wfreqs,1:24,squeeze(bsu_blfp(i+5,:,:)));shading flat
    xlim([10 100])
    xl = xlim();
    line(xl,[cur_probe cur_probe],'color','w','linewidth',2)
    caxis([0 0.3])
    xlabel('Frequency (Hz)','fontsize',16)
    ylabel('Channel Number','fontsize',16)
    title(sprintf('Cell%d Early',i+5),'fontsize',20)
    subplot(5,4,(i-1)*4+4)
    pcolor(wfreqs,1:24,squeeze(esu_elfp(i+5,:,:)));shading flat
    xlim([10 100])
    xl = xlim();
    line(xl,[cur_probe cur_probe],'color','w','linewidth',2)
    caxis([0 0.3])
    xlabel('Frequency (Hz)','fontsize',16)
    ylabel('Channel Number','fontsize',16)
    title(sprintf('Cell%d Late',i+5),'fontsize',20)
end   
fillPage(gcf,'margins',[0 0 0 0],'papersize',[30 30])

%%
for i = 1:10
    su_probe = Blocks{1}.suprobes(i);
    cur_kappai(i,:) = squeeze(ckapp_init(i,su_probe-1,:));
    cur_kappaf(i,:) = squeeze(ckapp_fin(i,su_probe-1,:));
end


%%
% clf
c = 6
% for c = 1:24
    if c <= 10
        cur_probe = Blocks{1}.suprobes(c);
    else
        cur_probe = Blocks{1}.muprobes(c-10);
    end
    % cur_max = max([max(max(ckapp_all(c,:,:))) max(max(ckapp_init(c,:,:))) ...
    %     max(max(ckapp_fin(c,:,:)))]);
    % cur_min = min([min(min(ckapp_all(c,:,:))) min(min(ckapp_init(c,:,:))) ...
    %     min(min(ckapp_fin(c,:,:)))]);
    cur_max = 0.5;
    cur_min = 0;
    
    % subplot(2,1,1)
    % pcolor(wfreqs,1:24,squeeze(ckapp_all(c,:,:)));shading flat
    % xl = xlim();
    % line(xl,[cur_probe cur_probe],'color','w','linewidth',2)
    % caxis([cur_min cur_max]);
    % title('Early spikes, Early LFP')
    % xlabel('Frequency (Hz)','fontsize',14)
    % ylabel('LFP Channel','fontsize',14)
    % colorbar
    subplot(2,1,1)
    pcolor(wfreqs,1:24,squeeze(ckapp_init(c,:,:)));shading flat
    xl = xlim();
    line(xl,[cur_probe cur_probe],'color','w','linewidth',2)
    caxis([cur_min cur_max]);
    title('Early spikes, Late LFP')
    xlabel('Frequency (Hz)','fontsize',14)
    ylabel('LFP Channel','fontsize',14)
    colorbar
    subplot(2,1,2)
    pcolor(wfreqs,1:24,squeeze(ckapp_fin(c,:,:)));shading flat
    xl = xlim();
    line(xl,[cur_probe cur_probe],'color','w','linewidth',2)
    caxis([cur_min cur_max]);
    title('Late spikes, Early LFP')
    xlabel('Frequency (Hz)','fontsize',14)
    ylabel('LFP Channel','fontsize',14)
    colorbar
    for i = 1:2
        subplot(2,1,i)
        set(gca,'xscale','log')
    end
%     if c <= 10
%         fname = sprintf('SU_%d_fullLFP_phaselock',c);
%     else
%         fname = sprintf('MUA_%d_fullLFP_phaselock',c-10);
%     end
%     fillPage(gcf,'margins',[0 0 0 0],'papersize',[15 20]);
%     print('-dpng',fname);
%     close
%     
% end
    
  c = 6
% for c = 1:24
    if c <= 10
        cur_probe = Blocks{1}.suprobes(c);
    else
        cur_probe = Blocks{1}.muprobes(c-10);
    end
    % cur_max = max([max(max(ckapp_all(c,:,:))) max(max(ckapp_init(c,:,:))) ...
    %     max(max(ckapp_fin(c,:,:)))]);
    % cur_min = min([min(min(ckapp_all(c,:,:))) min(min(ckapp_init(c,:,:))) ...
    %     min(min(ckapp_fin(c,:,:)))]);
    cur_max = 0.5;
    cur_min = 0;
    
    % subplot(2,1,1)
    % pcolor(wfreqs,1:24,squeeze(ckapp_all(c,:,:)));shading flat
    % xl = xlim();
    % line(xl,[cur_probe cur_probe],'color','w','linewidth',2)
    % caxis([cur_min cur_max]);
    % title('Early spikes, Early LFP')
    % xlabel('Frequency (Hz)','fontsize',14)
    % ylabel('LFP Channel','fontsize',14)
    % colorbar
    subplot(2,1,1)
    pcolor(wfreqs,1:24,squeeze(ckapp_init(c,:,:)));shading flat
    xl = xlim();
    line(xl,[cur_probe cur_probe],'color','w','linewidth',2)
    caxis([cur_min cur_max]);
    title('Early spikes, Late LFP')
    xlabel('Frequency (Hz)','fontsize',14)
    ylabel('LFP Channel','fontsize',14)
    colorbar
    subplot(2,1,2)
    pcolor(wfreqs,1:24,squeeze(ckapp_fin(c,:,:)));shading flat
    xl = xlim();
    line(xl,[cur_probe cur_probe],'color','w','linewidth',2)
    caxis([cur_min cur_max]);
    title('Late spikes, Early LFP')
    xlabel('Frequency (Hz)','fontsize',14)
    ylabel('LFP Channel','fontsize',14)
    colorbar
    for i = 1:2
        subplot(2,1,i)
        set(gca,'xscale','log')
    end
%     if c <= 10
%         fname = sprintf('SU_%d_fullLFP_phaselock',c);
%     else
%         fname = sprintf('MUA_%d_fullLFP_phaselock',c-10);
%     end
%     fillPage(gcf,'margins',[0 0 0 0],'papersize',[15 20]);
%     print('-dpng',fname);
%     close
%     
% end
  
%%
close all
c = 4
cur_probe = Blocks{1}.suprobes(c);
figure
pcolor(ind_bins_cents/Fsd,wfreqs,squeeze(ckapp_bin(c,cur_probe,:,:)));shading flat
 
figure
pcolor(ind_bins_cents/Fsd,wfreqs,squeeze(cmean_bin(c,cur_probe,:,:)));shading flat
   
%%
close all
cd /Users/James/James_scripts/bruce/modelfits
for c = 1:24
    if c <= 10
        cur_lfp = Blocks{1}.suprobes(c);
    else
        cur_lfp = Blocks{1}.muprobes(c-10);
    end
    cur_max = max([max(max(bsu_blfp(c,:,:))) max(max(bsu_elfp(c,:,:))) ...
        max(max(esu_blfp(c,:,:))) max(max(esu_elfp(c,:,:)))]);
    cur_min = min([min(min(bsu_blfp(c,:,:))) min(min(bsu_elfp(c,:,:))) ...
        min(min(esu_blfp(c,:,:))) min(min(esu_elfp(c,:,:)))]);
    
    subplot(2,3,1)
    pcolor(wfreqs,1:24,squeeze(bsu_blfp(c,:,:)));shading flat
    xl = xlim();
    line(xl,[cur_lfp cur_lfp],'color','w','linewidth',2)
    caxis([cur_min cur_max]);
    title('Early spikes, Early LFP','fontsize',16)
    xlabel('Frequency (Hz)','fontsize',16)
    ylabel('LFP Channel','fontsize',16)
    
    subplot(2,3,2)
    pcolor(wfreqs,1:24,squeeze(bsu_elfp(c,:,:)));shading flat
    xl = xlim();
    line(xl,[cur_lfp cur_lfp],'color','w','linewidth',2)
    caxis([cur_min cur_max]);
    title('Early spikes, Late LFP','fontsize',16)
    xlabel('Frequency (Hz)','fontsize',16)
    ylabel('LFP Channel','fontsize',16)
    
    subplot(2,3,3)
    pcolor(wfreqs,1:24,squeeze(ckapp_init(c,:,:)));shading flat
    xl = xlim();
    line(xl,[cur_lfp cur_lfp],'color','w','linewidth',2)
    caxis([0 0.5]);
    title('Early Spikes Phase-locking','fontsize',16)
    xlabel('Frequency (Hz)','fontsize',16)
    ylabel('LFP Channel','fontsize',16)    
    
    subplot(2,3,4)
    pcolor(wfreqs,1:24,squeeze(esu_blfp(c,:,:)));shading flat
    xl = xlim();
    line(xl,[cur_lfp cur_lfp],'color','w','linewidth',2)
    caxis([cur_min cur_max]);
    title('Late spikes, Early LFP','fontsize',16)
    xlabel('Frequency (Hz)','fontsize',16)
    ylabel('LFP Channel','fontsize',16)
    
    subplot(2,3,5)
    pcolor(wfreqs,1:24,squeeze(esu_elfp(c,:,:)));shading flat
    xl = xlim();
    line(xl,[cur_lfp cur_lfp],'color','w','linewidth',2)
    caxis([cur_min cur_max]);
    title('Late spikes, Late LFP','fontsize',16)
    xlabel('Frequency (Hz)','fontsize',16)
    ylabel('LFP Channel','fontsize',16)
    
    subplot(2,3,6)
    pcolor(wfreqs,1:24,squeeze(ckapp_fin(c,:,:)));shading flat
    xl = xlim();
    line(xl,[cur_lfp cur_lfp],'color','w','linewidth',2)
    caxis([0 0.5]);
    title('Late spikes phase-locking','fontsize',16)
    xlabel('Frequency (Hz)','fontsize',16)
    ylabel('LFP Channel','fontsize',16)

    for i = 1:6
        subplot(2,3,i)
        set(gca,'xscale','log')
        %colorbar
    end
    if c <= 10
        fname = sprintf('SU_%d_full_powkapp',c);
    else
        fname = sprintf('MUA_%d_full_powkapp',c-10);
    end
    fillPage(gcf,'margins',[0 0 0 0],'papersize',[26 15]);
    print('-dpng',fname);
    close
end

    
    
    
    





