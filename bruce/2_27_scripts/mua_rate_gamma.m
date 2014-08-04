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
[b,a] = butter(2,[1 200]/niqf);

scales = logspace(log10(3),log10(350),100);
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);

    
    %%
cd ~/Data/bruce/2_27_12
n_units = 24;

avg_lfp_amps = nan(length(all_fix_start_times),length(scales));
avgb_lfp_amps = nan(length(all_fix_start_times),length(scales));
avge_lfp_amps = nan(length(all_fix_start_times),length(scales));
net_spk_cnts = nan(length(all_fix_start_times),1);
netb_spk_cnts = nan(length(all_fix_start_times),1);
nete_spk_cnts = nan(length(all_fix_start_times),1);
spk_cnts = nan(length(all_fix_start_times),n_units);
beg_spk_cnts = nan(length(all_fix_start_times),n_units);
fin_spk_cnts = nan(length(all_fix_start_times),n_units);
used_pts = [];

ind_bins_edges = round(Fsd*(0:0.05:0.4));
all_spk_phases = cell(n_units,1);
init_spk_phases = cell(n_units,1);
fin_spk_phases = cell(n_units,1);
bin_spk_phases = cell(n_units,length(ind_bins_edges)-1);
for blockid = 1:3;
    
    fprintf('Processing block %d of %d...\n',blockid,4);
    
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
    lfp_num = 2;
    coefs = cwt(lfp_sampsd(:,lfp_num),scales,'cmor1-1');
    
    ampgram = abs(coefs);
    ampgram = bsxfun(@minus,ampgram,mean(ampgram,2));
    ampgram = bsxfun(@rdivide,ampgram,std(ampgram,[],2));
    ampgram = ampgram';
    
    phasegram = angle(coefs);
    phasegram = phasegram';
    
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
        
        avg_lfp_amps(temp_fix_set(i),:) = mean(ampgram(cur_range1,:));
        avgb_lfp_amps(temp_fix_set(i),:) = mean(ampgram(cur_range3,:));
        avge_lfp_amps(temp_fix_set(i),:) = mean(ampgram(cur_range2,:));
        
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
    
    for cellid = 1:24;
        if cellid <= 10
        cur_spktimes = Blocks{blockid}.spktimes{cellid};
        else
            cur_spktimes = Blocks{blockid}.mutimes{cellid-10};
        end
        cur_spk_inds = round(interp1(lfp_timed,1:length(lfp_timed),cur_spktimes));
        cur_all_spk = cur_spk_inds(ismember(cur_spk_inds,all_used_inds));
        cur_init_spk = cur_spk_inds(ismember(cur_spk_inds,all_init_inds));
        cur_fin_spk = cur_spk_inds(ismember(cur_spk_inds,all_fin_inds));
        all_spk_phases{cellid} = [all_spk_phases{cellid}; phasegram(cur_all_spk,:)];
        init_spk_phases{cellid} = [init_spk_phases{cellid}; phasegram(cur_init_spk,:)];
        fin_spk_phases{cellid} = [fin_spk_phases{cellid}; phasegram(cur_fin_spk,:)];
%         for ii = 1:length(ind_bins_edges)-1
%             cur_spk = cur_spk_inds(ismember(cur_spk_inds,all_inds_bin{ii}));
%             bin_spk_phases{cellid,ii} = [bin_spk_phases{cellid,ii}; phasegram(cur_spk,:)];
%         end
    end
    
end
%%
cd /Users/James/James_scripts/bruce/modelfits
save lfp_spike_calcs_ch2 all_* init_* fin_* *spk_cnts avg*_lfp* wfreqs Fsd *spk_phases
%%
used_set = find(~isnan(net_spk_cnts) & ~isnan(avge_lfp_amps(:,1)));

%% Analyze relationship between net spiking and LFP power
anet_alfp = corr(net_spk_cnts(used_set),avg_lfp_amps(used_set,:));
anet_blfp = corr(net_spk_cnts(used_set),avgb_lfp_amps(used_set,:));
anet_elfp = corr(net_spk_cnts(used_set),avge_lfp_amps(used_set,:));
bnet_alfp = corr(netb_spk_cnts(used_set),avg_lfp_amps(used_set,:));
bnet_blfp = corr(netb_spk_cnts(used_set),avgb_lfp_amps(used_set,:));
bnet_elfp = corr(netb_spk_cnts(used_set),avge_lfp_amps(used_set,:));
enet_alfp = corr(nete_spk_cnts(used_set),avg_lfp_amps(used_set,:));
enet_blfp = corr(nete_spk_cnts(used_set),avgb_lfp_amps(used_set,:));
enet_elfp = corr(nete_spk_cnts(used_set),avge_lfp_amps(used_set,:));

subplot(3,1,1)
plot(wfreqs,anet_alfp)
hold on
plot(wfreqs,anet_blfp,'r')
plot(wfreqs,anet_elfp,'k')
xlim(wfreqs([end 1]))
set(gca,'xscale','log')
title('Overall spiking','fontsize',14)
xlabel('Frequency (Hz)','fontsize',14)
ylabel('Correlation','fontsize',14)

subplot(3,1,2)
plot(wfreqs,bnet_alfp)
hold on
plot(wfreqs,bnet_blfp,'r')
plot(wfreqs,bnet_elfp,'k')
xlim(wfreqs([end 1]))
set(gca,'xscale','log')
title('Initial spiking','fontsize',14)
xlabel('Frequency (Hz)','fontsize',14)
ylabel('Correlation','fontsize',14)

subplot(3,1,3)
plot(wfreqs,enet_alfp)
hold on
plot(wfreqs,enet_blfp,'r')
plot(wfreqs,enet_elfp,'k')
xlim(wfreqs([end 1]))
set(gca,'xscale','log')
title('Late spiking','fontsize',14)
xlabel('Frequency (Hz)','fontsize',14)
ylabel('Correlation','fontsize',14)

%% Compute correlations between SU spk cnts and LFP power
for i = 1:24
    fprintf('Cell %d of %d\n',i,24);
    asu_alfp(i,:) = corr(spk_cnts(used_set,i),avg_lfp_amps(used_set,:));
    asu_blfp(i,:) = corr(spk_cnts(used_set,i),avgb_lfp_amps(used_set,:));
    asu_elfp(i,:) = corr(spk_cnts(used_set,i),avge_lfp_amps(used_set,:));
    bsu_alfp(i,:) = corr(beg_spk_cnts(used_set,i),avg_lfp_amps(used_set,:));
    bsu_blfp(i,:) = corr(beg_spk_cnts(used_set,i),avgb_lfp_amps(used_set,:));
    bsu_elfp(i,:) = corr(beg_spk_cnts(used_set,i),avge_lfp_amps(used_set,:));
    esu_alfp(i,:) = corr(fin_spk_cnts(used_set,i),avg_lfp_amps(used_set,:));
    esu_blfp(i,:) = corr(fin_spk_cnts(used_set,i),avgb_lfp_amps(used_set,:));
    esu_elfp(i,:) = corr(fin_spk_cnts(used_set,i),avge_lfp_amps(used_set,:));
end

%%
for c = 1:24;
    subplot(3,1,1)
    plot(wfreqs,asu_alfp(c,:))
    hold on
    plot(wfreqs,asu_blfp(c,:),'r')
    plot(wfreqs,asu_elfp(c,:),'k')
    xlim(wfreqs([end 1]))
    set(gca,'xscale','log')
    title('Overall spiking','fontsize',14)
    xlabel('Frequency (Hz)','fontsize',14)
    ylabel('Correlation','fontsize',14)
    
    subplot(3,1,2)
    plot(wfreqs,bsu_alfp(c,:))
    hold on
    plot(wfreqs,bsu_blfp(c,:),'r')
    plot(wfreqs,bsu_elfp(c,:),'k')
    xlim(wfreqs([end 1]))
    set(gca,'xscale','log')
    title('Initial spiking','fontsize',14)
    xlabel('Frequency (Hz)','fontsize',14)
    ylabel('Correlation','fontsize',14)
    
    subplot(3,1,3)
    plot(wfreqs,esu_alfp(c,:))
    hold on
    plot(wfreqs,esu_blfp(c,:),'r')
    plot(wfreqs,esu_elfp(c,:),'k')
    xlim(wfreqs([end 1]))
    set(gca,'xscale','log')
    title('Late spiking','fontsize',14)
    xlabel('Frequency (Hz)','fontsize',14)
    ylabel('Correlation','fontsize',14)
    
    fillPage(gcf,'margins',[0 0 0 0],'papersize',[10 20]);
    if c > 10
        fname = sprintf('MUA_%d_LFP_pow_corr',c-10);
    else
        fname = sprintf('SU_%d_LFP_pow_corr',c);
    end
    print('-dpng',fname);close
end
%% Compute relationships to Saccade amplitude
[sacamp_asus,sacamp_asus_p] = corr(spk_cnts(used_set,:),all_sac_amps(used_set));
[sacamp_bsus,sacamp_bsus_p] = corr(beg_spk_cnts(used_set,:),all_sac_amps(used_set));
[sacamp_esus,sacamp_esus_p] = corr(fin_spk_cnts(used_set,:),all_sac_amps(used_set));
[sacamp_anet,sacamp_anet_p] = corr(net_spk_cnts(used_set,:),all_sac_amps(used_set));
[sacamp_bnet,sacamp_bnet_p] = corr(netb_spk_cnts(used_set,:),all_sac_amps(used_set));
[sacamp_enet,sacamp_enet_p] = corr(nete_spk_cnts(used_set,:),all_sac_amps(used_set));

[sacamp_alfp,sacamp_lfp_p] = corr(all_sac_amps(used_set),avg_lfp_amps(used_set,:));
[sacamp_blfp,sacamp_blfp_p] = corr(all_sac_amps(used_set),avgb_lfp_amps(used_set,:));
[sacamp_elfp,sacamp_elfp_p] = corr(all_sac_amps(used_set),avge_lfp_amps(used_set,:));

figure
plot(wfreqs,sacamp_alfp)
hold on
plot(wfreqs,sacamp_blfp,'r')
plot(wfreqs,sacamp_elfp,'k')
set(gca,'xscale','log')
xlabel('Frequency (Hz)','fontsize',14)
ylabel('Correlation','fontsize',14)
xlim(wfreqs([end 1]))

%% ANALYZE SU PHASE_LOCKING

cmean_all = zeros(n_units,length(wfreqs));
ckapp_all = zeros(n_units,length(wfreqs));
cmean_init = zeros(n_units,length(wfreqs));
ckapp_init = zeros(n_units,length(wfreqs));
cmean_fin = zeros(n_units,length(wfreqs));
ckapp_fin = zeros(n_units,length(wfreqs));
cmean_bin = zeros(n_units,length(wfreqs),length(ind_bins_edges)-1);
ckapp_bin = zeros(n_units,length(wfreqs),length(ind_bins_edges)-1);
for cellid = 1:n_units
    fprintf('Cell %d of %d\n',cellid,n_units);
    for wbin = 1:length(wfreqs)
        %     fprintf('Freq %d of %d\n',wbin,length(wfreqs));
        cmean_all(cellid,wbin) = circ_mean(all_spk_phases{cellid}(:,wbin));
        ckapp_all(cellid,wbin) = circ_kappa(all_spk_phases{cellid}(:,wbin));
        cmean_init(cellid,wbin) = circ_mean(init_spk_phases{cellid}(:,wbin));
        ckapp_init(cellid,wbin) = circ_kappa(init_spk_phases{cellid}(:,wbin));
        cmean_fin(cellid,wbin) = circ_mean(fin_spk_phases{cellid}(:,wbin));
        ckapp_fin(cellid,wbin) = circ_kappa(fin_spk_phases{cellid}(:,wbin));
        %         cpval(cellid,wbin) = circ_otest(all_spk_phases{cellid}(:,wbin));
        
%         for ii = 1:length(ind_bins_edges)-1
%            cmean_bin(cellid,wbin,ii) = circ_mean(bin_spk_phases{cellid,ii}(:,wbin)); 
%            ckapp_bin(cellid,wbin,ii) = circ_kappa(bin_spk_phases{cellid,ii}(:,wbin)); 
%         end
    end
end

%%
su_probes = Blocks{1}.suprobes;
mu_probes = Blocks{1}.muprobes;
all_probes = [su_probes mu_probes];
[~,probe_ord] = sort(all_probes);

%%
clf
c = 19
plot(wfreqs,ckapp_all(c,:));
hold on
plot(wfreqs,ckapp_init(c,:),'r')
plot(wfreqs,ckapp_fin(c,:),'k')
% set(gca,'yscale','log')
set(gca,'xscale','log')
xlim([5 wfreqs(1)])
ylim([0 0.5])













