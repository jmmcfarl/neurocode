clear all
close all
cd ~/Data/bruce/G075/
load jbeG075Expts.mat

%%
sim_sac_blocks = [10 14 18 22 24 28 32 37 40 42 50] - 6;
% sim_sac_blocks = [18 28 40] - 6; %sim sac
% sim_sac_blocks = [10 22 32 42] - 6; %sim sac a
% sim_sac_blocks = [14 24 37 50] - 6; %sim sac b

repeat_inds = [1e3 2e3 3e3 2.01e5 2.02e5 2.03e5 4.01e5 4.02e5 4.03e5];
all_expt_id = [];
for bb = 1:length(sim_sac_blocks)
    n_trials(bb) = length(Expts{sim_sac_blocks(bb)}.Trials);
    trial_start_times{bb} = [Expts{sim_sac_blocks(bb)}.Trials(:).Start]/1e4;
    trial_stop_times{bb} = [Expts{sim_sac_blocks(bb)}.Trials(:).End]/1e4;
    trial_seof{bb} = [Expts{sim_sac_blocks(bb)}.Trials(:).seof];
    trial_completed{bb} = [Expts{sim_sac_blocks(bb)}.Trials(:).Result];
    trial_durs{bb} = trial_stop_times{bb} - trial_start_times{bb};
    all_expt_id = [all_expt_id sim_sac_blocks(bb)*ones(1,n_trials(bb))];
end
all_trial_start_times = cell2mat(trial_start_times);
all_trial_stop_times = cell2mat(trial_stop_times);
all_trial_seof = cell2mat(trial_seof);
all_trial_completed = cell2mat(trial_completed);
all_trial_durs = cell2mat(trial_durs);

%%
stim_fs = 1e4/117.5;
Fs = 3e4;
% dsf = 120;Fsd = Fs/dsf;
dsf = 80;Fsd = Fs/dsf;
niqf = Fsd/2;
[filt_b,filt_a] = butter(2,[1]/niqf,'high');
use_lfps = [1:24:96];
use_units = 1:96;
% use_lfps = [1 63];
backlag = round(0*Fsd);
forwardlag = round(4*Fsd);
lags = (-backlag:forwardlag);

scales = logspace(log10(10),log10(125),30);
scales = scales*60/dsf;
% scales = scales*2;
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);

desired_dt = 0.003;
%%
all_interp_lfps = [];
all_interp_phasegrams = [];
all_interp_ampgrams = [];
all_binned_spks = [];
all_expt_t = [];
all_stim_id = [];
all_t_since_start = [];
all_trial_id = [];

for bb = 1:length(sim_sac_blocks)
    fprintf('Analyzing block %d of %d\n',bb,length(sim_sac_blocks));
    load(sprintf('Expt%dClusterTimes.mat',sim_sac_blocks(bb)));
    
    filename = sprintf('Expt%dFullVmean.mat',sim_sac_blocks(bb));
    load(filename);
    
    Vmat = [];
    %     Vmatf = [];
    phasegrams = [];
    ampgrams = [];
    for ll = 1:length(use_lfps)
        fprintf('Electrode %d of %d\n',ll,length(use_lfps));
        filename = sprintf('Expt%d.p%dFullV.mat',sim_sac_blocks(bb),use_lfps(ll));
        load(filename);
        V = double(FullV.V);
        V = V + FullV.sumscale*sumv;
        V = V*FullV.intscale(1)/FullV.intscale(2);
        nparts = length(FullV.blklen);
        %         dV = [];
        dVf = [];
        cur_pt = 1;
        for pp = 1:nparts
            cur_range = cur_pt:(cur_pt + FullV.blklen(pp)-1);
            cur_range(cur_range > length(V)) = [];
            %             dV = [dV decimate(V(cur_range),dsf)];
            dVf = [dVf filtfilt(filt_b,filt_a,decimate(V(cur_range),dsf))]; %do some high-pass filtering
            cur_pt = cur_pt + FullV.blklen(pp);
        end
        Vmat(ll,:) = dVf;
        temp = cwt(dVf,scales,'cmor1-1');
        phasegrams(:,:,ll) = angle(temp)';
        ampgrams(:,:,ll) = abs(temp)';
    end
    t_ax = [];
    for pp = 1:nparts
        cur_t_ax = linspace(FullV.blkstart(pp),FullV.blkstart(pp)+FullV.blklen(pp)/Fs,FullV.blklen(pp));
        t_ax = [t_ax downsample(cur_t_ax,dsf)];
    end
    t_ax(size(Vmat,2)+1:end) = [];
    
    cur_use_trials = find(all_expt_id==sim_sac_blocks(bb) & all_trial_completed==1);
    
    trial_dur = floor(3.75/desired_dt);
    expt_t_axis = [];
    t_since_start = [];
    block_trial_id = [];
    block_stim_id = [];
    trial_cnt = 1;
    for tt = 1:length(cur_use_trials)
        cur_t_edges = all_trial_start_times(cur_use_trials(tt)):desired_dt:all_trial_stop_times(cur_use_trials(tt));
        cur_t_edges(trial_dur+2:end) = [];
        cl(tt) = length(cur_t_edges);
        binned_spks = nan(length(use_units),length(cur_t_edges)-1);
        for cc = 1:length(use_units)
            temp = histc(Clusters{use_units(cc)}.times,cur_t_edges);
            binned_spks(cc,:) = temp(1:end-1);
        end
        all_binned_spks = cat(1,all_binned_spks,binned_spks');
        block_trial_id = [block_trial_id; ones(length(cur_t_edges)-1,1)*trial_cnt];
        block_stim_id = [block_stim_id; ones(length(cur_t_edges)-1,1)*all_trial_seof(cur_use_trials(tt))];
        trial_cnt = trial_cnt + 1;
        
        expt_t_axis = [expt_t_axis cur_t_edges(1:end-1)+desired_dt/2];
        t_since_start = [t_since_start cur_t_edges(1:end-1)-cur_t_edges(1)+desired_dt/2];
    end
    
    unwr_phasegram = unwrap(phasegrams);
    interp_phasegrams = interp1(t_ax,unwr_phasegram,expt_t_axis);
    interp_phasegrams = mod(interp_phasegrams+pi,2*pi)-pi;
    interp_ampgrams = interp1(t_ax,ampgrams,expt_t_axis);
    interp_lfps = interp1(t_ax,Vmat',expt_t_axis);
    
    all_interp_phasegrams = cat(1,all_interp_phasegrams,interp_phasegrams);
    all_interp_ampgrams = cat(1,all_interp_ampgrams,interp_ampgrams);
    all_interp_lfps = cat(1,all_interp_lfps,interp_lfps);
    all_expt_t = cat(1,all_expt_t,expt_t_axis');
    all_stim_id = cat(1,all_stim_id,block_stim_id);
    all_t_since_start = cat(1,all_t_since_start,t_since_start');
end

use_set = [1 2 3 7 8 9];
use_inds = find(ismember(all_stim_id,repeat_inds(use_set)));
all_interp_n_ampgrams = bsxfun(@rdivide,all_interp_ampgrams,nanstd(all_interp_ampgrams(use_inds,:,:)));

rate_sm = round(.01/desired_dt);
for cc = 1:96
    all_sm_binned_spks(:,cc) = jmm_smooth_1d_cor(all_binned_spks(:,cc),rate_sm);
end
%%
load ~/Data/bruce/7_15_12/G029/ArrayConfig.mat
X_pos = ArrayConfig.X;
Y_pos = ArrayConfig.Y;
nearest_lfps = nan(length(use_lfps),1);
for ll = 1:96
    all_dists = sqrt((X_pos-X_pos(ll)).^2 + (Y_pos-Y_pos(ll)).^2);
    all_dists(ll) = inf;
    [~,best_loc] = min(all_dists(use_lfps));
    nearest_lfps(ll) = best_loc;
end

%%
% use_freqs = find(wfreqs <= 50);
use_freqs = find(wfreqs <= 50 & wfreqs >= 5);

load ./ns_40_phase_mod_v2
phase_ampkern = sqrt(phase_cfilt.^2+phase_sfilt.^2);
phase_phasekern = -atan2(phase_cfilt,phase_sfilt)+pi/2;
ampphase_ampkern = sqrt(ampphase_cfilt.^2+ampphase_sfilt.^2);
ampphase_phasekern = -atan2(ampphase_cfilt,ampphase_sfilt)+pi/2;
phase_filt = [phase_cfilt(:,use_freqs) phase_sfilt(:,use_freqs)];
ampphase_filt = [ampphase_cfilt(:,use_freqs) ampphase_sfilt(:,use_freqs)];

% load ./sim_sac_phase_mod_v2
load ./sim_sac_phase_mod_late_v2.mat
sphase_ampkern = sqrt(phase_cfilt.^2+phase_sfilt.^2);
sphase_phasekern = -atan2(phase_cfilt,phase_sfilt)+pi/2;
sampphase_ampkern = sqrt(ampphase_cfilt.^2+ampphase_sfilt.^2);
sampphase_phasekern = -atan2(ampphase_cfilt,ampphase_sfilt)+pi/2;
sphase_filt = [phase_cfilt(:,use_freqs) phase_sfilt(:,use_freqs)];
sampphase_filt = [ampphase_cfilt(:,use_freqs) ampphase_sfilt(:,use_freqs)];

load ./sim_sac_phase_mod_to

NT = size(all_interp_phasegrams,1);
phase_set = [reshape(cos(all_interp_phasegrams(:,use_freqs,:)),NT,length(use_freqs)*length(use_lfps)) reshape(sin(all_interp_phasegrams(:,use_freqs,:)),NT,length(use_freqs)*length(use_lfps))];
ampphase_set = [reshape(all_interp_n_ampgrams(:,use_freqs,:),NT,length(use_freqs)*length(use_lfps)).*reshape(cos(all_interp_phasegrams(:,use_freqs,:)),NT,length(use_freqs)*length(use_lfps)) ...
    reshape(all_interp_n_ampgrams(:,use_freqs,:),NT,length(use_freqs)*length(use_lfps)).*reshape(sin(all_interp_phasegrams(:,use_freqs,:)),NT,length(use_freqs)*length(use_lfps))];
% phase_elec_set = [repmat(1:length(use_lfps),1,length(use_freqs)) repmat(1:length(use_lfps),1,length(use_freqs))];
% phase_elec_set = phase_elec_set(:);
phase_elec_set = ones(length(use_freqs),1)*(1:length(use_lfps));
phase_elec_set = [phase_elec_set(:); phase_elec_set(:)];

for cur_unit = 1:96;
    fprintf('Unit %d of %d\n',cur_unit,96);
    cur_lfp = nearest_lfps(cur_unit);
    for ss = 1:length(repeat_inds)
        cur_inds = find(all_stim_id==repeat_inds(ss));
        bad_set = find(isnan(all_interp_ampgrams(cur_inds,1,1)));
        cur_inds(bad_set) = [];
        use_set = find(phase_elec_set==cur_lfp);
        
        cur_phasemod_out = phase_set(cur_inds,use_set)*phase_filt(cur_unit,:)';
        cur_ampphasemod_out = ampphase_set(cur_inds,use_set)*ampphase_filt(cur_unit,:)';
        cur_sphasemod_out = phase_set(cur_inds,use_set)*sphase_filt(cur_unit,:)';
        cur_sampphasemod_out = ampphase_set(cur_inds,use_set)*sampphase_filt(cur_unit,:)';
        n_trials(ss) = length(cur_phasemod_out)/trial_dur;
        
        phasemod_out{ss,cur_unit} = reshape(cur_phasemod_out,trial_dur,n_trials(ss));
        ampmod_out{ss,cur_unit} = reshape(cur_ampphasemod_out,trial_dur,n_trials(ss));
        sphasemod_out{ss,cur_unit} = reshape(cur_sphasemod_out,trial_dur,n_trials(ss));
        sampmod_out{ss,cur_unit} = reshape(cur_sampphasemod_out,trial_dur,n_trials(ss));
        binned_spk_mat{ss,cur_unit} = reshape(all_sm_binned_spks(cur_inds,cur_unit),trial_dur,n_trials(ss));
        lfp_mat{ss,cur_unit} = reshape(all_interp_lfps(cur_inds,cur_lfp),trial_dur,n_trials(ss));
    end
end

%%
stim_fs = 1e4/117.5;
stim_dur = 40/stim_fs;
stim_times = (0:40:320)/stim_fs;
lags = (1:trial_dur)*desired_dt;

flen_t = stim_dur;
tent_centers = [0:.009:flen_t];

tent_centers = round(tent_centers/desired_dt);
tbmat = construct_tent_bases(tent_centers,1);
[ntents,tblen] = size(tbmat);

tbmat = [zeros(ntents,tblen-1) tbmat];

stim_inds = 1+round(stim_times/desired_dt);
stim_inds(stim_inds > length(lags)) = [];
trial_inds = zeros(size(lags));
trial_inds(stim_inds) = 1;
trial_Tmat = zeros(length(lags),ntents);
for i = 1:ntents
    trial_Tmat(:,i) = conv(trial_inds,tbmat(i,:),'same');
end

%%
load ./expt3_repeat_preds
rep_dt = median(diff(rep_t));
rep_pred_r = rep_pred_r/rep_dt;
%% PLOT ACROSS-TRIAL AVERAGES
use_stims = [1 2 3 4 5 6 7 8 9];
stim_times = (0:40:320)/stim_fs;
tss = zeros(size(lags));
for ll = 1:length(lags)
    prev_stim = find(stim_times <= lags(ll),1,'last');
    tss(ll) = lags(ll)-stim_times(prev_stim);
end
% corr_type = 'spearman';
corr_type = 'pearson';

for uu = 1:length(use_stims)
    use_stim = use_stims(uu);
    for cc = 1:96
        
        tout = trial_Tmat*to_avg(cc,:)';
        psth=mean(binned_spk_mat{use_stim,cc}');
        
        %         tout_r2(uu,cc) = sum((psth-tout').^2)/length(psth);
        sout = squeeze(rep_pred_r(use_stim,cc,:))*desired_dt;
        sout = interp1(rep_t,sout,lags);
        %         ui = find(~isnan(sout));
        %         sout_r2(uu,cc) = sum((psth(ui)-sout(ui)).^2)/length(ui);
        
        tres = psth-tout';
        sres = psth-sout;
        sres(isnan(sres)) = nanmean(sres);
        
        time_r2(uu,cc) = 1-var(tres)/var(psth);
        stim_r2(uu,cc) = 1-var(sres)/var(psth);
        
        ss_phase_out = mean(sphasemod_out{use_stim,cc}');
        fh_phase_out = mean(phasemod_out{use_stim,cc}');
        ss_ampphase_out = mean(sampmod_out{use_stim,cc}');
        fh_ampphase_out = mean(ampmod_out{use_stim,cc}');
        
        uset = find(tss >= 0.15);
        %             uset = find(tss >= 0);
        
        [a,b] = corr(psth(uset)',ss_ampphase_out(uset)','type',corr_type);
        p_ss_c(uu,cc) = a;
        p_ss_p(uu,cc) = b;
        [a,b] = corr(psth(uset)',fh_ampphase_out(uset)','type',corr_type);
        p_fh_c(uu,cc) = a;
        p_fh_p(uu,cc) = b;
        
        [a,b] = corr(tres(uset)',ss_ampphase_out(uset)','type',corr_type);
        t_ss_c(uu,cc) = a;
        t_ss_p(uu,cc) = b;
        [a,b] = corr(tres(uset)',fh_ampphase_out(uset)','type',corr_type);
        t_fh_c(uu,cc) = a;
        t_fh_p(uu,cc) = b;
        
        [a,b] = corr(sres(uset)',ss_ampphase_out(uset)','type',corr_type);
        s_ss_c(uu,cc) = a;
        s_ss_p(uu,cc) = b;
        [a,b] = corr(sres(uset)',fh_ampphase_out(uset)','type',corr_type);
        s_fh_c(uu,cc) = a;
        s_fh_p(uu,cc) = b;
        
        [a,b] = corr(ss_ampphase_out(uset)',fh_ampphase_out(uset)','type',corr_type);
        ss_fh_c(uu,cc) = a;
        [a,b] = corr(ss_ampphase_out(uset)',ss_phase_out(uset)','type',corr_type);
        ss_ssp_c(uu,cc) = a;
        
        [a,b] = corr(ss_phase_out(uset)',fh_phase_out(uset)','type',corr_type);
        ssp_fhp_c(uu,cc) = a;
        %         [a,b] = corr(sres(uset)',ss_phase_out(uset)','type',corr_type);
        %         s_ssp_c(uu,cc) = a;
        %         s_ssp_p(uu,cc) = b;
        %         [a,b] = corr(sres(uset)',fh_phase_out(uset)','type',corr_type);
        %         s_fhp_c(uu,cc) = a;
        %         s_fhp_p(uu,cc) = b;
        
        
    end
end
%%
use_stim1 = 7;
use_stim2 = 8;
use_stim3 = 9;

stim_fs = 1e4/117.5;
stim_times = (40:40:320)/stim_fs;
lags = (1:trial_dur)*desired_dt;
close all
for cc = 1:96
    tout = trial_Tmat*to_avg(cc,:)';
    psth=mean(binned_spk_mat{use_stim1,cc}');
    
    subplot(3,1,1)
    shadedErrorBar(lags,mean(binned_spk_mat{use_stim1,cc}')*Fsd,std(binned_spk_mat{use_stim1,cc}'*Fsd)/sqrt(n_trials(use_stim1)));
    hold on
%     plot(lags,tout,'r','linewidth',2)
%     plot(rep_t,squeeze(rep_pred_r(use_stim1,cc,:))*desired_dt,'b')
    yl = ylim();
    for i = 1:length(stim_times)
        line(stim_times([i i]),yl,'color','k')
    end
    ylim([0 yl(2)])
    xlim([0 3.75])
    xlabel('Trial time (s)','fontsize',16)
    ylabel('Firing rate (Hz)','fontsize',16)
    
    psth=mean(binned_spk_mat{use_stim2,cc}');
    subplot(3,1,2)
        shadedErrorBar(lags,mean(binned_spk_mat{use_stim2,cc}')*Fsd,std(binned_spk_mat{use_stim2,cc}')*Fsd/sqrt(n_trials(use_stim2)));
    hold on
%     plot(lags,tout,'r','linewidth',2)
%     plot(rep_t,squeeze(rep_pred_r(use_stim2,cc,:))*desired_dt,'b')
    yl = ylim();
    for i = 1:length(stim_times)
        line(stim_times([i i]),yl,'color','k')
    end
    ylim([0 yl(2)])
    xlim([0 3.75])
    xlabel('Trial time (s)','fontsize',16)
    ylabel('Firing rate (Hz)','fontsize',16)

     psth=mean(binned_spk_mat{use_stim3,cc}');
      subplot(3,1,3)
     shadedErrorBar(lags,mean(binned_spk_mat{use_stim3,cc}')*Fsd,std(binned_spk_mat{use_stim3,cc}')*Fsd/sqrt(n_trials(use_stim3)));
    hold on
%     plot(lags,tout,'r','linewidth',2)
%     plot(rep_t,squeeze(rep_pred_r(use_stim3,cc,:))*desired_dt,'b')
    yl = ylim();
    for i = 1:length(stim_times)
        line(stim_times([i i]),yl,'color','k')
    end
    ylim([0 yl(2)])
    xlim([0 3.75])
    xlabel('Trial time (s)','fontsize',16)
    ylabel('Firing rate (Hz)','fontsize',16)

    gname = sprintf('Unit%d_PSTHs_simsacs',cc);
    fillPage(gcf,'Papersize',[12 15]);
    print(gname,'-dpng');
    
% %     pause
% clf
close all
end
%%
for cc = 1:96
    tout = trial_Tmat*to_avg(cc,:)';
    psth=mean(binned_spk_mat{use_stim,cc}');
    
    subplot(3,1,1)
    shadedErrorBar(lags,mean(binned_spk_mat{use_stim,cc}'),std(binned_spk_mat{use_stim,cc}')/sqrt(n_trials(use_stim)));
    hold on
    plot(lags,tout,'r','linewidth',2)
    plot(rep_t,squeeze(rep_pred_r(use_stim,cc,:))*desired_dt,'b')
    yl = ylim();
    for i = 1:length(stim_times)
        line(stim_times([i i]),yl,'color','k')
    end
    xlim([0 3.75])
    
    
    
    subplot(3,1,2)
    s1_out = mean(sampmod_out{use_stim,cc}');
    s2_out = mean(ampmod_out{use_stim,cc}');
    hold on
    %     shadedErrorBar(lags,mean(phasemod_out{use_stim,cc}'),std(phasemod_out{use_stim,cc}')/sqrt(n_trials(use_stim)),{'color','r'});
    shadedErrorBar(lags,mean(sampmod_out{use_stim,cc}'),std(sampmod_out{use_stim,cc}')/sqrt(n_trials(use_stim)),{'color','r'});
    shadedErrorBar(lags,mean(ampmod_out{use_stim,cc}'),std(ampmod_out{use_stim,cc}')/sqrt(n_trials(use_stim)),{'color','b'});
    yl = ylim();
    for i = 1:length(stim_times)
        line(stim_times([i i]),yl,'color','k')
    end
    xlim([0 3.75])
    
    ss_ampphase_out = mean(sampmod_out{use_stim,cc}');
    fh_ampphase_out = mean(ampmod_out{use_stim,cc}');
    
    subplot(3,1,3)
    %     shadedErrorBar(lags,mean(lfp_mat{use_stim,cc}'),std(lfp_mat{use_stim,cc}')/sqrt(n_trials(use_stim)),{'color','k'});
    
    psth = psth - mean(psth);
    tout = tout - mean(tout);
    plot(lags,zscore(psth-tout')','k')
    hold on
    plot(lags,zscore(ss_ampphase_out)-3,'r')
    
    plot(lags,zscore(fh_ampphase_out)-6,'b')
    
    yl = ylim();
    for i = 1:length(stim_times)
        line(stim_times([i i]),yl,'color','k')
    end
    xlim([0 3.75])
    ylim([-10 5])
    
    [t_ss_c(use_stim,cc)]
    [a,b] = corr(zscore(psth-tout')',zscore(ss_ampphase_out)','type','spearman')
    cc
    pause
    clf
end

%%

use_stim = 7;
for cc = 1:96
    
    all_tout(cc,:) = trial_Tmat*to_avg(cc,:)';
    
    all_psth(cc,:) =mean(binned_spk_mat{use_stim,cc}');
    
    all_s1_out(cc,:) = mean(sampmod_out{use_stim,cc}');
    all_s2_out(cc,:) = mean(ampmod_out{use_stim,cc}');
    all_lfp_out(cc,:) = mean(lfp_mat{use_stim,cc}');
    for bb = 2:7
        temp =stim_inds(bb):stim_inds(bb+1);
        temp(158:end) = [];
        stim_s1_out(cc,bb,:) = squeeze(all_s1_out(cc,temp));
        stim_lfp_out(cc,bb,:) = squeeze(all_lfp_out(cc,temp));
    end
end
norm_psth = bsxfun(@rdivide,all_psth,mean(all_psth,2));
use_stim = 8;
for cc = 1:96
    
    all_tout(cc,:) = trial_Tmat*to_avg(cc,:)';
    
    all_psth(cc,:) =mean(binned_spk_mat{use_stim,cc}');
    
    all_s1_out(cc,:) = mean(sampmod_out{use_stim,cc}');
    all_s2_out(cc,:) = mean(ampmod_out{use_stim,cc}');
    all_lfp_out(cc,:) = mean(lfp_mat{use_stim,cc}');
    for bb = 2:7
        temp =stim_inds(bb):stim_inds(bb+1);
        temp(158:end) = [];
        stim_s1_out2(cc,bb,:) = squeeze(all_s1_out(cc,temp));
        stim_lfp_out2(cc,bb,:) = squeeze(all_lfp_out(cc,temp));
    end
end
use_stim = 9;
for cc = 1:96
    
    all_tout(cc,:) = trial_Tmat*to_avg(cc,:)';
    
    all_psth(cc,:) =mean(binned_spk_mat{use_stim,cc}');
    
    all_s1_out(cc,:) = mean(sampmod_out{use_stim,cc}');
    all_s2_out(cc,:) = mean(ampmod_out{use_stim,cc}');
    all_lfp_out(cc,:) = mean(lfp_mat{use_stim,cc}');
    
    for bb = 2:7
        temp =stim_inds(bb):stim_inds(bb+1);
        temp(158:end) = [];
        stim_s1_out3(cc,bb,:) = squeeze(all_s1_out(cc,temp));
        stim_lfp_out3(cc,bb,:) = squeeze(all_lfp_out(cc,temp));
    end
end
full_stim_out = (stim_s1_out+stim_s1_out2+stim_s1_out3)/3;
full_stim_out = squeeze(mean(full_stim_out,2));
full_lfp_out = (stim_lfp_out+stim_lfp_out2+stim_lfp_out3)/3;
full_lfp_out = squeeze(mean(full_lfp_out,2));
%% PLOT ACROSS-TRIAL AVERAGES
use_stims = [4 5 6];
stim_fs = 1e4/117.5;
stim_times = (40:40:320)/stim_fs;
lags = (1:trial_dur)*desired_dt;


for cc = 1:96
    for i = 1:3
        %         subplot(3,1,i)
        %         shadedErrorBar(lags,mean(phasemod_out{use_stims(i),cc}'),std(phasemod_out{use_stims(i),cc}')/sqrt(n_trials(use_stims(i))),{'color','r'});
        %         hold on
        %         shadedErrorBar(lags,mean(ampmod_out{use_stims(i),cc}'),std(ampmod_out{use_stims(i),cc}')/sqrt(n_trials(use_stims(i))),{'color','b'});
        subplot(3,1,i)
        shadedErrorBar(lags,mean(exp(phasemod_out{use_stims(i),cc})'),std(exp(phasemod_out{use_stims(i),cc})')/sqrt(n_trials(use_stims(i))),{'color','r'});
        hold on
        shadedErrorBar(lags,mean(exp(ampmod_out{use_stims(i),cc})'),std(exp(ampmod_out{use_stims(i),cc})')/sqrt(n_trials(use_stims(i))),{'color','b'});
        yl = ylim();
        for i = 1:length(stim_times)
            line(stim_times([i i]),yl,'color','k')
        end
        xlim([0 3.75])
    end
    
    cc
    pause
    clf
end

%% PLOT ACROSS-TRIAL AVERAGES
use_stim = 7;
stim_fs = 1e4/117.5;
stim_times = (40:40:320)/stim_fs;
lags = (1:trial_dur)*desired_dt;


for cc = 1:96
    subplot(3,1,[1 2])
    imagesc(lags,1:n_trials(use_stim),ampmod_out{use_stim,cc}')
    
    subplot(3,1,3)
    shadedErrorBar(lags,mean(phasemod_out{use_stim,cc}'),std(phasemod_out{use_stim,cc}')/sqrt(n_trials(use_stim)),{'color','r'});
    hold on
    shadedErrorBar(lags,mean(ampmod_out{use_stim,cc}'),std(ampmod_out{use_stim,cc}')/sqrt(n_trials(use_stim)),{'color','b'});
    yl = ylim();
    for i = 1:length(stim_times)
        line(stim_times([i i]),yl,'color','k')
    end
    xlim([0 3.75])
    
    cc
    pause
    clf
end


%% COHERENCE BETWEEN PHASE MODELS AND REDISUALS

params.Fs = 1/desired_dt;
params.tapers = [6 11];
params.trialave = 1;
params.err = [2 .05];
params.fpass = [0 60];
win = 1;

use_stims = [7 8 9];
stim_times = (0:40:320)/stim_fs;
tss = zeros(size(lags));
for ll = 1:length(lags)
    prev_stim = find(stim_times <= lags(ll),1,'last');
    tss(ll) = lags(ll)-stim_times(prev_stim);
end
clear Coh Cerr
for cc = 1:96
    all_tres = [];
    all_sres = [];
    all_ss_pred = [];
    all_fh_pred = [];
    for uu = 1:length(use_stims)
        use_stim = use_stims(uu);
        
        tout = trial_Tmat*to_avg(cc,:)';
        psth=mean(binned_spk_mat{use_stim,cc}');
        
        sout = squeeze(rep_pred_r(use_stim,cc,:))*desired_dt;
        sout = interp1(rep_t,sout,lags);
        
        tres = psth-tout';
        sres = psth-sout;
        sres(isnan(sres)) = nanmean(sres);
        
        ss_ampphase_out = mean(sampmod_out{use_stim,cc}');
        fh_ampphase_out = mean(ampmod_out{use_stim,cc}');
        
        all_tres = [all_tres; tres(:)];
        all_sres = [all_sres; sres(:)];
        all_ss_pred = [all_ss_pred; ss_ampphase_out(:)];
        all_fh_pred = [all_fh_pred; fh_ampphase_out(:)];
        
    end
    
    all_tres = all_tres - mean(all_tres);
    all_sres = all_sres - mean(all_sres);
    all_ss_pred = all_ss_pred - mean(all_ss_pred);
    all_fh_pred = all_fh_pred - mean(all_fh_pred);
    
    [Coh(cc,:),phi,S12,S1,S2,f,confC,phistd,Cerr(cc,:,:)]=coherencysegc(all_sres,all_ss_pred,win,params);
end

%% FIT PHASE MODEL GLM AND COMPUTE R2
use_stims = [7 8 9];
stim_times = (0:40:320)/stim_fs;
tss = zeros(size(lags));
for ll = 1:length(lags)
    prev_stim = find(stim_times <= lags(ll),1,'last');
    tss(ll) = lags(ll)-stim_times(prev_stim);
end
for cc = 1:96
    all_ss_pred = [];
    all_fh_pred = [];
    all_ass_pred = [];
    all_afh_pred = [];
    all_psth = [];
    all_tout = [];
    all_sout = [];
    
    tout = log(trial_Tmat*to_avg(cc,:)');
    
    for uu = 1:length(use_stims)
        use_stim = use_stims(uu);
        cur_n_trials = size(sampmod_out{use_stim,cc},2);
        
        psth=mean(binned_spk_mat{use_stim,cc}');
        
        temp = squeeze(rep_pred_r(use_stim,cc,:))*desired_dt;
        temp = interp1(rep_t,temp,lags);
        temp(isnan(temp)) = nanmean(temp);
        sout(uu,:) = log(temp);
        
        ss_mean_out = mean(sampmod_out{use_stim,cc}');
        %         fh_mean_out = mean(ampmod_out{use_stim,cc}');
        for nn = 1:cur_n_trials
            
            cur_bs = binned_spk_mat{use_stim,cc}(:,nn);
            cur_ss_out = sampmod_out{use_stim,cc}(:,nn);
            %             cur_fh_out = ampmod_out{use_stim,cc}(:,nn);
            
            all_ss_pred = [all_ss_pred; cur_ss_out(:)];
            %             all_fh_pred = [all_fh_pred; cur_fh_out(:)];
            all_ass_pred = [all_ass_pred; ss_mean_out(:)];
            %             all_afh_pred = [all_afh_pred; fh_mean_out(:)];
            all_psth = [all_psth; cur_bs(:)];
            all_tout = [all_tout; tout(:)];
            all_sout = [all_sout; sout(uu,:)'];
            
        end
    end
    [tobeta_ss(cc,:)] = glmfit([all_tout(:)],all_psth(:),'poisson');
    [sobeta_ss(cc,:)] = glmfit([all_sout(:)],all_psth(:),'poisson');
    [sopabeta_ss(cc,:)] = glmfit([all_sout(:) all_ass_pred(:)],all_psth(:),'poisson');
    [tbeta_ss(cc,:)] = glmfit([all_ss_pred(:) all_tout(:)],all_psth(:),'poisson');
    [sbeta_ss(cc,:)] = glmfit([all_ss_pred(:) all_sout(:)],all_psth(:),'poisson');
    
    for uu = 1:length(use_stims)
        use_stim = use_stims(uu);
        cur_n_trials = size(sampmod_out{use_stim,cc},2);
        ss_mean_out = mean(sampmod_out{use_stim,cc}');
        %         fh_mean_out = mean(ampmod_out{use_stim,cc}');
        
        tpred_rate_mat = nan(size(binned_spk_mat{use_stim,cc}));
        spred_rate_mat = nan(size(binned_spk_mat{use_stim,cc}));
        for nn = 1:cur_n_trials
            cur_bs = binned_spk_mat{use_stim,cc}(:,nn);
            cur_ss_out = sampmod_out{use_stim,cc}(:,nn);
            
            cur_pred_out = glmval(tbeta_ss(cc,:)',[cur_ss_out(:) tout(:)],'log');
            tpred_rate_mat(:,nn) = cur_pred_out;
            cur_pred_out = glmval(sbeta_ss(cc,:)',[cur_ss_out(:) sout(uu,:)'],'log');
            spred_rate_mat(:,nn) = cur_pred_out;
        end
        cur_pred_out = glmval(tobeta_ss(cc,:)',[tout(:)],'log');
        topred_rate = cur_pred_out;
        cur_pred_out = glmval(sobeta_ss(cc,:)',[sout(uu,:)'],'log');
        sopred_rate = cur_pred_out;
        cur_pred_out = glmval(sopabeta_ss(cc,:)',[sout(uu,:)' ss_mean_out(:)],'log');
        sopapred_rate = cur_pred_out;
        
        psth=mean(binned_spk_mat{use_stim,cc}');
        
        res_psth = psth-topred_rate';
        tomod_r2(uu,cc) = 1-var(res_psth)/var(psth);
        res_psth = psth-sopred_rate';
        somod_r2(uu,cc) = 1-var(res_psth)/var(psth);
        res_psth = psth-sopapred_rate';
        sopamod_r2(uu,cc) = 1-var(res_psth)/var(psth);
        
        res_mat = binned_spk_mat{use_stim,cc} - tpred_rate_mat;
        res_psth = mean(res_mat,2);
        phasemod_r2(uu,cc) = 1-var(res_psth)/var(psth);
        
        res_mat = binned_spk_mat{use_stim,cc} - spred_rate_mat;
        res_psth = mean(res_mat,2);
        sphasemod_r2(uu,cc) = 1-var(res_psth)/var(psth);
    end
    
    
    
end

%% FIT PHASE MODEL GLM AND COMPUTE LL
use_stims = [1 2 3 7 8 9];
stim_times = (0:40:320)/stim_fs;
tss = zeros(size(lags));
for ll = 1:length(lags)
    prev_stim = find(stim_times <= lags(ll),1,'last');
    tss(ll) = lags(ll)-stim_times(prev_stim);
end
for cc = 1:96
    all_ss_pred = [];
    all_fh_pred = [];
    all_ass_pred = [];
    all_afh_pred = [];
    all_psth = [];
    all_tout = [];
    all_sout = [];
    
    tout = log(trial_Tmat*to_avg(cc,:)');
    
    for uu = 1:length(use_stims)
        use_stim = use_stims(uu);
        cur_n_trials = size(sampmod_out{use_stim,cc},2);
        
        psth=mean(binned_spk_mat{use_stim,cc}');
        
        temp = squeeze(rep_pred_r(use_stim,cc,:))*desired_dt;
        temp = interp1(rep_t,temp,lags);
        temp(isnan(temp)) = nanmean(temp);
        sout(uu,:) = log(temp);
        
        ss_mean_out = mean(sampmod_out{use_stim,cc}');
        %         fh_mean_out = mean(ampmod_out{use_stim,cc}');
        for nn = 1:cur_n_trials
            
            cur_bs = binned_spk_mat{use_stim,cc}(:,nn);
            cur_ss_out = sampmod_out{use_stim,cc}(:,nn);
            %             cur_fh_out = ampmod_out{use_stim,cc}(:,nn);
            
            all_ss_pred = [all_ss_pred; cur_ss_out(:)];
            %             all_fh_pred = [all_fh_pred; cur_fh_out(:)];
            all_ass_pred = [all_ass_pred; ss_mean_out(:)];
            %             all_afh_pred = [all_afh_pred; fh_mean_out(:)];
            all_psth = [all_psth; cur_bs(:)];
            all_tout = [all_tout; tout(:)];
            all_sout = [all_sout; sout(uu,:)'];
            
        end
    end
    [tobeta_ss(cc,:)] = glmfit([all_tout(:)],all_psth(:),'poisson');
    [sobeta_ss(cc,:)] = glmfit([all_sout(:)],all_psth(:),'poisson');
    [sopabeta_ss(cc,:)] = glmfit([all_sout(:) all_ass_pred(:)],all_psth(:),'poisson');
    [tbeta_ss(cc,:)] = glmfit([all_ss_pred(:) all_tout(:)],all_psth(:),'poisson');
    [sbeta_ss(cc,:)] = glmfit([all_ss_pred(:) all_sout(:)],all_psth(:),'poisson');
    avg_rate(cc) = mean(all_psth);
    
    all_to_pred = [];
    all_so_pred = [];
    all_sopa_pred = [];
    all_tp_pred = [];
    all_sp_pred = [];
    all_psth_pred = [];
    all_bs = [];
    all_tss = [];
    all_lags = [];
    for uu = 1:length(use_stims)
        use_stim = use_stims(uu);
        cur_n_trials = size(sampmod_out{use_stim,cc},2);
        ss_mean_out = mean(sampmod_out{use_stim,cc}');
        %         fh_mean_out = mean(ampmod_out{use_stim,cc}');
        
        tpred_rate_mat = nan(size(binned_spk_mat{use_stim,cc}));
        spred_rate_mat = nan(size(binned_spk_mat{use_stim,cc}));
        for nn = 1:cur_n_trials
            cur_bs = binned_spk_mat{use_stim,cc}(:,nn);
            cur_ss_out = sampmod_out{use_stim,cc}(:,nn);
            
            cur_pred_out = glmval(tbeta_ss(cc,:)',[cur_ss_out(:) tout(:)],'log');
            tpred_rate_mat(:,nn) = cur_pred_out;
            cur_pred_out = glmval(sbeta_ss(cc,:)',[cur_ss_out(:) sout(uu,:)'],'log');
            spred_rate_mat(:,nn) = cur_pred_out;
        end
        cur_pred_out = glmval(tobeta_ss(cc,:)',[tout(:)],'log');
        topred_rate = repmat(cur_pred_out,1,cur_n_trials);
        cur_pred_out = glmval(sobeta_ss(cc,:)',[sout(uu,:)'],'log');
        sopred_rate = repmat(cur_pred_out,1,cur_n_trials);
        cur_pred_out = glmval(sopabeta_ss(cc,:)',[sout(uu,:)' ss_mean_out(:)],'log');
        sopapred_rate = repmat(cur_pred_out,1,cur_n_trials);
        
        psth=mean(binned_spk_mat{use_stim,cc}');
        psth_rate = repmat(psth,1,cur_n_trials);
        
        rep_tss = repmat(tss(:),1,cur_n_trials);
        rep_lags = repmat(lags(:),1,cur_n_trials);
        
        all_to_pred = [all_to_pred; topred_rate(:)];
        all_so_pred = [all_so_pred; sopred_rate(:)];
        all_sopa_pred = [all_sopa_pred; sopapred_rate(:)];
        all_tp_pred = [all_tp_pred; tpred_rate_mat(:)];
        all_sp_pred = [all_sp_pred; spred_rate_mat(:)];
        all_psth_pred = [all_psth_pred; psth_rate(:)];
        all_bs = [all_bs; binned_spk_mat{use_stim,cc}(:)];
        all_tss = [all_tss; rep_tss(:)];
        all_lags = [all_lags; rep_lags(:)];
    end
    all_psth_pred(all_psth_pred < 1e-20) = 1e-20;
         
    null_pred_rate = ones(size(all_bs))*avg_rate(cc);
    
%     use_set = find(all_tss >= 0.15);
%      use_set = find(all_tss >= 0.15 & all_lags >= 0.45);
   
    use_set = find(all_tss >= 0);
    to_LL(cc) = -sum(all_bs(use_set).*log(all_to_pred(use_set)) - all_to_pred(use_set))/sum(all_bs(use_set));
    so_LL(cc) = -sum(all_bs(use_set).*log(all_so_pred(use_set)) - all_so_pred(use_set))/sum(all_bs(use_set));
    sopa_LL(cc) = -sum(all_bs(use_set).*log(all_sopa_pred(use_set)) - all_sopa_pred(use_set))/sum(all_bs(use_set));
    tp_LL(cc) = -sum(all_bs(use_set).*log(all_tp_pred(use_set)) - all_tp_pred(use_set))/sum(all_bs(use_set));
    sp_LL(cc) = -sum(all_bs(use_set).*log(all_sp_pred(use_set)) - all_sp_pred(use_set))/sum(all_bs(use_set));
    psth_LL(cc) = -sum(all_bs(use_set).*log(all_psth_pred(use_set)) - all_psth_pred(use_set))/sum(all_bs(use_set));
   null_LL(cc) = -sum(all_bs(use_set).*log(null_pred_rate(use_set)) - null_pred_rate(use_set))/sum(all_bs(use_set));

   use_set = find(all_tss >= 0.1);
   to_LL1(cc) = -sum(all_bs(use_set).*log(all_to_pred(use_set)) - all_to_pred(use_set))/sum(all_bs(use_set));
   so_LL1(cc) = -sum(all_bs(use_set).*log(all_so_pred(use_set)) - all_so_pred(use_set))/sum(all_bs(use_set));
   sopa_LL1(cc) = -sum(all_bs(use_set).*log(all_sopa_pred(use_set)) - all_sopa_pred(use_set))/sum(all_bs(use_set));
   tp_LL1(cc) = -sum(all_bs(use_set).*log(all_tp_pred(use_set)) - all_tp_pred(use_set))/sum(all_bs(use_set));
   sp_LL1(cc) = -sum(all_bs(use_set).*log(all_sp_pred(use_set)) - all_sp_pred(use_set))/sum(all_bs(use_set));
   psth_LL1(cc) = -sum(all_bs(use_set).*log(all_psth_pred(use_set)) - all_psth_pred(use_set))/sum(all_bs(use_set));
   null_LL1(cc) = -sum(all_bs(use_set).*log(null_pred_rate(use_set)) - null_pred_rate(use_set))/sum(all_bs(use_set));

   use_set = find(all_tss >= 0.2);
   to_LL2(cc) = -sum(all_bs(use_set).*log(all_to_pred(use_set)) - all_to_pred(use_set))/sum(all_bs(use_set));
   so_LL2(cc) = -sum(all_bs(use_set).*log(all_so_pred(use_set)) - all_so_pred(use_set))/sum(all_bs(use_set));
   sopa_LL2(cc) = -sum(all_bs(use_set).*log(all_sopa_pred(use_set)) - all_sopa_pred(use_set))/sum(all_bs(use_set));
   tp_LL2(cc) = -sum(all_bs(use_set).*log(all_tp_pred(use_set)) - all_tp_pred(use_set))/sum(all_bs(use_set));
   sp_LL2(cc) = -sum(all_bs(use_set).*log(all_sp_pred(use_set)) - all_sp_pred(use_set))/sum(all_bs(use_set));
   psth_LL2(cc) = -sum(all_bs(use_set).*log(all_psth_pred(use_set)) - all_psth_pred(use_set))/sum(all_bs(use_set));
   null_LL2(cc) = -sum(all_bs(use_set).*log(null_pred_rate(use_set)) - null_pred_rate(use_set))/sum(all_bs(use_set));
end

poss_LL_imp = null_LL-psth_LL;
to_LL_imp = null_LL-to_LL;
so_LL_imp = null_LL-so_LL;
sopa_LL_imp = null_LL-sopa_LL;
tp_LL_imp = null_LL-tp_LL;
sp_LL_imp = null_LL-sp_LL;

poss_LL_imp1 = null_LL1-psth_LL1;
to_LL_imp1 = null_LL1-to_LL1;
so_LL_imp1 = null_LL1-so_LL1;
sopa_LL_imp1 = null_LL1-sopa_LL1;
tp_LL_imp1 = null_LL1-tp_LL1;
sp_LL_imp1 = null_LL1-sp_LL1;

poss_LL_imp2 = null_LL2-psth_LL2;
to_LL_imp2 = null_LL2-to_LL2;
so_LL_imp2 = null_LL2-so_LL2;
sopa_LL_imp2 = null_LL2-sopa_LL2;
tp_LL_imp2 = null_LL2-tp_LL2;
sp_LL_imp2 = null_LL2-sp_LL2;

%%
load ./spk_stability_data
load(sprintf('Expt%dClusterTimesDetails.mat',1));

f1 = figure;
f2 = figure;
for cc = 1:96
    back_set = find(ClusterDetails{cc}.clst == 1);
    clust_set = find(ClusterDetails{cc}.clst == 2);
    figure(f1)
    plot(ClusterDetails{cc}.xy(:,1),ClusterDetails{cc}.xy(:,2),'.')
    hold on
    plot(ClusterDetails{cc}.xy(clust_set,1),ClusterDetails{cc}.xy(clust_set,2),'r.')
    
    figure(f2)
    plot(block_avg_rates(:,cc))
    
    mean(s_ss_c(7:9,cc))
    
    pause
    figure(f1);clf;figure(f2);clf;
end

