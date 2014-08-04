clear all
close all
cd ~/Data/bruce/G075/
load jbeG075Expts.mat

stim_fs = 1e4/117.5;
Fs = 3e4;
% dsf = 120;Fsd = Fs/dsf;
dsf = 80;Fsd = Fs/dsf;
niqf = Fsd/2;
[filt_b,filt_a] = butter(2,[1]/niqf,'high');
use_lfps = [1:8:96];
use_units = 1:96;
% use_lfps = [1 63];
backlag = round(0*Fsd);
forwardlag = round(4*Fsd);
lags = (-backlag:forwardlag);

scales = logspace(log10(8),log10(75),20);
scales = [scales 85 100 115 130 150 200 250];
scales = scales*60/dsf;
% scales = scales*2;
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);

%% COMPILE TRIAL-LEVEL DATA
% sim_sac_blocks = [10 14 18 22 24 28 32 37 40 42 50] - 6;
% sim_sac_blocks = [18 28 40] - 6; %sim sac
% sim_sac_blocks = [10 22 32 42] - 6; %sim sac a
sim_sac_blocks = [14 24 37 50] - 6; %sim sac b

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

%% COMPILE STIMULUS-LEVEL DATA
rel_start_times = (0:40:280)/stim_fs;
stim_dur = 40/stim_fs;
stims_per_trial = 8;
tot_n_trials = sum(n_trials);
all_stim_start_times = [];
all_stim_rel_num = [];
all_stim_trial_num = [];
all_stim_block_num = [];
for tt = 1:tot_n_trials
    cur_n_stims = floor(all_trial_durs(tt)/stim_dur);
    for st = 1:cur_n_stims
        all_stim_start_times = [all_stim_start_times; all_trial_start_times(tt)+rel_start_times(st)];
        all_stim_rel_num = [all_stim_rel_num; st];
        all_stim_trial_num = [all_stim_trial_num; tt];
        all_stim_block_num = [all_stim_block_num; all_expt_id(tt)];
    end
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
trial_cnt = 1;
desired_dt = 0.005;
all_interp_phasegrams = [];
all_interp_ampgrams = [];
all_binned_spks = [];
all_expt_t = [];
all_t_since_start = [];
all_trial_id = [];
for bb = 1:length(sim_sac_blocks)
    
    %%
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
        dV = [];
        dVf = [];
        cur_pt = 1;
        for pp = 1:nparts
            cur_range = cur_pt:(cur_pt + FullV.blklen(pp)-1);
            cur_range(cur_range > length(V)) = [];
            dV = [dV decimate(V(cur_range),dsf)];
            dVf = [dVf filtfilt(filt_b,filt_a,decimate(V(cur_range),dsf))]; %do some high-pass filtering
            cur_pt = cur_pt + FullV.blklen(pp);
        end
        Vmat(ll,:) = dV;
        %         Vmatf(ll,:) = dVf;
        temp = cwt(dV,scales,'cmor1-1');
        phasegrams(:,:,ll) = angle(temp)';
        ampgrams(:,:,ll) = abs(temp)';
    end
    t_ax = [];
    for pp = 1:nparts
        cur_t_ax = linspace(FullV.blkstart(pp),FullV.blkstart(pp)+FullV.blklen(pp)/Fs,FullV.blklen(pp));
        t_ax = [t_ax downsample(cur_t_ax,dsf)];
    end
    t_ax(size(Vmat,2)+1:end) = [];
    
    cur_stims = find(all_stim_block_num==sim_sac_blocks(bb));
    
    expt_t_axis = [];
    t_since_start = [];
    for tt = 1:length(cur_stims)
        cur_t_edges = all_stim_start_times(cur_stims(tt)):desired_dt:(all_stim_start_times(cur_stims(tt)) + stim_dur);
        binned_spks = nan(length(use_units),length(cur_t_edges)-1);
        for cc = 1:length(use_units)
            temp = histc(Clusters{use_units(cc)}.times,cur_t_edges);
            binned_spks(cc,:) = temp(1:end-1);
        end
        all_binned_spks = cat(1,all_binned_spks,binned_spks');
        all_trial_id = [all_trial_id; ones(length(cur_t_edges)-1,1)*trial_cnt];        
        trial_cnt = trial_cnt + 1;
        
        expt_t_axis = [expt_t_axis cur_t_edges(1:end-1)+desired_dt/2];
        t_since_start = [t_since_start cur_t_edges(1:end-1)-cur_t_edges(1)+desired_dt/2];
    end
    
    unwr_phasegram = unwrap(phasegrams);
    interp_phasegrams = interp1(t_ax,unwr_phasegram,expt_t_axis);
    interp_phasegrams = mod(interp_phasegrams+pi,2*pi)-pi;
    interp_ampgrams = interp1(t_ax,ampgrams,expt_t_axis);
    
    all_interp_phasegrams = cat(1,all_interp_phasegrams,interp_phasegrams);
    all_interp_ampgrams = cat(1,all_interp_ampgrams,interp_ampgrams);
    all_expt_t = cat(1,all_expt_t,expt_t_axis');
    all_t_since_start = cat(1,all_t_since_start,t_since_start');
end

%%
trial_set = unique(all_trial_id);
n_trials = length(trial_set);

xv_frac = 0.2;
n_xv_trials = round(n_trials*xv_frac);
xv_set = randperm(n_trials);
xv_set(n_xv_trials+1:end) = [];
xv_inds = find(ismember(all_trial_id,xv_set));
tr_inds = find(~ismember(all_trial_id,xv_set))';

bad_inds = find(isnan(all_interp_phasegrams(:,1,1)));
xv_inds(ismember(xv_inds,bad_inds)) = [];
tr_inds(ismember(tr_inds,bad_inds)) = [];

xv_inds_first = xv_inds(all_t_since_start(xv_inds) <= 0.2);
tr_inds_first = tr_inds(all_t_since_start(tr_inds) <= 0.2);
xv_inds_second = xv_inds(all_t_since_start(xv_inds) > 0.2 & all_t_since_start(xv_inds)  <= 0.4);
tr_inds_second = tr_inds(all_t_since_start(tr_inds) > 0.2 & all_t_since_start(tr_inds)  <= 0.4);

%%
flen_t = stim_dur;
tent_centers = [0:desired_dt:flen_t];

tent_centers = round(tent_centers/desired_dt);
tbmat = construct_tent_bases(tent_centers,1);
[ntents,tblen] = size(tbmat);

tbmat = [zeros(ntents,tblen-1) tbmat];

cur_trial_start_inds = find(all_t_since_start==desired_dt/2);
trial_inds = zeros(size(all_expt_t));
trial_inds(cur_trial_start_inds) = 1;
trial_Tmat = zeros(length(all_expt_t),ntents);
for i = 1:ntents
    trial_Tmat(:,i) = conv(trial_inds,tbmat(i,:),'same');
end

%%
close all
NL_type = 0; 

reg_params.dl1_ind = 3000;
reg_params.dl1_dep = 0;
reg_params.dl2_ind = 30000;
reg_params.dl2_dep = 0;
reg_params.dl2_freq_ind = 300;
reg_params.dl2_freq_dep = 0;
reg_params.dl2_time_ind = 0;
reg_params.dl2_time_dep = 0;
reg_params.l1_ind = 0;
reg_params.l1_dep = 0;
reg_params.is_phase = 1;

reg_params2.dl1_ind = 0;
reg_params2.dl1_dep = 0;
reg_params2.dl2_ind = 5000;
reg_params2.dl2_dep = 5000;
reg_params2.dl2_freq_ind = 0;
reg_params2.dl2_freq_dep = 0;
reg_params2.dl2_time_ind = 200;
reg_params2.dl2_time_dep = 0;
reg_params2.l1_ind = 0;
reg_params2.l1_dep = 0;
reg_params2.l1t_ind = 0;
reg_params2.l1t_dep = 0;
reg_params2.is_phase = 0;

silent = 1;
NT = length(all_expt_t);
used_inds = 1:NT;

nbins = 30;
pax = linspace(-pi,pi,nbins+1);
doub_pax = [pax(1:end-1) 2*pi+pax(1:end-1)];

new_phase_set = [reshape(cos(all_interp_phasegrams),NT,length(wfreqs)*length(use_lfps)) reshape(sin(all_interp_phasegrams),NT,length(wfreqs)*length(use_lfps))];
phase_elec_set = [repmat(1:length(use_lfps),1,length(wfreqs)) repmat(1:length(use_lfps),1,length(wfreqs))];
phase_elec_set = phase_elec_set(:);

for ll = 1:length(use_lfps)
    fprintf('Electrode %d of %d\n',ll,length(use_lfps));
    cur_units = find(nearest_lfps==ll);
    
    use_elecs = ll;
    use_set = find(ismember(phase_elec_set,use_elecs));
    Xmat = [new_phase_set(tr_inds,use_set)];

    Tmat = [trial_Tmat(tr_inds,:)];

    Xmat_xv = [new_phase_set(xv_inds,use_set)];
    Tmat_xv = [trial_Tmat(xv_inds,:)];

    XTmat = [Xmat Tmat];
    XTmat_xv = [Xmat_xv Tmat_xv];
    
    cur_alpha_phase = squeeze(all_interp_phasegrams(tr_inds,:,ll));
    Pmat = nan(length(tr_inds),length(wfreqs)*nbins);
    for ww = 1:length(wfreqs)
        cur_tb = tbrep(cur_alpha_phase(:,ww),pax);
        cur_tb(:,1) = cur_tb(:,1) + cur_tb(:,end);
        cur_tb(:,end) = [];
        cur_set = ((ww-1)*nbins+1):ww*nbins;
        Pmat(:,cur_set) = cur_tb;
    end
    cur_alpha_phase = squeeze(all_interp_phasegrams(xv_inds,:,ll));
    Pmat_xv = nan(length(xv_inds),length(wfreqs)*nbins);
    for ww = 1:length(wfreqs)
        cur_tb = tbrep(cur_alpha_phase(:,ww),pax);
        cur_tb(:,1) = cur_tb(:,1) + cur_tb(:,end);
        cur_tb(:,end) = [];
        cur_set = ((ww-1)*nbins+1):ww*nbins;
        Pmat_xv(:,cur_set) = cur_tb;
    end


%     Xmat_first = [new_phase_set(tr_inds_first,use_set) trial_Tmat(tr_inds_first,:)];
%     Xmat_second = [new_phase_set(tr_inds_second,use_set) trial_Tmat(tr_inds_second,:)];
%     Xmat_first_xv = [new_phase_set(xv_inds_first,use_set) trial_Tmat(xv_inds_first,:)];
%     Xmat_second_xv = [new_phase_set(xv_inds_second,use_set) trial_Tmat(xv_inds_second,:)];
%  
%     Tmat_first = [trial_Tmat(tr_inds_first,:)];
%     Tmat_second = [trial_Tmat(tr_inds_second,:)];
%     Tmat_first_xv = [trial_Tmat(xv_inds_first,:)];
%     Tmat_second_xv = [trial_Tmat(xv_inds_second,:)];

    for cc = 1:length(cur_units)
        fprintf('Cell %d of %d\n',cc,length(cur_units));
        
        Robs = all_binned_spks(tr_inds,cur_units(cc));
        tr_spkbns = convert_to_spikebins(Robs);
        Robs_first = all_binned_spks(tr_inds_first,cur_units(cc));
        Robs_second = all_binned_spks(tr_inds_second,cur_units(cc));
        
        stim_params = [0,0,ntents];
        klen = size(Tmat,2);
        K0 = zeros(klen+1,1);
        [fitp_to] = fit_GLM_phase_model(Tmat, Robs, K0, 1, stim_params, reg_params2,1,NL_type);
        stim_ind_to_filt(cur_units(cc),:) = fitp_to.k((1):(ntents));
        
%         stim_params = [length(wfreqs),length(use_elecs)];
%         klen = size(Xmat,2);
%         K0 = zeros(klen+1,1);
%         [fitp_sphase] = fit_GLM_phase_model(Xmat, Robs, K0,silent, stim_params, reg_params2,0,NL_type);
%         stim_sinphase_cfilt(cur_units(cc),:) = fitp_phase.k(1:length(use_elecs)*length(wfreqs));
%         stim_sinphase_sfilt(cur_units(cc),:) = fitp_phase.k((length(use_elecs)*length(wfreqs)+1):length(use_elecs)*length(wfreqs)*2);
% 
%         stim_params = [length(wfreqs),length(use_elecs), ntents];
%         klen = size(XTmat,2);
%         K0 = zeros(klen+1,1);
%         [fitp_fsphase] = fit_GLM_phase_model(XTmat, Robs, K0,silent, stim_params, reg_params2,1,NL_type);
%         stim_full_sinphase_cfilt(cur_units(cc),:) = fitp_fphase.k(1:length(use_elecs)*length(wfreqs));
%         stim_full_sinphase_sfilt(cur_units(cc),:) = fitp_fphase.k((length(use_elecs)*length(wfreqs)+1):length(use_elecs)*length(wfreqs)*2);
%         stim_full_tfilt(cur_units(cc),:) = fitp_fphase.k(length(use_elecs)*length(wfreqs)*2+1:end-1);
        
        stim_params = [nbins,length(wfreqs)];
        klen = size(Pmat,2);
        K0 = zeros(klen+1,1);
        [fitp_phase] = fit_GLM_phase_model(Pmat, Robs, K0,silent, stim_params, reg_params,0,NL_type);
        stim_ind_phase_pfilt(cur_units(cc),:) = fitp_phase.k(1:nbins*length(wfreqs));
% 
%         stim_params = [0,0,ntents];
%         klen = size(Tmat,2);
%         K0 = zeros(klen+1,1);
%         [fitp_first] = fit_GLM_phase_model(Tmat_first, Robs_first, K0, 1, stim_params, reg_params2,1,NL_type);
%         stim_first_to_filt(cur_units(cc),:) = fitp_first.k((1):(ntents));
%         [fitp_second] = fit_GLM_phase_model(Tmat_second, Robs_second, K0, 1, stim_params, reg_params2,1,NL_type);
%         stim_second_to_filt(cur_units(cc),:) = fitp_first.k((1):(ntents));
% 
%         stim_params = [length(wfreqs),length(use_elecs), ntents];
%         klen = size(Xmat_first,2);
%         K0 = zeros(klen+1,1);
%         [fitp_first_phase] = fit_GLM_phase_model(Xmat_first, Robs_first, K0,silent, stim_params, reg_params2,1,NL_type);
%         stim_first_cfilt(cur_units(cc),:) = fitp_first_phase.k(1:length(use_elecs)*length(wfreqs));
%         stim_first_sfilt(cur_units(cc),:) = fitp_first_phase.k((length(use_elecs)*length(wfreqs)+1):length(use_elecs)*length(wfreqs)*2);
%         stim_first_tfilt(cur_units(cc),:) = fitp_first_phase.k(length(use_elecs)*length(wfreqs)*2+1:end-1);
% 
%         [fitp_second_phase] = fit_GLM_phase_model(Xmat_second, Robs_second, K0,silent, stim_params, reg_params2,1,NL_type);
%         stim_second_cfilt(cur_units(cc),:) = fitp_second_phase.k(1:length(use_elecs)*length(wfreqs));
%         stim_second_sfilt(cur_units(cc),:) = fitp_second_phase.k((length(use_elecs)*length(wfreqs)+1):length(use_elecs)*length(wfreqs)*2);
%         stim_second_tfilt(cur_units(cc),:) = fitp_second_phase.k(length(use_elecs)*length(wfreqs)*2+1:end-1);
        
        
        xv_Robs = all_binned_spks(xv_inds,cur_units(cc));
        xv_Robs_first = all_binned_spks(xv_inds_first,cur_units(cc));
        xv_Robs_second = all_binned_spks(xv_inds_second,cur_units(cc));
        
        xv_to_pred_rate = Tmat_xv*fitp_to.k(1:end-1) + fitp_to.k(end);
        if NL_type == 0
            xv_to_pred_rate = log(1+exp(xv_to_pred_rate));
        else
            xv_to_pred_rate = exp(xv_to_pred_rate);
        end
        xv_to_LL(cur_units(cc)) = -sum(xv_Robs.*log(xv_to_pred_rate)-xv_to_pred_rate)/sum(xv_Robs);

        xv_phase_pred_rate = Pmat_xv*fitp_phase.k(1:end-1) + fitp_phase.k(end);
        if NL_type == 0
            xv_phase_pred_rate = log(1+exp(xv_phase_pred_rate));
        else
            xv_phase_pred_rate = exp(xv_phase_pred_rate);
        end
        xv_phase_LL(cur_units(cc)) = -sum(xv_Robs.*log(xv_phase_pred_rate)-xv_phase_pred_rate)/sum(xv_Robs);

        xv_phase_pred_rate = Xmat_xv*fitp_sphase.k(1:end-1) + fitp_sphase.k(end);
        if NL_type == 0
            xv_phase_pred_rate = log(1+exp(xv_phase_pred_rate));
        else
            xv_phase_pred_rate = exp(xv_phase_pred_rate);
        end
        xv_sphase_LL(cur_units(cc)) = -sum(xv_Robs.*log(xv_phase_pred_rate)-xv_phase_pred_rate)/sum(xv_Robs);

        xv_fphase_pred_rate = XTmat_xv*fitp_fsphase.k(1:end-1) + fitp_fsphase.k(end);
        if NL_type == 0
            xv_fphase_pred_rate = log(1+exp(xv_fphase_pred_rate));
        else
            xv_fphase_pred_rate = exp(xv_fphase_pred_rate);
        end
        xv_fsphase_LL(cur_units(cc)) = -sum(xv_Robs.*log(xv_fphase_pred_rate)-xv_fphase_pred_rate)/sum(xv_Robs);
        
        avg_rate = mean(Robs);
        null_pred = avg_rate*ones(size(xv_Robs));
        xv_null_LL(cur_units(cc)) = -sum(xv_Robs.*log(null_pred)-null_pred)/sum(xv_Robs);
        
%          xv_first_to_pred_rate = Tmat_first_xv*fitp_first.k(1:end-1) + fitp_first.k(end);
%         if NL_type == 0
%             xv_first_to_pred_rate = log(1+exp(xv_first_to_pred_rate));
%         else
%             xv_first_to_pred_rate = exp(xv_first_to_pred_rate);
%         end
%         xv_to_first_LL(cur_units(cc)) = -sum(xv_Robs_first.*log(xv_first_to_pred_rate)-xv_first_to_pred_rate)/sum(xv_Robs_first);
% 
%                  xv_second_to_pred_rate = Tmat_second_xv*fitp_second.k(1:end-1) + fitp_second.k(end);
%         if NL_type == 0
%             xv_second_to_pred_rate = log(1+exp(xv_second_to_pred_rate));
%         else
%             xv_second_to_pred_rate = exp(xv_second_to_pred_rate);
%         end
%         xv_to_second_LL(cur_units(cc)) = -sum(xv_Robs_second.*log(xv_second_to_pred_rate)-xv_second_to_pred_rate)/sum(xv_Robs_second);
% 
%         xv_first_pred_rate = Xmat_first_xv*fitp_first_phase.k(1:end-1) + fitp_first_phase.k(end);
%         if NL_type == 0
%             xv_first_pred_rate = log(1+exp(xv_first_pred_rate));
%         else
%             xv_first_pred_rate = exp(xv_first_pred_rate);
%         end
%         xv_first_phase_LL(cur_units(cc)) = -sum(xv_Robs_first.*log(xv_first_pred_rate)-xv_first_pred_rate)/sum(xv_Robs_first);
%         
%         xv_second_pred_rate = Xmat_second_xv*fitp_second_phase.k(1:end-1) + fitp_second_phase.k(end);
%         if NL_type == 0
%             xv_second_pred_rate = log(1+exp(xv_second_pred_rate));
%         else
%             xv_second_pred_rate = exp(xv_second_pred_rate);
%         end
%         xv_second_phase_LL(cur_units(cc)) = -sum(xv_Robs_second.*log(xv_second_pred_rate)-xv_second_pred_rate)/sum(xv_Robs_second);
% 
% 
%         avg_rate = mean(Robs_first);
%         null_pred = avg_rate*ones(size(xv_Robs_first));
%         xv_firstnull_LL(cur_units(cc)) = -sum(xv_Robs_first.*log(null_pred)-null_pred)/sum(xv_Robs_first);
%         
%         avg_rate = mean(Robs_second);
%         null_pred = avg_rate*ones(size(xv_Robs_second));
%         xv_secondnull_LL(cur_units(cc)) = -sum(xv_Robs_second.*log(null_pred)-null_pred)/sum(xv_Robs_second);
    end
end

%%
save simsac_phase_mod_data_first pax wfreqs nbins stim_ind* 
%%
stim_ind_ampkern = sqrt(stim_full_sinphase_cfilt.^2+stim_full_sinphase_sfilt.^2);
stim_ind_phasekern = -atan2(stim_full_sinphase_cfilt,stim_full_sinphase_sfilt)+pi/2;

stim_ind_ampkern_first = sqrt(stim_first_cfilt.^2+stim_first_sfilt.^2);
stim_ind_phasekern_first = -atan2(stim_first_cfilt,stim_first_sfilt)+pi/2;

stim_ind_ampkern_second = sqrt(stim_second_cfilt.^2+stim_second_sfilt.^2);
stim_ind_phasekern_second = -atan2(stim_second_cfilt,stim_second_sfilt)+pi/2;

