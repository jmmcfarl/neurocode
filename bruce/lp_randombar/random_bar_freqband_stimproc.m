clear all
close all

addpath('~/James_scripts/CircStat2011f/');

ExptNum = 232;
cd(['~/Data/bruce/M' num2str(ExptNum)]);
load ./random_bar_eyedata_ftime.mat

load(['lemM' num2str(ExptNum) 'Expts.mat']);
load ./bar_params.mat

if ExptNum == 235
    bar_expts(bar_expts==51) = []; %this has different set of stim positions
end
if ExptNum == 239
    bar_expts(bar_expts==40) = []; %this has different set of stim positions
end

%%
axis_or = Expts{bar_expts(1)}.Stimvals.or*pi/180; %in radians

Fs = 1000;
dsf = 2;
Fsd = Fs/dsf;

niqf = Fs/2;
[bb,aa] = butter(2,[1 100]/niqf);

[b_delta,a_delta] = butter(2,[3 6]/niqf);
[b_alpha,a_alpha] = butter(2,[8 15]/niqf);
[b_gamma,a_gamma] = butter(2,[30 60]/niqf);

new_dt = .002;
dsfrac = 5;

tbspace = 3;
nLags = 15*dsfrac/tbspace;
n_bar_pos = bar_params.n_bars;
n_lin_dims = length(bar_expts)-1;
stim_params = NIMcreate_stim_params([nLags n_bar_pos],new_dt,1,tbspace,n_lin_dims);
flen = nLags*tbspace;

% bar_axis = [bar_params.bar_axis (bar_params.bar_axis(end) + 0.125)]; %add extra bin edge for histing
bar_axis = bar_params.bar_axis;
%%
load ./CellList.mat
good_sus = find(all(CellList(bar_expts,:,1) > 0));
use_sus = 1:24;
%%
% load ./un_bar_pos.mat
% un_bar_pos(1:2) = [];

full_Xmat = [];
full_spkbinned = [];
full_exptvec = [];
full_exptvec_new = [];
full_trialvec_new = [];
full_taxis = [];
full_taxis_new = [];
full_old_t_inds = [];
full_bar_pos = [];
full_t_since_tstart = [];
full_used_inds = [];
full_used_inds_new = [];
for ee = 1:length(bar_expts);
    fprintf('Processing expt %d of %d\n',ee,length(bar_expts));
    fname = sprintf('Expt%dClusterTimes.mat',bar_expts(ee));
    load(fname);
    
    trial_durs{ee} = [Expts{bar_expts(ee)}.Trials(:).dur];
    n_trials(ee) = length(Expts{bar_expts(ee)}.Trials);
    all_t_axis = [];
    all_used_inds = [];
    all_t_axis_new = [];
    all_old_t_inds = [];
    all_used_inds_new = [];
    all_expt_vec = [];
    all_trial_vec = [];
    all_bar_Op = [];
    all_bar_Xmat = [];
    all_binned_spikes = [];
    all_used_inds = [];
    
    all_t_since_tstart = [];
    for tt = 1:n_trials(ee)
        
        cur_t_axis = [Expts{bar_expts(ee)}.Trials(tt).Start]/1e4;
        cur_bar_Op = [Expts{bar_expts(ee)}.Trials(tt).Op];
        
        cur_bar_Op_new = repmat(cur_bar_Op,1,dsfrac)';
        cur_bar_Op_new = cur_bar_Op_new(:);
        
        cur_t_edges = [cur_t_axis; Expts{bar_expts(ee)}.Trials(tt).End(end)/1e4];
        cur_t_edges_new = Expts{bar_expts(ee)}.Trials(tt).Start/1e4:new_dt:Expts{bar_expts(ee)}.Trials(tt).End(end)/1e4;
        cur_t_axis_new = 0.5*cur_t_edges_new(1:end-1) + 0.5*cur_t_edges_new(2:end);
        
        cur_bar_Op_new(length(cur_t_axis_new)+1:end) = [];
        
        cur_binned_spks = nan(length(cur_t_axis_new),length(use_sus));
        for cc = 1:length(use_sus)
            cur_hist = histc(Clusters{use_sus(cc)}.times,cur_t_edges_new);
            cur_binned_spks(:,cc) = cur_hist(1:end-1);
        end
                
        cur_bar_mat = zeros(length(cur_t_axis_new),n_bar_pos);
        for b = 1:n_bar_pos
            cur_set = find(cur_bar_Op_new==bar_axis(b));
%             cur_set = find(cur_bar_Op_new >= bar_axis(b) & cur_bar_Op_new < bar_axis(b+1));
            cur_bar_mat(cur_set,b) = 1;
        end
%         bar_Xmat = makeStimRows(cur_bar_mat,flen);
        bar_Xmat = create_time_embedding(cur_bar_mat,stim_params);
        if size(bar_Xmat,1) ~= size(cur_bar_mat,1)
            pause
        end
        
        cur_used_inds_new = ones(length(cur_t_axis_new),1);
        cur_used_inds_new(1:flen) = 0;
        
        cur_t_since_tstart = cur_t_axis_new - cur_t_axis(1);
        
        cur_used_inds_new(end-round(0.2/new_dt):end) = 0; %reduce phase estimation edge artifacts
        
        all_used_inds_new = [all_used_inds_new; cur_used_inds_new(:)];
        all_t_axis = [all_t_axis; cur_t_axis(:)];
        all_t_since_tstart = [all_t_since_tstart; cur_t_since_tstart(:)];
        all_t_axis_new = [all_t_axis_new; cur_t_axis_new(:)];
        all_binned_spikes = [all_binned_spikes; cur_binned_spks];
        all_expt_vec = [all_expt_vec; ones(length(cur_t_axis_new),1)*bar_expts(ee)];
        all_trial_vec = [all_trial_vec; ones(length(cur_t_axis_new),1)*tt];
        all_bar_Op = [all_bar_Op; cur_bar_Op_new(:)];
        all_bar_Xmat = [all_bar_Xmat; bar_Xmat];
    end
    
    cur_blink_start_times = all_eye_ts{ee}(all_blink_startinds{ee});
    cur_blink_stop_times = all_eye_ts{ee}(all_blink_stopinds{ee});
%     cur_blink_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),cur_blink_start_times));
%     cur_blink_stop_inds = round(interp1(all_t_axis,1:length(all_t_axis),cur_blink_stop_times));
%     cur_poss_blinks = find(~isnan(cur_blink_start_inds) & ~isnan(cur_blink_stop_inds));
%     
%     all_blink_inds = zeros(size(all_t_axis));
%     for i = 1:length(cur_poss_blinks)
%         all_blink_inds(cur_blink_start_inds(cur_poss_blinks(i)):cur_blink_stop_inds(cur_poss_blinks(i))) = 1;
%     end
    
    cur_blink_start_inds = round(interp1(all_t_axis_new,1:length(all_t_axis_new),cur_blink_start_times));
    cur_blink_stop_inds = round(interp1(all_t_axis_new,1:length(all_t_axis_new),cur_blink_stop_times));
    cur_poss_blinks = find(~isnan(cur_blink_start_inds) & ~isnan(cur_blink_stop_inds));
    all_blink_inds_new = zeros(size(all_t_axis_new));
    for i = 1:length(cur_poss_blinks)
        all_blink_inds_new(cur_blink_start_inds(cur_poss_blinks(i)):cur_blink_stop_inds(cur_poss_blinks(i))) = 1;
    end
    
    all_used_inds_new(all_blink_inds_new == 1) = 0;
    all_used_inds_new = logical(all_used_inds_new);
    
    full_Xmat = [full_Xmat; all_bar_Xmat];
    full_spkbinned = [full_spkbinned; all_binned_spikes];
    full_exptvec = [full_exptvec; ones(length(all_used_inds),1)*ee];
    full_exptvec_new = [full_exptvec_new; ones(length(all_used_inds_new),1)*ee];
    full_taxis = [full_taxis; all_t_axis];
    full_taxis_new = [full_taxis_new; all_t_axis_new];
    full_t_since_tstart = [full_t_since_tstart; all_t_since_tstart];
    full_trialvec_new = [full_trialvec_new; all_trial_vec];
    full_bar_pos = [full_bar_pos; all_bar_Op];
    
    full_used_inds_new = [full_used_inds_new; all_used_inds_new];
    full_used_inds = [full_used_inds; all_used_inds];
    
end

%%
full_lfps = [];
full_deltaphase = [];
full_alphaphase = [];
full_gammaphase = [];
for ee = 1:length(bar_expts);
    fprintf('Expt %d of %d\n',ee,length(bar_expts));
    fname = sprintf('lemM%dA.%d.lfp.mat',ExptNum,bar_expts(ee));
    load(fname);
    
    Fs = 1/LFP.Header.CRsamplerate;
    
    n_trials(ee) = length(LFP.Trials);
    %     lfp_trial_starts = [LFP.Trials(:).Start]/1e4;
    lfp_trial_starts = [LFP.Trials(:).ftime]/1e4;
    lfp_trial_ends = [LFP.Trials(:).End]/1e4;
    expt_lfp_t_axis = [];
    expt_lfps = [];
    expt_delta = [];
    expt_alpha = [];
    expt_gamma = [];
    for tt = 1:n_trials(ee)
        
        cur_npts = size(LFP.Trials(tt).LFP,1);
        cur_t_end(tt) = lfp_trial_starts(tt)+(cur_npts-1)/Fs;
        cur_t_axis = lfp_trial_starts(tt):1/Fs:cur_t_end(tt);
        
        if ~isempty(expt_lfp_t_axis)
            cur_sp = find(cur_t_axis > max(expt_lfp_t_axis),1,'first');
        else
            cur_sp = 1;
        end
        cur_t_axis = cur_t_axis(cur_sp:end);
        cur_LFP = [LFP.Trials(tt).LFP];
        cur_LFP = cur_LFP(cur_sp:end,:);
        cur_delta = filtfilt(b_delta,a_delta,cur_LFP);
        cur_alpha = filtfilt(b_alpha,a_alpha,cur_LFP);
        cur_gamma = filtfilt(b_gamma,a_gamma,cur_LFP);
        cur_LFP = filtfilt(bb,aa,cur_LFP);               
        
        expt_lfp_t_axis = [expt_lfp_t_axis; cur_t_axis(:)];
        expt_lfps = [expt_lfps; cur_LFP];
        expt_delta = [expt_delta; cur_delta];
        expt_alpha = [expt_alpha; cur_alpha];
        expt_gamma = [expt_gamma; cur_gamma];
    end
    
    cur_set = find(full_exptvec_new==ee);
    interp_lfps = interp1(expt_lfp_t_axis,expt_lfps,full_taxis_new(cur_set));
    interp_delta = interp1(expt_lfp_t_axis,expt_delta,full_taxis_new(cur_set));
    interp_alpha = interp1(expt_lfp_t_axis,expt_alpha,full_taxis_new(cur_set));
    interp_gamma = interp1(expt_lfp_t_axis,expt_gamma,full_taxis_new(cur_set));
    
    interp_gamma_phase = (angle(hilbert(interp_gamma)));
    interp_delta_phase = (angle(hilbert(interp_delta)));
    interp_alpha_phase = (angle(hilbert(interp_alpha)));

%     %compute CSD
%     vars.Fs = Fs;
%     vars.BrainBound = 1;
%     vars.ChanSep = 0.05;
%      vars.diam = 2; %0.5
%     CSD = PettersenCSD(expt_lfps','spline',vars)';
%     expt_lfps = CSD;
%


full_lfps = [full_lfps; interp_lfps];
    full_deltaphase = cat(1,full_deltaphase, interp_delta_phase);
    full_alphaphase = cat(1,full_alphaphase, interp_alpha_phase);
    full_gammaphase = cat(1,full_gammaphase, interp_gamma_phase);
end

%% Parse into training and XV sets
[c,ia,ic] = unique([full_exptvec_new full_trialvec_new],'rows');
n_trials = length(ia);

% rp = randperm(n_trials);
% rand_trial_vec = rp(ic);
% [~,ind_shuff] = sort(rand_trial_vec);

xv_frac = 0.;
n_xv_trials = round(n_trials*xv_frac);
xv_tset = randperm(n_trials);
xv_tset(n_xv_trials+1:end) = [];
tr_tset = find(~ismember(1:n_trials,xv_tset));
xv_inds = find(ismember(ic,xv_tset));
tr_inds = find(ismember(ic,tr_tset))';

xv_inds(full_used_inds_new(xv_inds)==0) = [];
tr_inds(full_used_inds_new(tr_inds) == 0) = [];

% %  Probably have to debug the trial-shuffling procedure after making
% changes (to include edge samples in initial data parsing)
% tr_inds_shuffle = ind_shuff(tr_inds);
% xv_inds_shuffle = ind_shuff(xv_inds);


%% Make indicator matrix for experiment number
Xexpt = zeros(length(full_taxis_new),length(bar_expts)-1);
for i = 1:length(bar_expts)-1
    cur_set = find(full_exptvec_new==i);
    Xexpt(cur_set,i) = 1;
end

%% Fit STRFs
reg_params = NIMcreate_reg_params('lambda_d2XT',100);
for cc = 1:24
    fprintf('Fitting cell %d\n',cc);
    Robs = full_spkbinned(tr_inds,cc);
    
    mod_signs = [1]; %determines whether input is exc or sup (doesn't matter in the linear case)
    NL_types = {'lin'}; %define subunit as linear
    silent = 1; %display model optimization iterations
    
    fit0(cc) = NIMinitialize_model(stim_params,mod_signs,NL_types,reg_params); %initialize NIM
    fit0(cc) = NIMfit_filters(fit0(cc),Robs,full_Xmat(tr_inds,:),Xexpt(tr_inds,:),[],silent); %fit stimulus filters
    
end


%% GET SACCADE TIMES
all_sac_inds = [];
all_microsac_inds = [];
all_firstsac_inds = [];
all_secondsac_inds = [];
% full_eye_ts = [];
full_eye_speed = [];
for ee = 1:length(bar_expts)
    expt_trial_starts = [Expts{bar_expts(ee)}.Trials(:).TrialStart]/1e4;
    expt_trial_end = [Expts{bar_expts(ee)}.Trials(:).TrueEnd]/1e4;
    
    sac_start_times = all_eye_ts{ee}(all_sac_startinds{ee});
    sac_stop_times = all_eye_ts{ee}(all_sac_stopinds{ee});
    sac_amps = all_sac_peakamps{ee};
    
    
    cur_sac_set = find(all_expt_num==ee);
    if length(cur_sac_set) ~= length(sac_start_times)
        error('saccade mismatch')
    end
    intrial_sac_set = find(~isnan(all_trial_num(cur_sac_set)));

    infirst_sac_set = find(ismember(cur_sac_set,all_isfirst));
    insecond_sac_set = find(ismember(cur_sac_set,all_issecond));
    big_sac_set = [infirst_sac_set; insecond_sac_set];
    micro_sac_set = find(sac_amps' < 60 & ~isnan(all_trial_num(cur_sac_set)));
    micro_sac_set(ismember(micro_sac_set,big_sac_set)) = [];

    sac_times = sac_start_times(intrial_sac_set);
    first_sac_times = sac_start_times(infirst_sac_set);
    second_sac_times = sac_start_times(insecond_sac_set);
    micro_sac_times = sac_start_times(micro_sac_set);
    
    cur_e_set = find(full_exptvec_new == ee);
    
        interp_eye_speed = interp1(all_eye_ts{ee},all_eye_speed{ee},full_taxis_new(cur_e_set));
    full_eye_speed = [full_eye_speed; interp_eye_speed];
        
    sac_inds = round(interp1(full_taxis_new(cur_e_set),1:length(cur_e_set),sac_times));
    sac_inds(isnan(sac_inds)) = [];
    micro_sac_inds = round(interp1(full_taxis_new(cur_e_set),1:length(cur_e_set),micro_sac_times));
    micro_sac_inds(isnan(micro_sac_inds)) = [];
    first_sac_inds = round(interp1(full_taxis_new(cur_e_set),1:length(cur_e_set),first_sac_times));
    first_sac_inds(isnan(first_sac_inds)) = [];
    second_sac_inds = round(interp1(full_taxis_new(cur_e_set),1:length(cur_e_set),second_sac_times));
    second_sac_inds(isnan(second_sac_inds)) = [];
    
    all_sac_inds = [all_sac_inds; cur_e_set(sac_inds(:))];
    all_microsac_inds = [all_microsac_inds; cur_e_set(micro_sac_inds(:))];
    all_firstsac_inds = [all_firstsac_inds; cur_e_set(first_sac_inds(:))];
    all_secondsac_inds = [all_secondsac_inds; cur_e_set(second_sac_inds(:))];
end
all_bigsac_inds = sort([all_firstsac_inds; all_secondsac_inds]);

%% Create saccade predictor mats
fNT = length(full_taxis_new);

tbspace = round(0.01/new_dt);
tent_filter = [(1:tbspace)/tbspace 1-(1:tbspace-1)/tbspace]/tbspace;
poss_lags = round(-0.2/new_dt):tbspace:(round(0.5/new_dt));
n_saclags = length(poss_lags);

all_bigsac_inds = sort([all_firstsac_inds; all_secondsac_inds]);
all_bigsac_indvec = zeros(size(full_taxis_new));
all_bigsac_indvec(all_bigsac_inds) = 1;

%apply tb filter
filtered = zeros(size(all_bigsac_indvec));
for i = 1:length(tent_filter)
    filtered = filtered + shift_mat_zpad(all_bigsac_indvec,i-tbspace,1)*tent_filter(i);
end
all_bigsac_indvec = filtered;

bigsac_Xmat = zeros(fNT, n_saclags);
for n = 1:n_saclags
    bigsac_Xmat(:,n) = shift_mat_zpad( all_bigsac_indvec, poss_lags(n), 1);
end


%%
clear temp_kern 
for cc = 1:24
    k_mat = reshape(fit0(cc).mods(1).filtK,nLags,n_bar_pos);
    temp_kern(cc,:) = sqrt(mean(k_mat.^2,2));
    [~,best_lag(cc)] = max(temp_kern(cc,:));
end
pref_lag = nLags - best_lag - 1;
[bps,fls] = meshgrid(1:n_bar_pos,1:nLags);
%%
close all

n_phasebins = 6;
phase_bin_edges = linspace(-pi,pi,n_phasebins+1);
phase_bin_cents = 0.5*phase_bin_edges(1:end-1) + 0.5*phase_bin_edges(2:end);

nbins = 20;
bin_edges = linspace(0,100,nbins+1);
bin_cents = 0.5*bin_edges(1:end-1) + 0.5*bin_edges(2:end);
clear cur_stas_g g_bin_cents marg_phase_dep *stas* gamma* alpha* delta*
for cc = 1:24
    cc
    cur_set = find(fls == best_lag(cc));
    cur_X = full_Xmat(tr_inds,cur_set);
    
    Robs = full_spkbinned(tr_inds,cc);
    ov_sta(cc,:) = Robs'*cur_X./sum(cur_X);
    
    gout = full_Xmat(tr_inds,:)*fit0(cc).mods(1).filtK;
    pp = prctile(gout,bin_edges);
    g_bin_cents(cc,:) = 0.5*pp(1:end-1) + 0.5*pp(2:end);

    gamma_stas(cc,:,:) = nan(n_phasebins,n_bar_pos);
    gamma_stas_g(cc,:,:) = nan(n_phasebins,nbins);
    for tt = 1:n_phasebins
        cur_inds = find(full_gammaphase(tr_inds,cc) >= phase_bin_edges(tt) & ...
            full_gammaphase(tr_inds,cc) < phase_bin_edges(tt+1));
        gamma_stas(cc,tt,:) = (Robs(cur_inds)'*cur_X(cur_inds,:))./sum(cur_X(cur_inds,:));
        gamma_marg_dep(cc,tt) = mean(Robs(cur_inds));
         for ii = 1:nbins
            cur_set = cur_inds(gout(cur_inds) >= pp(ii) & gout(cur_inds) < pp(ii+1));
            gamma_stas_g(cc,tt,ii) = mean(Robs(cur_set));
        end
   end
    
    alpha_stas(cc,:,:) = nan(n_phasebins,n_bar_pos);
    alpha_stas_g(cc,:,:) = nan(n_phasebins,nbins);
    for tt = 1:n_phasebins
        cur_inds = find(full_alphaphase(tr_inds,cc) >= phase_bin_edges(tt) & ...
            full_alphaphase(tr_inds,cc) < phase_bin_edges(tt+1));
        alpha_stas(cc,tt,:) = (Robs(cur_inds)'*cur_X(cur_inds,:))./sum(cur_X(cur_inds,:));
        alpha_marg_dep(cc,tt) = mean(Robs(cur_inds));
         for ii = 1:nbins
            cur_set = cur_inds(gout(cur_inds) >= pp(ii) & gout(cur_inds) < pp(ii+1));
            alpha_stas_g(cc,tt,ii) = mean(Robs(cur_set));
        end
   end
    
    
    delta_stas(cc,:,:) = nan(n_phasebins,n_bar_pos);
    delta_stas_g(cc,:,:) = nan(n_phasebins,nbins);
    for tt = 1:n_phasebins
        cur_inds = find(full_deltaphase(tr_inds,cc) >= phase_bin_edges(tt) & ...
            full_deltaphase(tr_inds,cc) < phase_bin_edges(tt+1));
        delta_stas(cc,tt,:) = (Robs(cur_inds)'*cur_X(cur_inds,:))./sum(cur_X(cur_inds,:));
        delta_marg_dep(cc,tt) = mean(Robs(cur_inds));
         for ii = 1:nbins
            cur_set = cur_inds(gout(cur_inds) >= pp(ii) & gout(cur_inds) < pp(ii+1));
            delta_stas_g(cc,tt,ii) = mean(Robs(cur_set));
        end
   end
end

delta_marg_dep = bsxfun(@rdivide,delta_marg_dep,mean(delta_marg_dep,2));
alpha_marg_dep = bsxfun(@rdivide,alpha_marg_dep,mean(alpha_marg_dep,2));
gamma_marg_dep = bsxfun(@rdivide,gamma_marg_dep,mean(gamma_marg_dep,2));

%% align preferred phases
for cc = 1:24
   [~,mind] = max(gamma_marg_dep(cc,:));
   align_gamma_stas(cc,:,:) = circshift(gamma_stas(cc,:,:),[0 -mind+round(n_phasebins/2) 0]);
   align_gamma_g(cc,:,:) = circshift(gamma_stas_g(cc,:,:),[0 -mind+round(n_phasebins/2) 0]);
   [~,mind] = max(alpha_marg_dep(cc,:));
   align_alpha_stas(cc,:,:) = circshift(alpha_stas(cc,:,:),[0 -mind+round(n_phasebins/2) 0]);
   align_alpha_g(cc,:,:) = circshift(alpha_stas_g(cc,:,:),[0 -mind+round(n_phasebins/2) 0]);
   [~,mind] = max(delta_marg_dep(cc,:));
   align_delta_stas(cc,:,:) = circshift(delta_stas(cc,:,:),[0 -mind+round(n_phasebins/2) 0]);
   align_delta_g(cc,:,:) = circshift(delta_stas_g(cc,:,:),[0 -mind+round(n_phasebins/2) 0]);
end
%%
sig = 0.3;
[XX,YY] = meshgrid(-5:5);
DD = XX.^2 + YY.^2;
gfun = exp(-DD/(2*sig^2));
gfun = gfun/sum(gfun(:));

gcolor = [0.2 0.8 0.2];

close all
for cc = 11:14
    cc
    cur_sta = squeeze(align_delta_stas(cc,:,:))/new_dt;
    cur_g = squeeze(align_delta_g(cc,:,:))/new_dt;
%     cur_sta = squeeze(align_alpha_stas(cc,:,:))/new_dt;
%     cur_g = squeeze(align_alpha_g(cc,:,:))/new_dt;
    subplot(3,3,1)
%     imagesc(cur_sta);
    imagesc(bar_axis,phase_bin_cents*180/pi,convn(cur_sta,gfun,'same'));
    xlabel('Bar pos (deg)');ylabel('Phase (deg)')
    title('Delta')
    subplot(3,3,2)
    plot(bar_axis,cur_sta(1,:));hold on
    plot(bar_axis,cur_sta(round(n_phasebins/2),:),'r');
    b = regress(cur_sta(round(n_phasebins/2),:)',cur_sta(1,:)');
    plot(bar_axis,cur_sta(1,:)*b,'k--');
    off = mean(cur_sta(round(n_phasebins/2),:) - cur_sta(1,:));
    plot(bar_axis,cur_sta(1,:)+off,'color',gcolor,'linestyle','--');
    xlabel('Bar pos (deg)');ylabel('Firing rate (Hz)')
    subplot(3,3,3)
    g1 = jmm_smooth_1d_cor(cur_g(1,:),0.75);
    gp = jmm_smooth_1d_cor(cur_g(round(n_phasebins/2),:),0.75);
    b = regress(gp',g1');
    off = mean(gp - g1);
    plot(g1);hold on
    plot(gp,'r');
    plot(g1*b,'k--');
    plot(g1+off,'color',gcolor,'linestyle','--');
    xlabel('Gen fun')
    ylabel('Firing rate (Hz)')
    
    cur_sta = squeeze(align_alpha_stas(cc,:,:))/new_dt;
    cur_g = squeeze(align_alpha_g(cc,:,:))/new_dt;
subplot(3,3,4)
%     imagesc(cur_sta);
    imagesc(bar_axis,phase_bin_cents*180/pi,convn(cur_sta,gfun,'same'));
    xlabel('Bar pos (deg)');ylabel('Phase (deg)')
    title('Alpha')
    subplot(3,3,5)
    plot(bar_axis,cur_sta(1,:));hold on
    plot(bar_axis,cur_sta(round(n_phasebins/2),:),'r');
    b = regress(cur_sta(round(n_phasebins/2),:)',cur_sta(1,:)');
    plot(bar_axis,cur_sta(1,:)*b,'k--');
    off = mean(cur_sta(round(n_phasebins/2),:) - cur_sta(1,:));
    plot(bar_axis,cur_sta(1,:)+off,'color',gcolor,'linestyle','--');
    xlabel('Bar pos (deg)');ylabel('Firing rate (Hz)')
    subplot(3,3,6)
    g1 = jmm_smooth_1d_cor(cur_g(1,:),0.75);
    gp = jmm_smooth_1d_cor(cur_g(round(n_phasebins/2),:),0.75);
    b = regress(gp',g1');
    off = mean(gp - g1);
    plot(g1);hold on
    plot(gp,'r');
    plot(g1*b,'k--');
    plot(g1+off,'color',gcolor,'linestyle','--');
    xlabel('Gen fun')
    ylabel('Firing rate (Hz)')

    cur_sta = squeeze(align_gamma_stas(cc,:,:))/new_dt;
    cur_g = squeeze(align_gamma_g(cc,:,:))/new_dt;
subplot(3,3,7)
%     imagesc(cur_sta);
    imagesc(bar_axis,phase_bin_cents*180/pi,convn(cur_sta,gfun,'same'));
    xlabel('Bar pos (deg)');ylabel('Phase (deg)')
    title('Gamma')
    subplot(3,3,8)
    plot(bar_axis,cur_sta(1,:));hold on
    plot(bar_axis,cur_sta(round(n_phasebins/2),:),'r');
    b = regress(cur_sta(round(n_phasebins/2),:)',cur_sta(1,:)');
    plot(bar_axis,cur_sta(1,:)*b,'k--');
    off = mean(cur_sta(round(n_phasebins/2),:) - cur_sta(1,:));
    plot(bar_axis,cur_sta(1,:)+off,'color',gcolor,'linestyle','--');
    xlabel('Bar pos (deg)');ylabel('Firing rate (Hz)')
    subplot(3,3,9)
    g1 = jmm_smooth_1d_cor(cur_g(1,:),0.75);
    gp = jmm_smooth_1d_cor(cur_g(round(n_phasebins/2),:),0.75);
    b = regress(gp',g1');
    off = mean(gp - g1);
    plot(g1);hold on
    plot(gp,'r');
    plot(g1*b,'k--');
    plot(g1+off,'color',gcolor,'linestyle','--');
    xlabel('Gen fun')
    ylabel('Firing rate (Hz)')


    pause
    clf
end

