clear all
% close all

ExptNum = 235;
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
dt = 0.01;
new_dt = 0.002;
usfrac = dt/new_dt;
tbspace = 3;
nLags = 15*usfrac/tbspace;
flen = nLags*tbspace;


n_bar_pos = bar_params.n_bars;
stim_params = NIMcreate_stim_params([nLags n_bar_pos],new_dt,1,tbspace,length(bar_expts)-1);
bar_axis = bar_params.bar_axis;
% bar_axis = [bar_params.bar_axis (bar_params.bar_axis(end) + 0.125)]; %add extra bin edge for histing

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
full_taxis = [];
full_trialvec = [];
full_bar_pos = [];
full_t_since_tstart = [];
full_used_inds = [];
for ee = 1:length(bar_expts);
    fprintf('Processing expt %d of %d\n',ee,length(bar_expts));
    fname = sprintf('Expt%dClusterTimes.mat',bar_expts(ee));
    load(fname);
    
    trial_durs{ee} = [Expts{bar_expts(ee)}.Trials(:).dur];
    n_trials(ee) = length(Expts{bar_expts(ee)}.Trials);
    all_t_axis = [];
    all_used_inds = [];
    all_old_t_inds = [];
    all_expt_vec = [];
    all_trial_vec = [];
    all_bar_Op = [];
    all_bar_Xmat = [];
    all_binned_spikes = [];
    all_used_inds = [];
    
    all_t_since_tstart = [];
    for tt = 1:n_trials(ee)
        
        cur_bar_Op = [Expts{bar_expts(ee)}.Trials(tt).Op];
        
        %up-sample bar positions
        cur_bar_Op = repmat(cur_bar_Op,1,usfrac)';
        cur_bar_Op = cur_bar_Op(:);
                
        cur_t_edges = Expts{bar_expts(ee)}.Trials(tt).Start/1e4:new_dt:Expts{bar_expts(ee)}.Trials(tt).End(end)/1e4;
        cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        
        cur_bar_Op(length(cur_t_axis)+1:end) = [];
               
        cur_binned_spks = nan(length(cur_t_axis),length(use_sus));
        for cc = 1:length(use_sus)
            cur_hist = histc(Clusters{use_sus(cc)}.times,cur_t_edges);
            cur_binned_spks(:,cc) = cur_hist(1:end-1);
        end
                
        cur_bar_mat = zeros(length(cur_t_axis),n_bar_pos);
        for b = 1:n_bar_pos
%             cur_set = find(cur_bar_Op >= bar_axis(b) & cur_bar_Op < bar_axis(b+1));
            cur_set = find(cur_bar_Op == bar_axis(b));
            cur_bar_mat(cur_set,b) = 1;
        end
        bar_Xmat = create_time_embedding(cur_bar_mat,stim_params);
        
        cur_used_inds = ones(length(cur_t_axis),1);
        cur_used_inds(1:flen) = 0;
        
        cur_t_since_tstart = cur_t_axis - cur_t_axis(1);
                
        all_used_inds = [all_used_inds; cur_used_inds(:)];
        all_t_axis = [all_t_axis; cur_t_axis(:)];
        all_t_since_tstart = [all_t_since_tstart; cur_t_since_tstart(:)];
        all_binned_spikes = [all_binned_spikes; cur_binned_spks];
        all_expt_vec = [all_expt_vec; ones(length(cur_t_axis),1)*bar_expts(ee)];
        all_trial_vec = [all_trial_vec; ones(length(cur_t_axis),1)*tt];
        all_bar_Op = [all_bar_Op; cur_bar_Op(:)];
        all_bar_Xmat = [all_bar_Xmat; bar_Xmat];
    end
    
%     % for eliminating blinks
%     cur_blink_start_times = all_eye_ts{ee}(all_blink_startinds{ee});
%     cur_blink_stop_times = all_eye_ts{ee}(all_blink_stopinds{ee});
%     cur_blink_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),cur_blink_start_times));
%     cur_blink_stop_inds = round(interp1(all_t_axis,1:length(all_t_axis),cur_blink_stop_times));
%     cur_poss_blinks = find(~isnan(cur_blink_start_inds) & ~isnan(cur_blink_stop_inds));
%     
%     all_blink_inds = zeros(size(all_t_axis));
%     for i = 1:length(cur_poss_blinks)
%         all_blink_inds(cur_blink_start_inds(cur_poss_blinks(i)):cur_blink_stop_inds(cur_poss_blinks(i))) = 1;
%     end
%     
%     cur_blink_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),cur_blink_start_times));
%     cur_blink_stop_inds = round(interp1(all_t_axis,1:length(all_t_axis),cur_blink_stop_times));
%     cur_poss_blinks = find(~isnan(cur_blink_start_inds) & ~isnan(cur_blink_stop_inds));
%     all_blink_inds = zeros(size(all_t_axis));
%     for i = 1:length(cur_poss_blinks)
%         all_blink_inds(cur_blink_start_inds(cur_poss_blinks(i)):cur_blink_stop_inds(cur_poss_blinks(i))) = 1;
%     end
%     all_used_inds(all_blink_inds == 1) = 0;
    
    % 
    full_Xmat = [full_Xmat; all_bar_Xmat];
    full_spkbinned = [full_spkbinned; all_binned_spikes];
    full_exptvec = [full_exptvec; ones(length(all_used_inds),1)*ee];
    full_taxis = [full_taxis; all_t_axis];
    full_t_since_tstart = [full_t_since_tstart; all_t_since_tstart];
    full_trialvec = [full_trialvec; all_trial_vec];
    full_bar_pos = [full_bar_pos; all_bar_Op];
    full_used_inds = [full_used_inds; all_used_inds];
    
end
full_used_inds = logical(full_used_inds);

%% Parse into training and XV sets
[c,ia,ic] = unique([full_exptvec full_trialvec],'rows');
n_trials = length(ia);

rp = randperm(n_trials);
rand_trial_vec = rp(ic);
[~,ind_shuff] = sort(rand_trial_vec);

xv_frac = 0.2;
n_xv_trials = round(n_trials*xv_frac);
xv_tset = randperm(n_trials);
xv_tset(n_xv_trials+1:end) = [];
tr_tset = find(~ismember(1:n_trials,xv_tset));
xv_inds = find(ismember(ic,xv_tset));
tr_inds = find(ismember(ic,tr_tset))';

xv_inds(full_used_inds(xv_inds)==0) = [];
tr_inds(full_used_inds(tr_inds) == 0) = [];

%% Make indicator matrix for experiment number
Xexpt = zeros(length(full_taxis),length(bar_expts)-1);
for i = 1:length(bar_expts)-1
    cur_set = find(full_exptvec==i);
    Xexpt(cur_set,i) = 1;
end

%% Fit STRFs
reg_params = NIMcreate_reg_params('lambda_d2XT',10,'XTmix',[1 5]);
null_params = NIMcreate_stim_params([1 size(Xexpt,2)],dt);

silent = 1;
for cc = 1:24
% for cc = 2
    fprintf('Fitting Lin model cell %d of %d\n',cc,24);
    Robs = full_spkbinned(tr_inds,cc);
    
    mod_signs = [1]; %determines whether input is exc or sup (doesn't matter in the linear case)
    NL_types = {'lin'}; %define subunit as linear
    
    fit0(cc) = NIMinitialize_model(stim_params,mod_signs,NL_types,reg_params); %initialize NIM
    fit0(cc) = NIMfit_filters(fit0(cc),Robs,full_Xmat(tr_inds,:),Xexpt(tr_inds,:),[],silent); %fit stimulus filters

%     mod_signs = [1 -1]; %determines whether input is exc or sup (doesn't matter in the linear case)
%     NL_types = {'threshlin','threshlin'}; %define subunit as linear
%     fit1(cc) = NIMinitialize_model(stim_params,mod_signs,NL_types,reg_params); %initialize NIM
%     fit1(cc) = NIMfit_filters(fit1(cc),Robs,full_Xmat(tr_inds,:),Xexpt(tr_inds,:),[],silent); %fit stimulus filters

    null_mod(cc) = NIMinitialize_model(null_params,1,{'lin'},reg_params); %initialize NIM
    null_mod(cc) = NIMfit_filters(null_mod(cc),Robs,Xexpt(tr_inds,:)); %fit stimulus filters

    if xv_frac > 0
        Robsxv = full_spkbinned(xv_inds,cc);
        null_xvLL(cc) = NIMmodel_eval(null_mod(cc),Robsxv,Xexpt(xv_inds,:));

        fit0_xvLL(cc) = NIMmodel_eval(fit0(cc),Robsxv,full_Xmat(xv_inds,:),Xexpt(xv_inds,:));
        fit0_xvimp(cc) = (fit0_xvLL(cc)-null_xvLL(cc))/log(2);
%         fit1_xvLL(cc) = NIMmodel_eval(fit1(cc),Robsxv,full_Xmat(xv_inds,:),Xexpt(xv_inds,:));
%         fit1_xvimp(cc) = (fit1_xvLL(cc)-null_xvLL(cc))/log(2);
% 

%         fprintf('fit0 XVimp: %.3f\n fit1 XVimp: %.3f\n',fit0_xvimp(cc),fit1_xvimp(cc));
        fprintf('fit0 XVimp: %.3f\n ',fit0_xvimp(cc));
    end
    
end

%%
save strf_fits stim_params *xvLL *xvimp null_mod fit0 reg_params

%%
for cc = 1:-6
NIMdisplay_model(fit0(cc));
cc
pause
close all
end

%%
clear temp_kern 
for cc = 1:24
    k_mat = reshape(fit0(cc).mods(1).filtK,stim_params.stim_dims(1),stim_params.stim_dims(2));
    temp_kern(cc,:) = sqrt(mean(k_mat.^2,2));
    [~,best_lag(cc)] = max(temp_kern(cc,:));
end
pref_lag = stim_params.stim_dims(1) - best_lag - 1;
[bps,fls] = meshgrid(1:stim_params.stim_dims(2),1:stim_params.stim_dims(1));

subplot(2,1,1)
imagesc((1:stim_params.stim_dims(1))*stim_params.dt*stim_params.tent_spacing,1:24,temp_kern)
xlabel('Time lag (s)','fontsize',16)
ylabel('Channel','fontsize',16)
title(sprintf('Temporal Kernels M%d',ExptNum));
subplot(2,1,2)
imagesc((1:stim_params.stim_dims(1))*stim_params.dt*stim_params.tent_spacing,1:24,bsxfun(@rdivide,temp_kern,max(temp_kern,[],2)))
xlabel('Time lag (s)','fontsize',16)
ylabel('Channel','fontsize',16)
title(sprintf('Temporal Kernels M%d',ExptNum));