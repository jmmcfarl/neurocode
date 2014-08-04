clear all
% close all
cd ~/Data/bruce/G082/
load jbeG082Expts.mat

load ~/Data/bruce/G082/stims/expt_data.mat
%%
stim_fs = 100; %in Hz
dt = 0.01;
flen = 15;
nPix = 25;
stim_params = NIMcreate_stim_params([flen 2*nPix],dt);

%%
use_expts = find(expt_binoc == 0 & expt_npix == 25);
orRC = false(length(use_expts),1);
for ii = 1:length(use_expts)
   if strcmp(Expts{use_expts(ii)}.Header.expname,'rls.orRC');
       orRC(ii) = true;
   end
end
use_expts(~orRC) = [];
%%
cd ~/Data/bruce/G082/

all_stim_times = [];
all_used_inds = [];
all_bar_mat = [];
all_binned_spks = [];
all_exptvec = [];
all_trialvec = [];
for ee = 1:length(use_expts);
    fprintf('Expt %d of %d\n',ee,length(use_expts));
    cur_expt = use_expts(ee);
    fname = sprintf('Expt%dClusterTimes.mat',cur_expt);
    load(fname);
    
    fname = sprintf('~/Data/bruce/G082/stims/Expt%d_stim',cur_expt);
    load(fname);
    
    trial_start_times = [Expts{cur_expt}.Trials(:).TrialStart]/1e4;
    trial_end_times = [Expts{cur_expt}.Trials(:).TrueEnd]/1e4;
    trial_durs = [Expts{cur_expt}.Trials(:).dur]/1e4;
    
    use_trials = find(trial_durs >= 0.5);
    n_trials = length(use_trials);
    for tt = 1:n_trials
        cur_stim_times = Expts{cur_expt}.Trials(use_trials(tt)).Start/1e4;
        cur_or = Expts{cur_expt}.Trials(use_trials(tt)).or;
        zero_or = find(cur_or==0);
        ninety_or = find(cur_or == 90);
        
        cur_stim_mat = zeros(length(cur_stim_times),nPix*2);
        cur_stim_mat(zero_or,1:nPix) = left_stim_mats{use_trials(tt)}(zero_or,:);
        cur_stim_mat(ninety_or,nPix+1:end) = left_stim_mats{use_trials(tt)}(ninety_or,:);
                
        cur_t_edges = [cur_stim_times; Expts{cur_expt}.Trials(use_trials(tt)).End(end)/1e4];
        cur_binned_spks = nan(length(cur_stim_times),96);
        for cc = 1:96
            cur_hist = histc(Clusters{cc}.times,cur_t_edges);
            cur_binned_spks(:,cc) = cur_hist(1:end-1);
%             temp = convert_to_spikebins(cur_hist(1:end-1));
        end
        
        bar_Xmat = create_time_embedding(cur_stim_mat,stim_params);
        
        cur_used_inds = ones(length(cur_stim_times),1);
        cur_used_inds(1:flen) = 0;
        
        all_stim_times = [all_stim_times; cur_stim_times];
        all_used_inds = [all_used_inds; cur_used_inds];
        all_bar_mat = [all_bar_mat; bar_Xmat];
        all_binned_spks = [all_binned_spks; cur_binned_spks];
        all_exptvec = [all_exptvec; ones(size(cur_stim_times))*ee];
        all_trialvec = [all_trialvec; ones(size(cur_stim_times))*tt];
    end
    
end

%%
[c,ia,ic] = unique([all_exptvec all_trialvec],'rows');
n_trials = length(ia);

xv_frac = 0.2;
n_xv_trials = round(n_trials*xv_frac);
xv_set = randperm(n_trials);
xv_set(n_xv_trials+1:end) = [];
xv_inds = find(ismember(ic,xv_set));
tr_inds = find(~ismember(ic,xv_set))';

tr_inds(all_used_inds(tr_inds) == 0) = [];
xv_inds(all_used_inds(xv_inds) == 0) = [];

Xexpt = zeros(length(all_stim_times),length(use_expts)-1);
for i = 1:length(use_expts)-1
    cur_set = find(all_exptvec==i);
    Xexpt(cur_set,i) = 1;
end

%%
l_d2XT = 300;
l_L1 = 20;
% ql_L1 = 0.002;
% ql_d2XT = 50;
ql_L1 = 0.001;
ql_d2XT = 30;

reg_params = NIMcreate_reg_params('lambda_d2XT',l_d2XT,'lambda_L1',l_L1,'temporal_boundaries','zero');
stim_params = NIMcreate_stim_params([flen 2*nPix],dt,1,1,length(use_expts)-1);
stim_params1 = NIMcreate_stim_params([flen nPix],dt,1,1,length(use_expts)-1);
null_params = NIMcreate_stim_params([1 size(Xexpt,2)],dt);
silent = 0;

for cc = [19]
    fprintf('Fitting Lin model cell %d of %d\n',cc,96);
    fprintf('Mean spike rate %.3f\n',mean(all_binned_spks(tr_inds,cc))/dt);
    Robs = all_binned_spks(tr_inds,cc);
    
    mod_signs = [1 1 1]; %determines whether input is exc or sup (doesn't matter in the linear case)
    NL_types = {'lin','quad','quad'}; %define subunit as linear
    
    fit0(cc) = NIMinitialize_model(stim_params,mod_signs,NL_types,reg_params); %initialize NIM
    
    fit0(cc) = NIMadjust_regularization(fit0(cc),[2:length(mod_signs)],'lambda_d2XT',ql_d2XT);
    fit0(cc) = NIMadjust_regularization(fit0(cc),[2:length(mod_signs)],'lambda_L1',ql_L1);
    temp = fit0(cc); temp.stim_params = stim_params1;
    L2_mats = create_L2_matrices(temp); %single-ori L2_mats
    full_ndims = prod(stim_params.stim_dims); half_ndims = prod(stim_params1.stim_dims);
    full_L2_mats.L2_d2XT = speye(full_ndims);
    full_L2_mats.L2_d2XT(1:half_ndims,1:half_ndims) = L2_mats.L2_d2XT;
    full_L2_mats.L2_d2XT((half_ndims+1):end,(half_ndims+1):end) = L2_mats.L2_d2XT;
    
    fit0(cc) = NIMfit_filters(fit0(cc),Robs,all_bar_mat(tr_inds,:),Xexpt(tr_inds,:),[],silent,[],full_L2_mats); %fit stimulus filters

    null_mod(cc) = NIMinitialize_model(null_params,1,{'lin'},reg_params); %initialize NIM
    null_mod(cc) = NIMfit_filters(null_mod(cc),Robs,Xexpt(tr_inds,:)); %fit stimulus filters

    if xv_frac > 0
        Robsxv = all_binned_spks(xv_inds,cc);
        xvLL(cc) = NIMmodel_eval(fit0(cc),Robsxv,all_bar_mat(xv_inds,:),Xexpt(xv_inds,:));
        null_xvLL(cc) = NIMmodel_eval(null_mod(cc),Robsxv,Xexpt(xv_inds,:));
        xv_imp(cc) = (xvLL(cc)-null_xvLL(cc))/log(2);
        fprintf('LL imp: %.3f\n',xv_imp(cc));
    end
    
end

%% 
cd /Users/James/James_scripts/bruce/G082
save dual_ori_models *xvLL xv_imp fit0 null_mod stim_params reg_params 