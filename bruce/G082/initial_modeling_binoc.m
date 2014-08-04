clear all
% close all

expt_num = 82;
%%
cd(sprintf('~/Data/bruce/G0%d',expt_num));
load(sprintf('jbeG0%dExpts.mat',expt_num));
load(sprintf('~/Data/bruce/G0%d/stims/expt_data.mat',expt_num));
%%
stim_fs = 100; %in Hz
dt = 0.01;
flen = 15;
if expt_num==82
    nPix = 25;
elseif expt_num == 85
    nPix = 36;
else
    error('No preset nPix');
end
stim_params = NIMcreate_stim_params([flen 2*nPix],dt);
stim_params1 = NIMcreate_stim_params([flen nPix],dt);

use_40pix = false(40,1);
use_40pix(8:32) = true;
%%
% use_expts = [45 46 48];
use_expts = find(expt_binoc == 1);
orRC = false(length(use_expts),1);
for ii = 1:length(use_expts)
   if strcmp(Expts{use_expts(ii)}.Header.expname,'rls.orRC');
       orRC(ii) = true;
   end
end
use_expts(~orRC) = [];
%%
cd(sprintf('~/Data/bruce/G0%d/',expt_num));

all_stim_times = [];
all_used_inds = [];
all_bar_mat = [];
all_Lbar_mat = [];
all_Rbar_mat = [];
all_binned_spks = [];
all_exptvec = [];
all_trialvec = [];
for ee = 1:length(use_expts);
    fprintf('Expt %d of %d\n',ee,length(use_expts));
    cur_expt = use_expts(ee);
    fname = sprintf('Expt%dClusterTimes.mat',cur_expt);
    load(fname);
    
    fname = sprintf('~/Data/bruce/G0%d/stims/Expt%d_stim',expt_num,cur_expt);
    load(fname);
    
    if expt_npix(cur_expt) == 40
        cur_use_pix = use_40pix;
    else
        cur_use_pix = true(nPix,1);
    end
    
    trial_start_times = [Expts{cur_expt}.Trials(:).TrialStart]/1e4;
    trial_end_times = [Expts{cur_expt}.Trials(:).TrueEnd]/1e4;
    trial_durs = [Expts{cur_expt}.Trials(:).dur]/1e4;
    
    use_trials = find(trial_durs >= 0.5);
    
    %% CHECK FOR STIMULUS MISMATCHES
    n_stim = nan(length(use_trials),1);
    n_frm = nan(length(use_trials),1);
    for tt = 1:length(use_trials)
        n_stim(tt) = length(Expts{cur_expt}.Trials(use_trials(tt)).Start);
        n_frm(tt) = size(left_stim_mats{use_trials(tt)},1);
    end
    misaligned = find(n_frm ~= n_stim);
    fprintf('Eliminating %d of %d misaligned trials\n',length(misaligned),length(use_trials));
    use_trials(misaligned) = [];
    %%
    
    n_trials = length(use_trials);
    for tt = 1:n_trials
        if mod(tt,5) == 0
        fprintf('Loading trial %d of %d\n',tt,n_trials);
        end
        cur_stim_times = Expts{cur_expt}.Trials(use_trials(tt)).Start/1e4;
        
        cur_stim_mat = zeros(length(cur_stim_times),nPix*2);
        cur_stim_mat(:,1:nPix) = left_stim_mats{use_trials(tt)}(:,cur_use_pix);
        cur_stim_mat(:,nPix+1:end) = right_stim_mats{use_trials(tt)}(:,cur_use_pix);
        cur_lstim_mat = left_stim_mats{use_trials(tt)}(:,cur_use_pix);
        cur_rstim_mat = right_stim_mats{use_trials(tt)}(:,cur_use_pix);
                
        cur_t_edges = [cur_stim_times; Expts{cur_expt}.Trials(use_trials(tt)).End(end)/1e4];
        cur_binned_spks = nan(length(cur_stim_times),96);
        for cc = 1:96
            cur_hist = histc(Clusters{cc}.times,cur_t_edges);
            cur_binned_spks(:,cc) = cur_hist(1:end-1);
        end
        
        bar_Xmat = create_time_embedding(cur_stim_mat,stim_params);
        bar_LXmat = create_time_embedding(cur_lstim_mat,stim_params1);
        bar_RXmat = create_time_embedding(cur_rstim_mat,stim_params1);
        
        cur_used_inds = ones(length(cur_stim_times),1);
        cur_used_inds(1:flen) = 0;
        
        all_stim_times = [all_stim_times; cur_stim_times];
        all_used_inds = [all_used_inds; cur_used_inds];
        all_bar_mat = [all_bar_mat; bar_Xmat];
        all_Lbar_mat = [all_Lbar_mat; bar_LXmat];
        all_Rbar_mat = [all_Rbar_mat; bar_RXmat];
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

% used_expt_set = 1:7;
% % used_expt_set(3:4) = [];
% % used_expt_set([1 2 5 6 7]) = [];
% tr_inds(~ismember(all_exptvec(tr_inds),used_expt_set)) = [];
% xv_inds(~ismember(all_exptvec(xv_inds),used_expt_set)) = [];


%%
reg_params = NIMcreate_reg_params('lambda_d2XT',300,'lambda_L1',20,'temporal_boundaries','zero');
stim_params = NIMcreate_stim_params([flen 2*nPix],dt,1,1,length(use_expts)-1);
stim_params1 = NIMcreate_stim_params([flen nPix],dt,1,1,length(use_expts)-1);
null_params = NIMcreate_stim_params([1 size(Xexpt,2)],dt);
% stim_params = NIMcreate_stim_params([flen 2*nPix],dt);
% stim_params1 = NIMcreate_stim_params([flen nPix],dt);
silent = 1;

all_inds = sort([tr_inds'; xv_inds]);
for cc = [1:96]
    fprintf('Fitting Lin model cell %d of %d\n',cc,96);
    fprintf('Mean spike rate %.3f\n',mean(all_binned_spks(tr_inds,cc))/dt);
    Robs = all_binned_spks(tr_inds,cc);
    
%     mod_signs = [1 1]; %determines whether input is exc or sup (doesn't matter in the linear case)
%     NL_types = {'lin','quad'}; %define subunit as linear
   mod_signs = [1]; %determines whether input is exc or sup (doesn't matter in the linear case)
    NL_types = {'lin'}; %define subunit as linear
    
    fit0(cc) = NIMinitialize_model(stim_params,mod_signs,NL_types,reg_params); %initialize NIM
%     fit0(cc) = NIMadjust_regularization(fit0(cc),2,'lambda_d2XT',50);
%     fit0(cc) = NIMadjust_regularization(fit0(cc),2,'lambda_L1',0.001);
    fit0(cc) = NIMfit_filters(fit0(cc),Robs,all_bar_mat(tr_inds,:),Xexpt(tr_inds,:),[],silent); %fit stimulus filters
%     fit0(cc) = NIMfit_filters(fit0(cc),Robs,all_bar_mat(tr_inds,:),[],[],silent); %fit stimulus filters

    fitL0(cc) = NIMinitialize_model(stim_params1,mod_signs,NL_types,reg_params); %initialize NIM
%     fitL0(cc) = NIMadjust_regularization(fitL0(cc),2,'lambda_d2XT',10);
%     fitL0(cc) = NIMadjust_regularization(fitL0(cc),2,'lambda_L1',0.001);
    fitL0(cc) = NIMfit_filters(fitL0(cc),Robs,all_Lbar_mat(tr_inds,:),Xexpt(tr_inds,:),[],silent); %fit stimulus filters
%     fitL0(cc) = NIMfit_filters(fitL0(cc),Robs,all_Lbar_mat(tr_inds,:),[],[],silent); %fit stimulus filters

%     fitR0(cc) = NIMinitialize_model(stim_params1,mod_signs,NL_types,reg_params); %initialize NIM
%     fitR0(cc) = NIMadjust_regularization(fitR0(cc),2,'lambda_d2XT',50);
%     fitR0(cc) = NIMadjust_regularization(fitR0(cc),2,'lambda_L1',0.001);
    fitR0(cc) = fitL0(cc); %initialize NIM
    fitR0(cc) = NIMfit_filters(fitR0(cc),Robs,all_Rbar_mat(tr_inds,:),Xexpt(tr_inds,:),[],silent); %fit stimulus filters
%     fitR0(cc) = NIMfit_filters(fitR0(cc),Robs,all_Rbar_mat(tr_inds,:),[],[],silent); %fit stimulus filters

    null_mod(cc) = NIMinitialize_model(null_params,1,{'lin'},[]); %initialize NIM
    null_mod(cc) = NIMfit_filters(null_mod(cc),Robs,Xexpt(tr_inds,:)); %fit stimulus filters
% avg_rate = mean(Robs);

    if xv_frac > 0
        Robsxv = all_binned_spks(xv_inds,cc);
        xvLL(cc) = NIMmodel_eval(fit0(cc),Robsxv,all_bar_mat(xv_inds,:),Xexpt(xv_inds,:));
        LxvLL(cc) = NIMmodel_eval(fitL0(cc),Robsxv,all_Lbar_mat(xv_inds,:),Xexpt(xv_inds,:));
        RxvLL(cc) = NIMmodel_eval(fitR0(cc),Robsxv,all_Rbar_mat(xv_inds,:),Xexpt(xv_inds,:));
%         xvLL(cc) = NIMmodel_eval(fit0(cc),Robsxv,all_bar_mat(xv_inds,:));
%         LxvLL(cc) = NIMmodel_eval(fitL0(cc),Robsxv,all_Lbar_mat(xv_inds,:));
%         RxvLL(cc) = NIMmodel_eval(fitR0(cc),Robsxv,all_Rbar_mat(xv_inds,:));

        null_xvLL(cc) = NIMmodel_eval(null_mod(cc),Robsxv,Xexpt(xv_inds,:));
% null_rate = ones(size(Robsxv))*avg_rate;
% null_xvLL(cc) = sum(Robsxv.*log(null_rate) - null_rate)/sum(Robsxv);

xv_imp(cc) = (xvLL(cc)-null_xvLL(cc))/log(2);
        Lxv_imp(cc) = (LxvLL(cc)-null_xvLL(cc))/log(2);
        Rxv_imp(cc) = (RxvLL(cc)-null_xvLL(cc))/log(2);
        fprintf('LL imp: %.3f\n',xv_imp(cc));
        fprintf('Left LL imp: %.3f   Right LL imp: %.3f\n',Lxv_imp(cc),Rxv_imp(cc));
    end
    
    for jj = 1:length(use_expts)
       cur_set = all_inds(all_exptvec(all_inds) == jj);
       cur_Robs = all_binned_spks(cur_set,cc);
        spLL(cc,jj) = NIMmodel_eval(fit0(cc),cur_Robs,all_bar_mat(cur_set,:),Xexpt(cur_set,:));
        LspLL(cc,jj) = NIMmodel_eval(fitL0(cc),cur_Robs,all_Lbar_mat(cur_set,:),Xexpt(cur_set,:));
        RspLL(cc,jj) = NIMmodel_eval(fitR0(cc),cur_Robs,all_Rbar_mat(cur_set,:),Xexpt(cur_set,:));
        null_spLL(cc,jj) = NIMmodel_eval(null_mod(cc),Robsxv,Xexpt(xv_inds,:));
    end
    
end

%%

