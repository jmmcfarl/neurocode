clear all
% close all

Expt_name = 'G087';
data_dir = ['~/Data/bruce/' Expt_name];

cd(data_dir);

load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

anal_dir = ['~/Analysis/bruce/' Expt_name];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end

anal_name = '_rph_gratings';
%%
stim_fs = 100; %in Hz
dt = 0.01;
flen = 15;
nPix = 6; %number of phases
Fr = 3;
stim_params = NIMcreate_stim_params([flen nPix],dt,Fr,1);
%%
include_expts = {'square.orXnph'};
for ii = 1:length(Expts)
    expt_names{ii} = Expts{ii}.Header.expname;
    expt_dds(ii) = Expts{ii}.Stimvals.dd;
    expt_bar_ori(ii) = Expts{ii}.Stimvals.or;
    expt_sac_dir(ii) = mod(Expts{ii}.Stimvals.Fa,180);
    expt_Fr(ii) = Expts{ii}.Stimvals.Fr;
    included_type(ii) = strcmp(expt_names{ii},include_expts);
end

use_expts = find(included_type);
%%
cd(data_dir);

trial_cnter = 0;

all_stim_times = [];
all_used_inds = false(0);
all_phase_Xmat = [];
all_binned_spks = [];
all_exptvec = [];
all_trialvec = [];
all_trial_ors = [];
for ee = 1:length(use_expts);
    fprintf('Expt %d of %d\n',ee,length(use_expts));
    cur_expt = use_expts(ee);
    fname = sprintf('Expt%dClusterTimes.mat',cur_expt);
    load(fname);
    
    online_fname = ['jbe' Expt_name '.online'];
    [mtr_seq,mtr_is_match] = get_mtr_seq(online_fname,Expts,cur_expt);
         
    trial_start_times = [Expts{cur_expt}.Trials(:).TrialStart]/1e4;
    trial_end_times = [Expts{cur_expt}.Trials(:).TrueEnd]/1e4;
    trial_durs = [Expts{cur_expt}.Trials(:).dur]/1e4;
    trial_ids = [Expts{cur_expt}.Trials(:).id];
 
    [un_ids,id_inds] = unique(trial_ids);       
    trial_durs = trial_durs(id_inds);
    trial_start_times = trial_start_times(id_inds);
    trial_end_times = trial_end_times(id_inds);
        
    trial_nph = [Expts{cur_expt}.Trials(id_inds).nph]; %when trial_nph == 6 it's random phase
    trial_ors = [Expts{cur_expt}.Trials(id_inds).or];
    use_trials = find(trial_durs >= 0.5 & trial_nph == 6); 
    mtr_seq = mtr_seq(use_trials);
    mtr_is_match = mtr_is_match(use_trials);
    if ~all(mtr_is_match==1)
        fprintf('WARNING, MISSING TRIAL MATCHES');
    end
    
    all_trial_ors = [all_trial_ors; trial_ors(use_trials)'];
    n_trials = length(use_trials);
    for tt = 1:n_trials
        cur_stim_times = Expts{cur_expt}.Trials(use_trials(tt)).Start/1e4;
        cur_phase_vals = mtr_seq{tt}(2:3:end);
        n_frames = length(cur_phase_vals)*Fr;
        if n_frames < (trial_durs(use_trials(tt))/(dt*Fr)+1)
            fprintf('Too few frames\n');
        else
            cur_t_edges = (trial_start_times(use_trials(tt)):dt:trial_end_times(use_trials(tt)))';
            trial_tlen(tt+trial_cnter) = length(cur_t_edges);
            cur_binned_spks = nan(length(cur_t_edges),96);
            for cc = 1:96
                cur_hist = histc(Clusters{cc}.times,cur_t_edges);
                cur_binned_spks(:,cc) = cur_hist;
                %             temp = convert_to_spikebins(cur_hist(1:end-1));
            end
            
            cur_phase_mat = zeros(length(cur_phase_vals),6);
            for ii = 1:6
                cur_set = find(cur_phase_vals == num2str(ii-1));
                cur_phase_mat(cur_set,ii) = 1;
            end
            phase_Xmat = create_time_embedding(cur_phase_mat,stim_params);
            trial_plen(tt+trial_cnter) = size(phase_Xmat,1);
            phase_Xmat(length(cur_t_edges)+1:end,:) = [];
            
            cur_used_inds = true(length(cur_t_edges),1);
            cur_used_inds(1:flen) = false;
            
            all_stim_times = [all_stim_times; cur_t_edges];
            all_used_inds = [all_used_inds; cur_used_inds];
            all_phase_Xmat = [all_phase_Xmat; phase_Xmat];
            all_binned_spks = [all_binned_spks; cur_binned_spks];
            all_exptvec = [all_exptvec; ones(size(cur_t_edges))*ee];
            all_trialvec = [all_trialvec; ones(size(cur_t_edges))*(tt + trial_cnter)];
        end
    end   
    trial_cnter = trial_cnter + n_trials;
end


%% CREATE XV TRIAL SET
% [c,ia,ic] = unique([all_exptvec all_trialvec],'rows');
[c,ia,ic] = unique(all_trialvec);
n_trials = length(ia);

xv_frac = 0.2;
n_xv_trials = round(n_trials*xv_frac);
xv_set = randperm(n_trials);
xv_set(n_xv_trials+1:end) = [];
xv_inds = find(ismember(ic,xv_set));
tr_inds = find(~ismember(ic,xv_set))';

tr_inds(~all_used_inds(tr_inds)) = [];
xv_inds(~all_used_inds(xv_inds)) = [];

use_or = 90;
cur_ori_inds = find(all_trial_ors(all_trialvec) == use_or);
tr_inds(~ismember(tr_inds,cur_ori_inds)) = [];
xv_inds(~ismember(xv_inds,cur_ori_inds)) = [];

Xexpt = zeros(length(all_stim_times),length(use_expts)-1);
for i = 1:length(use_expts)-1
    cur_set = find(all_exptvec==i);
    Xexpt(cur_set,i) = 1;
end

%%
l_d2XT = 100;
l_L1 = 0;
XT_mix = [0.02 1];
reg_params = NIMcreate_reg_params('lambda_d2XT',l_d2XT,'lambda_L1',l_L1,'XTmix',XT_mix);
stim_params = NIMcreate_stim_params([flen nPix],dt,Fr,1,length(use_expts)-1);
null_params = NIMcreate_stim_params([1 size(Xexpt,2)],dt);
silent = 1;
optim_params.progTol = 1e-7;
optim_params.optTol = 1e-4;

NL_types = {'quad','quad'}; %define subunit as linear

for cc = 1:96
%     
%     spikebins = convert_to_spikebins(all_binned_spks(tr_inds,cc));
%     spike_cond_stim = all_phase_Xmat(tr_inds(spikebins),:);
%     sta      = mean(spike_cond_stim) - mean(all_phase_Xmat(tr_inds,:));
%     sta = sta/norm(sta);
%     
%     npos = 3; nneg = 3;
%     proj_mat = sta'/(sta*sta')*sta;
%     stim_proj = all_phase_Xmat(tr_inds,:) - all_phase_Xmat(tr_inds,:)*proj_mat;
%     % stim_proj = stim_emb;
%     stvcv = cov(stim_proj(spikebins,:));  utvcv = cov(stim_proj);
%     [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
%     stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
%     
% %     plot(diag(evals),'o')
    
    fprintf('Fitting Lin model cell %d of %d\n',cc,96);
    fprintf('Mean spike rate %.3f\n',mean(all_binned_spks(tr_inds,cc))/dt);
    Robs = all_binned_spks(tr_inds,cc);
    null_mod(cc) = NIMinitialize_model(null_params,1,{'lin'},reg_params); %initialize NIM
    null_mod(cc) = NIMfit_filters(null_mod(cc),Robs,Xexpt(tr_inds,:)); %fit stimulus filters
    
    if xv_frac > 0
        Robsxv = all_binned_spks(xv_inds,cc);
        null_xvLL(cc) = NIMmodel_eval(null_mod(cc),Robsxv,Xexpt(xv_inds,:));
    end
   
    %% FIT LINEAR MODEL
    mod_signs = [1 1]; %determines whether input is exc or sup (doesn't matter in the linear case)
    fit0(cc) = NIMinitialize_model(stim_params,mod_signs,NL_types,reg_params); %initialize NIM
    fit0(cc) = NIMfit_filters(fit0(cc),Robs,all_phase_Xmat(tr_inds,:),Xexpt(tr_inds,:),[],silent,optim_params); %fit stimulus filters
    fit0(cc) = NIMadjust_regularization(fit0(cc) ,[1:length(mod_signs)],'lambda_L1',2);
    fit0(cc) = NIMfit_filters(fit0(cc),Robs,all_phase_Xmat(tr_inds,:),Xexpt(tr_inds,:),[],silent,optim_params); %fit stimulus filters

    if xv_frac > 0
        [xvLL(cc),~,xv_predrate] = NIMmodel_eval(fit0(cc),Robsxv,all_phase_Xmat(xv_inds,:),Xexpt(xv_inds,:));
        xv_imp(cc) = (xvLL(cc)-null_xvLL(cc))/log(2);
        fprintf('LL imp: %.3f\n',xv_imp(cc));
        %     [fh] = NIMcheck_nonstat(fit0(cc),Robsxv,all_bar_mat(xv_inds,:),Xexpt(xv_inds,:),25);
    end
% NIMdisplay_model(fit0(cc));
%     pause
%     close all
end