clear all
% close all

Expt_name = 'G086';
data_dir = ['~/Data/bruce/' Expt_name];
cd(data_dir);

load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

anal_dir = ['~/Analysis/bruce/' Expt_name];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end

anal_name = 'binoc_eyecorr';


%%
stim_fs = 100; %in Hz
dt = 0.01;
flen = 10;
tot_nPix = 36;
use_nPix = 24;
full_stim_params = NIMcreate_stim_params([flen 2*tot_nPix],dt);
Fr = 1;
bar_ori = 90;

%%
exclude_expts = {'rls.orXme'};
for ii = 1:length(Expts)
    expt_names{ii} = Expts{ii}.Header.expname;
    expt_dds(ii) = Expts{ii}.Stimvals.dd;
    expt_bar_ori(ii) = Expts{ii}.Stimvals.or;
    expt_sac_dir(ii) = mod(Expts{ii}.Stimvals.Fa,180);
    expt_Fr(ii) = Expts{ii}.Stimvals.Fr;
    excluded_type(ii) = strcmp(expt_names{ii},exclude_expts);
end

cur_expt_set = find(expt_binoc(:) == 1 & expt_bar_ori(:) == 90 & expt_sac_dir(:) == 90 & ...
    expt_dds(:) == 67 & expt_Fr(:) == 1 & expt_npix(:) == 36 & ~excluded_type(:));
% cur_expt_set = find(expt_binoc(:) == 1 & expt_bar_ori(:) == bar_ori & expt_sac_dir(:) == 90 & ...
%     expt_dds(:) == 12 & expt_Fr(:) == 1 & ~excluded_type(:));

%% COMPUTE TRIAL DATA
all_stim_times = [];
all_used_inds = false(0);
Xmat = [];
all_binned_spks = [];
all_exptvec = [];
all_trialvec = [];
all_tsince_start = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_spk_t = cell(96,1);
all_spk_xy = cell(96,1);
all_spk_clst = cell(96,1);
trial_cnt = 0;
for ee = 1:length(cur_expt_set);
    fprintf('Expt %d of %d\n',ee,length(cur_expt_set));
    cur_expt = cur_expt_set(ee);
    fname = sprintf('Expt%dClusterTimesDetails.mat',cur_expt);
    load(fname);
    %initialize spike index cell array
    for cc = 1:96
        cur_spk_used{cc} = false(length(ClusterDetails{cc}.t),1);
    end
    
    fname = sprintf('%s/stims/Expt%d_stim',data_dir,cur_expt);
    load(fname);
    
    trial_start_times = [Expts{cur_expt}.Trials(:).TrialStart]/1e4;
    trial_end_times = [Expts{cur_expt}.Trials(:).TrueEnd]/1e4;
    trial_durs = [Expts{cur_expt}.Trials(:).dur]/1e4;
    trial_ids = [Expts{cur_expt}.Trials(:).id];
    [un_ids,id_inds] = unique(trial_ids);
    
    buffer_pix = floor((expt_npix(cur_expt) - tot_nPix)/2);
    cur_use_pix = (1:tot_nPix) + buffer_pix;
    
    trial_durs = trial_durs(id_inds);
    trial_start_times = trial_start_times(id_inds);
    trial_end_times = trial_end_times(id_inds);
    
    all_trial_start_times = [all_trial_start_times trial_start_times];
    all_trial_end_times = [all_trial_end_times trial_end_times];
    
    use_trials = find(trial_durs >= 0.5);
    extra_trials = find(use_trials > length(left_stim_mats));
    if ~isempty(extra_trials)
        fprintf('Dropping %d trials with missing stim data\n',length(extra_trials));
        use_trials(extra_trials) = [];
    end
    n_trials = length(use_trials);
    for tt = 1:n_trials
        cur_stim_times = Expts{cur_expt}.Trials(use_trials(tt)).Start/1e4;
        
        if length(cur_stim_times) == 1 %for blocks where the trial stimulus time stamps arent recorded
            n_frames = size(left_stim_mats{use_trials(tt)},1);
            cur_stim_times = (cur_stim_times:dt*Fr:(cur_stim_times + (n_frames-1)*dt*Fr))';
            cur_t_edges = [cur_stim_times; cur_stim_times(end) + dt*Fr];
        else %for trials where we do have the stim time stamps recorded
            cur_t_edges = [cur_stim_times; Expts{cur_expt}.Trials(use_trials(tt)).End(end)/1e4];
        end
        
        
        cur_stim_mat = zeros(length(cur_stim_times),tot_nPix*2);
        cur_stim_mat(:,1:tot_nPix) = left_stim_mats{use_trials(tt)}(:,cur_use_pix);
        cur_stim_mat(:,tot_nPix+1:end) = right_stim_mats{use_trials(tt)}(:,cur_use_pix);
        cur_lstim_mat = left_stim_mats{use_trials(tt)}(:,cur_use_pix);
        cur_rstim_mat = right_stim_mats{use_trials(tt)}(:,cur_use_pix);
        
%         cur_binned_spks = nan(length(cur_stim_times),96);
%         for cc = 1:96
%             [cur_hist,cur_bins] = histc(ClusterDetails{cc}.t,cur_t_edges);
%             cur_binned_spks(:,cc) = cur_hist(1:end-1);
%             cur_spk_used{cc}(cur_bins(1:end-1) > 0) = true;
%         end
        
        bar_Xmat = create_time_embedding(cur_stim_mat,full_stim_params);
        
        cur_used_inds = true(length(cur_stim_times),1);
        cur_used_inds(1:flen) = false;
        
        all_stim_times = [all_stim_times; cur_stim_times];
        all_tsince_start = [all_tsince_start; cur_stim_times - trial_start_times(tt)];
        all_used_inds = [all_used_inds; cur_used_inds];
        Xmat = [Xmat; bar_Xmat];
        all_binned_spks = [all_binned_spks; cur_binned_spks];
        all_exptvec = [all_exptvec; ones(size(cur_stim_times))*ee];
        all_trialvec = [all_trialvec; ones(size(cur_stim_times))*(tt + trial_cnt)];
    end
    for cc = 1:96
        %rotate xy points so that they are all in registry (across expts)
        cur_angle = ClusterDetails{cc}.angle;
        if ee == 1
            base_ang(cc) = cur_angle;
        end
        cur_rot = cur_angle - base_ang(cc);
        rot_mat = [cos(cur_rot) -sin(cur_rot); sin(cur_rot) cos(cur_rot)];
%         cur_spk_xy = ClusterDetails{cc}.xy(cur_spk_used{cc},:)*rot_mat;
        
%         cur_spk_clst = ClusterDetails{cc}.clst(cur_spk_used{cc}); %store assigned cluster indices
%        all_spk_xy{cc} = cat(1,all_spk_xy{cc},cur_spk_xy);
%         all_spk_t{cc} = cat(1,all_spk_t{cc},cur_spk_t);
%         all_spk_clst{cc} = cat(1,all_spk_clst{cc},cur_spk_clst);
       all_spk_xy{cc} = cat(1,all_spk_xy{cc},ClusterDetails{cc}.xy*rot_mat);
        all_spk_t{cc} = cat(1,all_spk_t{cc},ClusterDetails{cc}.t');
         all_spk_clst{cc} = cat(1,all_spk_clst{cc},ClusterDetails{cc}.clst);
    end
    trial_cnt = trial_cnt + 1;
end
all_t_axis = all_stim_times;

%% normalize stimulus variance
Xmat = Xmat/std(reshape(Xmat(all_used_inds,:),1,[]));

%% SET UP USEABLE TRIALS
% [c,ia,ic] = unique([all_exptvec all_trialvec],'rows');
[c,ia,ic] = unique(all_trialvec);
n_trials = length(ia);
% bad_trials = unique(ic(bad_pts)); %trials with putative blinks

xv_frac = 0.2;
n_xv_trials = round(n_trials*xv_frac);
xv_set = randperm(n_trials);
xv_set(n_xv_trials+1:end) = [];
xv_inds = find(ismember(ic,xv_set));
tr_inds = find(~ismember(ic,xv_set))';

tr_inds(~all_used_inds(tr_inds)) = [];
xv_inds(~all_used_inds(xv_inds)) = [];
NT = length(tr_inds);

%% CREATE EVENT PREDICTORS FOR REAL AND SIM SACCADES (POOLING MICRO AND MACRO SACS)
Xexpt = zeros(length(all_stim_times),length(cur_expt_set)-1);
for i = 1:length(cur_expt_set)-1
    cur_set = find(all_exptvec==i);
    Xexpt(cur_set,i) = 1;
end

%%
cd(data_dir);
use_expt = cur_expt_set(1);
fname = sprintf('Expt%dClusterTimes.mat',use_expt);
load(fname);
clust_mark = nan(96,1);
clust_thresh = nan(96,1);
for cc = 1:96
    clust_thresh(cc) = Clusters{cc}.crit;
    if isfield(Clusters{cc},'marked')
        clust_mark(cc) = Clusters{cc}.marked;
    end
end
single_units = find(clust_mark == 2);

%% FIT OR LOAD INITIAL MODELS
refit_mods = 1;
save_dir = [anal_dir '/' anal_name];
if ~exist(save_dir,'dir')
    system(['mkdir ' save_dir]);
end
init_mods_fname = [save_dir '/init_mods_d67_threshvar' '.mat'];

full_klen = flen*2*use_nPix; half_klen = flen*use_nPix;
stim_params = NIMcreate_stim_params([flen 2*use_nPix],dt,1,1,length(cur_expt_set)-1);
stim_params1 = NIMcreate_stim_params([flen use_nPix],dt,1,1,length(cur_expt_set)-1);
optim_params.progTol = 1e-7;
optim_params.optTol = 1e-4;
silent = 1;
NL_types = {'lin', 'quad'};

%find which dimensions of full Xmat were using for model fitting
buffer_pix = floor((tot_nPix-use_nPix)/2);
use_pix = (1:use_nPix) + buffer_pix;
[PixPix,LagLag] = meshgrid(1:tot_nPix,1:flen);
PixPix = [PixPix PixPix];
PixPix = PixPix(:);
use_k = find(ismember(PixPix,use_pix));

n_prc_bins = 20;
thresh_prctiles = nan(96,n_prc_bins);

l_d2XT = 5000;
ql_d2XT = 1000;
XTmix = [1 1];

reg_params = NIMcreate_reg_params('lambda_d2XT',l_d2XT,'temporal_boundaries','zero');
null_params = NIMcreate_stim_params([1 size(Xexpt,2)],dt);

%create binocular d2XT Laplacian
et = ones(flen,1)*sqrt(XTmix(2));
ex = ones(use_nPix,1)*sqrt(XTmix(1));
D1t = spdiags([et -2*et et], [-1 0 1], flen, flen)';
D1x = spdiags([ex -2*ex ex], [-1 0 1], use_nPix, use_nPix)';
It = speye(flen);Ix = speye(use_nPix(1));
L2_d2XT = kron(Ix,D1t) + kron(D1x,It);

full_L2_mats.L2_d2XT = speye(full_klen);
full_L2_mats.L2_d2XT(1:half_klen,1:half_klen) = L2_d2XT;
full_L2_mats.L2_d2XT((half_klen+1):end,(half_klen+1):end) = L2_d2XT;

if refit_mods == 1 || ~exist(init_mods_fname,'file')
    fprintf('Computing initial models\n');
    
    %% first for MUA
    for cc = 1:96
        fprintf('Fitting MUA %d of %d\n',cc,96);
%         spk_bins = convert_to_spikebins(all_binned_spks(:,cc));
        % use_set = find(all_spk_xx{cc} >= use_thresh);
        
%         if length(spk_bins) ~= length(all_spk_clst{cc})
%             error('Problem with spike processing!');
%         end
        
        %for probes without a SU
        if ~ismember(cc,single_units)
            %             use_set = find(all_spk_clst{cc} == 2);
            med_thresh(cc) = median(all_spk_xy{cc}(:,1));
%             use_set = find(all_spk_xy{cc}(:,1) >= clust_thresh(cc));
            use_set = find(all_spk_xy{cc}(:,1) >= med_thresh(cc));
            poss_mua_set = 1:length(all_spk_clst{cc});
        else
            poss_mua_set = find(all_spk_clst{cc} ~= 2);
            med_thresh(cc) = median(all_spk_xy{cc}(poss_mua_set,1));
            use_set = poss_mua_set(all_spk_xy{cc}(poss_mua_set,1) > med_thresh(cc));
        end
%         use_binned_spks = hist(spk_bins(use_set),1:length(all_stim_times));
        cur_spike_times = all_spk_t{cc}(use_set);
        use_spike_bins = round(interp1(all_stim_times,1:length(all_stim_times),cur_spike_times));
        bad = find(isnan(use_spike_bins)); use_spike_bins(bad) = []; cur_spike_times(bad) = []; poss_mua_set(bad) = [];
        tdiffs = abs(all_stim_times(use_spike_bins) - cur_spike_times);
        too_far = find(tdiffs > dt/2);
        use_spike_bins(too_far) = [];
        use_binned_spks = convert_to_binned_spks(use_spike_bins,length(all_stim_times));
        
        fprintf('Mean spike rate %.3f\n',mean(use_binned_spks(tr_inds))/dt);
        Robs = use_binned_spks(tr_inds);
        n_spks(cc) = sum(Robs);
        mean_rate(cc) = mean(Robs)/dt;
        
        null_mod(cc) = NIMinitialize_model(null_params,1,{'lin'},reg_params); %initialize NIM
        null_mod(cc) = NIMfit_filters(null_mod(cc),Robs,Xexpt(tr_inds,:)); %fit stimulus filters
        
        if xv_frac > 0
            Robsxv = use_binned_spks(xv_inds);
            xv_nspks(cc) = sum(Robsxv);
            null_xvLL(cc) = NIMmodel_eval(null_mod(cc),Robsxv,Xexpt(xv_inds,:))*xv_nspks(cc);
        end
        
        mod_signs = 1;
        fit0(cc,1) = NIMinitialize_model(stim_params,mod_signs,NL_types,reg_params); %initialize NIM
        fit0(cc,1) = NIMfit_filters(fit0(cc,1),Robs,Xmat(tr_inds,use_k),Xexpt(tr_inds,:),[],silent,optim_params,full_L2_mats); %fit stimulus filters
        %         fit0(cc,1) = NIMfit_logexp_spkNL(fit0(cc,1),Robs,Xmat(tr_inds,:),Xexpt(tr_inds,:));
        if xv_frac > 0
            [xvLL(cc,1),~,xv_predrate] = NIMmodel_eval(fit0(cc,1),Robsxv,Xmat(xv_inds,use_k),Xexpt(xv_inds,:));
            xvLL(cc,1) = xvLL(cc,1)*xv_nspks(cc);
            xv_imp(cc,1) = (xvLL(cc,1)-null_xvLL(cc))/log(2);
            %         fprintf('LL imp: %.3f\n',xv_imp(cc,1));
        end
        
        mod_signs = [1 1];
        %add a filter
        rand_filt = randn(prod(stim_params.stim_dims),1)/prod(stim_params.stim_dims) * 1;
        fit0(cc,2) = NIMadd_NLinput(fit0(cc,1),'quad',1,rand_filt);
        %adjust regularization of new filter
        fit0(cc,2) = NIMadjust_regularization(fit0(cc,2),[length(mod_signs)],'lambda_d2XT',ql_d2XT);
        fit0(cc,2) = NIMfit_filters(fit0(cc,2),Robs,Xmat(tr_inds,use_k),Xexpt(tr_inds,:),[],silent,optim_params,full_L2_mats); %fit stimulus filters
        %         fit0(cc,2) = NIMfit_logexp_spkNL(fit0(cc,2),Robs,Xmat(tr_inds,use_k),Xexpt(tr_inds,:),0);
        
        if xv_frac > 0
            [xvLL(cc,2),~,xv_predrate] = NIMmodel_eval(fit0(cc,2),Robsxv,Xmat(xv_inds,use_k),Xexpt(xv_inds,:));
            xvLL(cc,2) = xvLL(cc,2)*xv_nspks(cc);
            xv_imp(cc,2) = (xvLL(cc,2)-null_xvLL(cc))/log(2);
            fprintf('LL imp 1filt: %.4f  2filt: %.4f\n',xv_imp(cc,1),xv_imp(cc,2));
        end
        
        [init_xvimp(cc),best_mod_index(cc)] = max(xv_imp(cc,:));
        init_mod(cc) = fit0(cc,best_mod_index(cc));
        init_mod(cc) = NIMfit_logexp_spkNL(init_mod(cc),Robs,Xmat(tr_inds,use_k),Xexpt(tr_inds,:),0);
        [temp_xvLL] = NIMmodel_eval(init_mod(cc),Robsxv,Xmat(xv_inds,use_k),Xexpt(xv_inds,:))*xv_nspks(cc);
        init_xvimp(cc) = (temp_xvLL - null_xvLL(cc))/log(2);
        
        thresh_prctiles(cc,:) = my_prctile(all_spk_xy{cc}(poss_mua_set,1),(100/n_prc_bins/2):(100/n_prc_bins):(100-100/n_prc_bins/2));
        for ii = 1:n_prc_bins
            use_set = poss_mua_set(all_spk_xy{cc}(poss_mua_set,1) >= thresh_prctiles(cc,ii));
            %             use_binned_spks = hist(spk_bins(use_set),1:length(all_stim_times));
            cur_spike_times = all_spk_t{cc}(use_set);
            use_spike_bins = round(interp1(all_stim_times,1:length(all_stim_times),cur_spike_times));
            bad = find(isnan(use_spike_bins)); use_spike_bins(bad) = []; cur_spike_times(bad) = [];
            tdiffs = abs(all_stim_times(use_spike_bins) - cur_spike_times);
            too_far = find(tdiffs > dt/2);
            use_spike_bins(too_far) = [];
            use_binned_spks = convert_to_binned_spks(use_spike_bins,length(all_stim_times));
            %             fprintf('Mean spike rate %.3f\n',mean(use_binned_spks(tr_inds))/dt);
            Robs = use_binned_spks(tr_inds);
            %             fitr(cc,ii) = NIMfit_filters(fit0(cc,best_mod_index(cc)),Robs,Xmat(tr_inds,use_k),Xexpt(tr_inds,:),[],silent); %fit stimulus filters
            fitr(cc,ii) = NIMfit_logexp_spkNL(fit0(cc,best_mod_index(cc)),Robs,Xmat(tr_inds,use_k),Xexpt(tr_inds,:)); %fit stimulus filters
            n_spksr(cc,ii) = sum(Robs);
            LLr(cc,ii) = NIMmodel_eval(fitr(cc,ii),Robs,Xmat(tr_inds,use_k),Xexpt(tr_inds,:))*n_spksr(cc,ii);
            mean_rater(cc,ii) = mean(Robs)/dt;
            
            null_modr(cc,ii) = NIMfit_filters(null_mod(cc),Robs,Xexpt(tr_inds,:)); %fit stimulus filters
            null_LLr(cc,ii) = NIMmodel_eval(null_modr(cc,ii),Robs,Xexpt(tr_inds,:))*n_spksr(cc,ii);
            
            if xv_frac > 0
                Robsxv = use_binned_spks(xv_inds);
                xv_n_spksr(cc,ii) = sum(Robsxv);
                xvLLr(cc,ii) = NIMmodel_eval(fitr(cc,ii),Robsxv,Xmat(xv_inds,use_k),Xexpt(xv_inds,:))*xv_n_spksr(cc,ii);
                null_xvLLr(cc,ii) = NIMmodel_eval(null_modr(cc,ii),Robsxv,Xexpt(xv_inds,:))*xv_n_spksr(cc,ii);
                xv_impr(cc,ii) = (xvLLr(cc,ii)-null_xvLLr(cc,ii))/log(2);
            end
        end
        
        [best_tvar_xvimp(cc),best_thresh(cc)] = max(xv_impr(cc,:));
        init_modr(cc) = fitr(cc,best_thresh(cc));
        
        fprintf('Base xvLL: %.4f at %.3f mean rate\n',init_xvimp(cc),mean_rate(cc));
        fprintf('Best thresh-var xvLL: %.4f at %.3f mean rate\n',best_tvar_xvimp(cc),mean_rater(cc,best_thresh(cc)));
        
        if xv_frac > 0
            save(init_mods_fname,'fit*','best_mod*','null*','xv_imp*','LL*','null_xvLL*','xvLL*','thresh_prctiles','*n_spks*');
        else
            save(init_mods_fname,'fit*','best_mod*','null*','LL*','thresh_prctiles','*n_spks*');
        end
    end
    
    %% NOW for SUA
    for cc = 1:length(single_units)
        fprintf('SU %d of %d\n',cc,length(single_units));
        
        use_set = find(all_spk_clst{single_units(cc)} == 2);
        cur_spike_times = all_spk_t{single_units(cc)}(use_set);
        use_spike_bins = round(interp1(all_stim_times,1:length(all_stim_times),cur_spike_times));
        bad = find(isnan(use_spike_bins)); use_spike_bins(bad) = []; cur_spike_times(bad) = [];
        tdiffs = abs(all_stim_times(use_spike_bins) - cur_spike_times);
        too_far = find(tdiffs > dt/2);
        use_spike_bins(too_far) = [];
        use_binned_spks = convert_to_binned_spks(use_spike_bins,length(all_stim_times));

%         spk_bins = convert_to_spikebins(all_binned_spks(:,single_units(cc)));
%         use_set = find(all_spk_clst{single_units(cc)} == 2);
%         use_binned_spks = hist(spk_bins(use_set),1:length(all_stim_times));
%         fprintf('Mean spike rate %.3f\n',mean(use_binned_spks(tr_inds))/dt);
        Robs = use_binned_spks(tr_inds);
        n_spks(96+cc) = sum(Robs);
        mean_rate(96+cc) = mean(Robs)/dt;
        
        null_mod(96+cc) = NIMinitialize_model(null_params,1,{'lin'},reg_params); %initialize NIM
        null_mod(96+cc) = NIMfit_filters(null_mod(96+cc),Robs,Xexpt(tr_inds,:)); %fit stimulus filters
        
        if xv_frac > 0
            Robsxv = use_binned_spks(xv_inds);
            xv_nspks(96+cc) = sum(Robsxv);
            null_xvLL(96+cc) = NIMmodel_eval(null_mod(96+cc),Robsxv,Xexpt(xv_inds,:))*xv_nspks(96+cc);
        end
        
        mod_signs = 1;
        fit0(96+cc,1) = NIMinitialize_model(stim_params,mod_signs,NL_types,reg_params); %initialize NIM
        fit0(96+cc,1) = NIMfit_filters(fit0(96+cc,1),Robs,Xmat(tr_inds,use_k),Xexpt(tr_inds,:),[],silent,optim_params,full_L2_mats); %fit stimulus filters
        if xv_frac > 0
            [xvLL(96+cc,1),~,xv_predrate] = NIMmodel_eval(fit0(96+cc,1),Robsxv,Xmat(xv_inds,use_k),Xexpt(xv_inds,:));
            xvLL(96+cc,1) = xvLL(96+cc,1)*xv_nspks(96+cc);
            xv_imp(96+cc,1) = (xvLL(96+cc,1)-null_xvLL(96+cc))/log(2);
        end
        
        mod_signs = [1 1];
        %add a filter
        rand_filt = randn(prod(stim_params.stim_dims),1)/prod(stim_params.stim_dims) * 1;
        fit0(96+cc,2) = NIMadd_NLinput(fit0(96+cc,1),'quad',1,rand_filt);
        %adjust regularization of new filter
        fit0(96+cc,2) = NIMadjust_regularization(fit0(96+cc,2),[length(mod_signs)],'lambda_d2XT',ql_d2XT);
        fit0(96+cc,2) = NIMfit_filters(fit0(96+cc,2),Robs,Xmat(tr_inds,use_k),Xexpt(tr_inds,:),[],silent,optim_params,full_L2_mats); %fit stimulus filters
        
        if xv_frac > 0
            [xvLL(96+cc,2),~,xv_predrate] = NIMmodel_eval(fit0(96+cc,2),Robsxv,Xmat(xv_inds,use_k),Xexpt(xv_inds,:));
            xvLL(96+cc,2) = xvLL(96+cc,2)*xv_nspks(96+cc);
            xv_imp(96+cc,2) = (xvLL(96+cc,2)-null_xvLL(96+cc))/log(2);
            fprintf('LL imp 1filt: %.4f  2filt: %.4f\n',xv_imp(96+cc,1),xv_imp(96+cc,2));
        end
        
        [init_xvimp(96+cc),best_mod_index(96+cc)] = max(xv_imp(96+cc,:));
        init_mod(96+cc) = fit0(96+cc,best_mod_index(96+cc));
        init_mod(96+cc) = NIMfit_logexp_spkNL(init_mod(96+cc),Robs,Xmat(tr_inds,use_k),Xexpt(tr_inds,:),0);
        [temp_xvLL] = NIMmodel_eval(init_mod(96+cc),Robsxv,Xmat(xv_inds,use_k),Xexpt(xv_inds,:))*xv_nspks(96+cc);
        init_xvimp(96+cc) = (temp_xvLL - null_xvLL(96+cc))/log(2);
                
        fprintf('Base xvLL: %.4f at %.3f mean rate\n',init_xvimp(96+cc),mean_rate(96+cc));
        
        if xv_frac > 0
            save(init_mods_fname,'fit*','best_mod*','null*','xv_imp*','LL*','null_xvLL*','xvLL*','thresh_prctiles','*n_spks*');
        else
            save(init_mods_fname,'fit*','best_mod*','null*','LL*','thresh_prctiles','*n_spks*');
        end
    end
    
else
    fprintf('Loading precomputed models\n');
    load(init_mods_fname);
end

%% INITIALIZE MODELS TO USE
thresh_xvimp = 20;
for cc = 1:96
    if xv_imp(cc,2) > xv_imp(cc,1)
        use_quad(cc) = 1;
    else
        use_quad(cc) = 0;
    end
    
    if best_tvar_xvimp(cc) > thresh_xvimp
        useable_unit(cc) = 1;
    else
        useable_unit(cc) = 0;
    end    
end

for cc = 1:length(single_units)
    if xv_imp(96+cc,2) > xv_imp(96+cc,1)
        use_quad(96+cc) = 1;
    else
        use_quad(96+cc) = 0;
    end
    useable_unit(96+cc) = 1;
end

%% REDEFINE USED TIME POINTS AND RE-TUNE MODELS TO FULL TR SET
xv_frac = 0;
xv_set = [];
xv_inds = [];
tr_inds = find(~ismember(ic,xv_set))';

tr_inds(~all_used_inds(tr_inds)) = [];
NT = length(tr_inds);

use_Xmat = Xmat(tr_inds,use_k);
all_use_binned_spks = nan(size(fit0,1),length(tr_inds));

silent = 1;
for cc = 1:96
    fprintf('MU %d of %d\n',cc,96);
%     spk_bins = convert_to_spikebins(all_binned_spks(:,cc));
    
    %for probes without a SU
    if ~ismember(cc,single_units)
        use_set = find(all_spk_xy{cc}(:,1) >= thresh_prctiles(cc,best_thresh(cc)));
        poss_mua_set = 1:length(all_spk_clst{cc});
    else
        poss_mua_set = find(all_spk_clst{cc} ~= 2);
        use_set = poss_mua_set(all_spk_xy{cc}(poss_mua_set,1) > thresh_prctiles(cc,best_thresh(cc)));
    end
    
    cur_spike_times = all_spk_t{cc}(use_set);
    use_spike_bins = round(interp1(all_stim_times,1:length(all_stim_times),cur_spike_times));
        bad = find(isnan(use_spike_bins)); use_spike_bins(bad) = []; cur_spike_times(bad) = []; poss_mua_set(bad) = [];
    tdiffs = abs(all_stim_times(use_spike_bins) - cur_spike_times);
    too_far = find(tdiffs > dt/2);
    use_spike_bins(too_far) = [];
    use_binned_spks = convert_to_binned_spks(use_spike_bins,length(all_stim_times));

%     use_binned_spks = hist(spk_bins(use_set),1:length(all_stim_times));

    fprintf('Mean spike rate %.3f\n',mean(use_binned_spks(tr_inds))/dt);
    Robs = use_binned_spks(tr_inds);
    all_use_binned_spks(cc,:) = Robs;
    n_spks(cc) = sum(Robs);
    init_mod(cc) = NIMfit_filters(fit0(cc,best_mod_index(cc)),Robs,use_Xmat,Xexpt(tr_inds,:),[],silent,optim_params,full_L2_mats); %fit stimulus filters
    init_mod_LL(cc) = init_mod(cc).LL_seq(end)*n_spks(cc);
    null_mod(cc) = NIMfit_filters(null_mod(cc),Robs,Xexpt(tr_inds,:)); %fit stimulus filters
    null_LL(cc) = null_mod(cc).LL_seq(end)*n_spks(cc);
end

for cc = 1:length(single_units)
%     spk_bins = convert_to_spikebins(all_binned_spks(:,single_units(cc)));
        use_set = find(all_spk_clst{single_units(cc)} == 2);
        
        cur_spike_times = all_spk_t{single_units(cc)}(use_set);
        use_spike_bins = round(interp1(all_stim_times,1:length(all_stim_times),cur_spike_times));
        bad = find(isnan(use_spike_bins)); use_spike_bins(bad) = []; cur_spike_times(bad) = [];
        tdiffs = abs(all_stim_times(use_spike_bins) - cur_spike_times);
        too_far = find(tdiffs > dt/2);
        use_spike_bins(too_far) = [];
        use_binned_spks = convert_to_binned_spks(use_spike_bins,length(all_stim_times));
        
%         use_binned_spks = hist(spk_bins(use_set),1:length(all_stim_times));
        fprintf('Mean spike rate %.3f\n',mean(use_binned_spks(tr_inds))/dt);
        Robs = use_binned_spks(tr_inds);
        all_use_binned_spks(96+cc,:) = Robs;
        n_spks(96+cc) = sum(Robs);
        init_mod(96+cc) = NIMfit_filters(fit0(96+cc,best_mod_index(96+cc)),Robs,use_Xmat,Xexpt(tr_inds,:),[],silent,optim_params,full_L2_mats); %fit stimulus filters
        init_mod_LL(96+cc) = init_mod(96+cc).LL_seq(end)*n_spks(96+cc);
        null_mod(96+cc) = NIMfit_filters(null_mod(96+cc),Robs,Xexpt(tr_inds,:)); %fit stimulus filters
        null_LL(96+cc) = null_mod(96+cc).LL_seq(end)*n_spks(96+cc);
end

init_mod_LL_imp = (init_mod_LL - null_LL)/log(2);

save(init_mods_fname,'fit*','best_mod*','null*','xv_imp*','LL*','null_xvLL*',...
    'xvLL*','thresh_prctiles','*n_spks*','init_*');


%% PARSE EYE-TRACKING DATA (IGNORING RIGHT EYE BECAUSE OF NOISE)
emfile = [data_dir '/jbe' Expt_name '.em.mat'];
load(emfile);

all_eye_vals = [];
all_eye_speed = [];
all_eye_ts = [];
all_eye_exptvec = [];
eye_smooth = 3;
for ee = 1:length(cur_expt_set);
    fprintf('Loading eye data for expt %d of %d\n',ee,length(cur_expt_set));
    cur_set = find(all_exptvec==ee);
    [eye_vals_interp,eye_ts_interp,eye_insample] = get_eye_tracking_data_new(Expt,all_t_axis(cur_set([1 end])));
    
    eye_dt = Expt.Header.CRrates(1);
    eye_fs = 1/eye_dt;
    lEyeXY = eye_vals_interp(:,1:2);
    rEyeXY = eye_vals_interp(:,3:4);
    
    %slight smoothing before computing speed
    sm_avg_eyepos = lEyeXY; eye_vel = lEyeXY; %initialization
    sm_avg_eyepos(:,1) = smooth(lEyeXY(:,1),eye_smooth);
    sm_avg_eyepos(:,2) = smooth(lEyeXY(:,2),eye_smooth);
    eye_vel(:,1) = [0; diff(sm_avg_eyepos(:,1))]/eye_dt;
    eye_vel(:,2) = [0; diff(sm_avg_eyepos(:,2))]/eye_dt;
    
    eye_speed = sqrt(eye_vel(:,1).^2+eye_vel(:,2).^2);
    
    all_eye_vals = [all_eye_vals; lEyeXY rEyeXY];
    all_eye_speed = [all_eye_speed; eye_speed];
    all_eye_ts = [all_eye_ts; eye_ts_interp'];
    all_eye_exptvec = [all_eye_exptvec; ee*ones(size(eye_speed))];
end

back_pts = 1 + find(diff(all_eye_ts) <= 0);
double_samples = [];
for i = 1:length(back_pts)
    next_forward = find(all_eye_ts > all_eye_ts(back_pts(i)-1),1,'first');
    double_samples = [double_samples back_pts(i):next_forward];
end
all_eye_ts(double_samples) = [];
all_eye_speed(double_samples) = [];
all_eye_vals(double_samples,:) = [];

interp_eye_vals = interp1(all_eye_ts,all_eye_vals,all_t_axis);

orth_eye_pos = interp_eye_vals(:,1)*cos(bar_ori*pi/180+pi/2) + interp_eye_vals(:,2)*sin(bar_ori*pi/180 + pi/2);
par_eye_pos = interp_eye_vals(:,1)*cos(bar_ori*pi/180) + interp_eye_vals(:,2)*sin(bar_ori*pi/180);
bad_pts = find(abs(orth_eye_pos) > 1); %cheap way of detecting blinks because the right eye signal was unreliable even for that at times

%% detect saccades
sac_thresh = 10;
peak_sig = [0; diff(sign(diff(all_eye_speed))); 0];
saccade_inds = find(peak_sig == -2 & all_eye_speed > sac_thresh);

peri_thresh = 3; %threshold eye speed for defining saccade boundary inds
thresh_cross_up = 1 + find(all_eye_speed(1:end-1) < peri_thresh & all_eye_speed(2:end) >= peri_thresh);
thresh_cross_down = 1 + find(all_eye_speed(1:end-1) >= peri_thresh & all_eye_speed(2:end) < peri_thresh);
sac_start_inds = nan(size(saccade_inds));
sac_stop_inds = nan(size(saccade_inds));
for ii = 1:length(saccade_inds)
    next_tc = find(thresh_cross_down > saccade_inds(ii),1,'first');
    if ~isempty(next_tc)
        sac_stop_inds(ii) = thresh_cross_down(next_tc);
    end
    prev_tc = find(thresh_cross_up < saccade_inds(ii),1,'last');
    if ~isempty(prev_tc)
        sac_start_inds(ii) = thresh_cross_up(prev_tc);
    end
    
end

%get rid of double-peaks
min_isi = 0.05; max_isi = Inf;
isis = [Inf; diff(sac_start_inds)]/eye_fs;
bad_isis = (isis < min_isi | isis > max_isi);
bad_sacs = find(isnan(sac_stop_inds) | isnan(sac_start_inds) | bad_isis);
saccade_inds(bad_sacs) = []; isis(bad_sacs) = []; sac_start_inds(bad_sacs) = []; sac_stop_inds(bad_sacs) = [];

saccade_times = all_eye_ts(saccade_inds);
sac_start_times = all_eye_ts(sac_start_inds);
sac_stop_times = all_eye_ts(sac_stop_inds);
sac_durs = sac_stop_times - sac_start_times;

sac_dbuff = round(0.005/eye_dt);
pre_inds = saccade_inds - sac_dbuff;
pre_inds(pre_inds < 1) = 1;
sac_pre_pos = all_eye_vals(pre_inds,:);
post_inds = saccade_inds + sac_dbuff;
post_inds(post_inds > length(all_eye_ts)) = length(all_eye_ts);
sac_post_pos = all_eye_vals(post_inds,:);

%use only left-eye signal here
sac_delta_pos = sac_post_pos(:,1:2) - sac_pre_pos(:,1:2);
sac_amps = sqrt(sum(sac_delta_pos.^2,2));
sac_dirs = atan2(sac_delta_pos(:,2),sac_delta_pos(:,1));

temp = ones(length(saccade_times),1);
saccades = struct('peak_time',mat2cell(saccade_times,temp),'start_time',mat2cell(sac_start_times,temp),...
    'stop_time',mat2cell(sac_stop_times,temp),'isi',mat2cell(isis,temp),...
    'duration',mat2cell(sac_durs,temp),'amplitude',mat2cell(sac_amps,temp),'direction',mat2cell(sac_dirs,temp),...
    'pre_Lx',mat2cell(sac_pre_pos(:,1),temp),'post_Lx',mat2cell(sac_post_pos(:,1),temp),...
    'pre_Ly',mat2cell(sac_pre_pos(:,2),temp),'post_Ly',mat2cell(sac_post_pos(:,2),temp),...
    'pre_Rx',mat2cell(sac_pre_pos(:,3),temp),'post_Rx',mat2cell(sac_post_pos(:,3),temp),...
    'pre_Ry',mat2cell(sac_pre_pos(:,4),temp),'post_Ry',mat2cell(sac_post_pos(:,4),temp));

% is_blink = find(sac_durs > 0.1);
% saccades(is_blink) = [];

sac_start_times = [saccades(:).start_time];
interp_sac_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_start_times));
bad = find(isnan(interp_sac_start_inds));
interp_sac_start_inds(isnan(interp_sac_start_inds)) = [];
saccades(bad) = []; sac_start_times(bad) = [];

%throw out sacs where the interpolated time is too far from the true time
interp_sac_start_times = all_t_axis(interp_sac_start_inds);
time_diff = abs(interp_sac_start_times(:) - sac_start_times(:));
bad = find(time_diff > dt);
saccades(bad) = [];
interp_sac_start_inds(bad) = [];

%% DEFINE POINTS TO ALLOW RAPID CHANGES (TRIAL STARTS AFTER SACCADES)
trial_start_inds = [1; find(diff(all_trialvec(tr_inds)) ~= 0) + 1];

used_saccade_inds = find(ismember(tr_inds,interp_sac_start_inds));
used_saccade_set = find(ismember(interp_sac_start_inds,tr_inds));

%push the effects of saccades forward in time
sac_shift = round(0.05/dt);

%define times when to resort to independent prior
use_prior = false(size(tr_inds));
use_prior(trial_start_inds) = true;
for i = 1:length(used_saccade_inds)
    next_trial = trial_start_inds(find(trial_start_inds > used_saccade_inds(i),1,'first'));
    %mark a point (forward in time) as a saccade if it occurs within the
    %same trial as the saccade
    if next_trial > used_saccade_inds(i) + sac_shift
        use_prior(used_saccade_inds(i) + sac_shift) = 1;
    end
end

%% INITIALIZE TRANSITION PRIORS FOR HMM
NT = length(tr_inds);
sp_dx = 0.05;
max_shift = 12;
dshift = 1;
shifts = -max_shift:dshift:max_shift;
n_shifts = length(shifts);
zero_frame = find(shifts==0);

%overall prior on shifts
eps_prior_sigma = 0.125; %0.125 start
leps_prior = -(shifts*sp_dx).^2/(2*eps_prior_sigma^2);
leps_prior = bsxfun(@minus,leps_prior,logsumexp(leps_prior)); %normalize
lA_tflip = repmat(leps_prior,n_shifts,1);

cdist = squareform(pdist(shifts'*sp_dx));
deps_sigma = 0.015; %.015
lA = -cdist.^2/(2*deps_sigma^2);
lA = bsxfun(@minus,lA,logsumexp(lA,2)); %normalize

%% ESTIMATE LL for each shift in each stimulus frame
tr_set = find(useable_unit);
n_tr_chs = length(tr_set);

%create filter banks for all units
lin_filt_bank = zeros(length(tr_set),full_klen);
quad_filt_bank = zeros(length(tr_set),full_klen);
lin_kerns = zeros(length(tr_set),stim_params.lin_dims);
mod_offsets = nan(1,length(tr_set));
for cc = 1:length(tr_set)
    cur_k = [init_mod(tr_set(cc)).mods(:).filtK];
    lin_filt_bank(cc,:) = cur_k(:,1);
    if use_quad(tr_set(cc))
        quad_filt_bank(cc,:) = cur_k(:,2);
    end
    mod_offsets(cc) = init_mod(tr_set(cc)).spk_NL_params(1);
    lin_kerns(cc,:) = init_mod(tr_set(cc)).kLin';
end
lin_filt_bank = reshape(lin_filt_bank',[stim_params.stim_dims(1:2) length(tr_set)]);
quad_filt_bank = reshape(quad_filt_bank',[stim_params.stim_dims(1:2) length(tr_set)]);

tot_left_eye_inds = false(2*tot_nPix,1); tot_left_eye_inds(1:tot_nPix) = true;
use_left_eye_inds = false(2*use_nPix,1); use_left_eye_inds(1:use_nPix) = true;
lin_Lfilts = lin_filt_bank(:,use_left_eye_inds,:);
lin_Rfilts = lin_filt_bank(:,~use_left_eye_inds,:);
quad_Lfilts = quad_filt_bank(:,use_left_eye_inds,:);
quad_Rfilts = quad_filt_bank(:,~use_left_eye_inds,:);

%add zero-padding to match stimulus
lin_Lfilts = cat(2,lin_Lfilts,zeros(flen,buffer_pix,length(tr_set)));
lin_Lfilts = cat(2,zeros(flen,buffer_pix,length(tr_set)),lin_Lfilts);
lin_Rfilts = cat(2,lin_Rfilts,zeros(flen,buffer_pix,length(tr_set)));
lin_Rfilts = cat(2,zeros(flen,buffer_pix,length(tr_set)),lin_Rfilts);
quad_Lfilts = cat(2,quad_Lfilts,zeros(flen,buffer_pix,length(tr_set)));
quad_Lfilts = cat(2,zeros(flen,buffer_pix,length(tr_set)),quad_Lfilts);
quad_Rfilts = cat(2,quad_Rfilts,zeros(flen,buffer_pix,length(tr_set)));
quad_Rfilts = cat(2,zeros(flen,buffer_pix,length(tr_set)),quad_Rfilts);


%indicator predictions
expt_out = Xexpt(tr_inds,:)*lin_kerns';

%precompute LL at all shifts for all units
LLs = nan(NT,length(tr_set),n_shifts);
for xx = 1:length(shifts)
    fprintf('Shift %d of %d\n',xx,n_shifts);
    cur_lin_shift_L = shift_matrix_Nd(lin_Lfilts,shifts(xx),2);
    cur_lin_shift_R = shift_matrix_Nd(lin_Rfilts,shifts(xx),2);
    cur_lin_shift = cat(2,cur_lin_shift_L,cur_lin_shift_R);
    cur_lin_shift = reshape(cur_lin_shift,size(Xmat,2),n_tr_chs);
    
    cur_quad_shift_L = shift_matrix_Nd(quad_Lfilts,shifts(xx),2);
    cur_quad_shift_R = shift_matrix_Nd(quad_Rfilts,shifts(xx),2);
    cur_quad_shift = cat(2,cur_quad_shift_L,cur_quad_shift_R);
    cur_quad_shift = reshape(cur_quad_shift,size(Xmat,2),n_tr_chs);
    
    %outputs of stimulus models at current X-matrix shift
    lin_out = Xmat(tr_inds,:)*cur_lin_shift;
    quad_out = (Xmat(tr_inds,:)*cur_quad_shift).^2;
    gfun = lin_out + quad_out;
    gfun = bsxfun(@plus,gfun,mod_offsets);
    
    %add contributions from extra lin kernels
    gfun = gfun + expt_out;
    
    too_large = gfun > 50;
    pred_rate = log(1+exp(gfun));
    pred_rate(too_large) = gfun(too_large);
    
    pred_rate(pred_rate < 1e-20) = 1e-20;
    
    LLs(:,:,xx) = all_use_binned_spks(tr_set,:)'.*log(pred_rate) - pred_rate;
end

%% DO INITIAL EYE INFERENCE AND MODEL RE-FITTING FOR XV SET

poss_xv_set = tr_set;

resh_X = reshape(Xmat(tr_inds,:)',[flen 2*tot_nPix NT]);
resh_X_sh = zeros(size(resh_X));
for XV = 1:length(poss_xv_set)
    xv_cell = poss_xv_set(XV);
    
    cur_tr_set = find(tr_set~=xv_cell);
    
    fprintf('Using XV cell %d,  %d of %d\n',xv_cell,XV,length(poss_xv_set));
    frame_LLs = squeeze(sum(LLs(:,cur_tr_set,:),2));
    
    %%
    lB = frame_LLs;
    lalpha=zeros(NT,n_shifts);
    lbeta = zeros(NT,n_shifts);
    lscale=zeros(NT,1); %initialize rescaling parameters
    %compute rescaled forward messages
    lalpha(1,:) = leps_prior + lB(1,:);
    lscale(1)=logsumexp(lalpha(1,:));
    lalpha(1,:) = lalpha(1,:) - lscale(1);
    for t=2:NT
        %         if mod(t,1000)==0
        %             fprintf('%d of %d\n',t,NT);
        %         end
        if use_prior(t)
            cur_lA = lA_tflip;
        else
            cur_lA = lA;
        end
        lalpha(t,:) = logmulexp(lalpha(t-1,:),cur_lA) + lB(t,:);
        lscale(t) = logsumexp(lalpha(t,:));
        lalpha(t,:)= lalpha(t,:) - lscale(t);
    end
    
    %compute rescaled backward messages
    lbeta(NT,:)=log(ones(1,n_shifts)) - lscale(NT);
    for t=NT-1:-1:1
        %         if mod(t,1000)==0
        %             fprintf('%d of %d\n',t,NT);
        %         end
        if use_prior(t+1)
            cur_lA = lA_tflip;
        else
            cur_lA = lA;
        end
        lf1 = lbeta(t+1,:) + lB(t+1,:);
        lbeta(t,:) = logmulexp(lf1,cur_lA') - lscale(t);
    end
    lgamma= lalpha + lbeta;
    lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
    
    gamma = exp(lgamma);
    post_mean(:,XV) = sum(bsxfun(@times,gamma,shifts),2);
    post_std(:,XV) = sqrt(sum(bsxfun(@times,gamma,shifts.^2),2) - post_mean(:,XV).^2);
    
    %% RECONSTRUCT MAP STIMULUS (OLD METHOD)
    [max_post,max_loc] = max(lgamma,[],2);
    
    shift_cor(:,XV) = shifts(max_loc);
    % shift_cor(:,XV) = post_mean;
    for ii = 1:NT
        %         d2 = dist_shift2d(resh_X(:,:,ii), -shift_cor(ii,XV),2,0);
        d2_left = shift_matrix_Nd(resh_X(:,tot_left_eye_inds,ii), -shift_cor(ii,XV),2);
        d2_right = shift_matrix_Nd(resh_X(:,~tot_left_eye_inds,ii), -shift_cor(ii,XV),2);
        resh_X_sh(:,:,ii) = cat(2,d2_left,d2_right);
    end
    X_sh = reshape(resh_X_sh,size(Xmat,2),NT)';
    X_sh = X_sh(:,use_k);
    
    %% REFIT MODELS
    %adjust regularization of new filter
    ref_mod{1}(xv_cell) = init_mod(xv_cell);
    ref_mod{1}(xv_cell) = NIMfit_filters(ref_mod{1}(xv_cell),all_use_binned_spks(xv_cell,:),X_sh,Xexpt(tr_inds,:),[],silent,optim_params,full_L2_mats); %fit stimulus filters
    
    ref_LL_imp(1,xv_cell) = (ref_mod{1}(xv_cell).LL_seq(end)*n_spks(xv_cell) - null_LL(xv_cell))/log(2);
    
    fprintf('Original: %.4f  New: %.4f\n',init_mod_LL_imp(xv_cell),ref_LL_imp(1,xv_cell));
    
end

%% if there are any units in the training set that weren't allowed to be XV units, then refit those, using a all-cell eye position
if any(~ismember(tr_set,poss_xv_set))
    
    fprintf('Using XV cell %d,  %d of %d\n',xv_cell,XV,length(poss_xv_set));
    frame_LLs = squeeze(sum(LLs,2));
    
    %%
    lB = frame_LLs;
    lalpha=zeros(NT,n_shifts);
    lbeta = zeros(NT,n_shifts);
    lscale=zeros(NT,1); %initialize rescaling parameters
    %compute rescaled forward messages
    lalpha(1,:) = leps_prior + lB(1,:);
    lscale(1)=logsumexp(lalpha(1,:));
    lalpha(1,:) = lalpha(1,:) - lscale(1);
    for t=2:NT
        if use_prior(t)
            cur_lA = lA_tflip;
        else
            cur_lA = lA;
        end
        lalpha(t,:) = logmulexp(lalpha(t-1,:),cur_lA) + lB(t,:);
        lscale(t) = logsumexp(lalpha(t,:));
        lalpha(t,:)= lalpha(t,:) - lscale(t);
    end
    
    %compute rescaled backward messages
    lbeta(NT,:)=log(ones(1,n_shifts)) - lscale(NT);
    for t=NT-1:-1:1
        if use_prior(t+1)
            cur_lA = lA_tflip;
        else
            cur_lA = lA;
        end
        lf1 = lbeta(t+1,:) + lB(t+1,:);
        lbeta(t,:) = logmulexp(lf1,cur_lA') - lscale(t);
    end
    lgamma= lalpha + lbeta;
    lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
    
    gamma = exp(lgamma);
    fpost_mean = sum(bsxfun(@times,gamma,shifts),2);
    fpost_std = sqrt(sum(bsxfun(@times,gamma,shifts.^2),2) - fpost_mean.^2);
    
    %% RECONSTRUCT MAP STIMULUS (OLD METHOD)
    [max_post,max_loc] = max(lgamma,[],2);
    
    fshift_cor = shifts(max_loc);
    for ii = 1:NT
        %         d2 = dist_shift2d(resh_X(:,:,ii), -shift_cor(ii,XV),2,0);
        d2_left = shift_matrix_Nd(resh_X(:,tot_left_eye_inds,ii), -fshift_cor(ii),2);
        d2_right = shift_matrix_Nd(resh_X(:,~tot_left_eye_inds,ii), -fshift_cor(ii),2);
        resh_X_sh(:,:,ii) = cat(2,d2_left,d2_right);
    end
    X_sh = reshape(resh_X_sh,size(Xmat,2),NT)';
    X_sh = X_sh(:,use_k);
    
    %% REFIT MODELS
    cur_fit_set = setdiff(tr_set,poss_xv_set);
    for cc = 1:length(cur_fit_set)
        fprintf('Refitting unit %d of %d in tr-only set\n',cc,length(cur_fit_set));
        cur_unit = cur_fit_set(cc);
        %adjust regularization of new filter
        ref_mod{1}(cur_unit) = init_mod(cur_unit);
        ref_mod{1}(cur_unit) = NIMfit_filters(ref_mod{1}(cur_unit),all_use_binned_spks(cur_unit,:),X_sh,Xexpt(tr_inds,:),[],silent,optim_params,full_L2_mats); %fit stimulus filters
        
        ref_LL_imp(1,cur_unit) = (ref_mod{1}(cur_unit).LL_seq(end)*n_spks(cur_unit) - null_LL(cur_unit))/log(2);
        
        fprintf('Original: %.4f  New: %.4f\n',init_mod_LL_imp(cur_unit),ref_LL_imp(1,cur_unit));
    end
else
    fshift_cor = mean(shift_cor);
    %jacknife error est
    fpost_std = (length(poss_xv_set)-1)*mean(var(post_std,[],2));
    fpost_mean = mean(post_mean,2);
end


%% PLOT IMPROVEMENT
mdl = LinearModel.fit(init_mod_LL_imp(tr_set)'./n_spks(tr_set)',ref_LL_imp(1,tr_set)'./n_spks(tr_set)');
xx = linspace(0,0.2,100);
[ypred,pred_errs] = predict(mdl,xx');
figure;hold on
plot(xx,ypred,'k','linewidth',2);plot(xx,pred_errs,'r--')
plot(init_mod_LL_imp(tr_set)./n_spks(tr_set),ref_LL_imp(1,tr_set)./n_spks(tr_set),'o')
line([0 0.5],[0 0.5])
% xlim([0 0.5]); ylim([0 0.5])

beta = mdl.Coefficients.Estimate;
fprintf('%.3f-fold improvement on xval\n',beta(2));

%% INITIALIZE FOR ITERATIVE EYE_CORRECTIONS
it_shift_cor{1} = shift_cor;
it_post_mean{1} = post_mean;
it_post_std{1} = post_std;
it_fpost_mean{1} = fpost_mean;
it_fpost_std{1} = fpost_std;
it_fshift_cor{1} = fshift_cor;

%% NOW ITERATE

n_ITER = 3;
for nn = 2:n_ITER+1
    
    fprintf('Eye correction iteration %d of %d\n',nn-1,n_ITER);
    
    %create filter banks for all units
    lin_filt_bank = zeros(length(tr_set),full_klen);
    quad_filt_bank = zeros(length(tr_set),full_klen);
    lin_kerns = zeros(length(tr_set),stim_params.lin_dims);
    mod_offsets = nan(1,length(tr_set));
    for cc = 1:length(tr_set)
        cur_k = [ref_mod{nn-1}(tr_set(cc)).mods(:).filtK];
        lin_filt_bank(cc,:) = cur_k(:,1);
        if use_quad(tr_set(cc))
            quad_filt_bank(cc,:) = cur_k(:,2);
        end
        mod_offsets(cc) = ref_mod{nn-1}(tr_set(cc)).spk_NL_params(1);
        lin_kerns(cc,:) = ref_mod{nn-1}(tr_set(cc)).kLin';
    end
    lin_filt_bank = reshape(lin_filt_bank',[stim_params.stim_dims(1:2) length(tr_set)]);
    quad_filt_bank = reshape(quad_filt_bank',[stim_params.stim_dims(1:2) length(tr_set)]);
    
    lin_Lfilts = lin_filt_bank(:,left_eye_inds,:);
    lin_Rfilts = lin_filt_bank(:,~left_eye_inds,:);
    quad_Lfilts = quad_filt_bank(:,left_eye_inds,:);
    quad_Rfilts = quad_filt_bank(:,~left_eye_inds,:);
    
    %indicator predictions
    expt_out = Xexpt(tr_inds,:)*lin_kerns';
    
    %precompute LL at all shifts for all units
    LLs = nan(NT,length(tr_set),n_shifts);
    for xx = 1:length(shifts)
        fprintf('Shift %d of %d\n',xx,n_shifts);
        cur_lin_shift_L = shift_matrix_Nd(lin_Lfilts,shifts(xx),2);
        cur_lin_shift_R = shift_matrix_Nd(lin_Rfilts,shifts(xx),2);
        cur_lin_shift = cat(2,cur_lin_shift_L,cur_lin_shift_R);
        cur_lin_shift = reshape(cur_lin_shift,full_klen,n_tr_chs);
        
        cur_quad_shift_L = shift_matrix_Nd(quad_Lfilts,shifts(xx),2);
        cur_quad_shift_R = shift_matrix_Nd(quad_Rfilts,shifts(xx),2);
        cur_quad_shift = cat(2,cur_quad_shift_L,cur_quad_shift_R);
        cur_quad_shift = reshape(cur_quad_shift,full_klen,n_tr_chs);
        
        %outputs of stimulus models at current X-matrix shift
        lin_out = use_Xmat*cur_lin_shift;
        quad_out = (use_Xmat*cur_quad_shift).^2;
        gfun = lin_out + quad_out;
        gfun = bsxfun(@plus,gfun,mod_offsets);
        
        %add contributions from extra lin kernels
        gfun = gfun + expt_out;
        
        too_large = gfun > 50;
        pred_rate = log(1+exp(gfun));
        pred_rate(too_large) = gfun(too_large);
        
        pred_rate(pred_rate < 1e-20) = 1e-20;
        
    LLs(:,:,xx) = all_use_binned_spks(tr_set,:)'.*log(pred_rate) - pred_rate;
    end
    
    for XV = 1:length(tr_set)
        xv_cell = tr_set(XV);
        
        cur_tr_inds = find(tr_set~=xv_cell);
        
        fprintf('Using XV cell %d,  %d of %d\n',xv_cell,XV,length(tr_set));
        frame_LLs = squeeze(sum(LLs(:,cur_tr_inds,:),2));
        
        %%
        lB = frame_LLs;
        lalpha=zeros(NT,n_shifts);
        lbeta = zeros(NT,n_shifts);
        lscale=zeros(NT,1); %initialize rescaling parameters
        %compute rescaled forward messages
        lalpha(1,:) = leps_prior + lB(1,:);
        lscale(1)=logsumexp(lalpha(1,:));
        lalpha(1,:) = lalpha(1,:) - lscale(1);
        for t=2:NT
            if use_prior(t)
                cur_lA = lA_tflip;
            else
                cur_lA = lA;
            end
            lalpha(t,:) = logmulexp(lalpha(t-1,:),cur_lA) + lB(t,:);
            lscale(t) = logsumexp(lalpha(t,:));
            lalpha(t,:)= lalpha(t,:) - lscale(t);
        end
        
        %compute rescaled backward messages
        lbeta(NT,:)=log(ones(1,n_shifts)) - lscale(NT);
        for t=NT-1:-1:1
            if use_prior(t+1)
                cur_lA = lA_tflip;
            else
                cur_lA = lA;
            end
            lf1 = lbeta(t+1,:) + lB(t+1,:);
            lbeta(t,:) = logmulexp(lf1,cur_lA') - lscale(t);
        end
        lgamma= lalpha + lbeta;
        lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
        
        gamma = exp(lgamma);
        post_mean(:,XV) = sum(bsxfun(@times,gamma,shifts),2);
        post_std(:,XV) = sqrt(sum(bsxfun(@times,gamma,shifts.^2),2) - post_mean(:,XV).^2);
        
        %% RECONSTRUCT MAP STIMULUS (OLD METHOD)
        [max_post,max_loc] = max(lgamma,[],2);
        
        shift_cor(:,XV) = shifts(max_loc);
        for ii = 1:NT
            %         d2 = dist_shift2d(resh_X(:,:,ii), -shift_cor(ii,XV),2,0);
            d2_left = shift_matrix_Nd(resh_X(:,left_eye_inds,ii), -shift_cor(ii,XV),2);
            d2_right = shift_matrix_Nd(resh_X(:,~left_eye_inds,ii), -shift_cor(ii,XV),2);
            resh_X_sh(:,:,ii) = cat(2,d2_left,d2_right);
        end
        X_sh = reshape(resh_X_sh,full_klen,NT)';
        
        %% REFIT MODELS
        %adjust regularization of new filter
        ref_mod{nn}(xv_cell) = ref_mod{nn-1}(xv_cell);
        %     ref_mod{1}(xv_cell) = NIMadjust_regularization(ref_mod{1}(xv_cell),[1 2],'lambda_L1',0);
        ref_mod{nn}(xv_cell) = NIMfit_filters(ref_mod{nn}(xv_cell),all_use_binned_spks(xv_cell,:),X_sh,Xexpt(tr_inds,:),[],silent,optim_params,full_L2_mats); %fit stimulus filters
        
        ref_LL_imp(nn,xv_cell) = (ref_mod{nn}(xv_cell).LL_seq(end)*n_spks(xv_cell) - null_LL(xv_cell))/log(2);
        
        fprintf('Original: %.4f  Previous: %.4f  New: %.4f\n',init_mod_LL_imp(xv_cell),ref_LL_imp(nn-1,xv_cell),ref_LL_imp(nn,xv_cell));
        
    end
    
    %% if there are any units in the training set that weren't allowed to be XV units, then refit those, using a all-cell eye position
    if any(~ismember(tr_set,poss_xv_set))
        
        fprintf('Using XV cell %d,  %d of %d\n',xv_cell,XV,length(poss_xv_set));
        frame_LLs = squeeze(sum(LLs,2));
        
        %%
        lB = frame_LLs;
        lalpha=zeros(NT,n_shifts);
        lbeta = zeros(NT,n_shifts);
        lscale=zeros(NT,1); %initialize rescaling parameters
        %compute rescaled forward messages
        lalpha(1,:) = leps_prior + lB(1,:);
        lscale(1)=logsumexp(lalpha(1,:));
        lalpha(1,:) = lalpha(1,:) - lscale(1);
        for t=2:NT
            if use_prior(t)
                cur_lA = lA_tflip;
            else
                cur_lA = lA;
            end
            lalpha(t,:) = logmulexp(lalpha(t-1,:),cur_lA) + lB(t,:);
            lscale(t) = logsumexp(lalpha(t,:));
            lalpha(t,:)= lalpha(t,:) - lscale(t);
        end
        
        %compute rescaled backward messages
        lbeta(NT,:)=log(ones(1,n_shifts)) - lscale(NT);
        for t=NT-1:-1:1
            if use_prior(t+1)
                cur_lA = lA_tflip;
            else
                cur_lA = lA;
            end
            lf1 = lbeta(t+1,:) + lB(t+1,:);
            lbeta(t,:) = logmulexp(lf1,cur_lA') - lscale(t);
        end
        lgamma= lalpha + lbeta;
        lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
        
        gamma = exp(lgamma);
        fpost_mean = sum(bsxfun(@times,gamma,shifts),2);
        fpost_std = sqrt(sum(bsxfun(@times,gamma,shifts.^2),2) - fpost_mean.^2);
        
        %% RECONSTRUCT MAP STIMULUS (OLD METHOD)
        [max_post,max_loc] = max(lgamma,[],2);
        
        fshift_cor = shifts(max_loc);
        for ii = 1:NT
            %         d2 = dist_shift2d(resh_X(:,:,ii), -shift_cor(ii,XV),2,0);
            d2_left = shift_matrix_Nd(resh_X(:,left_eye_inds,ii), -fshift_cor(ii),2);
            d2_right = shift_matrix_Nd(resh_X(:,~left_eye_inds,ii), -fshift_cor(ii),2);
            resh_X_sh(:,:,ii) = cat(2,d2_left,d2_right);
        end
        X_sh = reshape(resh_X_sh,full_klen,NT)';
        
        %% REFIT MODELS
        cur_fit_set = setdiff(tr_set,poss_xv_set);
        for cc = 1:length(cur_fit_set)
            fprintf('Refitting unit %d of %d in tr-only set\n',cc,length(cur_fit_set));
            cur_unit = cur_fit_set(cc);
            %adjust regularization of new filter
            ref_mod{nn}(cur_unit) = ref_mod{nn-1}(cur_unit);
            ref_mod{nn}(cur_unit) = NIMfit_filters(ref_mod{nn}(cur_unit),aall_use_binned_spks(cur_unit,:),X_sh,Xexpt(tr_inds,:),[],silent,optim_params,full_L2_mats); %fit stimulus filters
            ref_LL_imp(nn,cur_unit) = (ref_mod{nn}(cur_unit).LL_seq(end)*n_spks(cur_unit) - null_LL(cur_unit))/log(2);
            fprintf('Original: %.4f  Previous: %.4f  New: %.4f\n',init_mod_LL_imp(cur_unit),ref_LL_imp(nn-1,cur_unit),ref_LL_imp(nn,cur_unit));
        end
    else
        fshift_cor = mean(shift_cor);
        %jacknife error est
        fpost_std = (length(poss_xv_set)-1)*mean(var(post_std,[],2));
        fpost_mean = mean(post_mean,2);
    end
    
    %%
    it_post_mean{nn} = post_mean;
    it_post_std{nn} = post_std;
    it_shift_cor{nn} = shift_cor;
    it_fpost_mean{nn} = fpost_mean;
    it_fpost_std{nn} = fpost_std;
    it_fshift_cor{nn} = fshift_cor;
    
    %% PLOT PROGRESS
    mdl = LinearModel.fit(init_mod_LL_imp(tr_set)'./n_spks(tr_set)',ref_LL_imp(nn,tr_set)'./n_spks(tr_set)');
    xx = linspace(0,0.2,100);
    [ypred,pred_errs] = predict(mdl,xx');
    figure;hold on
    plot(xx,ypred,'k','linewidth',2);plot(xx,pred_errs,'r--')
    plot(init_mod_LL_imp(tr_set)./n_spks(tr_set),ref_LL_imp(nn,tr_set)./n_spks(tr_set),'o')
    plot(init_mod_LL_imp(tr_set)./n_spks(tr_set),ref_LL_imp(nn-1,tr_set)./n_spks(tr_set),'r*')
    line([0 0.5],[0 0.5])
    % xlim([0 0.5]); ylim([0 0.5])
    
    beta = mdl.Coefficients.Estimate;
    fprintf('%.3f-fold improvement on xval\n',beta(2));
    
    mdl = LinearModel.fit(ref_LL_imp(nn-1,tr_set)'./n_spks(tr_set)',ref_LL_imp(nn,tr_set)'./n_spks(tr_set)');
    beta = mdl.Coefficients.Estimate;
    fprintf('%.3f-fold previous improvement on xval\n',beta(2));
    
    %     xx = linspace(0,0.2,100);
    %     [ypred,pred_errs] = predict(mdl,xx');
    %     subplot(2,1,2)
    %     hold on
    %     plot(xx,ypred,'k','linewidth',2);plot(xx,pred_errs,'r--')
    %     plot(prev_imp,new_imp,'o')
    %     plot(prev_imp_nonxv,new_imp_nonxv,'k*')
    %     line([0 0.2],[0 0.2])
    %     xlim([0 0.25]); ylim([0 0.25])
    
    
    %     %% SAVE PROGRESS
    %     eye_times = all_stim_times(tr_inds);
    %     save full_eye_correct_v6_90deg *_LL it* tr_set non_xv_cells tr_inds eye_times
    
end


%%
measured_pre_Lx = [saccades(:).pre_Lx];
measured_post_Lx = [saccades(:).post_Lx];
measured_delta_pos = measured_post_Lx(used_saccade_set) - measured_pre_Lx(used_saccade_set);

for ii = 1:length(it_fpost_mean)
    inferred_pos = mean(it_fpost_mean{ii},2);
    inferred_pre_pos = inferred_pos(used_saccade_inds-1);
    inferred_post_pos = inferred_pos(used_saccade_inds + 6);
    inferred_delta_pos = (inferred_post_pos - inferred_pre_pos)*sp_dx;
    
    sac_delta_corr(ii) =  corr(measured_delta_pos(:),inferred_delta_pos(:),'type','spearman');
end