clear all
% close all
cd ~/Data/bruce/G081/
load jbeG081Expts.mat
addpath('~/James_scripts/bruce/G081/');
%%
stim_fs = 100; %in Hz
use_sus = 1:96;

%% USE ONLY GRAY BACKGROUND DATA
flen = 12;
beg_buffer = round(stim_fs*0.15); %don't use the first X data after start of a trial.

cur_expt_set = [64 65];

%% COMPUTE TRIAL DATA
all_stim_times = [];
all_rel_stimes = [];
all_rel_etimes = [];
all_phase = [];
all_Or = [];
all_used_inds = [];
all_binned_spks = [];
all_exptvec = [];
all_trialvec = [];
all_trial_start_times = [];
all_trial_end_times = [];
for ee = 1:length(cur_expt_set);
    fprintf('Expt %d of %d\n',ee,length(cur_expt_set));
    cur_expt = cur_expt_set(ee);
    fname = sprintf('Expt%dClusterTimes.mat',cur_expt);
    load(fname);
    
    n_trials = length(Expts{cur_expt}.Trials);
    trial_start_times = [Expts{cur_expt}.Trials(:).TrialStart]/1e4;
    trial_end_times = [Expts{cur_expt}.Trials(:).TrueEnd]/1e4;
    trial_durs = [Expts{cur_expt}.Trials(:).dur]/1e4;
    
    all_trial_start_times = [all_trial_start_times; trial_start_times(:)];
    all_trial_end_times = [all_trial_end_times; trial_end_times(:)];
    for tt = 1:n_trials
        cur_stim_times = Expts{cur_expt}.Trials(tt).Start/1e4;
        cur_Or = Expts{cur_expt}.Trials(tt).or;
        cur_phase = Expts{cur_expt}.Trials(tt).ph;
        
        cur_t_edges = [cur_stim_times; Expts{cur_expt}.Trials(tt).End(end)/1e4];
        cur_binned_spks = nan(length(cur_stim_times),length(use_sus));
        for cc = 1:length(use_sus)
            cur_hist = histc(Clusters{use_sus(cc)}.times,cur_t_edges);
            cur_binned_spks(:,cc) = cur_hist(1:end-1);
        end
        
        cur_used_inds = ones(length(cur_stim_times),1);
        cur_used_inds(1:flen) = 0;
        cur_used_inds(1:beg_buffer) = 0;
        
        all_stim_times = [all_stim_times; cur_stim_times];
        all_rel_stimes = [all_rel_stimes; cur_stim_times- trial_start_times(tt)];
        all_rel_etimes = [all_rel_etimes; trial_end_times(tt) - cur_stim_times];
        all_Or = [all_Or; cur_Or];
        all_phase = [all_phase; cur_phase'];
        all_used_inds = [all_used_inds; cur_used_inds];
        all_binned_spks = [all_binned_spks; cur_binned_spks];
        all_exptvec = [all_exptvec; ones(size(cur_stim_times))*ee];
        all_trialvec = [all_trialvec; ones(size(cur_stim_times))*tt];
    end
    
end
%%
un_bar_oris = unique(all_Or);
un_bar_oris(1) = [];
n_bar_oris = length(un_bar_oris);

cur_bar_mat = zeros(length(all_stim_times),n_bar_oris);
for bb = 1:n_bar_oris
    cur_set = find(all_Or==un_bar_oris(bb));
    cur_bar_mat(cur_set,bb) = 1;
end

all_bar_mat = makeStimRows(cur_bar_mat,flen);

%% IF YOU WANT TO DO CROSS-VALIDATION ON INITIAL MODEL FITS, OTHERWISE,
[c,ia,ic] = unique([all_exptvec all_trialvec],'rows');
n_trials = length(ia);

xv_frac = 0.2;
xv_set = randperm(n_trials);
xv_set(round(n_trials*xv_frac)+1:end) = [];

tr_set = setdiff(1:n_trials,xv_set);

tr_inds = find(ismember(ic,tr_set));
xv_inds = find(ismember(ic,xv_set));
tr_inds(all_used_inds(tr_inds) == 0) = [];
xv_inds(all_used_inds(xv_inds) == 0) = [];

tr_Xmat = all_bar_mat(tr_inds,:);
xv_Xmat = all_bar_mat(xv_inds,:);

%%
% SDIM = length(un_bar_oris);
% for cc = 1:96
%     fprintf('Fitting cell %d of %d\n',cc,96);
%     Robs_tr = all_binned_spks(tr_inds,cc);
%     tr_spkbns = convert_to_spikebins(all_binned_spks(tr_inds,cc));
%     
%     %compute STA
%     sta = mean(tr_Xmat(tr_spkbns,:)) - mean(tr_Xmat);
%     
%     %don't project out STA
%     proj_mat = sta'*inv(sta*sta')*sta;
%     stim_proj = tr_Xmat - tr_Xmat*proj_mat;
%     %     stim_proj = cur_tr_stimemb;
%     stvcv = cov(stim_proj(tr_spkbns,:));
%     utvcv = cov(stim_proj);
%     [evecs,evals] = eig(stvcv-utvcv);
%     evs = diag(evals);
%     STCbvs  = evecs;
%     sta = sta';
%     npos = 4; nneg = 4;
%     stcs_compareset  = evecs(:,[1:nneg,length(evs)-npos+1:end]);
%     stcs_compareset  = stcs_compareset(:,end:-1:1);
%     rstcs = fliplr(stcs_compareset); %reversed STC kernels (suppressive first)
%     used_stcs = [sta stcs_compareset(:,1:npos) rstcs(:,1:nneg)];
%     used_stcs = bsxfun(@rdivide,used_stcs,sqrt(sum(used_stcs.^2)));
%     
%     %FIT NULL MODEL
%     avg_rate = mean(Robs_tr);
%     trpred_rate = ones(1,length(tr_inds))*avg_rate;
%     null_LL(cc) = -sum(Robs_tr.*log(trpred_rate') - trpred_rate')/sum(Robs_tr)
%     
%     subplot(3,4,1)
%     imagesc(reshape(sta,flen,SDIM));
%     cax = caxis(); cm = max(abs(cax));caxis([-cm cm]);
%     subplot(3,4,2)
%     plot(evs(1:30),'.')
%     subplot(3,4,3)
%     plot(evs(end-30:end),'.')
%     for i = 1:8
%         subplot(3,4,4+i)
%         imagesc(reshape(used_stcs(:,i+1),flen,SDIM))
%         cax = caxis(); cm = max(abs(cax));caxis([-cm cm]);
%     end
%     pause
%     clf
% end

%%
SDIM = length(un_bar_oris);
klen = flen*SDIM;
stim_params.spatial_dims = 1;
stim_params.sdim = SDIM;
stim_params.flen = flen;

for cc = 1:96
    fprintf('Fitting cell %d of %d\n',cc,96);
    tr_spkbns = convert_to_spikebins(all_binned_spks(tr_inds,cc));
    
    clear defmod
    defmod.lambda_L1x = 20;
    defmod.lambda_L2x = 100;
   defmod.lambda_d2XT = 25;
    defmod.lambda_d2X = 0;
    defmod.lambda_d2T = 100;
    clear kern_types
    kern_types{1} = 'lin';
    init_kerns = 0.01*randn(klen,1+npq+nnq);
    init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
    quad = createGNM(init_kerns,[1 ones(1,npq) -1*ones(1,nnq)],kern_types,defmod,stim_params);
    ori_fit(cc) = fitGNM_filters(quad,tr_Xmat,tr_spkbns,'none',[],1e-4,1e-6,1);
    
    [ori_LL(cc)] = getLL_GNM(ori_fit(cc),tr_Xmat,tr_spkbns,'none');

    Robs = all_binned_spks(tr_inds,cc);
    avg_rate = mean(Robs);
    null_prate = ones(length(tr_inds),1)*avg_rate;
    null_LL(cc) = -sum(Robs.*log(null_prate) - null_prate)/sum(Robs);
    
    if xv_frac > 0
        xv_Robs = all_binned_spks(xv_inds,use_sus(cc));
        xv_spkbns = convert_to_spikebins(xv_Robs);
        null_prate = ones(length(xv_inds),1)*avg_rate;
        null_xvLL(cc) = -sum(xv_Robs.*log(null_prate) - null_prate)/sum(xv_Robs);
        
        ori_xvLL(cc) = getLL_GNM(ori_fit(cc),xv_Xmat,xv_spkbns,'none');
        fprintf('xv Imp: %.3f\n',null_xvLL(cc)-ori_xvLL(cc));
    end
        fprintf('tr Imp: %.3f\n',null_LL(cc)-ori_LL(cc));
        
        

        cur_k = get_k_mat(ori_fit(cc));
        cur_k = reshape(cur_k,flen,SDIM);
        t_var = var(cur_k,[],2);
        [~,best_slice] = max(t_var);
        ori_tuning(cc,:) = cur_k(best_slice,:);
end


[~,best_ori] = max(ori_tuning,[],2);


%%
save unit_ori_models best_ori ori_tuning ori_* *_xvLL SDIM un_bar_oris

%%
for cc = 1:96
            fprintf('Cell %d  xv Imp: %.3f\n',cc,null_xvLL(cc)-ori_xvLL(cc));

        cur_k = get_k_mat(quad_fit(cc));
        cur_k = reshape(cur_k,flen,SDIM);
        imagesc(un_bar_oris,1:flen,cur_k); 
        mm = max(abs(cur_k(:)));
        caxis([-0.9*mm 0.9*mm])
        pause
        clf
end







