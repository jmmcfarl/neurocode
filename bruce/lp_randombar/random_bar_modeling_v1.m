clear all
close all

cd ~/Data/bruce/2_27_12/M232
% load ./random_bar_eyedata.matxlim([-0. 0.4])
load ./random_bar_eyedata_ftime.mat

load ./lemM232Expts.mat
%%
rf_cent = [4.73 -2.2];
axis_or = 50*pi/180; %in radians

%% COMPUTE PERPINDICULAR EYE POSITION
for ee = 1:length(bar_expts)
    
    avg_eyepos = 0.5*all_eye_vals_cal{ee}(:,1:2) + 0.5*all_eye_vals_cal{ee}(:,3:4);
    
    fix_r = sqrt(sum(avg_eyepos.^2,2));
    fix_angle = atan2(avg_eyepos(:,2),avg_eyepos(:,1));
    angle_diff = fix_angle - axis_or;
    fix_perp{ee} = fix_r.*sin(angle_diff);
    
end

%%
load ./CellList.mat
good_sus = find(all(CellList(bar_expts,:,1) > 0));
good_cellnum = CellList(bar_expts(1),good_sus,1);
% use_sus = find(all(CellList(bar_expts,:,1) > 0));
use_sus = 1:24;
% used_expts{1} = 1:6;
% used_expts{2} = 1:5;
% used_expts{3} = 1:6;
% used_expts{4} = 
%%
load ./un_bar_pos.mat
un_bar_pos(1:2) = [];

flen = 20;
full_Xmat = [];
full_spkbinned = [];
full_exptvec = [];
full_trialvec = [];
full_taxis = [];
full_bar_pos = [];
full_bar_Xmat = [];
for ee = 1:length(bar_expts);
    fprintf('Processing expt %d of %d\n',ee,length(bar_expts));
    fname = sprintf('Expt%dClusterTimes.mat',bar_expts(ee));
    load(fname);
    
    trial_durs{ee} = [Expts{bar_expts(ee)}.Trials(:).dur];
    n_trials(ee) = length(Expts{bar_expts(ee)}.Trials);
    all_t_axis = [];
    all_used_inds = [];
    all_expt_vec = [];
    all_trial_vec = [];
    all_bar_Op = [];
    all_bar_e0 = [];
    all_bar_Xmat = [];
    all_binned_spikes = [];
    all_used_inds = [];
    for tt = 1:n_trials(ee)
        
        cur_t_axis = [Expts{bar_expts(ee)}.Trials(tt).Start]/1e4;
        cur_bar_Op = [Expts{bar_expts(ee)}.Trials(tt).Op];
%         cur_bar_e0 = [Expts{bar_expts(ee)}.Trials(tt).e0];
        
        cur_t_edges = [cur_t_axis; Expts{bar_expts(ee)}.Trials(tt).End(end)/1e4];
        cur_binned_spks = nan(length(cur_t_axis),length(use_sus));
        for cc = 1:length(use_sus)
            cur_hist = histc(Clusters{use_sus(cc)}.times,cur_t_edges);
            cur_binned_spks(:,cc) = cur_hist(1:end-1);
        end
        
%         cur_interp_eye_perp = interp1(all_eye_ts{ee},fix_perp{ee},cur_t_axis);
%         cur_bar_Op_cor = cur_bar_Op - cur_interp_eye_perp;
        
        cur_bar_mat = zeros(length(cur_t_axis),length(un_bar_pos));
        for bb = 1:length(un_bar_pos)
            cur_set = find(cur_bar_Op==un_bar_pos(bb));
            cur_bar_mat(cur_set,bb) = 1;
%             cur_bar_mat(cur_set,bb) = cur_bar_e0(cur_set);
        end
        bar_Xmat = makeStimRows(cur_bar_mat,flen);
        
        cur_used_inds = ones(length(cur_t_axis),1);
        cur_used_inds(1:flen) = 0;
        
        all_t_axis = [all_t_axis; cur_t_axis(:)];
        all_used_inds = [all_used_inds; cur_used_inds(:)];
        all_binned_spikes = [all_binned_spikes; cur_binned_spks];
        all_expt_vec = [all_expt_vec; ones(length(cur_t_axis),1)*bar_expts(ee)];
        all_trial_vec = [all_trial_vec; ones(length(cur_t_axis),1)*tt];
        all_bar_Op = [all_bar_Op; cur_bar_Op(:)];
%         all_bar_e0 = [all_bar_e0; cur_bar_e0(:)];
        all_bar_Xmat = [all_bar_Xmat; bar_Xmat];
    end
    
    cur_blink_start_times = all_eye_ts{ee}(all_blink_startinds{ee});
    cur_blink_stop_times = all_eye_ts{ee}(all_blink_stopinds{ee});
    cur_blink_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),cur_blink_start_times));
    cur_blink_stop_inds = round(interp1(all_t_axis,1:length(all_t_axis),cur_blink_stop_times));
    cur_poss_blinks = find(~isnan(cur_blink_start_inds) & ~isnan(cur_blink_stop_inds));
    
    all_blink_inds = zeros(size(all_t_axis));
    for i = 1:length(cur_poss_blinks)
        all_blink_inds(cur_blink_start_inds(cur_poss_blinks(i)):cur_blink_stop_inds(cur_poss_blinks(i))) = 1;
    end
    
    all_used_inds(all_blink_inds == 1) = 0;
    all_used_inds = logical(all_used_inds);    
    
    full_Xmat = [full_Xmat; all_bar_Xmat(all_used_inds,:)];
    full_spkbinned = [full_spkbinned; all_binned_spikes(all_used_inds,:)];
    full_exptvec = [full_exptvec; ones(sum(all_used_inds),1)*ee];
    full_taxis = [full_taxis; all_t_axis(all_used_inds)];
    full_trialvec = [full_trialvec; all_trial_vec(all_used_inds)];
    full_bar_pos = [full_bar_pos; all_bar_Op(all_used_inds)];
end

%%
% load ./random_bar_eyedata.mat
all_fix_inds = [];
all_first_inds = [];
all_second_inds = [];
all_small_inds = [];
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
    cur_fir_sec = [infirst_sac_set; insecond_sac_set];
    micro_sac_set = find(sac_amps' < 60 & ~isnan(all_trial_num(cur_sac_set)));
    micro_sac_set(ismember(micro_sac_set,cur_fir_sec)) = [];
    
%     fix_times = sac_stop_times(intrial_sac_set);
%     first_fix_times = sac_stop_times(infirst_sac_set);
%     second_fix_times = sac_stop_times(insecond_sac_set);
    fix_times = sac_start_times(intrial_sac_set);
    first_fix_times = sac_start_times(infirst_sac_set);
    second_fix_times = sac_start_times(insecond_sac_set);
    small_fix_times = sac_start_times(micro_sac_set);
    
    cur_e_set = find(full_exptvec == ee);
    
    fix_inds = round(interp1(full_taxis(cur_e_set),1:length(cur_e_set),fix_times));
    first_fix_inds = round(interp1(full_taxis(cur_e_set),1:length(cur_e_set),first_fix_times));
    second_fix_inds = round(interp1(full_taxis(cur_e_set),1:length(cur_e_set),second_fix_times));
    small_fix_inds = round(interp1(full_taxis(cur_e_set),1:length(cur_e_set),small_fix_times));

    fix_inds(isnan(fix_inds)) = [];
    first_fix_inds(isnan(first_fix_inds)) = [];
    second_fix_inds(isnan(second_fix_inds)) = [];
    small_fix_inds(isnan(small_fix_inds)) = [];
    
    all_fix_inds = [all_fix_inds; cur_e_set(fix_inds(:))];
    all_first_inds = [all_first_inds; cur_e_set(first_fix_inds(:))];
    all_second_inds = [all_second_inds; cur_e_set(second_fix_inds(:))];
    all_small_inds = [all_small_inds; cur_e_set(small_fix_inds(:))];
end

used_sac_inds = sort([all_first_inds; all_second_inds]);
% used_sac_inds = sort(all_fix_inds);
% small_set = setdiff(all_fix_inds,used_sac_inds);

dt = 0.01;
cur_dt = 0.01;
% flen_t = 0.5;
tent_centers = [0:cur_dt:0.8];
tent_centers = round(tent_centers/dt);
tbmat = construct_tent_bases(tent_centers,1);
[ntents,tblen] = size(tbmat);

shift = round(0.3/dt);
tbmat = [zeros(ntents,tblen-2*shift-1) tbmat];
tent_centers = tent_centers-shift;

trial_inds = zeros(size(full_taxis));
trial_inds(used_sac_inds) = 1;
trial_Tmat = zeros(length(full_taxis),ntents);
for i = 1:ntents
    trial_Tmat(:,i) = conv(trial_inds,tbmat(i,:),'same');
end

% trial_inds = zeros(size(full_taxis));
% trial_inds(all_small_inds) = 1;
% Tmat1 = zeros(length(full_taxis),ntents);
% for i = 1:ntents
%     Tmat1(:,i) = conv(trial_inds,tbmat(i,:),'same');
% end
% trial_inds = zeros(size(full_taxis));
% trial_inds(all_second_inds) = 1;
% Tmat2 = zeros(length(full_taxis),ntents);
% for i = 1:ntents
%     Tmat2(:,i) = conv(trial_inds,tbmat(i,:),'same');
% end

trial_inds = zeros(size(full_taxis));
trial_inds(all_first_inds) = 1;
Tmat1 = zeros(length(full_taxis),ntents);
for i = 1:ntents
    Tmat1(:,i) = conv(trial_inds,tbmat(i,:),'same');
end
trial_inds = zeros(size(full_taxis));
trial_inds(all_second_inds) = 1;
Tmat2 = zeros(length(full_taxis),ntents);
for i = 1:ntents
    Tmat2(:,i) = conv(trial_inds,tbmat(i,:),'same');
end

lin_X = zeros(length(full_taxis),length(bar_expts)-1);
for i = 1:length(bar_expts)-1
    cur_set = find(full_exptvec==i);
    lin_X(cur_set,i) = 1;
end

%%
maxlag = round(0.5/dt);
lags = -maxlag:maxlag;
rate_sm =1;
sac_trig_avgs = nan(24,length(lags));
fsac_trig_avgs = nan(24,length(lags));
ssac_trig_avgs = nan(24,length(lags));
msac_trig_avgs = nan(24,length(lags));
cur_used_sacs = used_sac_inds;
cur_used_sacs(cur_used_sacs < maxlag | cur_used_sacs > length(full_taxis) - maxlag) = [];
cur_fused_sacs = all_first_inds;
cur_fused_sacs(cur_fused_sacs < maxlag | cur_fused_sacs > length(full_taxis) - maxlag) = [];
cur_sused_sacs = all_second_inds;
cur_sused_sacs(cur_sused_sacs < maxlag | cur_sused_sacs > length(full_taxis) - maxlag) = [];
cur_used_msacs = all_small_inds;
cur_used_msacs(cur_used_msacs < maxlag | cur_used_msacs > length(full_taxis) - maxlag) = [];
for cc = 1:24
    cur_sm_rate = jmm_smooth_1d_cor(full_spkbinned(:,cc),rate_sm);
    sac_trig_avgs(cc,:) = get_event_trig_avg(cur_sm_rate,cur_used_sacs,maxlag,maxlag);
    fsac_trig_avgs(cc,:) = get_event_trig_avg(cur_sm_rate,cur_fused_sacs,maxlag,maxlag);
    ssac_trig_avgs(cc,:) = get_event_trig_avg(cur_sm_rate,cur_sused_sacs,maxlag,maxlag);
    msac_trig_avgs(cc,:) = get_event_trig_avg(cur_sm_rate,cur_used_msacs,maxlag,maxlag);
end

%%
[c,ia,ic] = unique([full_exptvec full_trialvec],'rows');
n_trials = length(ia);

xv_frac = 0.2;
n_xv_trials = round(n_trials*xv_frac);
xv_set = randperm(n_trials);
xv_set(n_xv_trials+1:end) = [];
xv_inds = find(ismember(ic,xv_set));
tr_inds = find(~ismember(ic,xv_inds))';


% Xmat = bsxfun(@minus,full_Xmat,mean(full_Xmat));
Xmat = full_Xmat;
%%
load ./random_bar_unit_models_ftime glm_fit

n_bar_pos = length(un_bar_pos);
stim_params.spatial_dims = 1;
stim_params.sdim = length(un_bar_pos);
stim_params.flen = flen;
klen = length(un_bar_pos)*flen;
clear defmod
defmod.lambda_L1x = 1;
defmod.lambda_d2XT = 20;
% defmod.lambda_d2X = 100;
defmod.labmda_d2T = 100;

NLtype = 0;
for cc = 1:length(use_sus)
    fprintf('Cell %d of %d\n',cc,length(use_sus));
    
    tr_spkbns = convert_to_spikebins(full_spkbinned(tr_inds,cc));
    cur_sta = mean(Xmat(tr_spkbns,:));
    
    kern_types{1} = 'lin';
    init_kerns = 0.01*randn(klen,1);
    init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
    init_linK = zeros(length(bar_expts)-1,1);
    glm = createGNM_v2(init_kerns,1,kern_types,init_linK,defmod,stim_params);
%     glm.spk_nl = 'exp';
    glm_fit(cc) = fitGNM_filters_v2(glm,Xmat(tr_inds,:),lin_X(tr_inds,:),tr_spkbns,'none',[],1e-4,1e-6);
    
    
    kern_types{1} = 'threshlin';
    kern_types{2} = 'threshlin';
    init_kerns = 0.01*randn(klen,2);
    init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
    init_linK = zeros(length(bar_expts)-1,1);
    gnm = createGNM_v2(init_kerns,[1 -1],kern_types,init_linK,defmod,stim_params);
%     glm.spk_nl = 'exp';
    gnm_fit(cc) = fitGNM_filters_v2(gnm,Xmat(tr_inds,:),lin_X(tr_inds,:),tr_spkbns,'none',[],1e-4,1e-6);
    
    
%     [~, ~, ~, ~, g] = getLL_GNM(glm_fit(cc),Xmat,cur_spkbns,'none');
%     fitGNM_spkNL(glm_fit(cc),g,cur_spkbns,1,[]);
    
    glm_kern = get_k_mat(glm_fit(cc));
    stim_out = Xmat*glm_kern;
%     
%     l2_pen = 50;
%     l1_pen = 0.1;
%     NLtype = 0;
% %    cur_Xmat = [trial_Tmat stim_out lin_X];
%    cur_Xmat = [trial_Tmat stim_out];
%     lamrange2 = [l2_pen 1 ntents 0];
%      llist = [l1_pen 1:ntents];
%    K0 = zeros(size(cur_Xmat,2)+1,1);
%     [fitp,grad] = GLMsolve_jmm(cur_Xmat(tr_inds,:), tr_spkbns, K0, 1, [], lamrange2,[], [], llist, [], NLtype);
%     sac_kern_ind(cc,:) = fitp.k(1:ntents);
% 
%     l2_pen = 50;
%     l1_pen = 0.1;
%     NLtype = 0;
% %     cur_Xmat2 = [Tmat1 stim_out lin_X];
%     cur_Xmat2 = [Tmat1 stim_out];
%     %     lamrange2 = [l2_pen 1 ntents 0; l2_pen ntents+1 2*ntents 0];
%     lamrange2 = [l2_pen 1 ntents 0];
%     llist = [l1_pen 1:ntents];
%     K0 = zeros(size(cur_Xmat2,2)+1,1);
%     [fitp2,grad] = GLMsolve_jmm(cur_Xmat2(tr_inds,:), tr_spkbns, K0, 1, [], lamrange2,[], [], llist, [], NLtype);
%     sac1_kern_ind(cc,:) = fitp2.k(1:ntents);
% 
%     l2_pen = 50;
%     l1_pen = 0.1;
%     NLtype = 0;
% %     cur_Xmat2 = [Tmat1 stim_out lin_X];
%     cur_Xmat2 = [Tmat2 stim_out];
%     %     lamrange2 = [l2_pen 1 ntents 0; l2_pen ntents+1 2*ntents 0];
%     lamrange2 = [l2_pen 1 ntents 0];
%     llist = [l1_pen 1:ntents];
%     K0 = zeros(size(cur_Xmat2,2)+1,1);
%     [fitp2,grad] = GLMsolve_jmm(cur_Xmat2(tr_inds,:), tr_spkbns, K0, 1, [], lamrange2,[], [], llist, [], NLtype);
%     sac2_kern_ind(cc,:) = fitp2.k(1:ntents);
%     %     %     sac_kern_ind2(cc,:) = fitp2.k((ntents+1):2*ntents);
% %     
% %     l2_pen = 400;
% %     l1_pen = 0.1;
% %     NLtype = 0;
% %     cur_Xmatg = [trial_Tmat bsxfun(@times,trial_Tmat,stim_out) stim_out lin_X];
% %     lamrange2 = [l2_pen 1 ntents 0; 4*l2_pen ntents+1 2*ntents 0];
% %     llist = [l1_pen 1:2*ntents];
% %     K0 = zeros(size(cur_Xmatg,2)+1,1);
% %     [fitpg,grad] = GLMsolve_jmm(cur_Xmatg(tr_inds,:), tr_spkbns, K0, 1, [], lamrange2,[], [], llist, [], NLtype);
% %     sac_kern_indi(cc,:) = fitpg.k(1:ntents);
% %     sac_kern_indg(cc,:) = fitpg.k((ntents+1):2*ntents);
% % 
%     l2_pen = 50;
%     l1_pen = 0.1;
%     NLtype = 0;
% %     cur_Xmat3 = [trial_Tmat Tmat1 stim_out lin_X];
%     cur_Xmat3 = [Tmat1 Tmat2 stim_out];
%     lamrange2 = [l2_pen 1 ntents 0; l2_pen ntents+1 2*ntents 0];
%     llist = [l1_pen 1:2*ntents];
%     K0 = zeros(size(cur_Xmat3,2)+1,1);
%     [fitp3,grad] = GLMsolve_jmm(cur_Xmat3(tr_inds,:), tr_spkbns, K0, 1, [], lamrange2,[], [], llist, [], NLtype);
%     sac_kern_ind1(cc,:) = fitp3.k(1:ntents);
%     sac_kern_ind2(cc,:) = fitp3.k((ntents+1):2*ntents);
%     
    cur_nMat = lin_X;
    K0 = zeros(size(cur_nMat,2)+1,1);
    [fitp_n,grad] = GLMsolve_jmm(cur_nMat(tr_inds,:), tr_spkbns, K0, 1, [], [],[], [], [], [], NLtype);
    
    xv_Robs = full_spkbinned(xv_inds,cc);
     xv_spkbns = convert_to_spikebins(full_spkbinned(xv_inds,cc));
   
    xv_pred_rate = stim_out(xv_inds) + lin_X(xv_inds,:)*glm_fit(cc).linK+ glm_fit(cc).spk_theta;
    xv_pred_rate = log(1+exp(xv_pred_rate));
%     xv_pred_rate = exp(xv_pred_rate);
% %     

[nll, pnll, lpen, gnm_prate] = getLL_GNM_v2(gnm_fit(cc),full_Xmat(xv_inds,:),lin_X(xv_inds,:),xv_spkbns,'none');

% %     xv_sac_pred_rate = cur_Xmat(xv_inds,:)*fitp.k(1:end-1) + fitp.k(end);
% %     xv_sac_pred_rate = log(1+exp(xv_sac_pred_rate));
% % %     xv_sac_pred_rate = exp(xv_sac_pred_rate);
% % % 
% %     xv_sac_pred_rate2 = cur_Xmat2(xv_inds,:)*fitp2.k(1:end-1) + fitp2.k(end);
% %     xv_sac_pred_rate2 = log(1+exp(xv_sac_pred_rate2));
% % 
% %     xv_sac_pred_rate3 = cur_Xmat3(xv_inds,:)*fitp3.k(1:end-1) + fitp3.k(end);
% %     xv_sac_pred_rate3 = log(1+exp(xv_sac_pred_rate3));
% % % 
% % %         xv_sac_pred_rateg = cur_Xmatg(xv_inds,:)*fitpg.k(1:end-1) + fitpg.k(end);
% % %     xv_sac_pred_rateg = log(1+exp(xv_sac_pred_rateg));
% % 
    xv_glm_LL(cc) = -sum(xv_Robs.*log(xv_pred_rate) - xv_pred_rate)/sum(xv_Robs);
    xv_gnm_LL(cc) = -sum(xv_Robs.*log(gnm_prate) - gnm_prate)/sum(xv_Robs);
% %     xv_sac_LL(cc) = -sum(xv_Robs.*log(xv_sac_pred_rate) - xv_sac_pred_rate)/sum(xv_Robs);
% %     xv_sac2_LL(cc) = -sum(xv_Robs.*log(xv_sac_pred_rate2) - xv_sac_pred_rate2)/sum(xv_Robs);
% %     xv_sac3_LL(cc) = -sum(xv_Robs.*log(xv_sac_pred_rate3) - xv_sac_pred_rate3)/sum(xv_Robs);
% % %     xv_sacg_LL(cc) = -sum(xv_Robs.*log(xv_sac_pred_rateg) - xv_sac_pred_rateg)/sum(xv_Robs);
% %     
% % %     avg_rate = mean(full_spkbinned(tr_inds,cc));
% % %     null_pred = avg_rate*ones(size(xv_Robs));
    null_pred = cur_nMat(xv_inds,:)*fitp_n.k(1:end-1) + fitp_n.k(end);
    null_pred = log(1+exp(null_pred));
    xv_null_LL(cc) = -sum(xv_Robs.*log(null_pred)-null_pred)/sum(xv_Robs);
% 
end

%%
xv_imp = xv_null_LL-xv_LL;
xv_sac_imp = xv_null_LL-xv_sac_LL;
xv_sac2_imp = xv_null_LL-xv_sac2_LL;
xv_sac3_imp = xv_null_LL-xv_sac3_LL;
% xv_sacg_imp = xv_null_LL-xv_sacg_LL;

%%
save random_bar_unit_models_ftime glm_fit xv_* sac_kern* sac_trig_avgs lags msac* dt
%%
close all
fl_ax = (0:-1:-flen)*dt-dt/2;
for cc = 1:24
    subplot(2,1,1)
    imagesc(un_bar_pos,fl_ax,flipud(reshape(get_k_mat(glm_fit(cc)),flen,n_bar_pos)));
    cm = caxis();m = max(abs(cm));caxis([-m m]);
    xlabel('Bar Position (deg)','fontsize',20,'fontname','arial')
    ylabel('Time lag (s)','fontsize',20,'fontname','arial')
    set(gca,'fontsize',16,'fontname','arial')
    subplot(2,1,2)
    plot(tent_centers*dt,sac_kern_ind(cc,:))
    hold on
%     plot(tent_centers*dt,msac_kern_ind(cc,:),'r');
    xlabel('Time since saccade (s)','fontsize',20,'fontname','arial')
    set(gca,'fontsize',16,'fontname','arial')
%     legend('Saccade','Fixational-saccade')
%     plot(tent_centers*dt,sac_kern_ind1(cc,:),'k')
%     plot(tent_centers*dt,sac_kern_ind2(cc,:),'g')
%     plot(tent_centers*dt,sac_kern_ind2(cc,:),'k')
%     xlim([-0.1 0.5])
     xlim([-0.3 0.5])
     box off
     
     fname = sprintf('Stim_sac_unit_ftime%d',cc);
     fillPage(gcf,'papersize',[10 18])
     print(fname,'-dpdf')
     close
   cc
%     pause
%     clf
end

%%
ugood_sus = 1:24;
ugood_cellnum = 1:24;
close all
fl_ax = (0:-1:-flen)*dt-dt/2;
for cc = 1:length(ugood_sus)
    subplot(4,2,[1 2 3 4])
    imagesc(un_bar_pos,fl_ax,flipud(reshape(get_k_mat(glm_fit(ugood_sus(cc))),flen,n_bar_pos)));
    cm = caxis();m = max(abs(cm));caxis([-m m]);
    xlabel('Bar Position (deg)','fontsize',20,'fontname','arial')
    ylabel('Time lag (s)','fontsize',20,'fontname','arial')
    set(gca,'fontsize',16,'fontname','arial')
    subplot(4,2,5)
    plot(tent_centers*dt,sac_kern_ind(ugood_sus(cc),:),'linewidth',1)
    xlabel('Time since saccade (s)','fontsize',20,'fontname','arial')
    ylabel('Saccade kernel (AU)','fontsize',20,'fontname','arial')
    set(gca,'fontsize',16,'fontname','arial')
     xlim([-0.3 0.5])
     box off
     grid
hold on
    subplot(4,2,6)
    plot(tent_centers*dt,msac_kern_ind(ugood_sus(cc),:),'r','linewidth',1)
% plot(tent_centers*dt,sac_kern_ind2(cc,:),'r','linewidth',1);
% plot(tent_centers*dt,msac_kern_ind(cc,:),'k','linewidth',1);
% plot(tent_centers*dt,sac_kern_ind(cc,:),'g','linewidth',1);
    xlabel('Time since saccade (s)','fontsize',20,'fontname','arial')
    ylabel('Saccade kernel (AU)','fontsize',20,'fontname','arial')
    set(gca,'fontsize',16,'fontname','arial')
     xlim([-0.3 0.5])
     box off
     grid on
     subplot(4,2,7)
     plot(lags*dt,sac_trig_avgs(ugood_sus(cc),:)/dt,'b','linewidth',1)
     hold on; 
     plot(lags*dt,fsac_trig_avgs(ugood_sus(cc),:)/dt,'g','linewidth',1);
     plot(lags*dt,ssac_trig_avgs(ugood_sus(cc),:)/dt,'k','linewidth',1);
     axis tight
xlim([-0.3 0.5])
     box off
     grid on
     yl = ylim();
         xlabel('Time since saccade (s)','fontsize',20,'fontname','arial')
         ylabel('Average rate (Hz)','fontsize',20,'fontname','arial')
    set(gca,'fontsize',16,'fontname','arial')
%           plot(lags*dt,msac_trig_avgs(ugood_sus(cc),:)/dt,'r','linewidth',1)
     subplot(4,2,8)
     plot(lags*dt,msac_trig_avgs(ugood_sus(cc),:)/dt,'r','linewidth',1)
xlim([-0.3 0.5])
     grid on
     box off
     
     subplot(4,2,7); yl = ylim(); ylim(yl); subplot(4,2,8); ylim(yl);
         xlabel('Time since saccade (s)','fontsize',20,'fontname','arial')
         ylabel('Average rate (Hz)','fontsize',20,'fontname','arial')
    set(gca,'fontsize',16,'fontname','arial')

    fname = sprintf('Stim_bothsac_unit%d',ugood_cellnum(cc));
     fillPage(gcf,'papersize',[16 20])
     print(fname,'-dpng')
     close
   cc
%     pause
%     clf
end

% for cc = 1:-22
%     cur_k = get_k_mat(gnm_fit(cc));
% 
%     subplot(3,1,1)
%     imagesc(reshape(cur_k(:,1),flen,n_bar_pos));
%     cx = caxis();m = max(abs(cx));caxis([-m m]);
%     subplot(3,1,2)
%     imagesc(reshape(cur_k(:,0),flen,n_bar_pos));
%     cx = caxis();m = max(abs(cx));caxis([-m m]);
%     subplot(3,1,3)
%     imagesc(reshape(cur_k(:,3),flen,n_bar_pos));
%      cx = caxis();m = max(abs(cx));caxis([-m m]);
%    
%    pause
%    clf
% end

%%
close all
for cc = 1:24
    cc
    subplot(2,1,1)
        plot(tent_centers*dt,sac_kern_ind((cc),:),'k','linewidth',1)
         hold on
       plot(tent_centers*dt,sac1_kern_ind((cc),:),'b','linewidth',1)
        plot(tent_centers*dt,sac2_kern_ind((cc),:),'r','linewidth',1)
xlim([-0.3 0.5])
   subplot(2,1,2)     
        plot(tent_centers*dt,sac_kern_ind1((cc),:),'b','linewidth',1)
         hold on
       plot(tent_centers*dt,sac_kern_ind2((cc),:),'r','linewidth',1)
xlim([-0.3 5.5])
        
        pause
        clf
end