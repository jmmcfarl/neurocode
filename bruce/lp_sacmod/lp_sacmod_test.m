clear all
% close all

dir_prefix = '~';
Expt_name = 'M266';

%% SET ANALYSIS PARAMETERS
stim_fs = 100; %in Hz
dt = 0.01;

min_trial_dur = 0.5;
beg_buffer = 0.3;

flen = 12;
use_nPix = 32;
% full_nPix = 36;
stim_params = NIMcreate_stim_params([flen use_nPix],dt);
Fr = 1;

n_probes = 24;

%% LOOP OVER EXPTS

data_dir = [dir_prefix '/Data/bruce/' Expt_name];
cd(data_dir);

if Expt_name(1) == 'G'
    load(sprintf('jbe%sExpts.mat',Expt_name));
elseif Expt_name(1) == 'M'
    load(sprintf('lem%sExpts.mat',Expt_name));
end

save_dir = [dir_prefix '/Analysis/bruce/' Expt_name '/sac_mod'];
if ~exist(save_dir,'dir')
    system(['mkdir ' save_dir]);
end

trial_dur = 4;

load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

%%
include_expts = {'rls.Fa', 'rls.FaXimi','rls.FaXFaXFs'};
for ii = 1:length(Expts)
    if ~isempty(Expts{ii})
        expt_names{ii} = Expts{ii}.Header.expname;
        expt_dds(ii) = Expts{ii}.Stimvals.dd;
        expt_bar_ori(ii) = Expts{ii}.Stimvals.or;
        expt_sac_dir(ii) = mod(Expts{ii}.Stimvals.Fa,180);
        expt_Fr(ii) = Expts{ii}.Stimvals.Fr;
        expt_imback(ii) = isfield(Expts{ii}.Trials,'imi');
        included_type(ii) = any(strcmp(expt_names{ii},include_expts));
    end
end
expt_has_ds(isnan(expt_has_ds)) = 0;
expt_has_ds(expt_has_ds == -1) = 0;
expt_binoc(isnan(expt_binoc)) = 0;
cur_block_set = find(included_type & ~expt_binoc' & expt_Fr == 1);

if strcmp(Expt_name,'M270')
    cur_block_set(cur_block_set == 5) = [];
    cur_block_set(cur_block_set == 19) = []; %this has a different value of pixel size than earlier blocks
end

sim_sac_expts = find(~expt_has_ds(cur_block_set));
imback_gs_expts = find(expt_has_ds(cur_block_set) & expt_imback(cur_block_set)');
grayback_gs_expts = find(expt_has_ds(cur_block_set) & ~expt_imback(cur_block_set)');

%%
if exist('./fin_aclust_data.mat','file')
    load('./fin_aclust_data.mat');
    [n_aclust_expts,n_aclust_probes] = size(autoclust);
else
    disp('No fin_aclust_data found.');
    autoclust = [];
    n_aclust_expts = 0; n_aclust_probes = 0;
end
aclust_probenums = [autoclust(cur_block_set(1),:).probenum];
aclust_CellList = aclust_CellList(cur_block_set,:);
su_probes = find(any(aclust_CellList > 0));
su_ids = unique(aclust_CellList(aclust_CellList > 0));

best_corresp_su_probe = nan(length(su_ids),1);
for ii = 1:length(su_ids)
    cnts = sum(aclust_CellList == su_ids(ii));
    [~,best_corresp_su_probe(ii)] = max(cnts);
end

%% COMPUTE TRIAL DATA
cd(data_dir);

fprintf('Computing prep data\n');
trial_cnt = 0;
cur_toffset = 0;

all_stim_times = [];
all_stim_mat = [];
all_t_axis = [];
all_t_bin_edges = [];
all_tsince_start = [];
all_blockvec = [];
all_trialvec = [];
all_trial_Se = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_bin_edge_pts = [];
all_spk_times = cell(n_probes,1);
all_clust_ids = cell(length(su_probes),1);
trial_toffset = zeros(length(cur_block_set),1);
for ee = 1:length(cur_block_set);
    if ismember(ee,grayback_gs_expts)
        fprintf('Block %d of %d; grayback GS, ori:%d\n',ee,length(cur_block_set),expt_bar_ori(cur_block_set(ee)));
    elseif ismember(ee,imback_gs_expts)
        fprintf('Block %d of %d; imback GS, ori:%d\n',ee,length(cur_block_set),expt_bar_ori(cur_block_set(ee)));
    elseif ismember(ee,sim_sac_expts)
        fprintf('Block %d of %d; SimSac, ori:%d\n',ee,length(cur_block_set),expt_bar_ori(cur_block_set(ee)));
    else
        fprintf('Block %d of %d;  UNMATCHED EXPT TYPE\n',ee,length(cur_block_set));
    end
    cur_block = cur_block_set(ee);
    fname = sprintf('Expt%dClusterTimesDetails.mat',cur_block);
    load(fname);
    for cc = 1:n_probes
        all_spk_times{cc} = cat(1,all_spk_times{cc},ClusterDetails{cc}.t' + cur_toffset);
        cur_ind = find(su_probes == cc);
        if ~isempty(cur_ind)
            aclust_probe = find(aclust_probenums == cc);
            all_clust_ids{cur_ind} = cat(1,all_clust_ids{cur_ind},autoclust(cur_block_set(ee),aclust_probe).idx(:));
        end
    end
    
    trial_start_times = [Expts{cur_block}.Trials(:).TrialStart]/1e4;
    trial_end_times = [Expts{cur_block}.Trials(:).TrueEnd]/1e4;
    trial_durs = [Expts{cur_block}.Trials(:).dur]/1e4;
    trial_ids = [Expts{cur_block}.Trials(:).id];
    
    [un_ids,id_inds] = unique(trial_ids);
    rpt_trials = false;
    if length(un_ids) < length(trial_ids)
        rpt_trials = true;
        fprintf('Warning, repeat trial inds detected!\n');
        use_trials = [];
    else
        use_trials = find(trial_durs >= min_trial_dur);
    end
    %         trial_durs = trial_durs(id_inds);
    %         trial_start_times = trial_start_times(id_inds);
    %         trial_end_times = trial_end_times(id_inds);
    
    
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)' + cur_toffset);
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)' + cur_toffset);
    
    trial_Se = [Expts{cur_block}.Trials(:).se];
    trial_Se = trial_Se(id_inds);
    all_trial_Se = cat(1,all_trial_Se,trial_Se(use_trials)');
    
    fname = sprintf('%s/stims/Expt%d_stim',data_dir,cur_block);
    load(fname);
    buffer_pix = floor((expt_npix(cur_block) - use_nPix)/2);
    cur_use_pix = (1:use_nPix) + buffer_pix;
    
    
    n_trials = length(use_trials);
    for tt = 1:n_trials
        cur_stim_times = Expts{cur_block}.Trials(use_trials(tt)).Start'/1e4;
        n_frames = size(left_stim_mats{use_trials(tt)},1);
        if length(cur_stim_times) == 1
            cur_stim_times = (cur_stim_times:dt*Fr:(cur_stim_times + (n_frames-1)*dt*Fr))';
            cur_t_edges = [cur_stim_times; cur_stim_times(end) + dt*Fr];
        end
        %             cur_t_edges = (trial_start_times(use_trials(tt)):dt:trial_end_times(use_trials(tt)))';
        cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        
        cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
        
        if ~any(isnan(left_stim_mats{use_trials(tt)}(:))) && n_frames > min_trial_dur/dt
            use_frames = min(length(cur_stim_times),n_frames);
            cur_stim_mat = double(left_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix));
            
            all_stim_mat = [all_stim_mat; cur_stim_mat];   
            all_stim_times = [all_stim_times; cur_stim_times' + cur_toffset];
            all_t_axis = [all_t_axis; cur_t_axis + cur_toffset];
            all_t_bin_edges = [all_t_bin_edges; cur_t_edges + cur_toffset];
            all_tsince_start = [all_tsince_start; cur_tsince_start];
            all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
            all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
            all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
            
        end
    end
    trial_cnt = trial_cnt + n_trials;
    trial_toffset(ee) = all_t_bin_edges(end);
    cur_toffset = trial_toffset(ee);
end

%%
% mahal_thresh = su_data.mah_thresh;
all_binned_mua = nan(length(all_t_axis),n_probes);
%for only-MU probes
for cc = 1:n_probes
    [cur_spkhist,cur_bin_inds] = histc(all_spk_times{cc},all_t_bin_edges);
    cur_spkhist(all_bin_edge_pts) = [];
    all_binned_mua(:,cc) = cur_spkhist;
end

%for SU probes
fprintf('Using %d SUs\n',length(su_ids));
all_binned_sua = nan(length(all_t_axis),length(su_ids));
su_used_blocks = false(length(cur_block_set),length(su_ids));
probe_with_su = false(length(cur_block_set),n_probes);
for ss = 1:length(su_ids)
    all_su_spk_times{ss} = [];
    cur_su_probes = find(any(aclust_CellList(:,su_probes) == su_ids(ss)));
    for pp = 1:length(cur_su_probes)
        all_su_inds = find(all_clust_ids{cur_su_probes(pp)} == 1);
        cur_suahist = histc(all_spk_times{su_probes(cur_su_probes(pp))}(all_su_inds),all_t_bin_edges);
        cur_suahist(all_bin_edge_pts) = [];
        
        spk_block_inds = round(interp1(all_t_axis,all_blockvec,all_spk_times{su_probes(cur_su_probes(pp))}));
        for ee = 1:length(cur_block_set);
            cur_block_inds = find(all_blockvec == ee);
            %if SU is isolated in this block
            if aclust_CellList(ee,su_probes(cur_su_probes(pp))) == su_ids(ss)
                su_used_blocks(ee,ss) = true;
                probe_with_su(ee,su_probes(cur_su_probes(pp))) = true;
                all_binned_sua(cur_block_inds,ss) = cur_suahist(cur_block_inds);
                cur_set = all_su_inds(spk_block_inds(all_su_inds) == ee);
                all_su_spk_times{ss} = cat(1,all_su_spk_times{ss},all_spk_times{su_probes(cur_su_probes(pp))}(cur_set));
            end
        end
    end
end
%% DEFINE DATA USED FOR ANALYSIS
used_inds = find(all_tsince_start >= beg_buffer);
all_Xmat = create_time_embedding(all_stim_mat,stim_params);

n_blocks = length(cur_block_set);
Xblock = zeros(length(all_stim_times),n_blocks);
for i = 1:n_blocks
    cur_set = find(all_blockvec==i);
    Xblock(cur_set,i) = 1;
end

%%
init_stim_params = NMMcreate_stim_params([flen use_nPix],dt);
init_stim_params(2) = NMMcreate_stim_params([n_blocks 1],dt);

block_L2 = 1;
null_stim_params = init_stim_params(2:end);
null_reg_params = NMMcreate_reg_params('lambda_L2',[block_L2]);

base_lambda_d2XT = 10;
base_lambda_L1 = 2;

n_squared_filts = 2;
mod_signs = [1 ones(1,n_squared_filts) 1];
NL_types = [{'lin'} repmat({'quad'},1,n_squared_filts) {'lin'}];
init_d2XT = [ones(n_squared_filts+1,1); 0;];
init_L2 = [zeros(n_squared_filts+1,1); block_L2];
init_reg_params = NMMcreate_reg_params('lambda_d2XT',init_d2XT,'lambda_L2',init_L2);
init_Xtargs = [ones(n_squared_filts+1,1); 2];
init_optim_p.optTol = 1e-5; init_optim_p.progTol = 1e-8;
silent = 0;

%% FOR SUA
for cc = 1:length(su_ids)
% cc = 1

fprintf('Computing trig avgs for SUA %d of %d\n',cc,length(su_ids));
cur_used_blocks = find(su_used_blocks(:,cc)); %blocks when NO SU
cur_poss_inds = find(ismember(all_blockvec,cur_used_blocks));
cur_used_inds = used_inds(ismember(used_inds,cur_poss_inds));

if ~isempty(cur_used_inds)
    
    Robs = all_binned_sua(cur_used_inds,cc);
    cur_X{1} = all_Xmat(cur_used_inds,:);
    cur_X{2} = Xblock(cur_used_inds,:);
    
    gqm1 = NMMinitialize_model(init_stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs);
    gqm1 = NMMfit_filters(gqm1,Robs,cur_X,[],[],silent,init_optim_p);
    [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(gqm1,Robs,cur_X);
    gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_d2XT',base_lambda_d2XT./var(gint)');
    gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_L1',base_lambda_L1./std(gint)');
    gqm1 = NMMfit_filters(gqm1,Robs,cur_X,[],[],silent);
    
    null_mod = NMMinitialize_model(null_stim_params,[1],{'lin'},null_reg_params,[1]);
    null_mod = NMMfit_filters(null_mod,Robs,cur_X(2:end));
    
    su_qmods(cc) = gqm1;
    su_qmod_LLimp(cc) = gqm1.LL_seq(end) - null_mod.LL_seq(end);
end
end

%% FOR MUA
for cc = 1:24
    
    fprintf('Computing trig avgs for MUA %d of %d\n',cc,24);
    cur_used_blocks = find(~probe_with_su(:,cc));
    cur_poss_inds = find(ismember(all_blockvec,cur_used_blocks));
    cur_used_inds = used_inds(ismember(used_inds,cur_poss_inds));
    
    if ~isempty(cur_used_inds)
        
        Robs = all_binned_mua(cur_used_inds,cc);
        cur_X{1} = all_Xmat(cur_used_inds,:);
        cur_X{2} = Xblock(cur_used_inds,:);
        
        gqm1 = NMMinitialize_model(init_stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs);
        gqm1 = NMMfit_filters(gqm1,Robs,cur_X,[],[],silent,init_optim_p);
        [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(gqm1,Robs,cur_X);
        gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_d2XT',base_lambda_d2XT./var(gint)');
        gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_L1',base_lambda_L1./std(gint)');
        gqm1 = NMMfit_filters(gqm1,Robs,cur_X,[],[],silent);
        
%         null_mod = NMMinitialize_model(null_stim_params,[1],{'lin'},null_reg_params,[1]);
%         null_mod = NMMfit_filters(null_mod,Robs,cur_X(2:end));
        
        qmods(cc) = gqm1;
%         qmod_LLimp(cc) = gqm1.LL_seq(end) - null_mod.LL_seq(end);
    end
end


%% FOR MUA
% % for cc = 1:24
%     cc = 23;
%     for bb = 1:n_blocks
%         fprintf('Computing trig avgs for MUA %d of %d\n',cc,24);
%         cur_used_blocks = bb;
%         cur_poss_inds = find(ismember(all_blockvec,cur_used_blocks));
%         cur_used_inds = used_inds(ismember(used_inds,cur_poss_inds));
%         
%         if ~isempty(cur_used_inds)
%             
%             Robs = all_binned_mua(cur_used_inds,cc);
%             cur_X{1} = all_Xmat(cur_used_inds,:);
%             cur_X{2} = Xblock(cur_used_inds,:);
%             
%             gqm1 = qmods(cc);
%             gqm1 = NMMfit_filters(gqm1,Robs,cur_X,[],[],silent);
%             
%             null_mod = NMMinitialize_model(null_stim_params,[1],{'lin'},null_reg_params,[1]);
%             null_mod = NMMfit_filters(null_mod,Robs,cur_X(2:end));
%             
%             qmods_bb(cc,bb) = gqm1;
%             qmod_bb_LLimp(cc,bb) = gqm1.LL_seq(end) - null_mod.LL_seq(end);
%         end
%     end
% % end
%%
cd(save_dir)
save_name = 'init_mod_fits';
save(save_name,'*qmod*');
load('full_sacmod_data');

%%
close all
f2 = figure;
for ss = 1:length(sua_data)
    cur_used_blocks = find(su_used_blocks(:,ss)); %blocks when NO SU
    cur_poss_inds = find(ismember(all_blockvec,cur_used_blocks));
    cur_used_inds = used_inds(ismember(used_inds,cur_poss_inds));
    Robs = all_binned_sua(cur_used_inds,ss);
    cur_X{1} = all_Xmat(cur_used_inds,:);
    cur_X{2} = Xblock(cur_used_inds,:);

    fig_handles = NMMdisplay_model(su_qmods(ss),cur_X,Robs,1);
   
   figure(f2); clf;hold on
%    plot(anal_params.lags*anal_params.dt,sua_data(ss).gsac_avg);
   subplot(3,1,1);hold on
   plot(anal_params.lags*anal_params.dt,sua_data(ss).gsac_im_avg,'r');
   plot(anal_params.lags*anal_params.dt,sua_data(ss).gsac_gray_avg,'k');
   subplot(3,1,2);hold on
   plot(anal_params.lags*anal_params.dt,sua_data(ss).msac_im_avg,'r');
   plot(anal_params.lags*anal_params.dt,sua_data(ss).msac_gray_avg,'k');
   subplot(3,1,3);hold on
   plot(anal_params.lags*anal_params.dt,sua_data(ss).sim_gsac_avg,'r');
   plot(anal_params.lags*anal_params.dt,sua_data(ss).sim_msac_avg,'k');
   axis tight; box off
%    legend('Guided, Im-back','Guided, Gray-back','Micro, Im-back','Micro, Gray-back')
    xlabel('Time since saccade onset (s)','fontsize',14)
    ylabel('Relative firing rate','fontsize',14)
   pause; figure(fig_handles.stim_filts); close;
end

%%
close all
f2 = figure;
for ss = 23:length(mua_data)
    ss
   fig_handles = NMMdisplay_model(qmods(ss),[],[],1);
   
   figure(f2); clf;hold on
   subplot(3,1,1);hold on
   plot(anal_params.lags*anal_params.dt,mua_data(ss).gsac_im_avg,'r');
   plot(anal_params.lags*anal_params.dt,mua_data(ss).gsac_gray_avg,'k');
   subplot(3,1,2);hold on
   plot(anal_params.lags*anal_params.dt,mua_data(ss).msac_im_avg,'r');
   plot(anal_params.lags*anal_params.dt,mua_data(ss).msac_gray_avg,'k');
   subplot(3,1,3);hold on
   plot(anal_params.lags*anal_params.dt,mua_data(ss).sim_gsac_avg,'r');
   plot(anal_params.lags*anal_params.dt,mua_data(ss).sim_msac_avg,'k');
    
   pause; figure(fig_handles.stim_filts); close;
end

%%
close all
for ss = 1:24
    filts = [qmods(ss).mods(1:3).filtK];
    filts = reshape(filts,[flen use_nPix 3]);
    spatial_envs = squeeze(std(filts));
    avg_sp_envs(ss,:) = mean(spatial_envs,2);
    temporal_envs = squeeze(std(filts,[],2));
    avg_sp_envs(ss,:) = mean(spatial_envs,2);
    avg_temp_envs(ss,:) = mean(temporal_envs,2);
    lin_temp_envs(ss,:) = squeeze(temporal_envs(:,1));
    [fit_params,fit_z] = fitGaussianCurve((1:use_nPix)',avg_sp_envs(ss,:)');
    rf_meanloc(ss) = fit_params(1);
    rf_std(ss) = fit_params(2);
end

figure
temp_ax = (0:(flen-1))*dt;
pix_ax = ((1:use_nPix) - use_nPix/2)*Expts{6}.Stimvals.dw;
x_off = cos(Expts{2}.Stimvals.or*pi/180)*Expts{2}.Stimvals.xo;
y_off = sin(Expts{2}.Stimvals.or*pi/180)*Expts{2}.Stimvals.yo;
imagesc(pix_ax,1:24,avg_sp_envs);
xlabel('Relative position (deg)','fontsize',14)
ylabel('Probe depth','fontsize',14)
fillPage(gcf,'papersize',[5 5]);

Navg_temp_envs = bsxfun(@rdivide,avg_temp_envs,max(avg_temp_envs,[],2));
Nlin_temp_envs = bsxfun(@rdivide,lin_temp_envs,max(lin_temp_envs,[],2));
figure
imagesc(temp_ax,1:24,avg_temp_envs);
% imagesc(temp_ax,1:24,lin_temp_envs);
xlabel('Time lag (s)','fontsize',14)
ylabel('Probe depth','fontsize',14)
fillPage(gcf,'papersize',[5 5]);

%% SU
close all
su_num = 2;

repeat_trials = find(all_trial_Se == 1e4);
n_rpts = length(repeat_trials)

rpt_taxis = (1:round(trial_dur)/dt)*dt;
psth = nan(n_rpts,round(trial_dur/dt));
y_cnt = 0; y_space = 3;
for ii = 1:n_rpts
    cur_inds = find(all_trialvec == repeat_trials(ii));
    cur_inds(size(psth,2)+1:end) = [];
    psth(ii,1:length(cur_inds)) = all_binned_sua(cur_inds,su_num);

    cur_spk_set = find(all_su_spk_times{su_num} > all_trial_start_times(repeat_trials(ii)) & ...
        all_su_spk_times{su_num} < all_trial_end_times(repeat_trials(ii)));
    spk_times = all_su_spk_times{su_num}(cur_spk_set) - all_trial_start_times(repeat_trials(ii));
    temp_hist(ii,:) = hist(spk_times,rpt_taxis);
    plot(spk_times,ones(size(spk_times))*y_cnt,'k.')
    hold on
    y_cnt = y_cnt + y_space;
end
cur_X{1} = all_Xmat(cur_inds,:);
cur_X{2} = ones(length(cur_inds),n_blocks)/n_blocks;
[LL, penLL, pred_rate] = NMMmodel_eval(su_qmods(su_num),[],cur_X);

plot(rpt_taxis,nanmean(psth)/dt,'r','linewidth',2);
hold on
plot(rpt_taxis,pred_rate/dt,'b','linewidth',2);

%%
for cc = 1:24
    Xt = [qmods(cc).mods(:).Xtarget];
    filts = reshape([qmods(cc).mods(find(Xt==1)).filtK],[flen use_nPix sum(Xt==1)]);
    spatial_vars(cc,:,:) = squeeze(var(filts,[],2));
end
%% MU
% close all
% mu_num = 1;
% 
% repeat_trials = find(all_trial_Se == 1e4);
% n_rpts = length(repeat_trials)
% 
% rpt_taxis = (1:round(trial_dur)/dt)*dt;
% psth = nan(n_rpts,round(trial_dur/dt));
% y_cnt = 0; y_space = 3;
% for ii = 1:n_rpts
%     cur_inds = find(all_trialvec == repeat_trials(ii));
%     cur_inds(size(psth,2)+1:end) = [];
%     psth(ii,1:length(cur_inds)) = all_binned_sua(cur_inds,su_num);
% 
%     cur_spk_set = find(all_su_spk_times{su_num} > all_trial_start_times(repeat_trials(ii)) & ...
%         all_su_spk_times{su_num} < all_trial_end_times(repeat_trials(ii)));
%     spk_times = all_su_spk_times{su_num}(cur_spk_set) - all_trial_start_times(repeat_trials(ii));
%     temp_hist(ii,:) = hist(spk_times,rpt_taxis);
%     plot(spk_times,ones(size(spk_times))*y_cnt,'k.')
%     hold on
%     y_cnt = y_cnt + y_space;
% end
% cur_X{1} = all_Xmat(cur_inds,:);
% cur_X{2} = ones(length(cur_inds),n_blocks)/n_blocks;
% [LL, penLL, pred_rate] = NMMmodel_eval(su_qmods(su_num),[],cur_X);
% 
% plot(rpt_taxis,nanmean(psth)/dt,'r','linewidth',2);
% hold on
% plot(rpt_taxis,pred_rate/dt,'b','linewidth',2);
% 
