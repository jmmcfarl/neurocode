clear all
% close all

Expt_name = 'G086';
data_dir = ['~/Data/bruce/' Expt_name];

cd(data_dir);

load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

%%
stim_fs = 100; %in Hz
dt = 0.01;
flen = 12;
use_nPix = 24;
stim_params = NIMcreate_stim_params([flen 2*nPix],dt);
Fr = 1;
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

% use_expts = find(expt_binoc(:) == 1 & expt_bar_ori(:) == 90 & expt_sac_dir(:) == 90 & ...
%     expt_dds(:) == 67 & expt_Fr(:) == 1 & expt_npix(:) == 36 & ~excluded_type(:));
use_expts = find(expt_binoc(:) == 1 & expt_bar_ori(:) == 90 & expt_sac_dir(:) == 90 & ...
    expt_dds(:) == 12 & expt_Fr(:) == 1 & ~excluded_type(:));

%%

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
    
    fname = sprintf('%s/stims/Expt%d_stim',data_dir,cur_expt);
    load(fname);
    
    trial_start_times = [Expts{cur_expt}.Trials(:).TrialStart]/1e4;
    trial_end_times = [Expts{cur_expt}.Trials(:).TrueEnd]/1e4;
    trial_durs = [Expts{cur_expt}.Trials(:).dur]/1e4;
    trial_ids = [Expts{cur_expt}.Trials(:).id];
    [un_ids,id_inds] = unique(trial_ids);
    
    buffer_pix = floor((expt_npix(cur_expt) - nPix)/2);
    cur_use_pix = (1:nPix) + buffer_pix;
    
    trial_durs = trial_durs(id_inds);
    trial_start_times = trial_start_times(id_inds);
    trial_end_times = trial_end_times(id_inds);
    
    use_trials = find(trial_durs >= 0.5);
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
        
        
        cur_stim_mat = zeros(length(cur_stim_times),nPix*2);
        cur_stim_mat(:,1:nPix) = left_stim_mats{use_trials(tt)}(:,cur_use_pix);
        cur_stim_mat(:,nPix+1:end) = right_stim_mats{use_trials(tt)}(:,cur_use_pix);
        cur_lstim_mat = left_stim_mats{use_trials(tt)}(:,cur_use_pix);
        cur_rstim_mat = right_stim_mats{use_trials(tt)}(:,cur_use_pix);
        
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


%% CREATE XV TRIAL SET
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

all_used_inds = sort([tr_inds(:); xv_inds(:)]);

%% normalize stimulus variance
all_bar_mat = all_bar_mat/var(reshape(all_bar_mat(all_used_inds,:),1,[]));

%%
fig_out_dir = ['~/Analysis/bruce/' Expt_name '/stc' anal_name];
if ~exist(fig_out_dir,'dir')
    system(['mkdir ' fig_out_dir]);
end

nneg = 3;
npos = 3;
ndims = size(all_bar_mat,2);
unit_sta = nan(96,ndims);
unit_stcs = nan(96,ndims,npos+nneg);
for cc = 1:96
    spikebins = convert_to_spikebins(all_binned_spks(tr_inds,cc));
    spike_cond_stim = all_bar_mat(tr_inds(spikebins),:);
    sta      = mean(spike_cond_stim) - mean(all_bar_mat(tr_inds,:));
    sta = sta/norm(sta);
    proj_mat = sta'/(sta*sta')*sta;
    stim_proj = all_bar_mat(tr_inds,:) - all_bar_mat(tr_inds,:)*proj_mat;
    % stim_proj = stim_emb;
    stvcv = cov(stim_proj(spikebins,:));  utvcv = cov(stim_proj);
    [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
    stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);

    fname = sprintf('%s/unit%d',fig_out_dir,cc);
    figure('visible','off');
    subplot(3,3,1)
    imagesc(reshape(sta,flen,2*nPix));
    ca = max(abs(sta)); caxis([-ca ca]);
    for ii = 1:3
        subplot(3,3,3+ii)
        imagesc(reshape(stcs(:,ii),flen,2*nPix));
        ca = max(abs(stcs(:,ii))); caxis([-ca ca]);
    end
    for ii = 1:3
        subplot(3,3,6+ii)
        imagesc(reshape(stcs(:,end-ii+1),flen,2*nPix));
        ca = max(abs(stcs(:,end-ii+1))); caxis([-ca ca]);
    end
    fillPage(gcf,'papersize',[10 7.5]);
    print(fname,'-dpdf');
    close all

    unit_sta(cc,:) = sta;
    unit_stcs(cc,:,:) = stcs;
end

sname = [anal_dir '/stc_data' anal_name];
save(sname,'unit_sta','unit_stcs');
%%
to_print = 1;
fig_out_dir = ['~/Analysis/bruce/' Expt_name '/quad_mods' anal_name];
if ~exist(fig_out_dir,'dir')
    system(['mkdir ' fig_out_dir]);
end


l_d2XT = 40000;
l_L1 = 150;
ql_L1 = 70;
ql_d2XT = 700;
% ql_L1 = 0.001;
% ql_d2XT = 50;

reg_params = NIMcreate_reg_params('lambda_d2XT',l_d2XT,'lambda_L1',l_L1,'temporal_boundaries','zero');
stim_params = NIMcreate_stim_params([flen 2*nPix],dt,1,1,length(use_expts)-1);
stim_params1 = NIMcreate_stim_params([flen nPix],dt,1,1,length(use_expts)-1);
null_params = NIMcreate_stim_params([1 size(Xexpt,2)],dt);
silent = 1;
optim_params.progTol = 1e-7;
optim_params.optTol = 1e-4;
optim_params.stoch_sampling = round(length(all_stim_times)/20);
max_nfilts = 4;
NL_types = {'lin','quad','quad','quad'}; %define subunit as linear

for cc = [1:96]
    fprintf('Fitting Lin model cell %d of %d\n',cc,96);
    fprintf('Mean spike rate %.3f\n',mean(all_binned_spks(tr_inds,cc))/dt);
    Robs = all_binned_spks(tr_inds,cc);
    null_mod(cc) = NIMinitialize_model(null_params,1,{'lin'},reg_params); %initialize NIM
    null_mod(cc) = NIMfit_filters(null_mod(cc),Robs,Xexpt(tr_inds,:)); %fit stimulus filters
    
    if xv_frac > 0
        Robsxv = all_binned_spks(xv_inds,cc);
        null_xvLL(cc) = NIMmodel_eval(null_mod(cc),Robsxv,Xexpt(xv_inds,:));
    end
    xv_imp{cc} = [];
    %% FIT LINEAR MODEL
    mod_signs = 1; %determines whether input is exc or sup (doesn't matter in the linear case)
    fit0{cc}(1) = NIMinitialize_model(stim_params,mod_signs,NL_types,reg_params); %initialize NIM
    temp = fit0{cc}(1); temp.stim_params = stim_params1;
    L2_mats = create_L2_matrices(temp); %single-ori L2_mats
    full_ndims = prod(stim_params.stim_dims); half_ndims = prod(stim_params1.stim_dims);
    full_L2_mats.L2_d2XT = speye(full_ndims);
    full_L2_mats.L2_d2XT(1:half_ndims,1:half_ndims) = L2_mats.L2_d2XT;
    full_L2_mats.L2_d2XT((half_ndims+1):end,(half_ndims+1):end) = L2_mats.L2_d2XT;
    fit0{cc}(1) = NIMfit_filters(fit0{cc}(1),Robs,all_bar_mat(tr_inds,:),Xexpt(tr_inds,:),[],silent,optim_params,full_L2_mats); %fit stimulus filters
    
    if xv_frac > 0
        [xvLL{cc}(1),~,xv_predrate] = NIMmodel_eval(fit0{cc}(1),Robsxv,all_bar_mat(xv_inds,:),Xexpt(xv_inds,:));
        xv_imp{cc}(1) = (xvLL{cc}(1)-null_xvLL(cc))/log(2);
        fprintf('LL imp: %.3f\n',xv_imp{cc}(1));
        %     [fh] = NIMcheck_nonstat(fit0(cc),Robsxv,all_bar_mat(xv_inds,:),Xexpt(xv_inds,:),25);
    end
    
    fname = sprintf('%s/unit%d_varfilts',fig_out_dir,cc);
    NIMdisplay_model(fit0{cc}(1),[],[],[],0);
    title(sprintf('LL imp %.3f',xv_imp{cc}(1)));
    fillPage(gcf,'papersize',[10 5*1]);
    print(gcf,'-dpsc','-painters',fname);
    close all
    
    %% NOW TRY ADDING QUADRATIC TERMS
    cur_imp = xv_imp{cc}(1);
    cur_nfilts = 2;
    while cur_imp > 0 && cur_nfilts <= max_nfilts
        fprintf('Fitting model with %d filters\n',cur_nfilts);
        
        mod_signs = ones(cur_nfilts); %determines whether input is exc or sup (doesn't matter in the linear case)
        
        %add a filter
        rand_filt = randn(prod(stim_params.stim_dims),1)/prod(stim_params.stim_dims) * 1;
        fit0{cc}(cur_nfilts) = NIMadd_NLinput(fit0{cc}(cur_nfilts-1),'quad',1,rand_filt);
        %adjust regularization of new filter
        fit0{cc}(cur_nfilts) = NIMadjust_regularization(fit0{cc}(cur_nfilts),[length(mod_signs)],'lambda_d2XT',ql_d2XT);
        fit0{cc}(cur_nfilts) = NIMadjust_regularization(fit0{cc}(cur_nfilts),[length(mod_signs)],'lambda_L1',0);
        
        fit0{cc}(cur_nfilts) = NIMfit_filters(fit0{cc}(cur_nfilts),Robs,all_bar_mat(tr_inds,:),Xexpt(tr_inds,:),[],silent,optim_params,full_L2_mats); %fit stimulus filters
        
        %add in L1 for new filter and re-fit
        fit0{cc}(cur_nfilts) = NIMadjust_regularization(fit0{cc}(cur_nfilts),[length(mod_signs)],'lambda_L1',ql_L1);
        fit0{cc}(cur_nfilts) = NIMfit_filters(fit0{cc}(cur_nfilts),Robs,all_bar_mat(tr_inds,:),Xexpt(tr_inds,:),[],silent,optim_params,full_L2_mats); %fit stimulus filters
        
        if xv_frac > 0
            [xvLL{cc}(cur_nfilts),~,xv_predrate] = NIMmodel_eval(fit0{cc}(cur_nfilts),Robsxv,all_bar_mat(xv_inds,:),Xexpt(xv_inds,:));
            xv_imp{cc}(cur_nfilts) = (xvLL{cc}(cur_nfilts)-null_xvLL(cc))/log(2);
            fprintf('LL imp: %.3f\n',xv_imp{cc}(cur_nfilts));
            %     [fh] = NIMcheck_nonstat(fit0(cc),Robsxv,all_bar_mat(xv_inds,:),Xexpt(xv_inds,:),25);
        end
        
        NIMdisplay_model(fit0{cc}(cur_nfilts),[],[],[],0);
        subplot(cur_nfilts,2,2)
        title(sprintf('LL imp %.3f',xv_imp{cc}(cur_nfilts)));
        fillPage(gcf,'papersize',[10 5*cur_nfilts]);
        print(gcf,'-append','-dpsc','-painters',fname);
        close all
        
        cur_imp = xv_imp{cc}(end) - xv_imp{cc}(end-1);
        cur_nfilts = cur_nfilts + 1;
        
    end
end

%%
sname = [anal_dir '/quad_models' anal_name];
save(sname,'fit0','null_mod','xvLL','xv_imp','null_xvLL');

%%
