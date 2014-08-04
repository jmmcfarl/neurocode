clear all
% close all
cd ~/Data/bruce/G081/
load jbeG081Expts.mat
addpath('~/James_scripts/bruce/G081/');
%%
stim_fs = 100; %in Hz
use_sus = 1:96;
Fs = 3e4;
dsf = 60;Fsd = Fs/dsf;
niqf = Fsd/2;
[filt_b,filt_a] = butter(2,[1.5 40]/niqf);
use_lfps = [1:8:96];

%% PARSE TRIAL DATA STRUCTURES
for i = 1:length(Expts)
    if strcmp(Expts{i}.Header.expname,'grating.OpXseRC') | strcmp(Expts{i}.Header.expname,'grating.OpRC')
        is_bar_expt(i) = 1;
    else
        is_bar_expt(i) = 0;
    end
    
    if strcmp(Expts{i}.Stimvals.Bs,'image')
        expt_image_back(i) = 1;
    else
        expt_image_back(i) = 0;
    end
    
    expt_sim_sacs(i) = Expts{i}.Stimvals.ijump;
    expt_bar_ori(i) = Expts{i}.Stimvals.or;
    
end
expt_bar_ori(expt_bar_ori == -45) = 135;

load ./all_un_bar_pos
n_bar_pos = size(all_un_bar_pos,1);

%% USE ONLY GRAY BACKGROUND DATA
flen = 25;
old_flen = 12;
beg_buffer = round(stim_fs*0.15);
bar_oris = [0];
un_bar_pos = all_un_bar_pos(:,1);

fprintf('Analyzing %d ori expts\n',bar_oris);

%expts with X deg bars and any back (including sim sacs)
cur_expt_set = find(is_bar_expt==1 & expt_bar_ori == bar_oris);

% %expts with X deg bars and any back (including sim sacs)
% cur_expt_set = find(is_bar_expt==1 & expt_bar_ori == bar_oris & expt_image_back == 0);

% %expts with X deg bars and any back (including sim sacs)
% cur_expt_set = find(is_bar_expt==1 & expt_bar_ori == bar_oris & expt_image_back == 1);

cur_expt_set(cur_expt_set<=8) = []; %this expt had continuous bar luminance
cur_expt_set(cur_expt_set == 13) = []; %problem with LFP data?
cur_expt_set(cur_expt_set > 60) = [];

%% COMPUTE TRIAL DATA
all_stim_times = [];
all_rel_stimes = [];
all_rel_etimes = [];
all_phase = [];
all_Op = [];
all_bar_mat = [];
all_used_inds = [];
all_or_used_inds = [];
all_binned_spks = [];
all_exptvec = [];
all_trialvec = [];
for ee = 1:length(cur_expt_set);
    fprintf('Expt %d of %d\n',ee,length(cur_expt_set));
    cur_expt = cur_expt_set(ee);
    fname = sprintf('Expt%dClusterTimes.mat',cur_expt);
    load(fname);
    
    n_trials = length(Expts{cur_expt}.Trials);
    trial_start_times = [Expts{cur_expt}.Trials(:).TrialStart]/1e4;
    trial_end_times = [Expts{cur_expt}.Trials(:).TrueEnd]/1e4;
    trial_durs = [Expts{cur_expt}.Trials(:).dur]/1e4;
    
    for tt = 1:n_trials
        cur_stim_times = Expts{cur_expt}.Trials(tt).Start/1e4;
        cur_Op = Expts{cur_expt}.Trials(tt).Op;
        cur_phase = Expts{cur_expt}.Trials(tt).ph;
        
        cur_t_edges = [cur_stim_times; Expts{cur_expt}.Trials(tt).End(end)/1e4];
        cur_binned_spks = nan(length(cur_stim_times),length(use_sus));
        for cc = 1:length(use_sus)
            cur_hist = histc(Clusters{use_sus(cc)}.times,cur_t_edges);
            cur_binned_spks(:,cc) = cur_hist(1:end-1);
        end
        
        cur_bar_mat = zeros(length(cur_stim_times),n_bar_pos);
        for bb = 1:n_bar_pos
            cur_set = find(cur_Op==un_bar_pos(bb));
            pset = cur_phase(cur_set) == 0;
            nset = cur_phase(cur_set) == pi;
            cur_bar_mat(cur_set(pset),bb) = 1;
            cur_bar_mat(cur_set(nset),bb) = -1;
            %                 cur_bar_mat(cur_set,bb) = 1;
        end
        
        bar_Xmat = makeStimRows(cur_bar_mat,flen);
        cur_used_inds = ones(length(cur_stim_times),1);
        cur_used_inds(1:flen) = 0;
        cur_used_inds(1:beg_buffer) = 0;
        
        cur_or_used_inds = ones(length(cur_stim_times),1);
        cur_or_used_inds(1:old_flen) = 0;
        cur_or_used_inds(1:beg_buffer) = 0;
        
        all_stim_times = [all_stim_times; cur_stim_times];
        all_rel_stimes = [all_rel_stimes; cur_stim_times- trial_start_times(tt)];
        all_rel_etimes = [all_rel_etimes; trial_end_times(tt) - cur_stim_times];
        all_Op = [all_Op; cur_Op];
        all_phase = [all_phase; cur_phase'];
        all_used_inds = [all_used_inds; cur_used_inds];
        all_or_used_inds = [all_or_used_inds; cur_or_used_inds];
        all_bar_mat = [all_bar_mat; bar_Xmat];
        all_binned_spks = [all_binned_spks; cur_binned_spks];
        all_exptvec = [all_exptvec; ones(size(cur_stim_times))*ee];
        all_trialvec = [all_trialvec; ones(size(cur_stim_times))*tt];
    end
    
end

%%
all_Vmat_interp = [];
for ee = 1:length(cur_expt_set);
    cur_expt = cur_expt_set(ee);
    
    Vmat = [];
    for ll = 1:length(use_lfps)
        fprintf('Electrode %d of %d\n',ll,length(use_lfps));
        filename = sprintf('Expt%d.p%dFullV.mat',cur_expt,use_lfps(ll));
        load(filename);
        V = double(FullV.V);
        %     V = V + FullV.sumscale*sumv;
        V = V*FullV.intscale(1)/FullV.intscale(2);
        nparts = length(FullV.blklen);
        dV = [];
        %splice together multiple blocks
        cur_pt = 1;
        for pp = 1:nparts
            cur_range = cur_pt:(cur_pt + FullV.blklen(pp)-1);
            cur_range(cur_range > length(V)) = [];
            curV = decimate(V(cur_range),dsf);
            curV = filtfilt(filt_b,filt_a,curV);
            dV = [dV curV];
            cur_pt = cur_pt + FullV.blklen(pp);
        end
        Vmat(:,ll) = dV;
    end
    
    t_ax = [];
    for pp = 1:nparts
        cur_t_ax = linspace(FullV.blkstart(pp),FullV.blkstart(pp)+FullV.blklen(pp)/Fs,FullV.blklen(pp));
        t_ax = [t_ax downsample(cur_t_ax,dsf)];
    end
    t_ax(size(Vmat,1)+1:end) = [];
    
    fprintf('LFP len: %d\n',range(t_ax));
    
    use_set = find(all_exptvec==ee);
    Vmat_interp = interp1(t_ax,Vmat,all_stim_times(use_set));
    all_Vmat_interp = [all_Vmat_interp; Vmat_interp];
end

%% IF YOU WANT TO DO CROSS-VALIDATION ON INITIAL MODEL FITS, OTHERWISE, 
[c,ia,ic] = unique([all_exptvec all_trialvec],'rows');
n_trials = length(ia);
% bad_trials = unique(ic(bad_pts)); %trials with putative blinks

tr_set = 1:n_trials;
tr_set = find(ismember(1:n_trials,tr_set));

cur_tr_inds = find(ismember(ic,tr_set));

cur_tr_inds(all_used_inds(cur_tr_inds) == 0) = [];
% Xmat = abs(all_bar_mat);
Xmat = all_bar_mat;

%% LOAD AND INTERPOLATE MODELS

load ./full_eye_correct_0deg.mat
use_tr_inds = find(ismember(eye_times,all_stim_times(cur_tr_inds)));
% use_tr_inds = find(ismember(tr_inds,cur_tr_inds));

%% double sample Xmat
bar_dx = 0.08;
up_fac = 2;
bar_dxd = bar_dx/up_fac;
SDIM = n_bar_pos;
NT = length(cur_tr_inds);
resh_X = reshape(Xmat(cur_tr_inds,:)',[flen SDIM NT]);

if up_fac > 1
fprintf('Up-sampling X-matrix\n');
dresh_X = [];
for i = 1:SDIM-1
    dresh_X = cat(2,dresh_X,resh_X(:,i,:));
    dresh_X = cat(2,dresh_X,0.5*resh_X(:,i,:) + 0.5*resh_X(:,i+1,:));
end
dresh_X = cat(2,dresh_X,resh_X(:,end,:));

    up_SDIM = up_fac*SDIM-1;
else
    up_SDIM = SDIM;
    dresh_X = resh_X;
end
X_z = reshape(dresh_X,flen*up_SDIM,NT)';


%% RECONSTRUCT MAP STIMULUS
resh_X = reshape(X_z',[flen up_SDIM NT]);
resh_X_sh = zeros(size(resh_X));
for ii = 1:NT
    d2 = dist_shift2d(resh_X(:,:,ii), -it_shift_cor{end}(use_tr_inds(ii)),2,0);
    resh_X_sh(:,:,ii) = d2;
end
X_sh = reshape(resh_X_sh,flen*up_SDIM,NT)';

%%
Yobs = zscore(all_Vmat_interp(cur_tr_inds,:));

full_X{1} = X_sh; full_X{1}(full_X{1} == -1) = 0;
full_X{2} = X_sh; full_X{2}(full_X{2} == 1) = 0;

%%
stim_params.spatial_dims = 1;
stim_params.sdim = up_SDIM;
stim_params.flen = flen;
klen = up_SDIM*flen;
clear defmod

for cc = 1:length(use_lfps)
    fprintf('LFP %d of %d\n',cc,length(use_lfps));
    yobs_var = var(Yobs(:,cc));
    
    defmod.lambda_L1x = 20;
    defmod.lambda_d2T = 25;
    defmod.lambda_d2X = 75;
    defmod.lambda_d2XT = 50;
 
    nmods = 2;
    kern_types{1} = 'lin'; kern_types{2} = 'lin';
    kern_stim_inds = [1 2];
    init_kerns = 0.01*randn(klen,nmods);
    glm = createGNM_sepstim(init_kerns,[1 1],kern_types,kern_stim_inds,defmod,stim_params,'gauss');
    glm_fit_cor_full(cc) = fitGNM_filters_sepstim(glm,full_X,Yobs(:,cc),'none',[],1e-4,1e-6,0);
    [glm_sse, ~, ~, prate] = getLL_GNM_sepstim(glm_fit_cor_full(cc),full_X,Yobs(:,cc),'none');
    glm_cor_full_r2(cc) = 1-glm_sse/yobs_var;
    lfp_full_residual(cc,:) = Yobs(:,cc) - prate;
  
%     glm_fit_cor_full(cc) = fitGNM_filters_sepstim(glm,full_X,Yobs(:,cc),'none',[],1e-4,1e-6,0);
%     [glm_sse, ~, ~, prate] = getLL_GNM_sepstim(glm_fit_cor_full(cc),full_X,Yobs(:,cc),'none');
%     glm_cor_full_r2(cc) = 1-glm_sse/yobs_var;

%     nmods = 4;
%     kern_types{1} = 'lin'; kern_types{2} = 'lin'; kern_types{3} = 'quad'; kern_types{4} = 'quad';
%     kern_stim_inds = [2 1 1 1];
%     init_kerns = 0.01*randn(klen,nmods);
%     glm = createGNM_sepstim(init_kerns,[1 1 1 -1],kern_types,kern_stim_inds,defmod,stim_params,'gauss');
%     lin_reg_scale = 5;
%     glm.mods(3).lambda_L1x = defmod.lambda_L1x/lin_reg_scale/4;
%     glm.mods(3).lambda_d2T = defmod.lambda_d2T/lin_reg_scale;
%     glm.mods(3).lambda_d2X = defmod.lambda_d2X/lin_reg_scale;
%     glm.mods(3).lambda_d2XT = defmod.lambda_d2XT/lin_reg_scale;
%     glm.mods(4).lambda_L1x = defmod.lambda_L1x/lin_reg_scale/4;
%     glm.mods(4).lambda_d2T = defmod.lambda_d2T/lin_reg_scale;
%     glm.mods(4).lambda_d2X = defmod.lambda_d2X/lin_reg_scale;
%     glm.mods(4).lambda_d2XT = defmod.lambda_d2XT/lin_reg_scale;
%     
%     glm_fit_cor_full(cc) = fitGNM_filters_sepstim(glm,full_X,Yobs(:,cc),'none',[],1e-4,1e-6,0);
%     [glm_sse, ~, ~, prate] = getLL_GNM_sepstim(glm_fit_cor_full(cc),full_X,Yobs(:,cc),'none');
%     glm_cor_full_r2(cc) = 1-glm_sse/yobs_var;

    
        nmods = 1;
    kern_types{1} = 'lin';
    init_kerns = 0.01*randn(klen,nmods);
    glm = createGNM(init_kerns,[1],kern_types,defmod,stim_params,'gauss');
   
    glm_fit_cor(cc) = fitGNM_filters(glm,abs(X_sh),Yobs(:,cc),'none',[],1e-4,1e-6,0);
    [glm_sse, ~, ~, prate] = getLL_GNM(glm_fit_cor(cc),abs(X_sh),Yobs(:,cc),'none');
    glm_cor_r2(cc) = 1-glm_sse/yobs_var;

    lfp_residual(cc,:) = Yobs(:,cc) - prate;
    
%     fprintf('R2: %.3f  R2_cor: %.3f\n',glm_cor_r2(cc),glm_cor_full_r2(cc));
end
model_t_axis = all_stim_times(cur_tr_inds);

%%
save lfp_models_linear_0deg glm_fit_cor glm_cor_r2 lfp_residual use_lfps model_t_axis cur_tr_inds glm_fit_cor_full glm_cor_full_r2 lfp_full_residual

%%
close all
for cc = 1:length(use_lfps)
    
    cur_filts = get_k_mat(glm_fit_cor_full(cc));
    ca = max(abs(cur_filts(:)));
    subplot(2,1,1)
    imagesc(reshape(cur_filts(:,1),flen,up_SDIM));colormap(jet);
    caxis([-ca ca])
    subplot(2,1,2)
    imagesc(reshape(cur_filts(:,2),flen,up_SDIM))
    caxis([-ca ca])
    
    glm_cor_full_r2(cc)
   pause
   clf
end

%%
blag = round(0.1*Fsd);
flag = round(0.3*Fsd);
cur_inds = find(X_sh(:,end-20*flen) == 1);

trig_avg = get_event_trig_avg(Yobs(:,1),cur_inds,blag,flag);
trig_avgr = get_event_trig_avg(lfp_full_residual(1,:),cur_inds,blag,flag);

