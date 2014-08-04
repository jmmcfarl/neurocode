clear all
% close all
addpath('~/James_scripts/bruce/eye_tracking/');

Expt_num = 29;
Expt_name = sprintf('G0%d',Expt_num);
data_dir = ['/media/NTlab_data1/Data/bruce/' Expt_name];
if ~exist(data_dir,'dir');
    system(['mkdir ' data_dir]);
end
cd(data_dir);

load(sprintf('%sExpts.mat',Expt_name)); %load in Expts struct

anal_dir = ['~/Analysis/bruce/' Expt_name];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end
cluster_dir = [anal_dir '/clustering/'];

n_probes = 96;
dt = 2*118e-4;
min_trial_dur = 1.75;
trial_dur = 171*dt/2;

%%
temp = cellfun(@(x) x.Header.expname,Expts,'uniformoutput',false);

cur_block_set = find(strcmp(temp,'image.serangeRC'));
cur_block_set = sort([cur_block_set find(strcmp(temp,'image.serangeXseRC'))]);
n_blocks = length(cur_block_set);

%% load overall su data
% LOAD REFCLUSTERS
cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];
fname = [cluster_dir '/final_cluster.mat'];
if exist(fname,'file')
    load(fname);
    SU_numbers = unique(SU_ID_mat(~isnan(SU_ID_mat)));
    for ii = 1:length(SU_numbers)
        SU_tot_nblocks = sum(SU_ID_mat(:) == SU_numbers(ii));
    end
    fprintf('%d SUs Clustered\n',length(SU_numbers));
    
else
    disp('No Cluster data found.');
end


%% COMPUTE TRIAL DATA
cd(data_dir);

fprintf('Computing prep data\n');
trial_cnt = 0;

all_t_axis = [];
all_t_bin_edges = [];
all_tsince_start = [];
all_blockvec = [];
all_trialvec = [];
all_Seeds = [];
all_xo = [];
all_yo = [];
all_trial_blocknums = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_bin_edge_pts = [];
all_spk_times = cell(n_probes,1);
all_clust_ids = cell(n_probes,1);

trial_toffset = zeros(length(cur_block_set),1);
cur_toffset = 0;
for ee = 1:n_blocks;
    fprintf('Block %d of %d\n',ee,n_blocks);
    cur_block = cur_block_set(ee);
    
    fname = [cluster_dir sprintf('/Block%d_Clusters.mat',cur_block)];
    load(fname,'Clusters');
    for cc = 1:n_probes
        all_spk_times{cc} = cat(1,all_spk_times{cc},Clusters{cc}.times + cur_toffset);
        %         all_spk_inds{cc} = cat(1,all_spk_inds{cc},Clusters{cc}.spk_inds + cur_spkind_offset);
        all_clust_ids{cc} = cat(1,all_clust_ids{cc},Clusters{cc}.spike_clusts);
    end
    
    trial_start_times = [Expts{cur_block}.Trials(:).TrialStart]/1e4;
    trial_end_times = [Expts{cur_block}.Trials(:).TrueEnd]/1e4;
    trial_durs = [Expts{cur_block}.Trials(:).dur]/1e4;
    trial_ids = [Expts{cur_block}.Trials(:).id];
    [un_ids,id_inds] = unique(trial_ids);
    rpt_trials = false;
    
    n_rpt_frames = nan(length(trial_start_times),1);
    if isfield(Expts{cur_block}.Trials,'rptframes')
        for tt = 1:length(trial_start_times)
            n_rpt_frames(tt) = length(Expts{cur_block}.Trials(tt).rptframes);
        end
    end
    
    use_trials = find(trial_durs >= min_trial_dur);
    
    dset = find(n_rpt_frames(use_trials) > 0);
    fprintf('Removing %d of %d trials with rpt frames\n',length(dset),length(use_trials));
    use_trials(dset) = [];
    
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)');
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)');
    all_trial_blocknums = cat(1,all_trial_blocknums,ones(length(use_trials),1)*ee);
    
    n_trials = length(use_trials);
    for tt = 1:n_trials
        cur_t_edges = [trial_start_times(use_trials(tt)):dt:trial_end_times(use_trials(tt))];
        cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        n_frames = length(cur_t_axis);
                
        cur_yo = downsample(Expts{cur_block}.Trials(use_trials(tt)).yo,2);
        cur_xo = Expts{cur_block}.Trials(use_trials(tt)).xo;
        cur_Seeds = Expts{cur_block}.Trials(use_trials(tt)).Seedseq;
        
        cur_yo(n_frames+1:end) = [];
        cur_xo(n_frames+1:end) = [];
        cur_Seeds(n_frames+1:end) = [];

        extra =length(cur_t_axis)- length(cur_xo);
        cur_t_axis(end-extra+1:end) = [];
        cur_t_edges(end-extra+1:end) = [];
        
        cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
        
        all_yo = [all_yo; cur_yo(:)];
        all_xo = [all_xo; cur_xo(:)];
        all_Seeds = [all_Seeds; cur_Seeds(:)];
        
        all_t_axis = [all_t_axis; cur_t_axis'];
        all_t_bin_edges = [all_t_bin_edges; cur_t_edges'];
        all_tsince_start = [all_tsince_start; cur_tsince_start'];
        all_blockvec = [all_blockvec; ones(length(cur_t_axis),1)*ee];
        all_trialvec = [all_trialvec; ones(length(cur_t_axis),1)*(tt + trial_cnt)];
        all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
        
        if length(all_t_axis) ~= length(all_Seeds)
            pause
        end
    end
    trial_cnt = trial_cnt + n_trials;
end


%% BIN SPIKES FOR MU AND SU
% for SU probes
fprintf('Using %d SUs\n',length(SU_numbers));
all_binned_sua = nan(length(all_t_axis),length(SU_numbers));
cur_used_blocks = find(ismember(SU_target_blocks,cur_block_set));
SU_ID_mat = SU_ID_mat(cur_used_blocks,:);
[CC,BB] = meshgrid(1:length(SU_clust_data),1:length(cur_used_blocks));

su_probes = nan(1,length(SU_numbers));
all_su_spk_times = cell(length(SU_numbers),1);
% all_su_spk_inds = cell(length(SU_numbers),1);
for ss = 1:length(SU_numbers)
    used_clust_set = unique(CC(SU_ID_mat==ss)); %set of clusters used to capture this SU
    SU_block_probes(ss,:) = nan(1,length(cur_block_set));
    cur_su_spk_times = [];
    %     cur_su_spk_inds = [];
    cur_blocks = [];
    for cc = 1:length(used_clust_set)
        cur_clust = used_clust_set(cc);
        cur_probe = SU_clust_data(cur_clust).probe_num;
        cur_clust_label = SU_clust_data(cur_clust).cluster_label;
        cur_blocks = [cur_blocks find(SU_ID_mat(:,cur_clust) == ss)];
        SU_block_probes(ss,cur_blocks) = cur_probe;
        
        all_su_inds = all_clust_ids{cur_probe} == cur_clust_label;
        cur_su_spk_times = all_spk_times{cur_probe}(all_su_inds);
        %         cur_su_spk_inds = all_spk_inds{cur_probe}(all_su_inds);
        spk_block_inds = round(interp1(all_t_axis,all_blockvec,cur_su_spk_times));
        cur_su_spk_times = cur_su_spk_times(ismember(spk_block_inds,cur_blocks));
        %         cur_su_spk_inds = cur_su_spk_inds(ismember(spk_block_inds,cur_blocks));
        
        all_su_spk_times{ss} = cat(1,all_su_spk_times{ss},cur_su_spk_times(:));
        %         all_su_spk_inds{ss} = cat(1,all_su_spk_inds{ss},cur_su_spk_inds(:));
    end
    
    cur_suahist = histc(all_su_spk_times{ss},all_t_bin_edges);
    cur_suahist(all_bin_edge_pts) = [];
    cur_id_set = ismember(all_blockvec,cur_blocks);
    all_binned_sua(cur_id_set,ss) = cur_suahist(cur_id_set);
    su_probes(ss) = mode(SU_block_probes(ss,~isnan(SU_block_probes(ss,:))));
end

all_binned_mua = nan(length(all_t_axis),n_probes);
%for only-MU probes
for cc = 1:n_probes
    cur_mua_inds = all_clust_ids{cc} == 1;
    %     cur_mua_inds = 1:length(all_spk_times{cc});;
    [cur_spkhist,cur_bin_inds] = histc(all_spk_times{cc}(cur_mua_inds),all_t_bin_edges);
    cur_spkhist(all_bin_edge_pts) = [];
    all_binned_mua(:,cc) = cur_spkhist;
end

%%

Nyp = 1024;
Nxp = 1280;
Ny = round(Nyp/dsfrac);
Nx = round(Nxp/dsfrac);
xax = linspace(-Nx/2,Nx/2,Nx)/Fsd; yax = linspace(-Ny/2,Ny/2,Ny)/Fsd;

rfx = 0.34; rfy = -0.43;

im_patch = [rfx-1 rfx+1; rfy-0.99 rfy+0.99];
im_patch_pix = Fsd*im_patch;
RF_patch = [rfx-0.6 rfx+0.6; rfy-0.6 rfy+0.6];
RF_patch_pix = Fsd*im_patch;

xpatch_inds = find(xax >= im_patch(1,1) & xax <= im_patch(1,2));
ypatch_inds = find(yax >= im_patch(2,1) & yax <= im_patch(2,2));

RF_xpatch_inds = find(xax >= RF_patch(1,1) & xax <= RF_patch(1,2));
RF_ypatch_inds = find(yax >= RF_patch(2,1) & yax <= RF_patch(2,2));

NT = length(all_Seeds);
cd ~/James_scripts/data_processing/Images/image_set_A

all_im_patches = nan(NT,length(ypatch_inds),length(xpatch_inds));
unique_images = unique(all_Seeds);
for ii = 1:length(unique_images)
    
    filename = sprintf('%.4d.png',unique_images(ii));
    IMAGEorg = imread(filename);
    IMAGEorg = double(IMAGEorg); % convert to double format
    IMAGEorg = IMAGEorg - 127.5; %shift range to [-127.5:127.5]
    IMAGE = flipud(IMAGEorg); %flip y
    
    IMAGE = imresize(IMAGE,1/dsfrac); %image downsampling
    
    cur_samp_set = find(all_Seeds == unique_images(ii));
    fprintf('Analyzing image %d, %d samps\n',ii,length(cur_samp_set));
    
    for j = 1:length(cur_samp_set)
        
        ypatch_inds_adj = round(ypatch_inds - all_yo(cur_samp_set(j))*Fsd);
        xpatch_inds_adj = round(xpatch_inds - all_xo(cur_samp_set(j))*Fsd);
        
        cur_patch = IMAGE(ypatch_inds_adj,xpatch_inds_adj);
        all_im_patches(cur_samp_set(j),:,:) = cur_patch;
    end
end

%%
norm_stim = bsxfun(@minus,all_im_patches,mean(all_im_patches));
norm_stim = bsxfun(@rdivide,norm_stim,std(norm_stim));

%%
[XXcoords,YYcoords] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));

RF_subset = find(XXcoords >= RF_patch(1,1) & XXcoords <= RF_patch(1,2) & ...
    YYcoords >= RF_patch(2,1) & YYcoords <= RF_patch(2,2));

flen = 5;
SDIM = length(RF_xpatch_inds);
stim_params = NMMcreate_stim_params([flen SDIM SDIM],dt);
X = create_time_embedding(norm_stim(:,RF_subset),stim_params);


%%
beg_buffer = 0.15; end_buffer = 0.05;
used_inds = find(all_tsince_start >= 0.15 & (trial_dur-all_tsince_start) >= end_buffer);
xv_frac = 0.15;

use_trials = unique(all_trialvec(used_inds));
nuse_trials = length(use_trials);

n_xv_trials = round(xv_frac*nuse_trials);
xv_trials = randperm(nuse_trials);
xv_trials(n_xv_trials+1:end) = [];
xv_trials = use_trials(xv_trials);
tr_trials = setdiff(use_trials,xv_trials);

tr_inds = used_inds(ismember(all_trialvec(used_inds),tr_trials));
xv_inds = used_inds(ismember(all_trialvec(used_inds),xv_trials));

full_inds = [tr_inds; xv_inds];

%%
base_lambda_d2XT = 400;
base_lambda_L1 = 10;

n_squared_filts = 1;
mod_signs = ones(1,n_squared_filts+1);
NL_types = {'lin' 'quad'};
init_d2XT = [10*ones(n_squared_filts+1,1)];
init_reg_params = NMMcreate_reg_params('lambda_d2XT',init_d2XT);
init_Xtargs = [ones(n_squared_filts+1,1)];
silent = 1;

n_cells = size(all_binned_sua,2);
for tcell = 1:n_cells;
% tcell = 3;
fprintf('Fitting model for Cell %d of %d\n',tcell,n_cells);
cur_tr_inds = tr_inds(~isnan(all_binned_sua(tr_inds,tcell)));
cur_xv_inds = xv_inds(~isnan(all_binned_sua(xv_inds,tcell)));

    gqm1 = NMMinitialize_model(stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs);
    gqm1 = NMMfit_filters(gqm1,all_binned_sua(cur_tr_inds,tcell),X(cur_tr_inds,:),[],[],silent);
    [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(gqm1,all_binned_sua(cur_tr_inds,tcell),X(cur_tr_inds,:));
    gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_d2XT',base_lambda_d2XT./var(gint)');
    gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_L1',base_lambda_L1./std(gint)');
    gqm1 = NMMfit_filters(gqm1,all_binned_sua(cur_tr_inds,tcell),X(cur_tr_inds,:),[],[],silent);
    init_mod_fit(tcell) = gqm1;
     
    [sua_xvLL(tcell), penLL, pred_rate, G, gint,~,sua_nullxvLL(tcell)] = NMMmodel_eval(gqm1,all_binned_sua(cur_xv_inds,tcell),X(cur_xv_inds,:));
    fprintf('LL imp: %.4f\n',sua_xvLL(tcell)-sua_nullxvLL(tcell));
end
%%
n_cells = size(all_binned_mua,2);
for tcell = 1:n_cells;
% tcell = 3;
fprintf('Fitting model for Cell %d of %d\n',tcell,n_cells);
cur_tr_inds = tr_inds(~isnan(all_binned_mua(tr_inds,tcell)));
cur_xv_inds = xv_inds(~isnan(all_binned_mua(xv_inds,tcell)));

    gqm1 = NMMinitialize_model(stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs);
    gqm1 = NMMfit_filters(gqm1,all_binned_mua(cur_tr_inds,tcell),X(cur_tr_inds,:),[],[],silent);
    [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(gqm1,all_binned_mua(cur_tr_inds,tcell),X(cur_tr_inds,:));
    gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_d2XT',base_lambda_d2XT./var(gint)');
    gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_L1',base_lambda_L1./std(gint)');
    gqm1 = NMMfit_filters(gqm1,all_binned_mua(cur_tr_inds,tcell),X(cur_tr_inds,:),[],[],silent);
    init_mua_fit(tcell) = gqm1;
     
    [mua_xvLL(tcell), penLL, pred_rate, G, gint,~,mua_nullxvLL(tcell)] = NMMmodel_eval(gqm1,all_binned_mua(cur_xv_inds,tcell),X(cur_xv_inds,:));
       fprintf('LL imp: %.4f\n',mua_xvLL(tcell)-mua_nullxvLL(tcell));
 
end