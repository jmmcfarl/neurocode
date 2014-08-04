%% Load Data
clear all;
addpath(genpath('~/Code/James_scripts'));
cd ~/Data/bruce/7_15_12
cd G029/

load ./CellList.mat
single_units = find(CellList(1,:,1) > 0);

cd ~/Data/bruce/7_15_12/G029/
load ./Expt3_newcompiled_data_d1p25
resh_all_stims = resh_all_stims/std(resh_all_stims(:));

load ./all_eyedata_expt3
Pix2Deg = 0.018837;
NT = length(full_stim_ids);

%% crop stimulus for the purpose of faster gabor function fitting
sdim = length(xpatch_inds);
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));

%% RECREATE TIME-EMBEDDED STIM MAT
flen = 8;

% PARSE DATA INTO FIXATIONS
diff_used_inds = [1; diff(used_inds)];
rel_fix_start_inds = [1; find(diff_used_inds > 1)];
rel_fix_stop_inds = [(find(diff_used_inds > 1)-1); NT];
n_trials = length(rel_fix_start_inds);


maxdur = sum(rel_fix_stop_inds-rel_fix_start_inds);
fullX = nan(maxdur,sdim^2*flen);
nfull_binned_spks = nan(maxdur,96);

nused_inds = [];
ntrial_inds = [];
cnt = 0;
for i = 1:n_trials
    fprintf('%d of %d\n',i,n_trials)
    cur_inds = rel_fix_start_inds(i):rel_fix_stop_inds(i);
    temp = makeStimRows(resh_all_stims(full_stim_ids(cur_inds),:),flen,1);
    fullX(cnt + (1:size(temp,1)),:) = temp;
    nfull_binned_spks(cnt + (1:size(temp,1)),:) = full_binned_spks(cur_inds(flen:end),:);
    nused_inds = [nused_inds; cur_inds(flen:end)'];
    ntrial_inds = [ntrial_inds; i*ones(length(cur_inds(flen:end)),1)];
    cnt = cnt + size(temp,1);
end
fullX(cnt+1:end,:) = [];
fullX = fullX/std(fullX(:));
nfull_binned_spks(cnt+1:end,:) = [];

full_binned_spks = nfull_binned_spks;

full_expt_vec = full_expt_vec(nused_inds);
full_trial_vec = full_trial_vec(nused_inds);
full_t = full_t(nused_inds);

% PARSE DATA INTO FIXATIONS
diff_used_inds = [1; diff(nused_inds)];
rel_fix_start_inds = [1; find(diff_used_inds > 1)];
rel_fix_stop_inds = [(find(diff_used_inds > 1)-1); NT];
n_trials = length(rel_fix_start_inds);

%%
NT = length(nused_inds);
xv_frac = 0.;
xv_num = round(xv_frac*n_trials);
xv_set = randperm(n_trials);
xv_set(xv_num+1:end) = [];
tr_set = setdiff(1:n_trials,xv_set);

xv_inds = [];
for ii = 1:xv_num
    cur_inds = find(ntrial_inds == xv_set(ii));
    xv_inds = [xv_inds; cur_inds];
end
tr_inds = setdiff(1:NT,xv_inds);

tot_dims = flen*sdim^2;

[un_expts,~,full_expt_inds] = unique(full_expt_vec);
n_un_expts = length(un_expts);
linX = zeros(NT,n_un_expts);
for i = 1:n_un_expts
    linX(full_expt_inds==i,i) = 1;
end

%% COMPUTE OUTPUTS OF GABOR MODELS BASED ON RF MAPPING DATA
load ./expt1_eyecor_d1p25_nosac_v2.mat gabor*
gabor_params = gabor_params_f{end};

load ./spatiotemp_gaborfits_expt1
old_init_gnm = init_gnm;

%initialize model
defmod.lambda_dT = 0;
defmod.lambda_d2X = 0;
defmod.lambda_L1x = 0;

stim_params.flen = flen;
stim_params.spatial_dims = 1;
stim_params.sdim = 2;
cur_ndims = 2;

stim_params2.flen = 1;
stim_params2.spatial_dims = 0;
stim_params2.sdim = 1;

clear gabor_emp*
gabor_emp1_filt = nan(96,tot_dims);
gabor_emp2_filt = nan(96,tot_dims);
for t = 1:96
    fprintf('Fitting cell %d of %d\n',t,96);
    gabor_emp1 = get_pgabor_mask_v2(XX,YY,gabor_params(t,1:6),0);
    gabor_emp2 = get_pgabor_mask_v2(XX,YY,gabor_params(t,1:6),pi/2);
    
    gabor_emp1_vec = [gabor_emp1(:)'; zeros(flen-1,sdim^2)];
    gabor_emp1_mat = makeStimRows(gabor_emp1_vec,flen);
    
    gabor_emp2_vec = [gabor_emp2(:)'; zeros(flen-1,sdim^2)];
    gabor_emp2_mat = makeStimRows(gabor_emp2_vec,flen);
    
    gabor_outs1 = fullX*gabor_emp1_mat';
    gabor_outs2 = fullX*gabor_emp2_mat';
    
        spikebins = convert_to_spikebins(full_binned_spks(:,t));
    tr_spikebins = convert_to_spikebins(full_binned_spks(tr_inds,t));
%     xv_spikebins = convert_to_spikebins(full_binned_spks(xv_inds,t));
    
    X = [gabor_outs1 gabor_outs2];
    [NT,klen] = size(X);
    k0 = zeros(klen,1);
    
    init_signs = ones(cur_ndims,1);
    init_kerns = 0.05*randn(klen,cur_ndims);
    for i = 1:cur_ndims
        kern_types{i} = 'quad';
    end
    init_linK = zeros(n_un_expts,1);
    temp_gnm = createGNM_v2(init_kerns,init_signs,kern_types,init_linK,defmod,stim_params);
    init_gnm(t) = fitGNM_filters_v2(temp_gnm,X(tr_inds,:),linX(tr_inds,:),tr_spikebins,'none',[],1e-4,1e-6,1);
%     xvLL_gnm(t) = getLL_GNM_v2(init_gnm(t),X(xv_inds,:),linX(xv_inds,:),xv_spikebins,'none');
    
    conv_mat = [gabor_emp1_mat' gabor_emp2_mat'];
    kmat = get_k_mat(init_gnm(t));
%     kmat = get_k_mat(old_init_gnm(t));
%     kouts = X*kmat;
%     g = sum(kouts.^2,2);
%     init_signs = 1;
%     init_kerns = 0;
%     kern_types{1} = 'lin';
%     temp_gnm = createGNM_v2(init_kerns,init_signs,kern_types,init_linK,defmod,stim_params2);
%     init_gnm2(t) = fitGNM_filters_v2(temp_gnm,g(tr_inds),linX(tr_inds,:),tr_spikebins,'none',[],1e-4,1e-6,0);
%     xvLL_gnm2(t) = getLL_GNM_v2(init_gnm2(t),g(xv_inds),linX(xv_inds,:),xv_spikebins,'none');
      
    k_pix = conv_mat*kmat;
    gabor_emp1_filt(t,:) = k_pix(:,1)';
    gabor_emp2_filt(t,:) = k_pix(:,2)';

    
        [nll, pnll, lpen, prate,g,g_int] = getLL_GNM_v2(init_gnm(t),X,linX,spikebins,'none');
    
    %     new_gint = zeros(NT,1);
    %     for i = 1:2
    %         new_gint = new_gint + g_int{i};
    %     end
    %     new_linX = [bsxfun(@times,linX,new_gint) linX];
    %     init_signs = [];
    %     init_kerns = [];
    %     kern_types = [];
    %     init_linK = zeros(2*n_un_expts,1);
    %     temp_gnm = createGNM_v2(init_kerns,init_signs,kern_types,init_linK,defmod,stim_params);
    %     temp_gnm = fitGNM_filters_v2(temp_gnm,X(tr_inds,:),new_linX(tr_inds,:),tr_spikebins,'none',[],1e-4,1e-6);
    %     [nll, pnll, lpen, prate2] = getLL_GNM_v2(temp_gnm,X,new_linX,spikebins,'none');
    
    
%         figure
%         plot(full_t,smooth(full_binned_spks(:,t),5),'.-')
%         hold on
%         plot(full_t,prate,'r')
%     % %     plot(prate2,'k')
%         pause
%         clf
    
    %
    %         % FIT SPATIAL MODEL
    %     cur_gabor_emp1 = gabor_emp1_mat(3,:);
    %     cur_gabor_emp2 = gabor_emp2_mat(3,:);
    %     cur_gabor_outs1 = fullX*cur_gabor_emp1';
    %     cur_gabor_outs2 = fullX*cur_gabor_emp2';
    %     cur_gabor_energy = cur_gabor_outs1.^2 + cur_gabor_outs2.^2;
    %
    %      X = cur_gabor_energy;
    %     [NT,klen] = size(X);
    %     k0 = zeros(klen,1);
    %     init_signs = 1;
    %     init_kerns = 0.05*randn(klen,1);
    %     kern_types{1} = 'lin';
    %     init_linK = zeros(n_un_expts,1);
    %     temp_gnm = createGNM_v2(init_kerns,init_signs,kern_types,init_linK,defmod,stim_params2);
    %     init_sgnm(t) = fitGNM_filters_v2(temp_gnm,X(tr_inds,:),linX(tr_inds,:),tr_spikebins,'none',[],1e-4,1e-6);
    %     xvLL_sgnm(t) = getLL_GNM_v2(init_sgnm(t),X(xv_inds,:),linX(xv_inds,:),xv_spikebins,'none');
    %
end

tv_theta_params = [init_gnm(:).linK];
ov_theta = [init_gnm(:).spk_theta];
tv_theta_params = bsxfun(@plus,tv_theta_params,ov_theta);


%% SET UP XV CELL SET
NSIG = 96;
xv_frac = 0.2;
tr_set = randperm(NSIG);
tr_set = tr_set(1:round(length(tr_set)*(1-xv_frac)));
xv_set = setdiff(1:NSIG,tr_set);
n_tr_cells = length(tr_set);

%% ESTIMATE LL for each shift in each stimulus frame
orig_eyedt = all_t(2)-all_t(1);
interp_eyespeeds = interp1(all_t,all_eyespeed,full_t);
avg_interp_eyespeed = dt/orig_eyedt*mean(interp_eyespeeds,2);

max_shift = 32;
dshift = 2;
x_shifts = -max_shift:dshift:max_shift;
y_shifts = -max_shift:dshift:max_shift;
[Xsh,Ysh] = meshgrid(x_shifts,y_shifts);
SH = [Xsh(:) Ysh(:)];
n_shifts = size(SH,1);

chunk_dur = 2;

%overall prior on shifts
eps_prior_sigma = 0.2; %0.2
leps_prior = -sum((SH/Fsd).^2,2)/(2*eps_prior_sigma^2);
leps_prior = bsxfun(@minus,leps_prior,logsumexp(leps_prior)); %normalize
lA_tflip = repmat(leps_prior',n_shifts,1);

%state transition matrix (includes a 'constant' prior)
min_deps_sigma = 0.01*chunk_dur;
max_deps_sigma = 0.4;
deps_offset = 0.35;
deps_slope = 15;
sac_deps_sigmas = (max_deps_sigma-min_deps_sigma)./(1+exp(-(avg_interp_eyespeed-deps_offset)*deps_slope))+min_deps_sigma;

cdist = squareform(pdist(SH/Fsd));

ov_lgamma = nan(NT,n_shifts);

n_fixs = length(rel_fix_start_inds);
n_iter = 1;
for it = 1:n_iter
    
    frame_LLs = zeros(NT,n_shifts);
    Robs = full_binned_spks(:,tr_set);
    
    gabor_filt_bank1 = reshape(gabor_emp1_filt(tr_set,:)',[flen sdim sdim n_tr_cells]);
    gabor_filt_bank2 = reshape(gabor_emp2_filt(tr_set,:)',[flen sdim sdim n_tr_cells]);
    shifted_gabor_bank1 = nan(sdim^2*flen,n_tr_cells);
    shifted_gabor_bank2 = nan(sdim^2*flen,n_tr_cells);
    
    shift_cnt = 1;
    for xx = 1:length(x_shifts)
        for yy = 1:length(y_shifts)
            fprintf('Shift %d of %d\n',shift_cnt,n_shifts);
            d2 = dist_shift4d(gabor_filt_bank1,x_shifts(xx),3);
            d2 = dist_shift4d(d2,y_shifts(yy),2);
            shifted_gabor_bank1 = reshape(d2,sdim^2*flen,n_tr_cells);
            d2 = dist_shift4d(gabor_filt_bank2,x_shifts(xx),3);
            d2 = dist_shift4d(d2,y_shifts(yy),2);
            shifted_gabor_bank2 = reshape(d2,sdim^2*flen,n_tr_cells);
            
            gabor_outs1 = fullX*shifted_gabor_bank1;
            gabor_outs2 = fullX*shifted_gabor_bank2;
            
            gfun = gabor_outs1.^2 + gabor_outs2.^2;
            gfun = bsxfun(@plus,gfun,tv_theta_params(full_expt_inds,tr_set));
            
            too_large = gfun > 50;
            pred_rate = log(1+exp(gfun));
            pred_rate(too_large) = gfun(too_large);
            pred_rate(pred_rate < 1e-20) = 1e-20;
            
            LLs = Robs.*log(pred_rate) - pred_rate;
            frame_LLs(:,shift_cnt) = sum(LLs,2);
            shift_cnt = shift_cnt + 1;
        end
    end
    
    %% HMM for inferring sequence of stimulus translations
    fprintf('Processing Chunked Data...\n');
    chunk_assignments = [];
    lB = [];
    chunk_labels = [];
    chunk_fix_nums = [];
    chunk_eyespeed = [];
    chunk_cnt = 0;
    for cf = 1:n_fixs
        cur_im_nums = rel_fix_start_inds(cf):rel_fix_stop_inds(cf);
        n_chunks = ceil(length(cur_im_nums)/chunk_dur);
        chunk_starts = (0:n_chunks-1)*chunk_dur + 1;
        chunk_stops = chunk_starts + chunk_dur-1;
        chunk_stops(chunk_stops > length(cur_im_nums)) = length(cur_im_nums);
        
        temp_chunk_assignments = nan(length(cur_im_nums),1);
        temp_lB = nan(n_chunks,n_shifts);
        temp_chunk_eyespeeds = nan(n_chunks,1);
        for i = 1:n_chunks
            temp_chunk_assignments(chunk_starts(i):chunk_stops(i)) = i;
            temp_lB(i,:) = sum(frame_LLs(cur_im_nums(chunk_starts(i):chunk_stops(i)),:),1);
            temp_chunk_eyespeeds(i) = mean(avg_interp_eyespeed(cur_im_nums(chunk_starts(i):chunk_stops(i))));
        end
        chunk_assignments = [chunk_assignments; temp_chunk_assignments+chunk_cnt];
        chunk_fix_nums = [chunk_fix_nums; cf*ones(n_chunks,1)];
        chunk_eyespeed = [chunk_eyespeed; temp_chunk_eyespeeds];
        lB = [lB; temp_lB];
        chunk_labels = [chunk_labels; 1; zeros(n_chunks-1,1)];
        chunk_cnt = chunk_cnt + n_chunks;
    end
    
    tot_n_chunks = size(lB,1);
    
    lalpha=zeros(tot_n_chunks,n_shifts);
    lbeta = zeros(tot_n_chunks,n_shifts);
    lscale=zeros(tot_n_chunks,1); %initialize rescaling parameters
    %compute rescaled forward messages
    lalpha(1,:) = leps_prior' + lB(1,:);
    lscale(1)=logsumexp(lalpha(1,:));
    lalpha(1,:) = lalpha(1,:) - lscale(1);
    for t=2:tot_n_chunks
        fprintf('%d of %d\n',t,tot_n_chunks);
        if chunk_labels(t)==0
            cur_lA = -cdist.^2/(2*sac_deps_sigmas(t-1)^2);
                        cur_lA = bsxfun(@plus,cur_lA,leps_prior'); %factor in constant prior
            cur_lA = bsxfun(@minus,cur_lA,logsumexp(cur_lA,2)); %normalize
        elseif chunk_labels(t)==1
            cur_lA = lA_tflip;
        end
        lalpha(t,:) = logmulexp(lalpha(t-1,:),cur_lA) + lB(t,:);
        lscale(t) = logsumexp(lalpha(t,:));
        lalpha(t,:)= lalpha(t,:) - lscale(t);
    end
    
    %compute rescaled backward messages
    lbeta(tot_n_chunks,:)=log(ones(1,n_shifts)) - lscale(tot_n_chunks);
    for t=tot_n_chunks-1:-1:1
        fprintf('%d\n',t);
        if chunk_labels(t+1)==0
            cur_lA = -cdist.^2/(2*sac_deps_sigmas(t)^2);
                        cur_lA = bsxfun(@plus,cur_lA,leps_prior'); %factor in constant prior
            cur_lA = bsxfun(@minus,cur_lA,logsumexp(cur_lA,2)); %normalize
        elseif chunk_labels(t+1)==1
            cur_lA = lA_tflip;
        end
        lf1 = lbeta(t+1,:) + lB(t+1,:);
        lbeta(t,:) = logmulexp(lf1,cur_lA') - lscale(t);
    end
    
    %compute posteriors over hidden states
    lgamma= lalpha + lbeta;
    lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
    
    ov_lgamma = lgamma(chunk_assignments,:);
    
    %% RECONSTRUCT MAP STIMULUS
    [max_post,max_loc] = max(ov_lgamma,[],2);
    x_cor{it} = SH(max_loc,1);
    y_cor{it} = SH(max_loc,2);
    
    %% RECONSTRUCT NEW STIMULUS MATRIX
    resh_X = reshape(fullX',[flen sdim sdim NT]);
    resh_X_sh = zeros(size(resh_X));
    for ii = 1:NT
        %     if mod(ii,100)==0 fprintf('%d of %d\n',ii,NT); end
        d2 = dist_shift3d(resh_X(:,:,:,ii), -x_cor{it}(ii), 3);
        d2 = dist_shift3d(d2,-y_cor{it}(ii),2);
        resh_X_sh(:,:,:,ii) = d2;
    end
    fullX_sh = reshape(resh_X_sh,sdim^2*flen,NT)';
    
    %%
    clear resh_X resh_X_sh
    clear gabor_emp*
    for t = 1:96
        fprintf('Fitting cell %d of %d\n',t,96);
        gabor_emp1 = get_pgabor_mask_v2(XX,YY,gabor_params(t,1:6),0);
        gabor_emp2 = get_pgabor_mask_v2(XX,YY,gabor_params(t,1:6),pi/2);
        
        gabor_emp1_vec = [gabor_emp1(:)'; zeros(flen-1,sdim^2)];
        gabor_emp1_mat = makeStimRows(gabor_emp1_vec,flen);
        
        gabor_emp2_vec = [gabor_emp2(:)'; zeros(flen-1,sdim^2)];
        gabor_emp2_mat = makeStimRows(gabor_emp2_vec,flen);
        
        gabor_outs1 = fullX_sh*gabor_emp1_mat';
        gabor_outs2 = fullX_sh*gabor_emp2_mat';
        
        spikebins = convert_to_spikebins(full_binned_spks(:,t));
        
        X = [gabor_outs1 gabor_outs2];
        [NT,klen] = size(X);
        k0 = zeros(klen,1);
        
        init_signs = ones(cur_ndims,1);
        init_kerns = 0.05*randn(klen,cur_ndims);
        for i = 1:cur_ndims
            kern_types{i} = 'quad';
        end
        init_linK = zeros(n_un_expts,1);
        temp_gnm = createGNM_v2(init_kerns,init_signs,kern_types,init_linK,defmod,stim_params);
        fit_gnm{it}(t) = fitGNM_filters_v2(temp_gnm,X,linX,spikebins,'none',[],1e-4,1e-6);
        
        conv_mat = [gabor_emp1_mat' gabor_emp2_mat'];
        kmat = get_k_mat(fit_gnm{it}(t));
        
        k_pix = conv_mat*kmat;
        gabor_emp1_filt(t,:) = k_pix(:,1)';
        gabor_emp2_filt(t,:) = k_pix(:,2)';
    end
    
    tv_theta_params = [fit_gnm{it}(:).linK];
    ov_theta = [fit_gnm{it}(:).spk_theta];
    tv_theta_params = bsxfun(@plus,tv_theta_params,ov_theta);
    
end

