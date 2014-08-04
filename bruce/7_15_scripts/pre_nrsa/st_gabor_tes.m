%% Load Data
clear all;
addpath(genpath('~/Code/James_scripts'));
cd ~/Data/bruce/7_15_12
cd G029/

load ./CellList.mat
single_units = find(CellList(1,:,1) > 0);

cd ~/Data/bruce/7_15_12/G029/
load ./Expt1_newcompiled_data_flen5_d1p25_new.mat
fullX = fullX/std(fullX(:));
flen = 5;

Pix2Deg = 0.018837;
[NT,klen] = size(fullX);

%% crop stimulus for the purpose of faster gabor function fitting
sdim = length(xpatch_inds);
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));


% PARSE DATA INTO FIXATIONS
diff_used_inds = diff(used_inds);
rel_fix_start_inds = [1; 1+find(diff_used_inds > 1)];
rel_fix_stop_inds = [(find(diff_used_inds > 1)); NT];
n_fixs = length(rel_fix_start_inds);

%% RECONSTRUCT NEW STIMULUS MATRIX
cur_full_t = full_t;
load ./gabor_tracking_varmeans 
y_cor_vals = y_cor{end};
x_cor_vals = x_cor{end};
load ./Expt1_newcompiled_data_fixeddelay_d1p25_new.mat full_t %load old t-axis
y_cor_interp = round(interp1(full_t,y_cor_vals,cur_full_t));
x_cor_interp = round(interp1(full_t,x_cor_vals,cur_full_t));
x_cor_interp(isnan(x_cor_interp)) = 0;
y_cor_interp(isnan(y_cor_interp)) = 0;
full_t = cur_full_t;

gabor_params = gabor_params_f{end};

resh_X = reshape(fullX',[flen sdim sdim NT]);
resh_X_sh = zeros(size(resh_X));
for ii = 1:NT
    %     if mod(ii,100)==0 fprintf('%d of %d\n',ii,NT); end
    d2 = dist_shift3d(resh_X(:,:,:,ii), -x_cor_interp(ii), 3);
    d2 = dist_shift3d(d2,-y_cor_interp(ii),2);
    resh_X_sh(:,:,:,ii) = d2;
end
fullX_sh = reshape(resh_X_sh,sdim^2*flen,NT)';

%%
tot_dims = flen*sdim^2;

[un_expts,~,full_expt_inds] = unique(full_expt_vec); 
n_un_expts = length(un_expts);
linX = zeros(NT,n_un_expts);
for i = 1:n_un_expts
    linX(full_expt_inds==i,i) = 1;
end

idier = [full_expt_vec full_trial_vec];
unique_ids = unique(idier,'rows');
n_trials = size(unique_ids,1);
xv_frac =0;
xv_num = round(xv_frac*n_trials);
xv_set = randperm(n_trials);
xv_set(xv_num+1:end) = [];
tr_set = setdiff(1:n_trials,xv_set);

xv_inds = [];
for ii = 1:xv_num
   cur_inds = find(full_expt_vec == unique_ids(xv_set(ii),1) & full_trial_vec == unique_ids(xv_set(ii),2));
   xv_inds = [xv_inds; cur_inds];
end
tr_inds = setdiff(1:NT,xv_inds);

%% COMPUTE OUTPUTS OF GABOR MODELS BASED ON RF MAPPING DATA

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
    
    gabor_outs1 = fullX_sh*gabor_emp1_mat';
    gabor_outs2 = fullX_sh*gabor_emp2_mat';
    
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

%     init_gnm2(t) = fitGNM_filters_v2(temp_gnm,X(tr_inds,:),linX(tr_inds,:),tr_spikebins,'none',[],1e-4,1e-6,0);
%     xvLL_gnm2(t) = getLL_GNM_v2(init_gnm2(t),X(xv_inds,:),linX(xv_inds,:),xv_spikebins,'none');

    conv_mat = [gabor_emp1_mat' gabor_emp2_mat'];
    kmat = get_k_mat(init_gnm(t));
    
    k_pix = conv_mat*kmat;
    gabor_emp1_filt(t,:) = k_pix(:,1)';
    gabor_emp2_filt(t,:) = k_pix(:,2)';
    
    
%     % FIT SPATIAL MODEL
%     cur_gabor_emp1 = gabor_emp1_mat(3,:);
%     cur_gabor_emp2 = gabor_emp2_mat(3,:);
%     cur_gabor_outs1 = fullX_sh*cur_gabor_emp1';
%     cur_gabor_outs2 = fullX_sh*cur_gabor_emp2';
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

load ./all_eyedata_expt1_v2
orig_eyedt = all_t(2)-all_t(1);
interp_eyespeeds = interp1(all_t,all_eyespeed,full_t);
avg_interp_eyespeed = dt/orig_eyedt*mean(interp_eyespeeds,2);

max_shift = 22;
dshift = 1;
x_shifts = -max_shift:dshift:max_shift;
y_shifts = -max_shift:dshift:max_shift;
[Xsh,Ysh] = meshgrid(x_shifts,y_shifts);
SH = [Xsh(:) Ysh(:)];
n_shifts = size(SH,1);

chunk_dur = 1;

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

n_iter = 2;
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
%             cur_lA = bsxfun(@plus,cur_lA,leps_prior'); %factor in constant prior
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
%             cur_lA = bsxfun(@plus,cur_lA,leps_prior'); %factor in constant prior
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

%%
init_LL = [init_gnm(:).LL];
for it = 1:n_iter
    fit_LL(it,:) = [fit_gnm{it}(:).LL];
end
delta_LL = bsxfun(@minus,fit_LL,init_LL);

figure
errorbar(1:n_iter,mean(delta_LL(:,tr_set),2),std(delta_LL(:,tr_set),[],2)/sqrt(length(tr_set)),'ko-')
yl = ylim();
% ylim([yl(1) 0])
hold on
errorbar(1:n_iter,mean(delta_LL(:,xv_set),2),std(delta_LL(:,xv_set),[],2)/sqrt(length(xv_set)),'ro-')
xlim([0 n_iter+1])

%%
% save spatiotemp_eyetrack_v2 x_cor y_cor full_t fit_gnm init_gnm






%%
fsdim = sdim^2;
nxax = xax(xpatch_inds);
nyax = yax(ypatch_inds);
for t = 1:96
    fprintf('Fitting GEM: Cell %d of %d\n',t,96);
    gabor_emp1 = get_pgabor_mask_v2(XX,YY,gabor_params(t,1:6),0);
    gabor_emp2 = get_pgabor_mask_v2(XX,YY,gabor_params(t,1:6),pi/2);
    
    gabor_emp1_vec = [gabor_emp1(:)'; zeros(flen-1,sdim^2)];
    gabor_emp1_mat = makeStimRows(gabor_emp1_vec,flen);
    
    gabor_emp2_vec = [gabor_emp2(:)'; zeros(flen-1,sdim^2)];
    gabor_emp2_mat = makeStimRows(gabor_emp2_vec,flen);
    
    conv_mat = [gabor_emp1_mat' gabor_emp2_mat'];
    ref_kmat = get_k_mat(init_gnm(t));
%     old_kmat = get_k_mat(old_init_gnm(t));
    ref_k_pix = conv_mat*ref_kmat;
%     old_k_pix = conv_mat*old_kmat;
    
    avg_best_theta = gabor_params(t,3)+pi/2;
    avg_best_x = gabor_params(t,1);
    avg_best_y = gabor_params(t,2);
    
    n_filts = size(ref_k_pix,2);
    ref_k_norm = zeros(n_filts,1);
    for ii = 1:n_filts
        ref_k_norm(ii) = norm(ref_k_pix(:,ii));
    end
    [~,ord] = sort(ref_k_norm,'descend');
    ref_k_pix = ref_k_pix(:,ord);
%     old_k_pix = old_k_pix(:,ord);
    ref_k_norm = ref_k_norm(ord);
    [~,ref_best_slice_ids] = get2dMaxSlices(ref_k_pix,flen,sdim,1:(sdim^2));
%     [~,old_best_slice_ids] = get2dMaxSlices(old_k_pix,flen,sdim,1:(sdim^2));
    
    figure
    for ii = 1:n_filts
        cur_filtmat = reshape(ref_k_pix(:,ii),flen,fsdim);
        
        filt_projprof = project_2drf(cur_filtmat(:),avg_best_theta,sdim);
        
        subplot(n_filts,2,(ii-1)*2+1)
        
        imagesc(nxax,nyax,reshape(cur_filtmat(ref_best_slice_ids(ii),:),sdim,sdim));colormap(gray(256));
        zmax(ii) = max(abs(cur_filtmat(ref_best_slice_ids(ii),:))); caxis([-zmax(ii) zmax(ii)]); %colorbar
        title(sprintf('Norm: %.3f',ref_k_norm(ii)))
        set(gca,'ydir','normal')
        hold on
        plot(avg_best_x,avg_best_y,'wo','linewidth',2)
        cur_linex = linspace(nxax(1),nxax(end),50);
        cur_liney = tan(avg_best_theta)*(cur_linex - avg_best_x)+avg_best_y;
        plot(cur_linex,cur_liney,'color','w')
        xlim(nxax([1 end]));ylim(nyax([1 end]));
        
        subplot(n_filts,2,(ii-1)*2+2)
        imagesc(1:sdim,-(0:(flen-1))*dt,filt_projprof);
        title(sprintf('Norm: %.3f',ref_k_norm(ii)))
        xm(ii) = max(abs(filt_projprof(:)));
        caxis([-xm(ii) xm(ii)]/1.25);
        title('Space-time projection')
        
    end
    set(gcf,'Position',[200 200 800 600])
%     figure
%     for ii = 1:n_filts
%         cur_filtmat = reshape(old_k_pix(:,ii),flen,fsdim);
%         
%         filt_projprof = project_2drf(cur_filtmat(:),avg_best_theta,sdim);
%         
%         subplot(n_filts,2,(ii-1)*2+1)
%         imagesc(nxax,nyax,reshape(cur_filtmat(old_best_slice_ids(ii),:),sdim,sdim));colormap(gray(256));
%         caxis([-zmax(ii) zmax(ii)]); %colorbar
%         title(sprintf('Norm: %.3f',ref_k_norm(ii)))
%         set(gca,'ydir','normal')
%         hold on
%         plot(avg_best_x,avg_best_y,'wo','linewidth',2)
%         cur_linex = linspace(nxax(1),nxax(end),50);
%         cur_liney = tan(avg_best_theta)*(cur_linex - avg_best_x)+avg_best_y;
%         plot(cur_linex,cur_liney,'color','w')
%         xlim(nxax([1 end]));ylim(nyax([1 end]));
%         
%         subplot(n_filts,2,(ii-1)*2+2)
%         imagesc(1:sdim,-(0:(flen-1))*dt,filt_projprof);
%         title(sprintf('Norm: %.3f',ref_k_norm(ii)))
% %         xm = max(abs(filt_projprof(:)));
%         caxis([-xm(ii) xm(ii)]/1.25);
%         title('Space-time projection')
%         
%     end
%     set(gcf,'Position',[1000 800 800 600])
    
%     disp('LL')
%     fprintf('Ref: %.3f\n',ref_LL(t));
%     fprintf('Old: %.3f\n',old_LL(t));
%     
    pause
    close all
end

%%
load ./all_eyedata_expt1_v2
Fs_t = 1/(all_t(2)-all_t(1));
lcf = 1/20;
[fb,fa] = butter(2,lcf/(Fs_t/2),'high');
eyepos_f = filtfilt(fb,fa,all_eyepos);

interp_eyepos_f = interp1(all_t,eyepos_f,full_t);


