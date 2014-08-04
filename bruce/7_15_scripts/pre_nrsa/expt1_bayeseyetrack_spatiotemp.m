%% Load Data
clear all;
addpath(genpath('~/Code/James_scripts'));
cd ~/Data/bruce/7_15_12
cd G029/

load ./CellList.mat
single_units = find(CellList(1,:,1) > 0);

cd ~/Data/bruce/7_15_12/G029/
load ./Expt1_newcompiled_data_flen6_d1p25_new.mat
fullX = fullX/std(fullX(:));

load ./eye_calibration_data

Pix2Deg = 0.018837;
[NT,klen] = size(fullX);

%% crop stimulus for the purpose of faster gabor function fitting
sdim = length(xpatch_inds);
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));

% SET UP XV CELL SET
NSIG = 96;
xv_frac = 0.2;
tr_set = randperm(NSIG);
tr_set = tr_set(1:round(length(tr_set)*(1-xv_frac)));
xv_set = setdiff(1:NSIG,tr_set);
n_tr_cells = length(tr_set);

% PARSE DATA INTO FIXATIONS
diff_used_inds = diff(used_inds);
rel_fix_start_inds = [1; 1+find(diff_used_inds > 1)];
rel_fix_stop_inds = [(find(diff_used_inds > 1)); NT];
n_fixs = length(rel_fix_start_inds);

%% COMPUTE OUTPUTS OF GABOR MODELS BASED ON RF MAPPING DATA
flen = 6;
tot_dims = flen*sdim^2;

%initialize model
defmod.lambda_dT = 0;
defmod.lambda_d2X = 0;
defmod.lambda_L1x = 0;

stim_params.flen = flen;
stim_params.spatial_dims = 1;
stim_params.sdim = 2;
cur_ndims = 2;

load ./expt1_eyecor_d1p25_nosac_v4 gabor*
gabor_params = gabor_params_f{end};
% load ./spatiotemporal_mods

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
    
    X = [gabor_outs1 gabor_outs2];
    [NT,klen] = size(X);
    k0 = zeros(klen,1);
    
    init_signs = ones(cur_ndims,1);
    init_kerns = 0.05*randn(klen,cur_ndims);
    for i = 1:cur_ndims
        kern_types{i} = 'quad';
    end
    temp_gnm = createGNM(init_kerns,init_signs,kern_types,defmod,stim_params);
    init_gnm(t) = fitGNM_filters(temp_gnm,X,spikebins,'none',[],1e-4,1e-6);
    init_gnm_spk_theta(t) = init_gnm(t).spk_theta;
    
    conv_mat = [gabor_emp1_mat' gabor_emp2_mat'];
    kmat = get_k_mat(init_gnm(t));
    
    k_pix = conv_mat*kmat;
    gabor_emp1_filt(t,:) = k_pix(:,1)';
    gabor_emp2_filt(t,:) = k_pix(:,2)';
end
gnm_spk_theta = [init_gnm(:).spk_theta];


%% ESTIMATE LL for each shift in each stimulus frame
n_iter = 5;
for it = 1:n_iter
    
    max_shift = 20;
    dshift = 1;
    x_shifts = -max_shift:dshift:max_shift;
    y_shifts = -max_shift:dshift:max_shift;
    [Xsh,Ysh] = meshgrid(x_shifts,y_shifts);
    SH = [Xsh(:) Ysh(:)];
    n_shifts = size(SH,1);
    
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
            gfun = bsxfun(@plus,gfun,gnm_spk_theta(tr_set));
            
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
    chunk_dur = 1;
    
    %overall prior on shifts
    eps_prior_sigma = 0.3; %0.2
    leps_prior = -sum((SH/Fsd).^2,2)/(2*eps_prior_sigma^2);
    leps_prior = bsxfun(@minus,leps_prior,logsumexp(leps_prior)); %normalize
    
    %state transition matrix (includes a 'constant' prior)
    deps_sigma = 0.01*chunk_dur; %0.06
    cdist = squareform(pdist(SH/Fsd));
    lA = -cdist.^2/(2*deps_sigma^2);
    lA = bsxfun(@plus,lA,leps_prior'); %factor in constant prior
    lA = bsxfun(@minus,lA,logsumexp(lA,2)); %normalize
    
    ov_lgamma = nan(NT,n_shifts);
    
    for cf = 1:n_fixs
        fprintf('Fixation %d of %d\n',cf,n_fixs);
        cur_im_nums = rel_fix_start_inds(cf):rel_fix_stop_inds(cf);
        
        n_chunks = ceil(length(cur_im_nums)/chunk_dur);
        
        chunk_starts = (0:n_chunks-1)*chunk_dur + 1;
        chunk_stops = chunk_starts + chunk_dur;
        chunk_stops(chunk_stops > length(cur_im_nums)) = length(cur_im_nums);
        
        chunk_assignments = nan(length(cur_im_nums),1);
        lB = nan(n_chunks,n_shifts);
        for i = 1:n_chunks
            chunk_assignments(chunk_starts(i):chunk_stops(i)) = i;
            lB(i,:) = sum(frame_LLs(cur_im_nums(chunk_starts(i):chunk_stops(i)),:));
        end
        
        lalpha=zeros(n_chunks,n_shifts);
        lbeta = zeros(n_chunks,n_shifts);
        lscale=zeros(n_chunks,1); %initialize rescaling parameters
        %compute rescaled forward messages
        lalpha(1,:) = leps_prior' + lB(1,:);
        lscale(1)=logsumexp(lalpha(1,:));
        lalpha(1,:) = lalpha(1,:) - lscale(1);
        for t=2:n_chunks
            lalpha(t,:) = logmulexp(lalpha(t-1,:),lA) + lB(t,:);
            lscale(t) = logsumexp(lalpha(t,:));
            lalpha(t,:)= lalpha(t,:) - lscale(t);
        end
        
        %compute rescaled backward messages
        lbeta(n_chunks,:)=log(ones(1,n_shifts)) - lscale(n_chunks);
        for t=n_chunks-1:-1:1
            lf1 = lbeta(t+1,:) + lB(t+1,:);
            lbeta(t,:) = logmulexp(lf1,lA') - lscale(t);
        end
        
        %compute posteriors over hidden states
        lgamma= lalpha + lbeta;
        lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
        
        ov_lgamma(cur_im_nums,:) = lgamma(chunk_assignments,:);
        
    end
    
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
        temp_gnm = createGNM(init_kerns,init_signs,kern_types,defmod,stim_params);
        fit_gnm{it}(t) = fitGNM_filters(temp_gnm,X,spikebins,'none',[],1e-4,1e-6);
        
        conv_mat = [gabor_emp1_mat' gabor_emp2_mat'];
        kmat = get_k_mat(fit_gnm{it}(t));
        
        k_pix = conv_mat*kmat;
        gabor_emp1_filt(t,:) = k_pix(:,1)';
        gabor_emp2_filt(t,:) = k_pix(:,2)';
    end
    
    gnm_spk_theta = [fit_gnm{it}(:).spk_theta];
    
    
end

%%
init_LL = [init_gnm(:).LL];
for it = 1:3
    fit_LL(it,:) = [fit_gnm{it}(:).LL];
end
delta_LL = bsxfun(@minus,fit_LL,init_LL);

figure
plot(1:3,mean(delta_LL(:,tr_set),2),'ko-')
yl = ylim();
ylim([yl(1) 0])
hold on
plot(1:3,mean(delta_LL(:,xv_set),2),'ro-')


%%
save spatiotemp_eyetrack_v2 x_cor y_cor full_t fit_gnm init_gnm






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
    ref_kmat = get_k_mat(fit_gnm{3}(t));
    old_kmat = get_k_mat(init_gnm(t));
    ref_k_pix = conv_mat*ref_kmat;
    old_k_pix = conv_mat*old_kmat;
    
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
    old_k_pix = old_k_pix(:,ord);
    ref_k_norm = ref_k_norm(ord);
    [~,ref_best_slice_ids] = get2dMaxSlices(ref_k_pix,flen,sdim,1:(sdim^2));
    [~,old_best_slice_ids] = get2dMaxSlices(old_k_pix,flen,sdim,1:(sdim^2));
    
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
    set(gcf,'Position',[200 800 800 600])
    figure
    for ii = 1:n_filts
        cur_filtmat = reshape(old_k_pix(:,ii),flen,fsdim);
        
        filt_projprof = project_2drf(cur_filtmat(:),avg_best_theta,sdim);
        
        subplot(n_filts,2,(ii-1)*2+1)
        imagesc(nxax,nyax,reshape(cur_filtmat(old_best_slice_ids(ii),:),sdim,sdim));colormap(gray(256));
        caxis([-zmax(ii) zmax(ii)]); %colorbar
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
%         xm = max(abs(filt_projprof(:)));
        caxis([-xm(ii) xm(ii)]/1.25);
        title('Space-time projection')
        
    end
    set(gcf,'Position',[1000 800 800 600])
    
%     disp('LL')
%     fprintf('Ref: %.3f\n',ref_LL(t));
%     fprintf('Old: %.3f\n',old_LL(t));
    
    pause
    close all
end

%%
load ./all_eyedata_expt1_v2
Fs_t = 1/(all_t(2)-all_t(1));
lcf = 1/10;
[fb,fa] = butter(2,lcf/(Fs_t/2),'high');
eyepos_f = filtfilt(fb,fa,all_eyepos);

interp_eyepos_f = interp1(all_t,eyepos_f,full_t);


