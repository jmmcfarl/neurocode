%% Load Data
clear all;
addpath(genpath('~/Code/James_scripts'));
cd ~/Data/bruce/7_15_12
cd G029/

load ./CellList.mat
single_units = find(CellList(1,:,1) > 0);

cd ~/Data/bruce/7_15_12/G029/
load ./Expt3_newcompiled_data_d1p5
resh_all_stims = resh_all_stims/std(resh_all_stims(:));

load ./all_eyedata_expt3
Pix2Deg = 0.018837;
NT = length(full_stim_ids);

%% crop stimulus for the purpose of faster gabor function fitting
sdim = length(xpatch_inds);
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));

% PARSE DATA INTO FIXATIONS
diff_used_inds = [1; diff(used_inds)];
rel_fix_start_inds = [1; find(diff_used_inds > 1)];
rel_fix_stop_inds = [(find(diff_used_inds > 1)-1); NT];
n_trials = length(rel_fix_start_inds);

%% RECREATE TIME-EMBEDDED STIM MAT
flen = 6;
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

%%
NT = length(nused_inds);
xv_frac = 0.2;
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

%% COMPUTE OUTPUTS OF GABOR MODELS BASED ON RF MAPPING DATA
load ./expt1_eyecor_d1p25_nosac_v2.mat gabor*
gabor_params = gabor_params_f{end};

load ./spatiotemporal_mods

clear gabor_emp*
%initialize model

defmod.lambda_dT = 0;
defmod.lambda_d2X = 0;
defmod.lambda_L1x = 0;

defmod2.lambda_L2x = 50;
defmod2.lambda_dT = 0;
defmod2.lambda_d2X = 0;
defmod2.lambda_L1x = 0;

spatial_stim_params.flen = flen;
spatial_stim_params.spatial_dims = 0;
spatial_stim_params.sdim = 1;

stim_params.flen = flen;
stim_params.spatial_dims = 1;
stim_params.sdim = 2;

for t = 1:96
    fprintf('Cell %d of %d\n',t,96);
    tr_spikebins = convert_to_spikebins(full_binned_spks(tr_inds,t));
    xv_spikebins = convert_to_spikebins(full_binned_spks(xv_inds,t));

    gabor_emp1 = get_pgabor_mask_v2(XX,YY,gabor_params(t,1:6),0);
    gabor_emp2 = get_pgabor_mask_v2(XX,YY,gabor_params(t,1:6),pi/2);
    gabor_emp1_vec = [gabor_emp1(:)'; zeros(flen-1,sdim^2)];
    gabor_emp1_mat = makeStimRows(gabor_emp1_vec,flen);
    
    gabor_emp2_vec = [gabor_emp2(:)'; zeros(flen-1,sdim^2)];
    gabor_emp2_mat = makeStimRows(gabor_emp2_vec,flen);
    
    gabor_outs1 = fullX*gabor_emp1_mat';
    gabor_outs2 = fullX*gabor_emp2_mat';
        
%     X2 = gabor_outs1.^2 + gabor_outs2.^2;
%     [NT,klen] = size(X2);
%     k0 = zeros(klen,1);
%     cur_ndims = 1;
%     init_signs = ones(cur_ndims,1);
%     init_kerns = 0.05*randn(klen,cur_ndims);
%     for i = 1:cur_ndims
%         kern_types{i} = 'lin';
%     end
%     spatial_gnm(t) = createGNM(init_kerns,init_signs,kern_types,defmod,spatial_stim_params);
%     spatial_gnm(t) = fitGNM_filters(spatial_gnm(t),X2(tr_inds,:),tr_spikebins,'none',[],1e-4,1e-6);
%     spatial_xvLL(t) = getLL_GNM(spatial_gnm(t),X2(xv_inds,:),xv_spikebins,'none');    
        
    X = [gabor_outs1 gabor_outs2];
    [NT,klen] = size(X);
    k0 = zeros(klen,1);
    
    cur_ndims = 2;
    init_signs = ones(cur_ndims,1);
    init_kerns = 0.05*randn(klen,cur_ndims);
    for i = 1:cur_ndims
        kern_types{i} = 'quad';
    end
    st_gnm(t) = createGNM(init_kerns,init_signs,kern_types,defmod,stim_params);
    st_gnm(t) = fitGNM_filters(st_gnm(t),X(tr_inds,:),tr_spikebins,'none',[],1e-4,1e-6);
    st_xvLL(t) = getLL_GNM(st_gnm(t),X(xv_inds,:),xv_spikebins,'none');    
       
    st_gnm2(t) = createGNM(init_kerns,init_signs,kern_types,defmod2,stim_params);
    st_gnm2(t) = fitGNM_filters(st_gnm2(t),X(tr_inds,:),tr_spikebins,'none',[],1e-4,1e-6);
    st_xvLL2(t) = getLL_GNM(st_gnm2(t),X(xv_inds,:),xv_spikebins,'none');    

    [nll, pnll, lpen, prate, g, int_g] = getLL_GNM(gnm(t),X(tr_inds,:),tr_spikebins,'none');
    old_gnm(t) = fitGNM_spkNL(gnm(t),g,tr_spikebins,0);
    old_xvLL(t) = getLL_GNM(old_gnm(t),X(xv_inds,:),xv_spikebins,'none');    
 
%     [nll, pnll, lpen, prate_st] = getLL_GNM(st_gnm(t),X(tr_inds,:),tr_spikebins,'none');
%     [nll, pnll, lpen, prate_old] = getLL_GNM(old_gnm(t),X(tr_inds,:),tr_spikebins,'none');
%     [nll, pnll, lpen, prate_old2] = getLL_GNM(old_gnm2(t),X(tr_inds,:),tr_spikebins,'none');
    
    conv_mat = [gabor_emp1_mat' gabor_emp2_mat'];
    kmat = get_k_mat(old_gnm(t));
    k_pix = conv_mat*kmat;
    gabor_emp1_filt(t,:) = k_pix(:,1)';
    gabor_emp2_filt(t,:) = k_pix(:,2)';
end

%%
% SET UP XV CELL SET
NSIG = 96;
xv_frac = 0.2;
tr_set = randperm(NSIG);
tr_set = tr_set(1:round(length(tr_set)*(1-xv_frac)));
xv_set = setdiff(1:NSIG,tr_set);
n_tr_cells = length(tr_set);

init_gnm_spk_theta = [old_gnm(:).spk_theta];
%% ESTIMATE LL for each shift in each stimulus frame
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
        gfun = bsxfun(@plus,gfun,init_gnm_spk_theta(tr_set));
        
        too_large = gfun > 50;
        pred_rate = log(1+exp(gfun));
        pred_rate(too_large) = gfun(too_large);
        pred_rate(pred_rate < 1e-20) = 1e-20;
        
        LLs = Robs.*log(pred_rate) - pred_rate;
        frame_LLs(:,shift_cnt) = sum(LLs,2);
        shift_cnt = shift_cnt + 1;
    end
end


%%
load ./spatiotemporal_mods
fsdim = sdim^2;
nxax = xax(xpatch_inds);
nyax = yax(ypatch_inds);
for t = single_units
    fprintf('Fitting GEM: Cell %d of %d\n',t,96);
    gabor_emp1 = get_pgabor_mask_v2(XX,YY,gabor_params(t,1:6),0);
    gabor_emp2 = get_pgabor_mask_v2(XX,YY,gabor_params(t,1:6),pi/2);
    
    gabor_emp1_vec = [gabor_emp1(:)'; zeros(flen-1,sdim^2)];
    gabor_emp1_mat = makeStimRows(gabor_emp1_vec,flen);
    
    gabor_emp2_vec = [gabor_emp2(:)'; zeros(flen-1,sdim^2)];
    gabor_emp2_mat = makeStimRows(gabor_emp2_vec,flen);
    
    gabor_outs1 = fullX*gabor_emp1_mat';
    gabor_outs2 = fullX*gabor_emp2_mat';
    X = [gabor_outs1 gabor_outs2];
    X2 = gabor_outs1.^2 + gabor_outs2.^2;
    spikebins = convert_to_spikebins(full_binned_spks(:,t));

    conv_mat = [gabor_emp1_mat' gabor_emp2_mat'];
    ref_kmat = get_k_mat(st_gnm(t));
    old_kmat = get_k_mat(gnm(t));
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
    
    [nll, pnll, lpen, prate_st] = getLL_GNM(st_gnm(t),X,spikebins,'none');
    [nll, pnll, lpen, prate_old] = getLL_GNM(old_gnm(t),X,spikebins,'none');
    [nll, pnll, lpen, prate_spatial] = getLL_GNM(spatial_gnm(t),X2,spikebins,'none');
    
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
        caxis([-xm(ii) xm(ii)]/1.5);
        title('Space-time projection')
        
    end
    set(gcf,'Position',[200 200 800 600])
    figure
    for ii = 1:n_filts
        cur_filtmat = reshape(old_k_pix(:,ii),flen,fsdim);
        
        filt_projprof = project_2drf(cur_filtmat(:),avg_best_theta,sdim);
        
        subplot(n_filts,2,(ii-1)*2+1)
        imagesc(nxax,nyax,reshape(cur_filtmat(old_best_slice_ids(ii),:),sdim,sdim));colormap(gray(256));
        zmax(ii) = max(abs(cur_filtmat(old_best_slice_ids(ii),:))); caxis([-zmax(ii) zmax(ii)]); %colorbar
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
        caxis([-xm(ii) xm(ii)]/1.5);
        title('Space-time projection')
        
    end
    set(gcf,'Position',[1000 200 800 600])
    
%     disp('LL')
%     fprintf('Ref: %.3f\n',ref_LL(t));
%     fprintf('Old: %.3f\n',old_LL(t));
    
    figure
    plot(prate_st)
    hold on
    plot(prate_old,'r')
    plot(prate_spatial,'k')
%     plot(smooth(full_binned_spks(:,t),3),'g')
    
    pause
    close all
end

%%
fullX = resh_all_stims(full_stim_ids,:);
gabor_outs1 = fullX*gabor_emp1_filt';
gabor_outs2 = fullX*gabor_emp2_filt';
gabor_outs = sqrt(gabor_outs1.^2 + gabor_outs2.^2);

%%
% use_inds = find(full_xo_vec ~= 0);
use_inds = 1:length(full_xo_vec);
for t = 1:96
    [beta(t,:),dev(t),stats(t)] = glmfit(gabor_outs(use_inds,t),full_binned_spks(use_inds,t),'poisson');
end

%%
resp_win = 0.1;
n_resp_bins = ceil(resp_win/dt);
stim_start_inds = [1 1+find(diff(full_stim_ids) ~= 0)];
n_stims = length(stim_start_inds);
stim_spks = nan(n_stims,96);
stim_ids = nan(n_stims,1);
for i = 1:n_stims
    cur_inds = stim_start_inds(i):(stim_start_inds(i)+n_resp_bins);
    stim_spks(i,:) = sum(full_binned_spks(cur_inds,:));
    stim_ids(i) = full_stim_ids(stim_start_inds(i));
end

gabor_outs1 = fullX(stim_start_inds,:)*gabor_emp1_filt';
gabor_outs2 = fullX(stim_start_inds,:)*gabor_emp2_filt';
gabor_outs = sqrt(gabor_outs1.^2 + gabor_outs2.^2);

flen = 25;
for t = 1:96
    X = makeStimRows(gabor_outs(:,t),flen,0);
    [beta(t,:),dev(t),stats(t)] = glmfit(gabor_outs(use_inds,t),stim_spks(use_inds,t),'poisson');
end

%%
stim_impulse = zeros(length(full_stim_ids),96);
stim_impulse(stim_start_inds,:) = gabor_outs;
flen = 25;
sm_pen = 100;
clear beta
for t = 1:96
    X = makeStimRows(stim_impulse(:,t),flen,0);
    %    [beta(t,:),dev(t),stats(t)] = glmfit(X,full_binned_spks(:,t),'poisson');
    spikebins = convert_to_spikebins(full_binned_spks(:,t));
    [fitp,grad] = GLMsolve_jmm(X, spikebins, [0; 0], 1, [], [], [sm_pen 1 flen], [], [], [], 0);
    beta(t,:) = fitp.k(1:end-1);
end
