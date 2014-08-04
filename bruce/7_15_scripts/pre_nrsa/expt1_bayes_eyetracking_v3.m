%% Load Data
clear all;
addpath(genpath('~/Code/James_scripts'));
cd ~/Data/bruce/7_15_12
cd G029/
% load ./eye_calibration_data
% load ./G029Expts.mat
% cd G034/

load ./CellList.mat
single_units = find(CellList(1,:,1) > 0);

cd ~/Data/bruce/7_15_12/G029/
load ./Expt1_newcompiled_data_fixedlag_d1p25_nosac.mat
fullX = fullX/std(fullX(:));

Pix2Deg = 0.018837;
[NT,klen] = size(fullX);

%% crop stimulus for the purpose of faster gabor function fitting
new_RF_patch = [-0.11 0.9; -0.9 0.1]; %location of RFs in degrees [x1 x2;y1 y2]
[curXX,curYY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));
curXX = curXX(:);
curYY = curYY(:);
new_crop = find(curXX >= new_RF_patch(1,1) & curXX <= new_RF_patch(1,2) & ...
    curYY >= new_RF_patch(2,1) & curYY <= new_RF_patch(2,2));
fullX_cropped = fullX(:,new_crop);

xpatch_inds_cropped = find(xax >= new_RF_patch(1,1) & xax <= new_RF_patch(1,2));
ypatch_inds_cropped = find(yax >= new_RF_patch(2,1) & yax <= new_RF_patch(2,2));

sdim = length(xpatch_inds);
[XXc,YYc] = meshgrid(xax(xpatch_inds_cropped),yax(ypatch_inds_cropped));
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));

cd ~/Data/bruce/7_15_12/
cd G029/
load ./grating_mu_data
gr_oris = unique_oris;
[~,sf_inds] = max(avg_mu_sf_profile,[],2);
pref_sfs = unique_sfs(sf_inds);
load ./oned_fixation_fits_v3.mat

% SET UP XV CELL SET
NSIG = 96;
xv_frac = 0.2;
tr_set = randperm(NSIG);
tr_set = tr_set(1:round(length(tr_set)*(1-xv_frac)));
xv_set = setdiff(1:NSIG,tr_set);
n_tr_cells = length(tr_set);

% PARSE DATA INTO FIXATIONS

diff_used_inds = diff(used_inds);
rel_fix_start_inds = [1; find(diff_used_inds > 1)];
rel_fix_stop_inds = [(find(diff_used_inds > 1)-1); NT];
n_fixs = length(rel_fix_start_inds);

%% COMPUTE OUTPUTS OF GABOR MODELS BASED ON RF MAPPING DATA

gab_priors(1).type = 'gauss';
gab_priors(1).theta(2) = 0.25; %prior std

gab_priors(2).type = 'gauss';
gab_priors(2).theta(2) = 0.25;

gab_priors(4).type = 'gam';
gab_priors(4).theta(1) = 8; %shape
gab_priors(4).theta(2) = 0.03; %scale

gab_priors(5).type = 'gam';
gab_priors(5).theta(1) = 8; %shape
gab_priors(5).theta(2) = 0.011; %scale

gab_priors(6).type = 'gam';
gab_priors(6).theta(1) = 2; %shape
gab_priors(6).theta(2) = 2; %scale

LB = [-0.1 -0.8 0 0.1 0.02 0.2 0 -Inf];
UB = [0.8 0.1 pi 0.4 0.25 6 Inf Inf];
hold_const = [0 0 1 0 0 0 0 0];
all_const = [1 1 1 1 1 1 0 0];

for t = 1:96
    fprintf('Fitting GEM: Cell %d of %d\n',t,96);
    
    gab_priors(1).theta(1) = mean_x(t);
    gab_priors(2).theta(1) = mean_y(t);

    init_params(1) = mean_x(t); %x0
    init_params(2) = mean_y(t); %y0
    init_params(3) = degtorad(pref_mu_ori(t)); %theta
    init_params(4) = 1/6; %lambda
    init_params(5) = 0.4*init_params(4); %sigma
    init_params(6) = 2; %eccentricity of gaussian ellipse
    init_params(7) = 0; %weight of quad term
    init_params(8) = 0; %const offset
    
    [gabor_params_f{1}(t,:),LL(1,t)] = fit_gabor_energy_mod(XXc,YYc,init_params,fullX_cropped,full_binned_spks(:,t),hold_const,LB,UB,gab_priors);
    %     [gabor_params_n(t,:),LLn(t)] = fit_gabor_energy_mod(XXc,YYc,init_params,fullX_cropped,full_binned_spks(:,t),all_const,LB,UB);
    
    gabor_emp1 = get_pgabor_mask_v2(XX,YY,gabor_params_f{1}(t,1:6),0);
    gabor_emp2 = get_pgabor_mask_v2(XX,YY,gabor_params_f{1}(t,1:6),pi/2);
    
    gabor_emp1_filt(t,:) = gabor_emp1(:);
    gabor_emp2_filt(t,:) = gabor_emp2(:);
    
end

%% SAVE THESE INITIAL GABOR MODELS
cd ~/Data/bruce/7_15_12/G029/
save gabor_initfits_d1p25_nodrift_nosac gabor_params_* LL gabor_emp*_filt gab_priors

%%
n_iter = 3;
for it = 1:n_iter
    
    %% ESTIMATE LL for each shift in each stimulus frame
    max_shift = 27;
    dshift = 1;
    x_shifts = -max_shift:dshift:max_shift;
    y_shifts = -max_shift:dshift:max_shift;
    [Xsh,Ysh] = meshgrid(x_shifts,y_shifts);
    SH = [Xsh(:) Ysh(:)];
    n_shifts = size(SH,1);
    
    frame_LLs = zeros(NT,n_shifts);
    Robs = full_binned_spks(:,tr_set);
    
    gabor_filt_bank1 = reshape(gabor_emp1_filt(tr_set,:)',[sdim sdim n_tr_cells]);
    gabor_filt_bank2 = reshape(gabor_emp2_filt(tr_set,:)',[sdim sdim n_tr_cells]);
    shifted_gabor_bank1 = nan(sdim^2,n_tr_cells);
    shifted_gabor_bank2 = nan(sdim^2,n_tr_cells);
    
    shift_cnt = 1;
    for xx = 1:length(x_shifts)
        for yy = 1:length(y_shifts)
            fprintf('Shift %d of %d\n',shift_cnt,n_shifts);
            d2 = dist_shift3d(gabor_filt_bank1,x_shifts(xx),2);
            d2 = dist_shift3d(d2,y_shifts(yy),1);
            shifted_gabor_bank1 = reshape(d2,sdim^2,n_tr_cells);
            d2 = dist_shift3d(gabor_filt_bank2,x_shifts(xx),2);
            d2 = dist_shift3d(d2,y_shifts(yy),1);
            shifted_gabor_bank2 = reshape(d2,sdim^2,n_tr_cells);
            
            gabor_outs1 = fullX*shifted_gabor_bank1;
            gabor_outs2 = fullX*shifted_gabor_bank2;
            energy_out = sqrt(gabor_outs1.^2 + gabor_outs2.^2);
            
            gfun = bsxfun(@times,energy_out,gabor_params_f{it}(tr_set,7)');
            gfun = bsxfun(@plus,gfun,gabor_params_f{it}(tr_set,8)');
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
    chunk_dur = 4;
    
    %overall prior on shifts
    eps_prior_sigma = 0.15; %0.2
    leps_prior = -sum((SH/Fsd).^2,2)/(2*eps_prior_sigma^2);
    leps_prior = bsxfun(@minus,leps_prior,logsumexp(leps_prior)); %normalize
    % cur_cent = [ov_eps(cur_fix,1) ov_eps(cur_fix,2)];
    % cdist = sum((bsxfun(@minus,SH,cur_cent)/Fsd).^2,2);
    % leps_prior = -cdist/(2*eps_prior_sigma^2);
    % leps_prior = bsxfun(@minus,leps_prior,logsumexp(leps_prior));
 
    eps_loose_prior_sigma = 0.4; %0.2
    leps_loose_prior = -sum((SH/Fsd).^2,2)/(2*eps_loose_prior_sigma^2);
    leps_loose_prior = bsxfun(@minus,leps_loose_prior,logsumexp(leps_loose_prior)); %normalize
    
    lA_tflip = repmat(leps_prior',n_shifts,1);
    
    %sig = (1-exp(-d/sc))*(max-min)+min
    deps_sigma_sac_min = 0.05;
    deps_sigma_sac_max = 0.4;
    dist_scale = 0.15;    
    
    %state transition matrix (includes a 'constant' prior)
    deps_sigma = 0.03; %0.06
    cdist = squareform(pdist(SH/Fsd));
    lA = -cdist.^2/(2*deps_sigma^2);
    lA = bsxfun(@plus,lA,leps_loose_prior'); %factor in constant prior
    lA = bsxfun(@minus,lA,logsumexp(lA,2)); %normalize
    
    ov_lgamma = nan(NT,n_shifts);
    
    fprintf('Processing Chunked Data...\n');
    chunk_assignments = [];
    lB = [];
    chunk_labels = [];
    chunk_fix_nums = [];
    chunk_cnt = 0;
    for cf = 1:n_fixs
        cur_im_nums = rel_fix_start_inds(cf):rel_fix_stop_inds(cf);
        
        n_chunks = ceil(length(cur_im_nums)/chunk_dur);
        
        chunk_starts = (0:n_chunks-1)*chunk_dur + 1;
        chunk_stops = chunk_starts + chunk_dur;
        chunk_stops(chunk_stops > length(cur_im_nums)) = length(cur_im_nums);
        
        temp_chunk_assignments = nan(length(cur_im_nums),1);
        temp_lB = nan(n_chunks,n_shifts);
        for i = 1:n_chunks
            temp_chunk_assignments(chunk_starts(i):chunk_stops(i)) = i;
            temp_lB(i,:) = sum(frame_LLs(cur_im_nums(chunk_starts(i):chunk_stops(i)),:));
        end
        chunk_assignments = [chunk_assignments; temp_chunk_assignments+chunk_cnt];
        chunk_fix_nums = [chunk_fix_nums; cf*ones(n_chunks,1)];
        lB = [lB; temp_lB];
        if isnan(eye_dx(cf))
            chunk_labels = [chunk_labels; 1; zeros(n_chunks-1,1)];
        else
            chunk_labels = [chunk_labels; 2; zeros(n_chunks-1,1)];
        end
        chunk_cnt = chunk_cnt + n_chunks;
    end
    
    tot_n_chunks = size(lB,1);
    avg_eye_dx = mean(eye_dx,2);
    avg_eye_dy = mean(eye_dy,2);
    eye_diff_d = sqrt((eye_dx(:,1)-eye_dx(:,2)).^2 + (eye_dy(:,1)-eye_dy(:,2)).^2);
    
    sac_deps_sigmas = (1-exp(-eye_diff_d/dist_scale))*(deps_sigma_sac_max-deps_sigma_sac_min) + deps_sigma_sac_min;

    
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
            cur_lA = lA;
        elseif chunk_labels(t)==1
            cur_lA = lA_tflip;
        else
            cdist = pdist2(SH/Fsd,bsxfun(@minus,SH/Fsd,[avg_eye_dx(chunk_fix_nums(t)) avg_eye_dy(chunk_fix_nums(t))]));
            cur_lA = -cdist.^2/(2*sac_deps_sigmas(chunk_fix_nums(t))^2);
            cur_lA = bsxfun(@plus,cur_lA,leps_loose_prior'); %factor in constant prior
            cur_lA = bsxfun(@minus,cur_lA,logsumexp(cur_lA,2)); %normalize
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
            cur_lA = lA;
        elseif chunk_labels(t+1)==1
            cur_lA = lA_tflip;
        else
            cdist = pdist2(SH/Fsd,bsxfun(@minus,SH/Fsd,[avg_eye_dx(chunk_fix_nums(t+1)) avg_eye_dy(chunk_fix_nums(t+1))]));
            cur_lA = -cdist.^2/(2*sac_deps_sigmas(chunk_fix_nums(t+1))^2^2);
            cur_lA = bsxfun(@plus,cur_lA,leps_loose_prior'); %factor in constant prior
            cur_lA = bsxfun(@minus,cur_lA,logsumexp(cur_lA,2)); %normalize
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
    
%     %% TO VISUALIZE EYE TRACKING
%     close all
%     figure(1)
%     set(gcf,'Position',[500 800 800 800])
%             set(gca,'ydir','normal');
%     figure(2)
%     subplot(2,1,1)
%     plot(x_cor{it}/Fsd)
%     hold on
%     plot(x_drifts_l,'k')
%     plot(x_drifts_r,'g')
%     ylim([-0.25 0.25])
%     subplot(2,1,2)
%     plot(y_cor{it}/Fsd)
%     hold on
%     plot(y_drifts_l,'k')
%     plot(y_drifts_r,'g')
%     ylim([-0.25 0.25])
%     set(gcf,'Position',[1400 800 400 800])
%     for cf = 1:n_fixs
%         fprintf('Fixation %d of %d\n',cf,n_fixs);
%         cur_im_nums = rel_fix_start_inds(cf):rel_fix_stop_inds(cf);
%         
%         n_chunks = ceil(length(cur_im_nums)/chunk_dur);
%         
%         chunk_starts = (0:n_chunks-1)*chunk_dur + 1;
%         chunk_stops = chunk_starts + chunk_dur;
%         chunk_stops(chunk_stops > length(cur_im_nums)) = length(cur_im_nums);
%         for tt = 1:length(chunk_starts)
%             cur = cur_im_nums(chunk_starts(tt));
%             figure(1)
%             imagesc(x_shifts/Fsd,y_shifts/Fsd,reshape(ov_lgamma(cur,:),length(y_shifts),length(x_shifts)));
%             figure(2)
%             subplot(2,1,1)
%             xlim([cur-20 cur+20])
%             subplot(2,1,2)
%             xlim([cur-20 cur+20])
%             pause(0.2)
%         end
%         pause(1)
%     end
    
    %% RECONSTRUCT NEW STIMULUS MATRIX
    resh_X = reshape(fullX',[sdim sdim NT]);
    resh_X_sh = zeros(size(resh_X));
    for ii = 1:NT
        %     if mod(ii,100)==0 fprintf('%d of %d\n',ii,NT); end
        d2 = dist_shift2d(resh_X(:,:,ii), -x_cor{it}(ii), 2,0);
        d2 = dist_shift2d(d2,-y_cor{it}(ii),1,0);
        resh_X_sh(:,:,ii) = d2;
    end
    
    %% REFIT GEM PARAMETERS
    fullX_sh = reshape(resh_X_sh,sdim^2,NT)';
    fullX_cropped = fullX_sh(:,new_crop);
    
%     %if you want to test the xv set only 
%     for t = xv_set
%         fprintf('Fitting GEM: Cell %d of %d\n',t,length(xv_set));
%         gab_priors(1).theta(1) = gabor_params_f{it}(t,1);
%         gab_priors(2).theta(1) = gabor_params_f{it}(t,2);
%         
%         init_params = gabor_params_f{it}(t,:);
%         [gabor_params_test(t,:),LL_test(t)] = fit_gabor_energy_mod(XXc,YYc,init_params,fullX_cropped,full_binned_spks(:,t),hold_const,LB,UB,gab_priors);
%     end
    
    for t = 1:96
        fprintf('Fitting GEM: Cell %d of %d\n',t,96);
        gab_priors(1).theta(1) = gabor_params_f{it}(t,1);
        gab_priors(2).theta(1) = gabor_params_f{it}(t,2);
        
        init_params = gabor_params_f{it}(t,:);
        [gabor_params_f{it+1}(t,:),LL(it+1,t)] = fit_gabor_energy_mod(XXc,YYc,init_params,fullX_cropped,full_binned_spks(:,t),hold_const,LB,UB,gab_priors);
        
        gabor_emp1 = get_pgabor_mask_v2(XX,YY,gabor_params_f{it+1}(t,1:6),0);
        gabor_emp2 = get_pgabor_mask_v2(XX,YY,gabor_params_f{it+1}(t,1:6),pi/2);
        
        gabor_emp1_filt(t,:) = gabor_emp1(:);
        gabor_emp2_filt(t,:) = gabor_emp2(:);
    end
end

%%
save expt1_eyecor_d1p25_nosac_useDeps gabor_* LL xv_set tr_set x_cor y_cor ov_lgamma

%%
Fs_t = 1/dt;
lcf = 1/10;
[fb,fa] = butter(2,lcf/(Fs_t/2),'high');

x_drifts_l = nan(n_fixs,1);
y_drifts_l = nan(n_fixs,1);
x_drifts_r = nan(n_fixs,1);
y_drifts_r = nan(n_fixs,1);
for i = 1:n_fixs
    cur_set = rel_fix_start_inds(i):rel_fix_stop_inds(i);
    medx_pos = median(full_eyepos(cur_set,1));
    medy_pos = median(full_eyepos(cur_set,2));
    x_drifts_l(cur_set) = full_eyepos(cur_set,1) - medx_pos;
    y_drifts_l(cur_set) = full_eyepos(cur_set,2) - medy_pos;
%     x_drifts_l(cur_set) = full_eyepos(cur_set,1);
%     y_drifts_l(cur_set) = full_eyepos(cur_set,2);
    medx_pos = median(full_eyepos(cur_set,3));
    medy_pos = median(full_eyepos(cur_set,4));
    x_drifts_r(cur_set) = full_eyepos(cur_set,3) - medx_pos;
    y_drifts_r(cur_set) = full_eyepos(cur_set,4) - medy_pos;
%     x_drifts_r(cur_set) = full_eyepos(cur_set,3) ;
%     y_drifts_r(cur_set) = full_eyepos(cur_set,4);
end

params.Fs = Fs_t;
params.tapers = [8 15];
params.err = [2 .05];
params.trialave =1;
win = 100;
[xl_C,~,~,~,~,f,confCxl] = coherencysegc(x_drifts_l,x_cor{end},win,params);
[xr_C,~,~,~,~,f,confCxr] = coherencysegc(x_drifts_r,x_cor{end},win,params);
[yl_C,~,~,~,~,f,confCyl] = coherencysegc(y_drifts_l,y_cor{end},win,params);
[yr_C,~,~,~,~,f,confCyr] = coherencysegc(y_drifts_r,y_cor{end},win,params);
% [xl_C,~,~,~,~,f,confCxl,~] = coherencyc(x_drifts_l,x_cor{end},params);
% [xr_C,~,~,~,~,f,confCxr,~] = coherencyc(x_drifts_r,x_cor{end},params);
% [yl_C,~,~,~,~,f,confCyl,~] = coherencyc(y_drifts_l,y_cor{end},params);
% [yr_C,~,~,~,~,f,confCyr,~] = coherencyc(y_drifts_r,y_cor{end},params);
% 
figure
subplot(2,1,1)
plot(f,xl_C)
hold on
plot(f,xr_C,'r')
subplot(2,1,2)
plot(f,yl_C)
hold on
plot(f,yr_C,'r')


x_drifts_r = filtfilt(fb,fa,x_drifts_r);
x_drifts_l = filtfilt(fb,fa,x_drifts_l);
y_drifts_r = filtfilt(fb,fa,y_drifts_r);
y_drifts_l = filtfilt(fb,fa,y_drifts_l);

yc_l = corr(y_drifts_l,y_cor{end},'type','spearman')
yc_r = corr(y_drifts_r,y_cor{end},'type','spearman')

xc_l = corr(x_drifts_l,x_cor{end},'type','spearman')
xc_r = corr(x_drifts_r,x_cor{end},'type','spearman')

%%
