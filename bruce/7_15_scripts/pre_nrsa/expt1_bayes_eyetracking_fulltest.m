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
load ./Expt1_newcompiled_data_fixedlag_d1p25_full.mat
fullX = fullX/std(fullX(:));

Pix2Deg = 0.018837;
[NT,klen] = size(fullX);

%% crop stimulus for the purpose of faster gabor function fitting
new_RF_patch = [-0.11 0.9; -0.9 0.1]; %location of RFs in degrees [x1 x2;y1 y2]
[curXX,curYY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));
curXX = curXX(:);
curYY = curYY(:);

xpatch_inds_cropped = find(xax >= new_RF_patch(1,1) & xax <= new_RF_patch(1,2));
ypatch_inds_cropped = find(yax >= new_RF_patch(2,1) & yax <= new_RF_patch(2,2));

% if length(ypatch_inds_cropped) > length(xpatch_inds_cropped)
%     ypatch_inds_cropped(end) = [];
% elseif length(ypatch_inds_cropped) < length(xpatch_inds_cropped)
%     xpatch_inds_cropped(end) = [];
% end

new_crop = find(curXX >= new_RF_patch(1,1) & curXX <= new_RF_patch(1,2) & ...
    curYY >= new_RF_patch(2,1) & curYY <= new_RF_patch(2,2));
fullX_cropped = fullX(:,new_crop);

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

% PARSE DATA INTO FIXATIONS
diff_used_inds = [1; diff(used_inds)];
rel_fix_start_inds = [1; find(diff_used_inds > 1)];
rel_fix_stop_inds = [(find(diff_used_inds > 1)-1); NT];
n_fixs = length(rel_fix_start_inds);


%%
load ./expt1_eyecor_d1p25_nosac_v2 gabor* *set
n_tr_cells = length(tr_set);
for t = 1:96
    gabor_emp1 = get_pgabor_mask_v2(XX,YY,gabor_params_f{end}(t,1:6),0);
    gabor_emp2 = get_pgabor_mask_v2(XX,YY,gabor_params_f{end}(t,1:6),pi/2);
    
    gabor_emp1_filt(t,:) = gabor_emp1(:);
    gabor_emp2_filt(t,:) = gabor_emp2(:);
end
    
    %% ESTIMATE LL for each shift in each stimulus frame
    max_shift = 25;
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
            
            gfun = bsxfun(@times,energy_out,gabor_params_f{end}(tr_set,7)');
            gfun = bsxfun(@plus,gfun,gabor_params_f{end}(tr_set,8)');
            too_large = gfun > 50;
            pred_rate = log(1+exp(gfun));
            pred_rate(too_large) = gfun(too_large);
            pred_rate(pred_rate < 1e-20) = 1e-20;
            
            LLs = Robs.*log(pred_rate) - pred_rate;
            frame_LLs(:,shift_cnt) = sum(LLs,2);
            shift_cnt = shift_cnt + 1;
        end
    end
    
    
    % HMM for inferring sequence of stimulus translations
    chunk_dur = 2;
    
    %overall prior on shifts
    eps_prior_sigma = 0.3; %0.2
    leps_prior = -sum((SH/Fsd).^2,2)/(2*eps_prior_sigma^2);
    leps_prior = bsxfun(@minus,leps_prior,logsumexp(leps_prior)); %normalize
    % cur_cent = [ov_eps(cur_fix,1) ov_eps(cur_fix,2)];
    % cdist = sum((bsxfun(@minus,SH,cur_cent)/Fsd).^2,2);
    % leps_prior = -cdist/(2*eps_prior_sigma^2);
    % leps_prior = bsxfun(@minus,leps_prior,logsumexp(leps_prior));
    
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
    
    % RECONSTRUCT MAP STIMULUS
    [max_post,max_loc] = max(ov_lgamma,[],2);
    x_cor_new = SH(max_loc,1);
    y_cor_new = SH(max_loc,2);
    

%%
Fs_t = 1/dt;
lcf = 1/10;
[fb,fa] = butter(2,lcf/(Fs_t/2),'high');

x_sig_l = nan(NT,1);
y_sig_l = nan(NT,1);
x_sig_r = nan(NT,1);
y_sig_r = nan(NT,1);
for i = 1:n_fixs
    cur_set = rel_fix_start_inds(i):rel_fix_stop_inds(i);
    x_sig_l(cur_set) = full_eyepos(cur_set,1);
    y_sig_l(cur_set) = full_eyepos(cur_set,2);
    x_sig_r(cur_set) = full_eyepos(cur_set,3) ;
    y_sig_r(cur_set) = full_eyepos(cur_set,4);
end

inferred_dx = nan(n_fixs,1);
inferred_dy = nan(n_fixs,1);
for i = 2:n_fixs
    inferred_dx(i) = x_cor_new(rel_fix_start_inds(i)) - x_cor_new(rel_fix_stop_inds(i-1));
    inferred_dy(i) = y_cor_new(rel_fix_start_inds(i)) - y_cor_new(rel_fix_stop_inds(i-1));
end

% params.Fs = Fs_t;
% params.tapers = [8 15];
% params.err = [2 .05];
% params.trialave =1;
% win = 100;
% [xl_C,~,~,~,~,f,confCxl] = coherencysegc(x_sig_l,x_cor{end},win,params);
% [xr_C,~,~,~,~,f,confCxr] = coherencysegc(x_sig_r,x_cor{end},win,params);
% [yl_C,~,~,~,~,f,confCyl] = coherencysegc(y_sig_l,y_cor{end},win,params);
% [yr_C,~,~,~,~,f,confCyr] = coherencysegc(y_sig_r,y_cor{end},win,params);
% % [xl_C,~,~,~,~,f,confCxl,~] = coherencyc(x_sig_l,x_cor{end},params);
% % [xr_C,~,~,~,~,f,confCxr,~] = coherencyc(x_sig_r,x_cor{end},params);
% % [yl_C,~,~,~,~,f,confCyl,~] = coherencyc(y_sig_l,y_cor{end},params);
% % [yr_C,~,~,~,~,f,confCyr,~] = coherencyc(y_sig_r,y_cor{end},params);
% % 
% sm = 50;
% figure
% subplot(2,1,1)
% plot(f,smooth(xl_C,sm))
% hold on
% plot(f,smooth(xr_C,sm),'r')
% line(f([1 end]),[confCxl confCxl],'color','k')
% xlim([0 20])
% subplot(2,1,2)
% plot(f,smooth(yl_C,sm))
% hold on
% plot(f,smooth(yr_C,sm),'r')
% line(f([1 end]),[confCxl confCxl],'color','k')
% xlim([0 20])

% 
x_fsig_r = filtfilt(fb,fa,x_sig_r);
x_fsig_l = filtfilt(fb,fa,x_sig_l);
y_fsig_r = filtfilt(fb,fa,y_sig_r);
y_fsig_l = filtfilt(fb,fa,y_sig_l);

yc_l = corr(y_fsig_l,y_cor_new,'type','spearman')
yc_r = corr(y_fsig_r,y_cor_new,'type','spearman')

xc_l = corr(x_fsig_l,x_cor_new,'type','spearman')
xc_r = corr(x_fsig_r,x_cor_new,'type','spearman')

sacs = find(~isnan(eye_dx(:,1)));
nsacs = setdiff(2:n_fixs,sacs);
%%
close all
x_drift_l = [0; diff(x_sig_l)]/dt;
x_drift_r = [0; diff(x_sig_r)]/dt;
y_drift_l = [0; diff(y_sig_l)]/dt;
y_drift_r = [0; diff(y_sig_r)]/dt;

x_drift_l(rel_fix_start_inds) = nan;
x_drift_r(rel_fix_start_inds) = nan;
y_drift_l(rel_fix_start_inds) = nan;
y_drift_r(rel_fix_start_inds) = nan;

figure
subplot(2,2,1)
hist(x_drift_l,500)
xlim([-1.25 1.25])
yl = ylim();
line([0 0],yl,'color','r')
title('Left X','fontsize',16)
xlabel('Drift Velocity (deg/s)','fontsize',16)
subplot(2,2,2)
hist(x_drift_r,500)
yl = ylim();
line([0 0],yl,'color','r')
xlim([-1.25 1.25])
title('Right X','fontsize',16)
xlabel('Drift Velocity (deg/s)','fontsize',16)
subplot(2,2,3)
hist(y_drift_l,500)
yl = ylim();
line([0 0],yl,'color','r')
xlim([-1.25 1.25])
title('Left Y','fontsize',16)
xlabel('Drift Velocity (deg/s)','fontsize',16)
subplot(2,2,4)
hist(y_drift_r,500)
yl = ylim();
line([0 0],yl,'color','r')
xlim([-1.25 1.25])
title('Right Y','fontsize',16)
xlabel('Drift Velocity (deg/s)','fontsize',16)

figure
subplot(2,2,1)
hist(eye_dx(:,1),30)
xlim([-0.8 0.8])
yl = ylim();
line([0 0],yl,'color','r')
title('Left X','fontsize',16)
xlabel('Saccade Amplitude (deg)','fontsize',16)
subplot(2,2,2)
hist(eye_dx(:,2),30)
xlim([-0.8 0.8])
yl = ylim();
line([0 0],yl,'color','r')
title('Right X','fontsize',16)
xlabel('Saccade Amplitude (deg)','fontsize',16)
subplot(2,2,3)
hist(eye_dy(:,1),30)
xlim([-0.8 0.8])
yl = ylim();
line([0 0],yl,'color','r')
title('Left Y','fontsize',16)
xlabel('Saccade Amplitude (deg)','fontsize',16)
subplot(2,2,4)
hist(eye_dy(:,2),30)
xlim([-0.8 0.8])
yl = ylim();
line([0 0],yl,'color','r')
title('Right Y','fontsize',16)
xlabel('Saccade Amplitude (deg)','fontsize',16)

figure
subplot(2,1,1)
plot(eye_dx(:,1),eye_dx(:,2),'o')
subplot(2,1,2)
plot(eye_dy(:,1),eye_dy(:,2),'ro')

%%
%%
load ~/Data/bruce/7_15_12/G029/ArrayConfig.mat
X_pos = ArrayConfig.X;
Y_pos = ArrayConfig.Y;
