%% Load Data
% clear all;
addpath(genpath('~/Code/James_scripts'));
cd ~/Data/bruce/7_15_12
cd G029/
% load ./eye_calibration_data
% load ./G029Expts.mat
% cd G034/

load ./CellList.mat
single_units = find(CellList(1,:,1) > 0);

load ./Expt1_compiled_data_fixedlag.mat

% fullX = bsxfun(@minus,fullX,mean(fullX));
% fullX = bsxfun(@rdivide,fullX,std(fullX));

Pix2Deg = 0.018837;
%% COMPUTE OUTPUTS OF GABOR MODELS BASED ON RF MAPPING DATA
sdim = length(xpatch_inds);
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));

cd ~/Data/bruce/7_15_12/
cd G029/
load ./grating_mu_data
gr_oris = unique_oris;
[~,sf_inds] = max(avg_mu_sf_profile,[],2);
pref_sfs = unique_sfs(sf_inds);
load ./oned_fixation_fits_v2.mat

clear gabor_emp*
for t = 1:96
    init_params(1) = mean_x(t); %x0
    init_params(2) = mean_y(t); %y0
    init_params(3) = degtorad(pref_mu_ori(t)); %theta
    init_params(4) = 1/6; %lambda
    %     init_params(4) = 1/pref_sfs(t);
    init_params(5) = 0.3*init_params(4); %sigma
    init_params(6) = 1; %eccentricity of gaussian ellipse
    gabor_emp1 = get_pgabor_mask_v2(XX,YY,init_params,0);
    gabor_emp2 = get_pgabor_mask_v2(XX,YY,init_params,pi/2);
    
    gabor_emp1_filt(t,:) = gabor_emp1(:);
    gabor_emp2_filt(t,:) = gabor_emp2(:);
end

[NT,klen] = size(fullX);

gabor_emp_X1 = fullX*gabor_emp1_filt';
gabor_emp_X2 = fullX*gabor_emp2_filt';
gabor_emp_X = sqrt(gabor_emp_X1.^2 + gabor_emp_X2.^2);
% gabor_emp_X = zscore(gabor_emp_X);

%% DIVIDE INTO TRAINING AND XV SETS
xv_frac = 0;
NT = size(fullX,1);
xv_NT = round(xv_frac*NT);
% xv_samp = randperm(NT);
xv_samp = 1:NT;
xv_samp(xv_NT+1:end) = [];
tr_samp = setdiff(1:NT,xv_samp);

%% FIT AND VALIDATE GE MODELS

mult_factor = 2;

for t = 1:96;
    fprintf('SU %d of %d\n',t,96);
    cur_binned_spks = full_binned_spks(:,t);
    unique_spk_cnts = unique(cur_binned_spks(tr_samp));
    tr_spikebins = [];
    for i = 2:length(unique_spk_cnts)
        cur_set = find(cur_binned_spks(tr_samp) == unique_spk_cnts(i));
        tr_spikebins = [tr_spikebins; repmat(cur_set,unique_spk_cnts(i),1)];
    end
    tr_spikebins = sort(tr_spikebins);
    
    [beta(t,:),dev(t),stats(t)] = glmfit(squeeze(gabor_emp_X(tr_samp,t)),cur_binned_spks(tr_samp),'poisson');
    p_vals(t,:) = stats.p;
    pred_r = squeeze(gabor_emp_X(tr_samp,t))*beta(t,2) + beta(t,1);
    pred_r = exp(pred_r);    
    LL(t) = -sum(cur_binned_spks(tr_samp).*log(pred_r) - pred_r)/sum(cur_binned_spks(tr_samp));
%     pred_r = squeeze(gabor_emp_X(xv_samp,t))*beta(t,2) + beta(t,1);
%     pred_r = exp(pred_r);
%     xvLL(t) = -sum(cur_binned_spks(xv_samp).*log(pred_r) - pred_r)/sum(cur_binned_spks(xv_samp));
    
    avg_rate(t) = mean(cur_binned_spks(tr_samp));
    null_LL(t) = -sum(bsxfun(@minus,cur_binned_spks(tr_samp)*log(avg_rate(t)),avg_rate(t)))/sum(cur_binned_spks(tr_samp));
%     null_xv5LL(t) = -sum(bsxfun(@minus,cur_binned_spks(xv_samp)*log(avg_rate(t)),avg_rate(t)))/sum(cur_binned_spks(xv_samp));
                
    [fitp,grad] = GLMsolve_jmm( gabor_emp_X(tr_samp,t), tr_spikebins, [0; 0], 1, [], [], [], [], [], [], 0);
    lexp_beta(t,:) = fitp.k;
    [fitp,grad] = GLMsolve_jmm( gabor_emp_X(tr_samp,t), tr_spikebins, [mult_factor*lexp_beta(t,1); 0], 1, [], [], [], [], [], [1], 0);
    lexp_beta(t,:) = fitp.k;
    
end


pvals = [stats(:).p];
sig_fits = find(pvals(2,:) < .05);

%% PARSE DATA INTO FIXATIONS
min_fix_dur = 0.1;
trial_flips = 1 + find(diff(full_trial_vec) ~= 0);
sac_inds = find(full_insac==1);
all_break_inds = unique([trial_flips; sac_inds]);
fix_start_inds = [1; all_break_inds+1];
fix_stop_inds = [all_break_inds-1; NT];     
fix_durs = (fix_stop_inds-fix_start_inds)*dt;
used_fixs = find(fix_durs > min_fix_dur);
fprintf('Using %d of %d fixations\n',length(used_fixs),length(fix_durs));

%% GENERATE SIMULATED DATA
fullX_image = reshape(fullX,NT,sdim,sdim);
fullX_sh = fullX;

max_delta = 0.4;
drift_sigma = 0.1;
sim_dx = round(Fsd*(rand(length(used_fixs),1)*2*max_delta - max_delta));
sim_dy = round(Fsd*(rand(length(used_fixs),1)*2*max_delta - max_delta));

for fix = 1:length(used_fixs)
    fprintf('Fix %d of %d\n',fix,length(used_fixs));
    cur_ims = fix_start_inds(used_fixs(fix)):fix_stop_inds(used_fixs(fix));
    
    for ii = 1:length(cur_ims)
        d2 = dist_shift2d(squeeze(fullX_image(cur_ims(ii),:,:)), sim_dx(fix)+round(randn*drift_sigma*Fsd), 2,0);
        d2 = dist_shift2d(d2,sim_dy(fix)+round(randn*drift_sigma*Fsd),1,0);
        fullX_sh(cur_ims(ii),:) = d2(:);
    end
end

gabor_emp_X1 = fullX_sh*gabor_emp1_filt';
gabor_emp_X2 = fullX_sh*gabor_emp2_filt';
gabor_emp_X = sqrt(gabor_emp_X1.^2 + gabor_emp_X2.^2);

sim_spks = zeros(96,NT);
for t = 1:96;
    fprintf('SU %d of %d\n',t,96);
    
    %     pred_r = squeeze(gabor_emp_X(:,t))*beta(t,2) + beta(t,1);
    %     pred_r = exp(pred_r);
    
    pred_r = squeeze(gabor_emp_X(:,t))*lexp_beta(t,1) + lexp_beta(t,2);
    pred_r = log(1+exp(pred_r));
    
    sim_spks(t,:) = poissrnd(pred_r);
end

%%
full_eyepos = bsxfun(@minus,full_eyepos,mean(full_eyepos));

clear max_*
disp = 1;
tr_set = find(pvals(2,:) < 1e-3);

max_shift = 30;
x_shifts = -max_shift:max_shift;
y_shifts = -max_shift:max_shift;

for fix = 81:length(used_fixs)
    tot_LL = zeros(length(y_shifts),length(x_shifts));
    tot_spks = 0;
    cur_ims = fix_start_inds(used_fixs(fix)):fix_stop_inds(used_fixs(fix));
    
    for ii = 1:length(cur_ims);
        fprintf('Image %d of %d\n',ii,length(cur_ims));
        cur_im = squeeze(fullX_image(cur_ims(ii),:,:));
        
        shifted_ims = nan(length(x_shifts)*length(y_shifts),sdim^2);
        x_shift_vec = nan(length(x_shifts)*length(y_shifts),1);
        y_shift_vec = nan(length(x_shifts)*length(y_shifts),1);
        cnt = 1;
        for xx = 1:length(x_shifts)
            for yy = 1:length(y_shifts)
                d2 = dist_shift2d(cur_im, x_shifts(xx), 2,0);
                d2 = dist_shift2d(d2,y_shifts(yy),1,0);
                shifted_ims(cnt,:) = d2(:);
                x_shift_vec(cnt) = x_shifts(xx);
                y_shift_vec(cnt) = y_shifts(yy);
                cnt = cnt + 1;
            end
        end
        gabor_outs1 = gabor_emp1_filt*shifted_ims';
        gabor_outs2 = gabor_emp2_filt*shifted_ims';
        gabor_outs = gabor_outs1.^2 + gabor_outs2.^2;
        gabor_outs = sqrt(gabor_outs);
        
        %         pred_Rs = bsxfun(@times,gabor_outs,beta(:,2));
        %         pred_Rs = bsxfun(@plus,pred_Rs,beta(:,1));
        %         pred_Rs = exp(pred_Rs);
        
        pred_Rs = bsxfun(@times,gabor_outs,lexp_beta(:,1));
        pred_Rs = bsxfun(@plus,pred_Rs,lexp_beta(:,2));
        pred_Rs = log(1+exp(pred_Rs));
        
        LLs = bsxfun(@times,log(pred_Rs),sim_spks(:,cur_ims(ii)));
        LLs = LLs - pred_Rs;
        net_LLs = sum(LLs(tr_set,:));
        net_Ls = exp(net_LLs);
        
        tot_spks = tot_spks + sum(sim_spks(:,cur_ims(ii)));
        
        [max_LLs(cur_ims(ii)),max_LL_loc(cur_ims(ii))] = max(net_LLs);
        %     fullX_shifted(cur_im_num,:) = shifted_ims(max_LL_loc(cur_im_num),:);
        
        tot_LL = tot_LL + reshape(net_LLs,length(y_shifts),length(x_shifts));
        
    end
    
    if disp == 1
        
        tot_Ls = exp(tot_LL);
        figure
        imagesc(x_shifts/Fsd,y_shifts/Fsd,reshape(tot_LL/tot_spks,length(y_shifts),length(x_shifts)));
        colorbar
        set(gcf,'Position',[200 200 500 400])
        title(sprintf('dx=%.4f  dy= %.4f',sim_dx(fix)/Fsd,sim_dy(fix)/Fsd));
        hold on
        plot(sim_dx(fix)/Fsd,sim_dy(fix)/Fsd,'wo','linewidth',4)
        set(gca,'ydir','normal')
%         caxis([-1.55 -1.5])
        
%         figure
%         plot(cur_ims,x_shift_vec(max_LL_loc(cur_ims))/Fsd,'.-')
%         hold on
%         plot(cur_ims,y_shift_vec(max_LL_loc(cur_ims))/Fsd,'r.-')
%          set(gcf,'Position',[700 200 500 400])      

        pause
        close all
    end
    
end

