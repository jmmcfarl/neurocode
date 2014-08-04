%% Load Data
% clear all;
addpath(genpath('~/Code/James_scripts'));
cd ~/Data/bruce/7_15_12
cd G029/
% cd G034/

% load ./Expt1_compiled_data_ms.mat
% load ./Expt1_compiled_data_yinv.mat
load ./CellList.mat
single_units = find(CellList(1,:,1) > 0);
%
% load ./Expt1_compiled_data_yinv_xinv.mat
% fullX_xy = fullX;
% load ./Expt1_compiled_data_yinv.mat
% fullX_y = fullX;
% load ./Expt1_compiled_data_xinv.mat
% fullX_x = fullX;
load ./Expt1_compiled_data.mat

% fullX = bsxfun(@minus,fullX,mean(fullX));
% fullX = bsxfun(@rdivide,fullX,std(fullX));
% fullX_x = bsxfun(@minus,fullX_x,mean(fullX_x));
% fullX_x = bsxfun(@rdivide,fullX_x,std(fullX_x));
% fullX_y = bsxfun(@minus,fullX_y,mean(fullX_y));
% fullX_y = bsxfun(@rdivide,fullX_y,std(fullX_y));
% fullX_xy = bsxfun(@minus,fullX_xy,mean(fullX_xy));
% fullX_xy = bsxfun(@rdivide,fullX_xy,std(fullX_xy));

Pix2Deg = 0.018837;
dsfrac = 1.5;
Fsd = 1/Pix2Deg/dsfrac;
Nyp = 1024;
Nxp = 1280;
Ny = round(Nyp/dsfrac);
Nx = round(Nxp/dsfrac);
% RF_patch = [0.1 0.7; -0.7 -0.1]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch = [0 0.8; -0.8 -0.]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch_pix = Fsd*RF_patch;
RF_patch_cent = mean(RF_patch,2);
RF_patch_width = diff(RF_patch,[],2);
xax = linspace(-Nx/2,Nx/2,Nx)/Fsd; yax = linspace(-Ny/2,Ny/2,Ny)/Fsd;
xpatch_inds = find(xax >= RF_patch(1,1) & xax <= RF_patch(1,2));
ypatch_inds = find(yax >= RF_patch(2,1) & yax <= RF_patch(2,2));

%%
sdim = length(xpatch_inds);
% bandwidth = 0.3;
% ar = 1.25;
% kern_len = 21-1;
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));

cd ~/Data/bruce/7_15_12/
cd G029/
load ./grating_mu_data
gr_oris = unique_oris;
load ./oned_fixation_fits_v2.mat

clear gabor_emp*
for t = 1:96
    init_params(1) = mean_x(t); %x0
    init_params(2) = mean_y(t); %y0
%         init_params(1) = rand*0.6+0.1; %x0
%         init_params(2) = rand*0.6-0.8; %y0
    init_params(3) = degtorad(pref_mu_ori(t)); %theta
%         init_params(3) = degtorad(pref_mu_ori(t))+rand*2*pi; %theta
    init_params(4) = 1/6; %lambda
    init_params(5) = 0.3*init_params(4); %sigma
    init_params(6) = 1; %eccentricity of gaussian ellipse
    gabor_emp1 = get_pgabor_mask_v2(XX,YY,init_params,0);
    gabor_emp2 = get_pgabor_mask_v2(XX,YY,init_params,pi/2);
    gabor_emp1 = repmat(gabor_emp1(:),[1 flen])';
    gabor_emp2 = repmat(gabor_emp2(:),[1 flen])';
    
    for i = 1:flen       
        temp_kern = zeros(1,flen);
        temp_kern(i) = 1;
        gabor_emp1n = bsxfun(@times,gabor_emp1,temp_kern');
        gabor_emp2n = bsxfun(@times,gabor_emp2,temp_kern');
        
        gabor_emp1_filt(t,i,:) = gabor_emp1n(:);
        gabor_emp2_filt(t,i,:) = gabor_emp2n(:);
    end
end

gabor_emp1_filt = reshape(gabor_emp1_filt,96*flen,size(fullX,2));
gabor_emp2_filt = reshape(gabor_emp2_filt,96*flen,size(fullX,2));

%%
% gabor_out_X1 = fullX*gabor_filts';
% gabor_out_X2 = fullX*gabor_filts2';
% gabor_out_X = sqrt(gabor_out_X1.^2 + gabor_out_X2.^2);
% gabor_out_X = zscore(gabor_out_X);
% n_gabor_bases = size(gabor_out_X,2);

gabor_emp_X1 = fullX*gabor_emp1_filt';
gabor_emp_X2 = fullX*gabor_emp2_filt';
gabor_emp_X = sqrt(gabor_emp_X1.^2 + gabor_emp_X2.^2);
gabor_emp_X = zscore(gabor_emp_X);

gabor_emp_X = reshape(gabor_emp_X',[96 flen size(gabor_emp_X,1)]);
gabor_emp_X = permute(gabor_emp_X,[3 1 2]);

% gabor_emp_X1 = fullX_x*gabor_emp1_filt';
% gabor_emp_X2 = fullX_x*gabor_emp2_filt';
% gabor_emp_X_x = sqrt(gabor_emp_X1.^2 + gabor_emp_X2.^2);
% gabor_emp_X_x = zscore(gabor_emp_X_x);
% 
% gabor_emp_X1 = fullX_y*gabor_emp1_filt';
% gabor_emp_X2 = fullX_y*gabor_emp2_filt';
% gabor_emp_X_y = sqrt(gabor_emp_X1.^2 + gabor_emp_X2.^2);
% gabor_emp_X_y = zscore(gabor_emp_X_y);
% 
% gabor_emp_X1 = fullX_xy*gabor_emp1_filt';
% gabor_emp_X2 = fullX_xy*gabor_emp2_filt';
% gabor_emp_X_xy = sqrt(gabor_emp_X1.^2 + gabor_emp_X2.^2);
% gabor_emp_X_xy = zscore(gabor_emp_X_xy);

% sdim = 0;
% stim_params.spatial_dims = 0;
% stim_params.sdim = sdim;
% stim_params.flen = n_gabor_bases;
% [stimlen,k_len] = size(gabor_out_X);

%%
% for t = 1:15;
%     fprintf('SU %d of %d\n',t,15);
%     cur_binned_spks = full_binned_spks(:,t);
%     unique_spk_cnts = unique(cur_binned_spks);
%     spikebins = [];
%     for i = 2:length(unique_spk_cnts)
%         cur_set = find(cur_binned_spks == unique_spk_cnts(i));
%         spikebins = [spikebins; repmat(cur_set,unique_spk_cnts(i),1)];
%     end
%     spikebins = sort(spikebins);
%
%     gabor_sta(t,:) = mean(gabor_out_X(spikebins,:));
% end

%%
xv_frac = 0.2;
NT = size(fullX,1);
xv_NT = round(xv_frac*NT);
% xv_samp = randperm(NT);
xv_samp = 1:NT;
xv_samp(xv_NT+1:end) = [];
tr_samp = setdiff(1:NT,xv_samp);

%%

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
    
    %     for or = 1:length(orientations)
    %         cur_set = find(oris == orientations(or));
    %         gabor_sta(t,or,:) = mean(gabor_out_X(spikebins,cur_set));
    %     end
    
    [beta(t,:),dev,stats] = glmfit(squeeze(gabor_emp_X(tr_samp,t,:)),cur_binned_spks(tr_samp),'poisson');
    p_vals(t,:) = stats.p;
    beta(t,p_vals(t,:) > 0.05) = 0;
    pred_r = squeeze(gabor_emp_X(tr_samp,t,:))*beta(t,2:end)' + beta(t,1);
    pred_r = exp(pred_r);
    LL(t) = -sum(cur_binned_spks(tr_samp).*log(pred_r) - pred_r)/sum(cur_binned_spks(tr_samp));
    pred_r = squeeze(gabor_emp_X(xv_samp,t,:))*beta(t,2:end)' + beta(t,1);
    pred_r = exp(pred_r);
    xvLL(t) = -sum(cur_binned_spks(xv_samp).*log(pred_r) - pred_r)/sum(cur_binned_spks(xv_samp));
    
    all_pred_r = squeeze(gabor_emp_X(:,t,:))*beta(t,2:end)' + beta(t,1);
    all_pred_r = exp(all_pred_r);
    all_LL(t,:) = cur_binned_spks.*log(all_pred_r) - all_pred_r;
    
%     [bestk(t,:),LL_sm(t)] = fitGLM_weights_gabor(squeeze(gabor_emp_X(tr_samp,t,:)),tr_spikebins,beta(t,2:end)',1);
    
%     [beta_x(t,:),dev,stats] = glmfit(gabor_emp_X_x(tr_samp,t),cur_binned_spks(tr_samp),'poisson');
%     p_vals_x(t,:) = stats.p;
%     pred_r = exp(beta(t,2)*gabor_emp_X_x(tr_samp,t)+beta_x(t,1));
%     LL_x(t) = -sum(cur_binned_spks(tr_samp).*log(pred_r) - pred_r)/sum(cur_binned_spks(tr_samp));
%     pred_r = exp(beta(t,2)*gabor_emp_X_x(xv_samp,t)+beta_x(t,1));
%     xvLL_x(t) = -sum(cur_binned_spks(xv_samp).*log(pred_r) - pred_r)/sum(cur_binned_spks(xv_samp));
%     
%     [beta_y(t,:),dev,stats] = glmfit(gabor_emp_X_y(tr_samp,t),cur_binned_spks(tr_samp),'poisson');
%     p_vals_y(t,:) = stats.p;
%     pred_r = exp(beta(t,2)*gabor_emp_X_y(tr_samp,t)+beta_y(t,1));
%     LL_y(t) = -sum(cur_binned_spks(tr_samp).*log(pred_r) - pred_r)/sum(cur_binned_spks(tr_samp));
%     pred_r = exp(beta(t,2)*gabor_emp_X_y(xv_samp,t)+beta_y(t,1));
%     xvLL_y(t) = -sum(cur_binned_spks(xv_samp).*log(pred_r) - pred_r)/sum(cur_binned_spks(xv_samp));
%     
%     [beta_xy(t,:),dev,stats] = glmfit(gabor_emp_X_xy(tr_samp,t),cur_binned_spks(tr_samp),'poisson');
%     p_vals_xy(t,:) = stats.p;
%     pred_r = exp(beta(t,2)*gabor_emp_X_xy(tr_samp,t)+beta_xy(t,1));
%     LL_xy(t) = -sum(cur_binned_spks(tr_samp).*log(pred_r) - pred_r)/sum(cur_binned_spks(tr_samp));
%     pred_r = exp(beta(t,2)*gabor_emp_X_xy(xv_samp,t)+beta_xy(t,1));
%     xvLL_xy(t) = -sum(cur_binned_spks(xv_samp).*log(pred_r) - pred_r)/sum(cur_binned_spks(xv_samp));
    
    avg_rate(t) = mean(cur_binned_spks(tr_samp));
    null_LL(t) = -sum(bsxfun(@minus,cur_binned_spks(tr_samp)*log(avg_rate(t)),avg_rate(t)))/sum(cur_binned_spks(tr_samp));
    null_xvLL(t) = -sum(bsxfun(@minus,cur_binned_spks(xv_samp)*log(avg_rate(t)),avg_rate(t)))/sum(cur_binned_spks(xv_samp));
    
end

%%
