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

% load ./Expt1_compiled_data_fixedlag.mat
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
    pred_r = squeeze(gabor_emp_X(xv_samp,t))*beta(t,2) + beta(t,1);
    pred_r = exp(pred_r);
    xvLL(t) = -sum(cur_binned_spks(xv_samp).*log(pred_r) - pred_r)/sum(cur_binned_spks(xv_samp));
    
    avg_rate(t) = mean(cur_binned_spks(tr_samp));
    null_LL(t) = -sum(bsxfun(@minus,cur_binned_spks(tr_samp)*log(avg_rate(t)),avg_rate(t)))/sum(cur_binned_spks(tr_samp));
    null_xvLL(t) = -sum(bsxfun(@minus,cur_binned_spks(xv_samp)*log(avg_rate(t)),avg_rate(t)))/sum(cur_binned_spks(xv_samp));
    
end
pvals = [stats(:).p];
sig_fits = find(pvals(2,:) < .05);

%%
% cd ~/Data/bruce/7_15_12/G029
% load ./grating_glm_fits.mat
% grating_beta = beta;

pvals = [stats(:).p];
sig_fits = find(pvals(2,:) < 0.01);

mult_fact = 3;
% beta(:,2) = beta(:,2)*mult_fact;
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
            
    [fitp,grad] = GLMsolve_jmm( gabor_emp_X(tr_samp,t), tr_spikebins, [0; 0], 1, [], [], [], [], [], [], 0);
    beta(t,:) = fitp.k;
    [fitp,grad] = GLMsolve_jmm( gabor_emp_X(tr_samp,t), tr_spikebins, [mult_fact*beta(t,1); 0], 1, [], [], [], [], [], [1], 0);
    beta_m(t,:) = fitp.k;
    pred_rate = gabor_emp_X(:,t)*beta(t,1) + beta(t,2);
    pred_rate = log(1 + exp(pred_rate));
    min_pred_rate(t) = min(pred_rate);
    max_pred_rate(t) = max(pred_rate);
    pred_rate95(t) = prctile(pred_rate,95);
    
end


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

%%
full_eyepos = bsxfun(@minus,full_eyepos,mean(full_eyepos));

fullX_image = reshape(fullX,NT,sdim,sdim);

clear max_*
disp = 1;
tr_set = find(pvals(2,:) < 0.01);

max_shift = 30;
dshift = 2;
x_shifts = -max_shift:dshift:max_shift;
y_shifts = -max_shift:dshift:max_shift;

for fix = 300:length(used_fixs)
    tot_LL = zeros(length(y_shifts),length(x_shifts));
    tot_LL_m = zeros(length(y_shifts),length(x_shifts));
    cur_ims = fix_start_inds(used_fixs(fix)):fix_stop_inds(used_fixs(fix));
    
    tot_spks = 0;
    for ii = 1:length(cur_ims);
        fprintf('Image %d of %d\n',ii,length(cur_ims));
        cur_im = reshape(fullX(cur_ims(ii),:),sdim,sdim);
        
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
        
        pred_Rs = bsxfun(@times,gabor_outs,beta(:,1));
        pred_Rs = bsxfun(@plus,pred_Rs,beta(:,2));
%         pred_Rs = exp(pred_Rs);
        pred_Rs = log(1+exp(pred_Rs));
        LLs = bsxfun(@times,log(pred_Rs),full_binned_spks(cur_ims(ii),:)');
        LLs = LLs - pred_Rs;
        net_LLs = sum(LLs(tr_set,:));
%         net_Ls = exp(net_LLs);

        pred_Rs = bsxfun(@times,gabor_outs,beta_m(:,1));
        pred_Rs = bsxfun(@plus,pred_Rs,beta_m(:,2));
        %         pred_Rs = exp(pred_Rs);
        pred_Rs = log(1+exp(pred_Rs));
        LLs = bsxfun(@times,log(pred_Rs),full_binned_spks(cur_ims(ii),:)');
        LLs = LLs - pred_Rs;
        net_LLs_m = sum(LLs(tr_set,:));
%         net_Ls = exp(net_LLs);
        
        [max_LLs(cur_ims(ii)),max_LL_loc(cur_ims(ii))] = max(net_LLs);
        %     fullX_shifted(cur_im_num,:) = shifted_ims(max_LL_loc(cur_im_num),:);
        
        tot_LL = tot_LL + reshape(net_LLs,length(y_shifts),length(x_shifts));
        tot_LL_m = tot_LL_m + reshape(net_LLs_m,length(y_shifts),length(x_shifts));
        
        tot_spks = tot_spks + sum(full_binned_spks(cur_ims(ii),tr_set),2);
        
    end
    
    if disp == 1
        
        tot_Ls = exp(tot_LL);
        figure
        imagesc(x_shifts/Fsd,y_shifts/Fsd,reshape(tot_LL/tot_spks,length(y_shifts),length(x_shifts)));
        colorbar
        set(gcf,'Position',[200 200 500 400])
        
        figure
        imagesc(x_shifts/Fsd,y_shifts/Fsd,reshape(tot_LL_m/tot_spks,length(y_shifts),length(x_shifts)));
        colorbar
        set(gcf,'Position',[700 200 500 400])

        figure
        plot(cur_ims,x_shift_vec(max_LL_loc(cur_ims))/Fsd,'.-')
        hold on
        plot(cur_ims,y_shift_vec(max_LL_loc(cur_ims))/Fsd,'r.-')
         set(gcf,'Position',[1200 200 500 400])
       
%         figure
%         subplot(2,1,1)
%         plot(cur_ims,bsxfun(@minus,full_eyepos(cur_ims,[1 3]),mean(full_eyepos(cur_ims,[1 3]))),'o-')
%         axis tight
%         subplot(2,1,2)
%         plot(cur_ims,bsxfun(@minus,full_eyepos(cur_ims,[2 4]),mean(full_eyepos(cur_ims,[2 4]))),'o-')
%         axis tight
%         set(gcf,'Position',[1200 200 500 800])

        pause
        close all
    end
    
end



%%
% full_eyepos = bsxfun(@minus,full_eyepos,mean(full_eyepos));
%
% clear max_*
% disp = 0;
% % tr_set = randperm(length(sig_fits));
% % tr_set = tr_set(1:round(length(tr_set)*0.8));
% % xv_set = setdiff(1:length(sig_fits),tr_set);
%
% % tr_set = sig_fits(tr_set);
% % xv_set = sig_fits(xv_set);
%
% tr_set = find(pvals(2,:) < 1e-6);
%
% max_shift = 25;
% x_shifts = -max_shift:max_shift;
% y_shifts = -max_shift:max_shift;
% tot_LL = zeros(length(y_shifts),length(x_shifts));
%
%
% % fullX_shifted = nan(size(fullX));
% % for cur_im_num = 1:NT;
% for cur_im_num = 25345:25380;
%     fprintf('Image %d of %d\n',cur_im_num,NT);
%
%
%     cur_im = reshape(fullX(cur_im_num,:),sdim,sdim);
%     if disp == 1
%         figure
%         imagesc(xax(xpatch_inds),yax(ypatch_inds),cur_im);
%         set(gca,'ydir','normal');colormap(gray)
%         rectangle('Position',[0.1 -0.7 0.6 0.6],'edgecolor','w','linewidth',2)
%         set(gcf,'Position',[200 200 500 400])
%     end
%
%     shifted_ims = nan(length(x_shifts)*length(y_shifts),sdim^2);
%     x_shift_vec = nan(length(x_shifts)*length(y_shifts),1);
%     y_shift_vec = nan(length(x_shifts)*length(y_shifts),1);
%     cnt = 1;
%     for xx = 1:length(x_shifts)
%         for yy = 1:length(y_shifts)
%             d2 = dist_shift2d(cur_im, x_shifts(xx), 2,0);
%             d2 = dist_shift2d(d2,y_shifts(yy),1,0);
%             shifted_ims(cnt,:) = d2(:);
%             x_shift_vec(cnt) = x_shifts(xx);
%             y_shift_vec(cnt) = y_shifts(yy);
%             cnt = cnt + 1;
%         end
%     end
%
%     gabor_outs1 = gabor_emp1_filt*shifted_ims';
%     gabor_outs2 = gabor_emp2_filt*shifted_ims';
%     gabor_outs = gabor_outs1.^2 + gabor_outs2.^2;
%     gabor_outs % tr_set = find(pvals(2,:) < 1e-6);
%
% max_shift = 25;
% x_shifts = -max_shift:max_shift;
% y_shifts = -max_shift:max_shift;
% tot_LL = zeros(length(y_shifts),length(x_shifts));
= sqrt(gabor_outs);
%
%     pred_Rs = bsxfun(@times,gabor_outs,beta(:,2));
%     pred_Rs = bsxfun(@plus,pred_Rs,beta(:,1));
%     pred_Rs = exp(pred_Rs);
%     LLs = bsxfun(@times,log(pred_Rs),full_binned_spks(cur_im_num,:)');
%     LLs = LLs - pred_Rs;
%     net_LLs = sum(LLs(tr_set,:));
%     net_Ls = exp(net_LLs);
%
%     [max_LLs(cur_im_num),max_LL_loc(cur_im_num)] = max(net_LLs);
% %     fullX_shifted(cur_im_num,:) = shifted_ims(max_LL_loc(cur_im_num),:);
%
% tot_LL = tot_LL + reshape(net_LLs,length(y_shifts),length(x_shifts));
%
%     if disp == 1
%
%         figure
%         subplot(2,1,1)
%         imagesc(x_shifts/Fsd,y_shifts/Fsd,reshape(net_Ls,length(y_shifts),length(x_shifts)));
%         colorbar
%         subplot(2,1,2)
%         imagesc(x_shifts/Fsd,y_shifts/Fsd,reshape(net_LLs,length(y_shifts),length(x_shifts)));
%         colorbar
%         set(gcf,'Position',[700 200 500 800])
%
%         maxlag = 50;
%         cur_ax = (cur_im_num-maxlag):(cur_im_num+maxlag);
%         cur_ax(cur_ax < 1) = [];
%         trial_diff = [0; diff(full_trial_vec(cur_ax))];
%         trial_changes = find(trial_diff > 0);
%        figure
%         subplot(2,1,1)
%         plot(cur_ax-cur_im_num,bsxfun(@minus,full_eyepos(cur_ax,[1 3]),mean(full_eyepos(cur_ax,[1 3]))),'o-')
%         axis tight
%         yl = ylim();
%         line([-2 -2],yl)
%         line([0 0],yl,'color','r')
%         for i = 1:length(trial_changes)
%         line(cur_ax([trial_changes(i) trial_changes(i)])-cur_im_num,yl,'color','g')
%         end
%         subplot(2,1,2)
%         plot(cur_ax-cur_im_num,bsxfun(@minus,full_eyepos(cur_ax,[2 4]),mean(full_eyepos(cur_ax,[2 4]))),'o-')
%         axis tight
%         yl = ylim();
%         line([-2 -2],yl)
%         line([0 0],yl,'color','r')
%         for i = 1:length(trial_changes)
%         line(cur_ax([trial_changes(i) trial_changes(i)])-cur_im_num,yl,'color','g')
%         end
%         set(gcf,'Position',[1200 200 500 800])
%         pause
%         close all
%     end
%
% end
%
% all_xloc = x_shift_vec(max_LL_loc);
% all_yloc = y_shift_vec(max_LL_loc);

