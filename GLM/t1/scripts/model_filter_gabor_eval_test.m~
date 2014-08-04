%%
clear all
close all
cd '/Users/James/James_scripts/GLM/t1/allcell_fits/'
addpath('~/James_scripts/GLM/functions/')
addpath('~/James_scripts/GLM/t1/functions/')
addpath('~/James_scripts/general_functions/')

%% USABLE FILTERS, JUDGED BY EYE
used_filts = [...
    1 0 0 1 0;
    1 0 0 0 0;
    1 1 0 0 1;
    1 1 0 1 0;
    1 0 1 1 0;
    1 1 0 0 0;
    0 1 1 0 0;
    1 0 0 0 0;
    1 0 0 0 0;
    1 0 0 1 0;
    1 0 0 0 0;
    1 1 1 0 0;
    1 1 1 1 0;
    1 0 0 0 0;
    1 0 0 0 1;
    1 0 1 0 0;
    1 0 0 0 0;
    1 0 0 0 0;
    1 1 0 0 0;
    1 1 0 0 0;
    1 0 0 0 0;
    0 1 0 0 0;
    1 0 1 0 0;
    1 0 0 0 0;
    1 0 0 0 0;
    1 1 1 1 1
    1 0 0 0 1];

cell_mat = repmat((1:26)',1,5);

%% CONNECTED PAIRS INSTANTEOUS
conn_pairs = [1 3 -1;1 4 1; 1 7 -1; 1 20 1; 1 22 1; 1 24 -1; 1 25 1; 2 16 1;
    2 17 1; 2 18 -1; 2 20 -1; 3 4 -1; 3 7 1; 4 7 1; 4 22 1; 4 4 -1; 4 25 -1; 
    5 15 1; 5 21 1; 6 20 1; 6 22 1; 7 25 -1; 8 9 -1; 8 21 -1; 8 26 -1; 9 21 -1;
    9 26 -1; 10 14 1; 10 19 1; 10 23 -1; 11 22 1; 14 19 1; 14 23 -1; 11 22 1; 
    14 19 1; 14 23 -1; 16 20 -1; 17 18 -1; 18 20 -1; 20 22 1; 24 25 1]; 
connection_mat = zeros(26,26);
for i = 1:size(conn_pairs,1)
   connection_mat(conn_pairs(i,1),conn_pairs(i,2)) = conn_pairs(i,3);
   connection_mat(conn_pairs(i,2),conn_pairs(i,1)) = conn_pairs(i,3);
end

%% Cycle through all cells and fit Gabors to best slices or overall power
%% profiles.
for cc = 1:26
    fprintf('Cell %d\n',cc);
    eval(sprintf('load cell%d_500wdims_5mods_350dX_40L1.mat',cc));
    
    fsdim = full_glm.mods(1).fsdim;
    sdim = sqrt(fsdim);
    flen = length(full_glm.mods(1).pix)/fsdim;
    pids = 1:fsdim;
    
    nmods = length(full_glm.mods);
    w(cc,:) = arrayfun(@(x) x.w,full_glm.mods);
    nlx = full_glm.mods(1).nlx;
    
    pix_mat = get_pix_mat(full_glm);
    powprofs(cc,:,:)  = get2dPowerProfs(pix_mat,flen,sdim,pids);
    [best_slices(cc,:,:),best_slice_ids] = get2dMaxSlices(pix_mat,flen,sdim,pids);
    
    weighted_slices = squeeze(best_slices(cc,:,:)).*repmat(w(cc,:),length(pids),1);
    %     mod_structure(cc,:) = var(weighted_slices);
    weighted_powprofs = squeeze(powprofs(cc,:,:)).*repmat(w(cc,:),length(pids),1);
    smooth_powprofs = smoothks(weighted_powprofs,8,8,pids,sdim);
    mod_structure(cc,:) = var(weighted_powprofs);
    f1 =figure('visible','off');
    for i = 1:nmods
        fprintf('Filter %d of %d\n',i,nmods);
        
        nl_y_mat(cc,i,:) = full_glm.mods(i).nly;
%         cur_kern = squeeze(best_slices(:,:,i))';
        cur_kern = powprofs(cc,:,i)';
        [gabor_fit{cc,i},gabor_fitvals(:,cc,i),gauss_fitvals(:,cc,i)] = james_2d_gaborfit(cur_kern,sdim,0);
        
        if ~isnan(gabor_fit{cc,i}.x0) & w(cc,i) > 0
            cur_r2(cc,i) = (1-gabor_fit{cc,i}.ssq_err/gabor_fit{cc,i}.ssq);
            if cur_r2(cc,i) > 0.
                cur_mu = [gabor_fit{cc,i}.x0 gabor_fit{cc,i}.y0];
                cur_sig = gabor_fit{cc,i}.sigma;
                cur_theta = gabor_fit{cc,i}.theta;
                cur_gamma = gabor_fit{cc,i}.gamma;
                ra = 2*cur_sig;
                rb = 2*cur_sig/cur_gamma;
                
                subplot(nmods,3,3*(i-1)+1); colormap(jet);
                xm = max(abs(cur_kern(:)));
                pcolor(1:sdim,1:sdim,reshape(cur_kern,sdim,sdim));shading flat;caxis([-xm xm]);
                subplot(nmods,3,3*(i-1)+2);
                pcolor(1:sdim,1:sdim,reshape(gabor_fitvals(:,cc,i),sdim,sdim));shading flat;caxis([-xm xm])
                subplot(nmods,3,3*(i-1)+3);
                pcolor(1:sdim,1:sdim,reshape(gauss_fitvals(:,cc,i),sdim,sdim));shading flat;
            end
        end
    end
    cur_used_filts = find(used_filts(cc,:));
    rel_weights = w(cc,cur_used_filts);
    rel_weights = rel_weights/sum(rel_weights);
    avg_powprofs(cc,:) = (squeeze(powprofs(cc,:,cur_used_filts))*rel_weights')';
    avg_slices(cc,:) = (squeeze(best_slices(cc,:,cur_used_filts))*rel_weights')';
    
    [avg_gabor_fit{cc},avg_gabor_fitvals(:,cc),avg_gauss_fitvals(:,cc)] = james_2d_gaborfit(avg_powprofs(cc,:)',sdim,0);
    
    set(f1,'PaperUnits','centimeters');
    set(f1, 'PaperSize', [25 8*nmods]);
    set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))])
    pname = sprintf('cell%d_gab_fit',cc);
    print('-dpng',pname); close
end

save gabor_fitdata mod_structure powprofs best_slices gabor_fit w nl_y_mat nlx avg_* gabor_* gauss_*

%% Cycle through all cells and fit Gabors at each time slice
for cc = 1:27
% cc = 13;
    fprintf('Cell %d\n',cc);
    eval(sprintf('load cell%d_500wdims_5mods_350dX_40L1.mat',cc));
    
    fsdim = full_glm.mods(1).fsdim;
    sdim = sqrt(fsdim);
    flen = length(full_glm.mods(1).pix)/fsdim;
    pids = 1:fsdim;
    
    nmods = length(full_glm.mods);
    w(cc,:) = arrayfun(@(x) x.w,full_glm.mods);
    nlx = full_glm.mods(1).nlx;
    
    pix_mat = get_pix_mat(full_glm);
    [best_slices(cc,:,:),best_slice_ids] = get2dMaxSlices(pix_mat,flen,sdim,pids);
    
        f1 =figure('visible','off');
%     f1 =figure;
    for i = 1:nmods
        fprintf('Filter %d of %d\n',i,nmods);
        nl_y_mat(cc,i,:) = full_glm.mods(i).nly;
        cur_kern = squeeze(best_slices(cc,:,i))';
        cur_gabor_fit = james_2d_gaborfit(cur_kern,sdim,0);
        if ~isnan(cur_gabor_fit.x0) & w(cc,i) > 0
            cur_r2(cc,i) = (1-cur_gabor_fit.ssq_err/cur_gabor_fit.ssq);
            if cur_r2(cc,i) > 0.25 %if the best slice is OK use this as our initial guess and go through all slices
                filt_mat = reshape(pix_mat(:,i),6,fsdim)';
                xm = max(abs(filt_mat(:)));
                for ss = 1:6
                    cur_kern = filt_mat(:,ss);
                    [gabor_fit_slice{cc,i,ss},gabor_fitvals_slice(:,cc,i,ss)] = james_2d_gaborfit(cur_kern,sdim,0,cur_gabor_fit);
                    
                    subplot(2*nmods,8,16*(i-1)+ss); colormap(jet);
                    imagesc(1:sdim,1:sdim,reshape(cur_kern,sdim,sdim));caxis([-xm xm]);
                    subplot(2*nmods,8,16*(i-1)+8+ss); colormap(jet);
                    imagesc(1:sdim,1:sdim,reshape(gabor_fitvals_slice(:,cc,i,ss),sdim,sdim));caxis([-xm xm])
                    
                end
                gabor_filtmat = squeeze(gabor_fitvals_slice(:,cc,i,:));
                gabor_filtmat = gabor_filtmat'; filt_mat = filt_mat';
                filt_projprof = project_2drf(filt_mat(:),gabor_fit_slice{cc,i,best_slice_ids(i)}.theta,sdim);
                gabor_projprof = project_2drf(gabor_filtmat(:),gabor_fit_slice{cc,i,best_slice_ids(i)}.theta,sdim);
                filt_projprof_orth = project_2drf(filt_mat(:),gabor_fit_slice{cc,i,best_slice_ids(i)}.theta+pi/2,sdim);
                gabor_projprof_orth = project_2drf(gabor_filtmat(:),gabor_fit_slice{cc,i,best_slice_ids(i)}.theta+pi/2,sdim);
                xm = max(abs(filt_projprof(:)));
                subplot(2*nmods,8,16*(i-1)+7);
                imagesc(filt_projprof); caxis([-xm xm]);
                subplot(2*nmods,8,16*(i-1)+15)
                imagesc(gabor_projprof); caxis([-xm xm]);
                subplot(2*nmods,8,16*(i-1)+8);
                imagesc(filt_projprof_orth); caxis([-xm xm]);
                subplot(2*nmods,8,16*(i-1)+16)
                imagesc(gabor_projprof_orth); caxis([-xm xm]);
            end
        end
    end
    
    set(f1,'PaperUnits','centimeters');
    set(f1, 'PaperSize', [50 8*2*nmods]);
    set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))])
    pname = sprintf('cell%d_gab_fit_slices',cc);
    print('-dpng',pname); close
end

save gabor_fitdata_slices best_slices gabor_fit_slice w nl_y_mat nlx gabor_* 

%% OUTPUT SUMMARY FIGURES OF USED FILTERS
cd '/Users/James/James_scripts/GLM/t1/allcell_fits/'

for cc = 1:27
%     cc=13
    fprintf('Cell %d\n',cc);
    eval(sprintf('load cell%d_500wdims_5mods_350dX_40L1.mat',cc));
    
    fsdim = full_glm.mods(1).fsdim;
    sdim = sqrt(fsdim);
    flen = length(full_glm.mods(1).pix)/fsdim;
    pids = 1:fsdim;
    
    nmods = length(full_glm.mods);
    w(cc,:) = arrayfun(@(x) x.w,full_glm.mods);
    nlx = full_glm.mods(1).nlx;
    
    pix_mat = get_pix_mat(full_glm);
    powprofs(cc,:,:)  = get2dPowerProfs(pix_mat,flen,sdim,pids);
%     f1 = plot_2dmod_summary_withfit(nlx,squeeze(nl_y_mat(cc,:,:)),squeeze(best_slices(cc,:,:)),...
%         squeeze(powprofs(cc,:,:)),squeeze(gabor_fitvals(:,cc,:)),gabor_fit(cc,:),used_filts(cc,:));
    f1 = plot_2dmod_summary_withfit_timecx(nlx,squeeze(nl_y_mat(cc,:,:)),pix_mat,squeeze(powprofs(cc,:,:)),...
        squeeze(gabor_fitvals_slice(:,cc,:,:)),squeeze(gabor_fit_slice(cc,:,:)),used_filts(cc,:));
    colormap(jet);
    set(f1,'PaperUnits','centimeters');
    set(f1, 'PaperSize', [40 6*sum(used_filts(cc,:))]);
    set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))])
    pname = sprintf('cell%d_filtsum_v2',cc);
    print('-dpdf',pname); close
end

%% LOAD XCORR MAT DATA FOR THIS DATASET
cd ~/James_scripts/GLM/t1/
load ./xcorr_mat_7576_usedconds

%% CREATE SUMMARY PLOTS OF INSTANTANEOUS CELL-CELL RELATIONSHIPS
close all
% target_cell = 1;
% pos_r_cells = [3 4 20 22 25]; n_prs = length(pos_r_cells);
% neg_r_cells = [7 24]; n_nrs = length(neg_r_cells);
% target_cell = 2;
% pos_r_cells = [16 17]; n_prs = length(pos_r_cells);
% neg_r_cells = [18 20]; n_nrs = length(neg_r_cells);
% target_cell = 4;
% pos_r_cells = [7 22]; n_prs = length(pos_r_cells);
% neg_r_cells = [24 25]; n_nrs = length(neg_r_cells);
% target_cell = 10;
% pos_r_cells = [14 19]; n_prs = length(pos_r_cells);
% neg_r_cells = [23]; n_nrs = length(neg_r_cells);
target_cell = 3;
pos_r_cells = [1 7]; n_prs = length(pos_r_cells);
neg_r_cells = [4]; n_nrs = length(neg_r_cells);
% target_cell = 6;
% pos_r_cells = [20 22]; n_prs = length(pos_r_cells);
% neg_r_cells = []; n_nrs = length(neg_r_cells);


green = [0.2 0.7 0.2];

f1 = figure; hold on
for i = 1:n_prs
    if pos_r_cells(i) > target_cell
        plot(t_ax,xcorr_mat{target_cell,pos_r_cells(i)},'color',green,'linewidth',2)
        xlim([-.005 .005])
    else
        plot(t_ax,fliplr(xcorr_mat{pos_r_cells(i),target_cell}),'color',green,'linewidth',2)
        xlim([-.005 .005])
    end
end
for i = 1:n_nrs
    if neg_r_cells(i) > target_cell
        plot(t_ax,xcorr_mat{target_cell,neg_r_cells(i)},'r','linewidth',2)
        xlim([-.005 .005])
    else
        plot(t_ax,fliplr(xcorr_mat{neg_r_cells(i),target_cell}),'r','linewidth',2)
        xlim([-.005 .005])
    end
end
yl = ylim();
line([0 0],yl,'color','k')
xlabel('Time lag (s)','fontsize',14)
ylabel('Correlation','fontsize',14)

nmods = 5;
wnorm = bsxfun(@rdivide,w,sum(w,2));

powprofs_perm = permute(powprofs,[1 3 2]);
powprofs_mat = reshape(powprofs_perm,26*nmods,1024);
used_filts_vec = logical(reshape(used_filts,26*nmods,1));
cell_vec = reshape(cell_mat,26*nmods,1);
gabor_vec = reshape(gabor_fit,26*nmods,1);
w_vec = reshape(wnorm,26*nmods,1);
x0_vec = cellfun(@(x) x.x0,gabor_vec);
y0_vec = cellfun(@(x) x.y0,gabor_vec);
sigma_vec = cellfun(@(x) x.sigma,gabor_vec);
theta_vec = cellfun(@(x) x.theta,gabor_vec);
gamma_vec = cellfun(@(x) x.gamma,gabor_vec);
lambda_vec = cellfun(@(x) x.lambda,gabor_vec);

% define average stats for all good filters for each cell 
for i = 1:26
   cur_filts = find(cell_vec == i & used_filts_vec);
   utf = cur_filts(lambda_vec(cur_filts) < 20);
   if ~isempty(utf)
       cell_avg_theta(i) = 180/pi*circ_mean(mod(theta_vec(utf),pi));
   else
       cell_avg_theta(i) = nan;
   end
   fprintf('Cell%d\n',i)
   cell_avg_gamma(i) = mean(gamma_vec(cur_filts));
   cell_avg_sigma(i) = mean(sigma_vec(cur_filts));
   cell_avg_lambda(i) = mean(lambda_vec(cur_filts));
end

ra = 1.5*sigma_vec;
rb = 1.5*sigma_vec./gamma_vec;

f2 = figure; hold on
all_filts = find(used_filts_vec);
%if you want ellipses to have thickness reflecting model weights
% for i = 1:length(all_filts)
%     h=ellipse(ra(all_filts(i)),rb(all_filts(i)),theta_vec(all_filts(i)),x0_vec(all_filts(i)),...
%         y0_vec(all_filts(i)),'k',300,w_vec(all_filts(i))*0.5);
% end
% for j = 1:length(reciever_cells)
%     c2_filts = find(cell_vec == reciever_cells(j) & used_filts_vec);
%     for i = 1:length(c2_filts)
%         h=ellipse(ra(c2_filts(i)),rb(c2_filts(i)),theta_vec(c2_filts(i)),x0_vec(c2_filts(i)),...
%             y0_vec(c2_filts(i)),cmap(j,:),300,w_vec(c2_filts(i))*4);
%     end
%     target_filts = find(cell_vec == target_cell & used_filts_vec);
%     for i = 1:length(target_filts)
%         h=ellipse(ra(target_filts(i)),rb(target_filts(i)),theta_vec(target_filts(i)),x0_vec(target_filts(i)),...
%             y0_vec(target_filts(i)),[0.5 0.5 0.5],300,w_vec(target_filts(i))*4);
%     end
% end
for i = 1:length(all_filts)
    h=ellipse(ra(all_filts(i)),rb(all_filts(i)),theta_vec(all_filts(i)),x0_vec(all_filts(i)),...
        y0_vec(all_filts(i)),[0.8 0.8 0.8],300,0.1);
end
for j = 1:n_prs
    c2_filts = find(cell_vec == pos_r_cells(j) & used_filts_vec);
    for i = 1:length(c2_filts)
        h=ellipse(ra(c2_filts(i)),rb(c2_filts(i)),theta_vec(c2_filts(i)),x0_vec(c2_filts(i)),...
            y0_vec(c2_filts(i)),green,300,1.5);
    end
end
for j = 1:n_nrs
    c2_filts = find(cell_vec == neg_r_cells(j) & used_filts_vec);
    for i = 1:length(c2_filts)
        h=ellipse(ra(c2_filts(i)),rb(c2_filts(i)),theta_vec(c2_filts(i)),x0_vec(c2_filts(i)),...
            y0_vec(c2_filts(i)),'r',300,1.5);
    end
end
target_filts = find(cell_vec == target_cell & used_filts_vec);
for i = 1:length(target_filts)
    h=ellipse(ra(target_filts(i)),rb(target_filts(i)),theta_vec(target_filts(i)),x0_vec(target_filts(i)),...
        y0_vec(target_filts(i)),'b',300,1.5);
end
xlim([1 32]); ylim([1 32]);

%% CREATE SUMMARY PLOTS OF INSTANTANEOUS CELL-CELL RELATIONSHIPS
close all
% target_cell = 12;
% r_cells = [10 14 19]; n_prs = length(r_cells);
% target_cell = 3;
% r_cells = [1 7 4]; n_prs = length(r_cells);
% target_cell = 4;
% r_cells = [7 22 24 25]; n_prs = length(r_cells);
target_cell = 2;
r_cells = [16 17 18 20]; n_prs = length(r_cells);

f1 = figure; hold on
cmap = colormap(jet(n_prs));
cmap(end,:) = [1 0 0];
for i = 1:n_prs
    if r_cells(i) > target_cell
        plot(t_ax,xcorr_mat{target_cell,r_cells(i)},'color',cmap(i,:),'linewidth',3)
        xlim([-.005 .005])
    else
        plot(t_ax,fliplr(xcorr_mat{r_cells(i),target_cell}),'color',cmap(i,:),'linewidth',3)
        xlim([-.005 .005])
    end
end
yl = ylim();
line([0 0],yl,'color','k')
xlabel('Time lag (s)','fontsize',14)
ylabel('Correlation','fontsize',14)

nmods = 5;
wnorm = bsxfun(@rdivide,w,sum(w,2));

powprofs_perm = permute(powprofs,[1 3 2]);
powprofs_mat = reshape(powprofs_perm,26*nmods,1024);
used_filts_vec = logical(reshape(used_filts,26*nmods,1));
cell_vec = reshape(cell_mat,26*nmods,1);
gabor_vec = reshape(gabor_fit,26*nmods,1);
w_vec = reshape(wnorm,26*nmods,1);
x0_vec = cellfun(@(x) x.x0,gabor_vec);
y0_vec = cellfun(@(x) x.y0,gabor_vec);
sigma_vec = cellfun(@(x) x.sigma,gabor_vec);
theta_vec = cellfun(@(x) x.theta,gabor_vec);
gamma_vec = cellfun(@(x) x.gamma,gabor_vec);
lambda_vec = cellfun(@(x) x.lambda,gabor_vec);

ra = 1.5*sigma_vec;
rb = 1.5*sigma_vec./gamma_vec;

f2 = figure; hold on
all_filts = find(used_filts_vec);
for i = 1:length(all_filts)
    h=ellipse(ra(all_filts(i)),rb(all_filts(i)),theta_vec(all_filts(i)),x0_vec(all_filts(i)),...
        y0_vec(all_filts(i)),[0.8 0.8 0.8],300,0.5);
end
for j = 1:n_prs
    c2_filts = find(cell_vec == r_cells(j) & used_filts_vec);
    for i = 1:length(c2_filts)
        h=ellipse(ra(c2_filts(i)),rb(c2_filts(i)),theta_vec(c2_filts(i)),x0_vec(c2_filts(i)),...
            y0_vec(c2_filts(i)),cmap(j,:),300,2);
    end
end
target_filts = find(cell_vec == target_cell & used_filts_vec);
for i = 1:length(target_filts)
    h=ellipse(ra(target_filts(i)),rb(target_filts(i)),theta_vec(target_filts(i)),x0_vec(target_filts(i)),...
        y0_vec(target_filts(i)),'k',300,2);
end
xlim([1 32]); ylim([1 32]);

%%
nmods = 5;
wnorm = bsxfun(@rdivide,w,sum(w,2));

powprofs_perm = permute(powprofs,[1 3 2]);
powprofs_mat = reshape(powprofs_perm,26*nmods,1024);
used_filts_vec = logical(reshape(used_filts,26*nmods,1));
cell_vec = reshape(cell_mat,26*nmods,1);
gabor_vec = reshape(gabor_fit,26*nmods,1);
w_vec = reshape(wnorm,26*nmods,1);
x0_vec = cellfun(@(x) x.x0,gabor_vec);
y0_vec = cellfun(@(x) x.y0,gabor_vec);
sigma_vec = cellfun(@(x) x.sigma,gabor_vec);
theta_vec = cellfun(@(x) x.theta,gabor_vec);
gamma_vec = cellfun(@(x) x.gamma,gabor_vec);
lambda_vec = cellfun(@(x) x.lambda,gabor_vec);

lambda_vec(lambda_vec > 20) = nan; %these are just gaussians

lambda_diff_mat = sqrt(bsxfun(@minus,lambda_vec,lambda_vec').^2);
orientation_diff_mat = sqrt(bsxfun(@minus,theta_vec,theta_vec').^2); 
center_diff_mat = sqrt(bsxfun(@minus,x0_vec,x0_vec').^2 + bsxfun(@minus,y0_vec,y0_vec')); 

not_connected = find(connection_mat==0);
pos_connection = find(connection_mat==1);
neg_connection = find(connection_mat==-1);

figure
plot(center_diff_mat(not_connected),lambda_diff_mat(not_connected),'k.')
hold on
plot(center_diff_mat(pos_connection),lambda_diff_mat(pos_connection),'bo')
plot(center_diff_mat(neg_connection),lambda_diff_mat(neg_connection),'ro')

%%
% powprofs_perm = permute(powprofs,[1 3 2]);
% powprofs_mat = reshape(powprofs_perm,26*nmods,1024);
% used_filts_vec = logical(reshape(used_filts,26*nmods,1));
% cell_vec = reshape(cell_mat,26*nmods,1);
% gabor_vec = reshape(gabor_fit,26*nmods,1);
% w_vec = reshape(w,26*nmods,1);
%
% x0_vec = cellfun(@(x) x.x0,gabor_vec);
% y0_vec = cellfun(@(x) x.y0,gabor_vec);
% sigma_vec = cellfun(@(x) x.sigma,gabor_vec);
% theta_vec = cellfun(@(x) x.theta,gabor_vec);
% gamma_vec = cellfun(@(x) x.gamma,gabor_vec);
% lambda_vec = cellfun(@(x) x.lambda,gabor_vec);
%
% c1 = 1;
% c2 = 13;
% c1_filts = find(cell_vec == c1 & used_filts_vec);
% c2_filts = find(cell_vec == c2 & used_filts_vec);
% all_filts = find(used_filts_vec);
% figure(1);
% for i = 1:length(c1_filts)
%     subplot(length(c1_filts),1,i)
%     imagesc(reshape(powprofs_mat(c1_filts(i),:)',32,32))
% end
% figure(2);
% for i = 1:length(c2_filts)
%     subplot(length(c2_filts),1,i)
%     imagesc(reshape(powprofs_mat(c2_filts(i),:)',32,32))
% end
%
% ra = 1*sigma_vec;
% rb = 1*sigma_vec./gamma_vec;
%
% figure(3);hold on
% for i = 1:length(all_filts)
%     h=ellipse(ra(all_filts(i)),rb(all_filts(i)),theta_vec(all_filts(i)),x0_vec(all_filts(i)),...
%         y0_vec(all_filts(i)),'k',300,w_vec(all_filts(i))*0.1);
% end
% for i = 1:length(c1_filts)
%     h=ellipse(ra(c1_filts(i)),rb(c1_filts(i)),theta_vec(c1_filts(i)),x0_vec(c1_filts(i)),...
%         y0_vec(c1_filts(i)),'r',300,w_vec(c1_filts(i))*4);
% end
% for i = 1:length(c2_filts)
%     h=ellipse(ra(c2_filts(i)),rb(c2_filts(i)),theta_vec(c2_filts(i)),x0_vec(c2_filts(i)),...
%         y0_vec(c2_filts(i)),'b',300,w_vec(c2_filts(i))*4);
% end
% xlim([1 sdim]); ylim([1 sdim]);


%% For 2-cell comparisons
% c1 = 12;
% c2 = 26;
%
% cd '/Users/James/James_scripts/GLM/2d/allcell_fits/'
% eval(sprintf('load cell%d_500wdims_5mods_350dX_40L1.mat',c1));
%
% fsdim = full_glm.mods(1).fsdim;
% sdim = sqrt(fsdim);
% flen = length(full_glm.mods(1).pix)/fsdim;
% pids = 1:fsdim;
%
% nmods = length(full_glm.mods);
% w = arrayfun(@(x) x.w,full_glm.mods);
% pix_mat = get_pix_mat(full_glm);
% powprofs  = get2dPowerProfs(pix_mat,flen,sdim,pids);
% weighted_powprofs = powprofs.*repmat(w,length(pids),1);
% smooth_powprofs = smoothks(weighted_powprofs,8,8,pids,sdim);
% mod_structure = var(weighted_powprofs);
%
% cmap = colormap(jet(nmods));
%
%
% cur_r2 = nan(2,nmods);
% figure(1); hold on
% for i = 1:nmods
%     fprintf('Filter %d of %d\n',i,nmods);
%     cur_kern = powprofs(:,i);
%     [gabor_fit{c1,i},cur_gabor_fitvals,cur_gauss_fitvals] = james_2d_gaborfit(cur_kern,sdim,0);
%
%     if ~isnan(gabor_fit{c1,i}.x0) & w(i) > 0
%         cur_r2(c1,i) = (1-gabor_fit{c1,i}.ssq_err/gabor_fit{c1,i}.ssq);
%         if cur_r2(c1,i) > 0.5
%             cur_mu = [gabor_fit{c1,i}.x0 gabor_fit{c1,i}.y0];
%             cur_sig = gabor_fit{c1,i}.sigma;
%             cur_theta = gabor_fit{c1,i}.theta;
%             cur_gamma = gabor_fit{c1,i}.gamma;
%             %     pcolor(1:sdim,1:sdim,reshape(gauss_fitvals{c1}(i,:),sdim,sdim));shading flat
%             %     hold on
%             ra = 2*cur_sig;
%             rb = 2*cur_sig/cur_gamma;
%             figure(1);
%             h=ellipse(ra,rb,cur_theta,cur_mu(1),cur_mu(2),'k',300,w(i)*2);
%             xlim([1 sdim]); ylim([1 sdim]);
%
%             figure(2); subplot(nmods,3,3*(i-1)+1); colormap(jet);
%             pcolor(1:sdim,1:sdim,reshape(cur_kern,sdim,sdim));shading flat
%             subplot(nmods,3,3*(i-1)+2);
%             pcolor(1:sdim,1:sdim,reshape(cur_gabor_fitvals,sdim,sdim));shading flat;
%             subplot(nmods,3,3*(i-1)+3);
%             pcolor(1:sdim,1:sdim,reshape(cur_gauss_fitvals,sdim,sdim));shading flat;
%         end
%     end
% end
% full_glm1 = full_glm;
%
% %%
% cd '/Users/James/James_scripts/GLM/2d/allcell_fits/'
% eval(sprintf('load cell%d_500wdims_5mods_350dX_40L1.mat',c2));
%
% fsdim = full_glm.mods(1).fsdim;
% sdim = sqrt(fsdim);
% flen = length(full_glm.mods(1).pix)/fsdim;
% pids = 1:fsdim;
%
% nmods = length(full_glm.mods);
% w2 = arrayfun(@(x) x.w,full_glm.mods);
% pix_mat2 = get_pix_mat(full_glm);
% powprofs2  = get2dPowerProfs(pix_mat2,flen,sdim,pids);
% weighted_powprofs2 = powprofs2.*repmat(w2,length(pids),1);
% smooth_powprofs = smoothks(weighted_powprofs2,8,8,pids,sdim);
% mod_structure2 = var(weighted_powprofs2);
%
% cmap = colormap(jet(nmods));
%
% figure(1); hold on
% for i = 1:nmods
%     fprintf('Filter %d of %d\n',i,nmods);
%     cur_kern = powprofs2(:,i);
%     [gabor_fit{c2,i},cur_gabor_fitvals,cur_gauss_fitvals] = james_2d_gaborfit(cur_kern,sdim,0);
%
%     if ~isnan(gabor_fit{c2,i}.x0) & w2(i) > 0
%         cur_r2(c2,i) = (1-gabor_fit{c2,i}.ssq_err/gabor_fit{c2,i}.ssq);
%         if cur_r2(c2,i) > 0.5
%             cur_mu = [gabor_fit{c2,i}.x0 gabor_fit{c2,i}.y0];
%             cur_sig = gabor_fit{c2,i}.sigma;
%             cur_theta = gabor_fit{c2,i}.theta;
%             cur_gamma = gabor_fit{c2,i}.gamma;
%             %     pcolor(1:sdim,1:sdim,reshape(gauss_fitvals{c2}(i,:),sdim,sdim));shading flat
%             %     hold on
%             ra = 2*cur_sig;
%             rb = 2*cur_sig/cur_gamma;
%             figure(1);
%             h=ellipse(ra,rb,cur_theta,cur_mu(1),cur_mu(2),'r',300,w2(i)*2);
%             xlim([1 sdim]); ylim([1 sdim]);
%
%             figure(3); subplot(nmods,3,3*(i-1)+1); colormap(jet);
%             pcolor(1:sdim,1:sdim,reshape(cur_kern,sdim,sdim));shading flat
%             subplot(nmods,3,3*(i-1)+2);
%             pcolor(1:sdim,1:sdim,reshape(cur_gabor_fitvals,sdim,sdim));shading flat;
%             subplot(nmods,3,3*(i-1)+3);
%             pcolor(1:sdim,1:sdim,reshape(cur_gauss_fitvals,sdim,sdim));shading flat;
%         end
%     end
% end
%
%
% % plot2d_mod(full_glm1); colormap(gray);
% % plot2d_mod(full_glm); colormap(gray);
% %
