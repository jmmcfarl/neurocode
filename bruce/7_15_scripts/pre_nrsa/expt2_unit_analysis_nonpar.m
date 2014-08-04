%% Load Data
clear all;
addpath(genpath('~/Code/James_scripts'));
% cd ~/Data/bruce/7_15_12/G029/
cd ~/Data/bruce/7_15_12/G034/

load ./CellList.mat
single_units = find(CellList(1,:,1) > 0);

load ./Expt2_compiled_windata_d1p5.mat
resh_all_stims = resh_all_stims/std(resh_all_stims(:));

fix_win_dur = 0.15;

Pix2Deg = 0.018837;
NT = length(full_fix_ids);
sdim = length(xpatch_inds);
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));

%%
n_fixs = length(full_fix_wends);
fix_binned_spks = nan(n_fixs,96);
fix_expt_num = nan(n_fixs,1);
for i = 1:n_fixs
    cur_inds = full_fix_starts(i):full_fix_wends(i);
    fix_binned_spks(i,:) = sum(full_binned_spks(cur_inds,:));
    
    fix_expt_num(i) = unique(full_expt_vec(cur_inds));
end
[un_expts,~,full_expt_inds] = unique(fix_expt_num);
n_un_expts = length(un_expts);
linX = zeros(n_fixs,n_un_expts-1);
for i = 1:n_un_expts-1
    linX(fix_expt_num==i,i) = 1;
end

%%
% load ./expt1_eyecor_d1p25_nosac_v2.mat gabor*
load ./gabor_tracking_varmeans.mat gabor*
gabor_params = gabor_params_f{end};
clear gabor_*filt
for t = 1:96
    
    gabor_emp1 = get_pgabor_mask_v2(XX,YY,gabor_params(t,1:6),0);
    gabor_emp2 = get_pgabor_mask_v2(XX,YY,gabor_params(t,1:6),pi/2);
    
    gabor_emp1_filt(t,:) = gabor_emp1(:);
    gabor_emp2_filt(t,:) = gabor_emp2(:);
end


%% RECONSTRUCT MAP STIMULUS
load ./fixation_based_corrections_34
it = length(x_cor);
resh_X = reshape(resh_all_stims',[sdim sdim n_fixs]);
resh_X_sh = zeros(size(resh_X));
for ii = 1:n_fixs
    d2 = dist_shift2d(resh_X(:,:,ii), -x_cor{it}(ii), 2,0);
    d2 = dist_shift2d(d2,-y_cor{it}(ii),1,0);
    resh_X_sh(:,:,ii) = d2;
end
fullX_sh = reshape(resh_X_sh,sdim^2,n_fixs)';

all_gabor_out1 = fullX_sh*gabor_emp1_filt';
all_gabor_out2 = fullX_sh*gabor_emp2_filt';
spatial_mod_out = sqrt(all_gabor_out1.^2 + all_gabor_out2.^2);

%%
maxlag = 0.4;
tax_inds = 0:(round(maxlag/dtd)-1);
trial_binned_spks = nan(n_fixs,length(tax_inds),96);
for i = 1:n_fixs
   cur_inds = full_fix_starts(i)+tax_inds;
   cur_inds(cur_inds > full_fix_ends(i)) = [];
   cl = length(cur_inds);
   trial_binned_spks(i,1:cl,:) = full_binned_spks(cur_inds,:);
end
trial_binned_spks = trial_binned_spks/dtd;
%%
n_p_bins = 5;
trial_prc_avgs = nan(length(tax_inds),n_p_bins,96);
prctile_bin_edges = linspace(0,100,n_p_bins+1);
for cc = 1:96
    cur_bin_edges = prctile(spatial_mod_out(:,cc),prctile_bin_edges);
   for b = 1:n_p_bins
       cur_fixs = find(spatial_mod_out(:,cc) >= cur_bin_edges(b) & spatial_mod_out(:,cc) < cur_bin_edges(b+1));
       trial_prc_avgs(:,b,cc) = nanmean(trial_binned_spks(cur_fixs,:,cc));
   end    
end

%%
load ./expt2_unit_tempmods_v2.mat

stim_outs = norminv((prctile_bin_edges(1:end-1)+prctile_bin_edges(2:end))/2/100);
cmap = jet(n_p_bins);
for i = 1:96
    
    e2_offset = mean(fblock_kern(i,:),2) + fov_const(i);
    e2_preds = repmat(fstim_dep_kern(i,:),n_p_bins,1);
    e2_preds = bsxfun(@times,e2_preds,stim_outs');
    e2_preds = bsxfun(@plus,e2_preds,fstim_ind_kern(i,:));
    e2_preds = log(1+exp(e2_preds + e2_offset));
    
    subplot(2,1,1)
    for b = 1:n_p_bins
        plot(tax_inds*dtd,smooth(squeeze(trial_prc_avgs(:,b,i)),3),'color',cmap(b,:))
        hold on
    end
    grid on
    xlim([0 0.5])
    subplot(2,1,2)
    plot(tent_centers,e2_preds'/dtd);
    title('Expt 2')
    grid on
    xlim([0 0.5])
    
    pause
    clf
end
