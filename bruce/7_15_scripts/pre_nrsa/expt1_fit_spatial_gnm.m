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

% load ./Expt1_compiled_data_fixedlag.mat
% load ./Expt1_compiled_data_fixedlag.mat
cd ~/Data/bruce/7_15_12/G029/
% load ./Expt1_compiled_data_fixedlag_driftcor
% load ./Expt1_compiled_data_fixedlag_driftcor_leftonly
% load ./Expt1_compiled_data_fixedlag_driftcor_leftonly
% load ./Expt1_newcompiled_data_fixedlag_leftonly_d1p5.mat
% load ./Expt1_newcompiled_data_fixedlag_withdrift_leftonly_d1p5.mat
load ./Expt1_newcompiled_data_fixedlag_leftonly_d1p5.mat
fullX_or = fullX/std(fullX(:));
used_inds_or = used_inds;
load ./Expt1_newcompiled_data_fixedlag_withdrift_leftonly_d1p5.mat
fullX = fullX/std(fullX(:));

Pix2Deg = 0.018837;
[NT,k_len] = size(fullX);

%%
new_RF_patch = [-0.11 0.9; -0.9 0.1]; %location of RFs in degrees [x1 x2;y1 y2]
[curXX,curYY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));
curXX = curXX(:);
curYY = curYY(:);
new_crop = find(curXX >= new_RF_patch(1,1) & curXX <= new_RF_patch(1,2) & ...
    curYY >= new_RF_patch(2,1) & curYY <= new_RF_patch(2,2));
fullX = fullX(:,new_crop);
fullX_or = fullX_or(:,new_crop);

xpatch_inds = find(xax >= new_RF_patch(1,1) & xax <= new_RF_patch(1,2));
ypatch_inds = find(yax >= new_RF_patch(2,1) & yax <= new_RF_patch(2,2));

%% PARSE DATA INTO FIXATIONS
% min_fix_dur = 0.1;
% trial_flips = 1 + find(diff(full_trial_vec) ~= 0);
% sac_inds = find(full_insac==1);
% all_break_inds = unique([trial_flips; sac_inds]);
% fix_start_inds = [1; all_break_inds+1];
% fix_stop_inds = [all_break_inds-1; NT];
% fix_durs = (fix_stop_inds-fix_start_inds)*dt;
% used_fixs = find(fix_durs > min_fix_dur);
% fprintf('Using %d of %d fixations\n',length(used_fixs),length(fix_durs));
% 
% used_ims = [];
% for i = 1:length(used_fixs)
%     used_ims = [used_ims fix_start_inds(used_fixs(i)):fix_stop_inds(used_fixs(i))];
% end

%%
flen = 1;
sdim = sqrt(size(fullX,2)/flen);
stim_params.spatial_dims = 2;
stim_params.sdim = sdim;
stim_params.flen = flen;
[stimlen,k_len] = size(fullX);

%initialize model
defmod.lambda_dT = 0;
defmod.lambda_d2X = 1500;
defmod.lambda_L1x = 30;

cd ~/James_scripts/bruce/7_15_scripts/gnm_fits/  
for cur_cell = 1:96;
    fprintf('Cell %d\n',cur_cell)
    spikebins = convert_to_spikebins(full_binned_spks(used_ims,cur_cell));
    
    cur_ndims = 4;
    init_signs = ones(cur_ndims,1);
    init_kerns = 0.01*randn(k_len,cur_ndims);
    for i = 1:cur_ndims
        kern_types{i} = 'threshlin';
    end
    gnm(cur_cell) = createGNM(init_kerns,init_signs,kern_types,defmod,stim_params);
    gnm(cur_cell) = fitGNM_filters(gnm(cur_cell),fullX(used_ims,:),spikebins,'none',[],1e-4,1e-6);
    
    cur_k = get_k_mat(gnm(cur_cell));
    bad_mods = find(max(abs(cur_k)) == 0);
    gnm(cur_cell).mods(bad_mods) = [];
    plot2d_mod(gnm(cur_cell))
    fname = sprintf('GNM_C%d',cur_cell);
    fillPage(gcf,'Papersize',[6 10]);
    print(fname,'-dpng')
    close all
end

%%
flen = 1;
sdim = sqrt(size(fullX,2)/flen);
stim_params.spatial_dims = 2;
stim_params.sdim = sdim;
stim_params.flen = flen;
[stimlen,k_len] = size(fullX);

%initialize model
defmod.lambda_dT = 0;
defmod.lambda_d2X = 1000;
defmod.lambda_L1x =  30;

cd ~/Data/bruce/7_15_12/
cd G029/
load ./grating_mu_data
gr_oris = unique_oris;
[~,sf_inds] = max(avg_mu_sf_profile,[],2);
pref_sfs = unique_sfs(sf_inds);
load ./oned_fixation_fits_v3.mat
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));

cd ~/James_scripts/bruce/7_15_scripts/gnm_fits/  
for cur_cell = 1:96;
    fprintf('Cell %d\n',cur_cell)
    spikebins = convert_to_spikebins(full_binned_spks(:,cur_cell));
    
    cur_ndims = 4;
    init_signs = ones(cur_ndims,1);
    init_kerns = 0.01*randn(k_len,cur_ndims);
%     kern_types{1} = 'lin';
    for i = 1:cur_ndims
        kern_types{i} = 'quad';
    end
    quad_mod(cur_cell) = createGNM(init_kerns,init_signs,kern_types,defmod,stim_params);
    quad_mod(cur_cell) = fitGNM_filters(quad_mod(cur_cell),fullX,spikebins,'none',[],1e-4,1e-6);
    
    cur_k = get_k_mat(quad_mod(cur_cell));

    init_params(1) = mean_x(cur_cell); %x0
    init_params(2) = mean_y(cur_cell); %y0
    init_params(3) = degtorad(pref_mu_ori(cur_cell)); %theta
    init_params(4) = 1/pref_sfs(cur_cell); %lambda
    init_params(5) = 0.5*init_params(4); %sigma
    init_params(6) = 1; %eccentricity of gaussian ellipse
    gabor_emp1 = get_pgabor_mask_v2(XX,YY,init_params,0);

    figure('visible','off');
    subplot(4,2,[1 3])
    imagesc(reshape(cur_k,sdim,sdim));
    set(gca,'ydir','normal'); colormap(gray)
    subplot(4,2,[2 4])
    imagesc(gabor_emp1);
    set(gca,'ydir','normal'); colormap(gray)
    subplot(4,2,5)
    plot(unique_oris,avg_mu_ori_profile(cur_cell,:))
    axis tight
    subplot(4,2,6)
    plot(unique_sfs,avg_mu_sf_profile(cur_cell,:))
    axis tight
    subplot(4,2,7)
    plot(parallel_profile(cur_cell,:));
    axis tight
    subplot(4,2,8)
    plot(orth_profile(cur_cell,:),'r')
    axis tight        
    
    fname = sprintf('QUAD2_C%d',cur_cell);
    fillPage(gcf,'Papersize',[10 8]);
    print(fname,'-dpng')
    close all
end

%%
flen = 1;
sdim = sqrt(size(fullX,2)/flen);
stim_params.spatial_dims = 2;
stim_params.sdim = sdim;
stim_params.flen = flen;
[stimlen,k_len] = size(fullX);

%initialize model
defmod.lambda_dT = 0;
defmod.lambda_d2X = 1000;
defmod.lambda_L1x =  30;

cur_cell = 69;

spikebins = convert_to_spikebins(full_binned_spks(:,cur_cell));
cur_ndims = 1;
init_signs = ones(cur_ndims,1);
init_kerns = 0.01*randn(k_len,cur_ndims);
%     kern_types{1} = 'lin';
for i = 1:cur_ndims
    kern_types{i} = 'quad';
end
quad_mod(cur_cell) = createGNM(init_kerns,init_signs,kern_types,defmod,stim_params);
quad_mod(cur_cell) = fitGNM_filters(quad_mod(cur_cell),fullX,spikebins,'none',[],1e-4,1e-6);

uu = find(ismember(used_inds_or,used_inds));
spikebins = convert_to_spikebins(full_binned_spks(:,cur_cell));
cur_ndims = 1;
init_signs = ones(cur_ndims,1);
init_kerns = 0.01*randn(k_len,cur_ndims);
%     kern_types{1} = 'lin';
for i = 1:cur_ndims
    kern_types{i} = 'quad';
end
quad_mod_or(cur_cell) = createGNM(init_kerns,init_signs,kern_types,defmod,stim_params);
quad_mod_or(cur_cell) = fitGNM_filters(quad_mod_or(cur_cell),fullX_or(uu,:),spikebins,'none',[],1e-4,1e-6);

