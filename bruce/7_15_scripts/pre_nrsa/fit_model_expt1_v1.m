%% Load Data
clear all;
addpath(genpath('~/Code/James_scripts'));
cd ~/Data/bruce/7_15_12
cd G029/
% cd G034/

load ./Expt1_compiled_data.mat
% load ./Expt1_compiled_data_ms.mat
% load ./Expt1_compiled_data_yinv_xinv.mat
% load ./Expt1_compiled_data_xinv
% load ./Expt1_compiled_data_xyinv

sdim = sqrt(size(fullX,2)/flen);
stim_params.spatial_dims = 2;
stim_params.sdim = sdim;
stim_params.flen = flen;
[stimlen,k_len] = size(fullX);

fullX = bsxfun(@minus,fullX,mean(fullX));

fullX =  fullX/std(fullX(:));

load ./CellList.mat
Expt_nu = [1 6 16 17 20 25 28]; %these are the grating expts
% Expt_nu = [1 2 17 18 19 20 23 24]; %these are the grating expts
single_units = find(CellList(Expt_nu(1),:,1) > 0);

%%
% for t = 1:15
t = single_units(13);
    fprintf('SU %d of %d\n',t,15);
    cur_binned_spks = full_binned_spks(:,t);
    unique_spk_cnts = unique(cur_binned_spks);
    spikebins = [];
    for i = 2:length(unique_spk_cnts)
        cur_set = find(cur_binned_spks == unique_spk_cnts(i));
        spikebins = [spikebins; repmat(cur_set,unique_spk_cnts(i),1)];
    end
    spikebins = sort(spikebins);
    
    % TRY FITTING SEQUENCE OF STC MODELS WTIH INCREASING DIMENSIONALITY
    %initialize model
    defmod.lambda_dT = 0;
    defmod.lambda_d2X = 500;
    defmod.lambda_L1x = 0;
    
    cur_ndims = 1;
    init_signs = ones(cur_ndims,1);
    init_kerns = 0.001*randn(k_len,cur_ndims);
    kern_types{1} = 'quad';
%     for i = 2:cur_ndims
%         kern_types{i} = 'quad';
%     end
    quad_mod(t) = createGNM(init_kerns,init_signs,kern_types,defmod,stim_params);
    quad_mod(t) = fitGNM_filters(quad_mod(t),fullX,spikebins,'none',50,1e-3,1e-5);
    
% end