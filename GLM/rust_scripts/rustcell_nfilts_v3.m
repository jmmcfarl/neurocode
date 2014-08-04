clear all;
close all;

% addpath('~/Timm/rust/SparseFilterSelection/')
% addpath('~/Timm/MatlabRepository/')
addpath(genpath('~/James_Scripts'))

pars = load('~/Data/rust/infos/StandardParametersRust.mat');
defmod = pars.defmod;
defmod.h(1:end-1) = [];%restrict PSC term to delta function
flen = pars.flen;
hilen = length(defmod.h);
foff = flen + pars.hilen;

datdir = '~/Data/rust/stcbar/DATA/';
cfiles = dir([datdir,'*stc.mat']); ncells = length(cfiles);
fnames = arrayfun(@(x)x.name,cfiles,'UniformOutput',0);
uids   = cellfun(@(x)[x(2:3),'-',x(6:7)],fnames,'UniformOutput',0);
cd ~/Data/blanche/rec_75/matlabdata/
load stdparsRec75.mat

cd ~/James_scripts/stc_sparse_test_figs/stac_allcells
load ./used_stcims.mat


ncells = length(uids);
nn = 5; flen = 14;

cd ~/James_scripts/stc_sparse_test_figs/revised_modfits2/

fprintf('ANALYZING CELL %d OF %d\n\n',nn,ncells);

%% load data
eval(['load ',['~/Data/rust/stcbar/DATA/',fnames{nn}]]);
spikebins = convert_to_spikebins(spikes_per_frm(:));
% create XV data
[stimlen,sdim] = size(stim);
nparts = 60;
partlen = floor(stimlen/nparts);
nfold = 5;
partsperfold = nparts/nfold;
nxvparts = nparts/nfold;
shift_size = (nfold/nxvparts);

%boundaries of parts
pbounds  = [(0:nparts-1)*partlen+1;(1:nparts)*partlen]';

%compute average firing rate within each chunk
for i = 1:nparts
    part_inds{i} = pbounds(i,1):pbounds(i,2);
    cur_spkbs = find(ismember(part_inds{i},spikebins));
    part_rate(i) = length(cur_spkbs)/length(part_inds{i});
end
[~,part_rate_orders] = sort(part_rate);

%fold group assignments of parts
group_parts1 = [];
group_parts2 = [];
for i = 1:partsperfold
    cur_part_range = (i-1)*nfold + (1:nfold);
    group_assignments = circshift(1:nfold,[0 round(shift_size*(i-1))]);
    group_parts1 = [group_parts1; cur_part_range(group_assignments)];
    group_assignments = circshift(1:nfold,[0 1+round(shift_size*(i-1))]);
    group_parts2 = [group_parts2; cur_part_range(group_assignments)];
end
group_parts1 = part_rate_orders(group_parts1);
group_parts2 = part_rate_orders(group_parts2);
for i = 1:nfold
    group_avg_rate1(i) = mean(part_rate(group_parts1(:,i)));
    group_avg_rate2(i) = mean(part_rate(group_parts2(:,i)));
end

clear tr_inds xv_inds1 xv_inds2
%create Nfold different sets of XV and TR data
for i = 1:nfold
    xv_inds1{i} = [];
    xv_inds2{i} = [];
    tr_inds{i} = [];
    
    xv_inds1{i} = [];
    xv_inds2{i} = [];
    for j = 1:partsperfold
        xv_inds1{i} = [xv_inds1{i} part_inds{group_parts1(j,i)}];
        xv_inds2{i} = [xv_inds2{i} part_inds{group_parts2(j,i)}];
    end
    comb_xv_parts = [group_parts1(:,i); group_parts2(:,i)];
    tr_set = setdiff(1:nparts,comb_xv_parts);
    
    tr_inds{i} = [];
    for j = 1:length(tr_set)
        tr_inds{i} = [tr_inds{i} part_inds{tr_set(j)}];
    end
    
end

for i = 1:nfold
    tr_spkbns{i} = find(ismember(tr_inds{i},spikebins));
    xv_spkbns1{i} = find(ismember(xv_inds1{i},spikebins));
    xv_spkbns2{i} = find(ismember(xv_inds2{i},spikebins));
    tr_rate(i) = length(tr_spkbns{i})/length(tr_inds{i});
    xv_rate1(i) = length(xv_spkbns1{i})/length(xv_inds1{i});
    xv_rate2(i) = length(xv_spkbns2{i})/length(xv_inds2{i});
end

%%
klen = flen*sdim;
stim_params.spatial_dims = 1;
stim_params.sdim = sdim;
stim_params.flen = flen;
cd ~/James_scripts/GLM/rust_scripts/
% load ./rust_gnm_filt_v4.mat
% for xv = 1:nfold
for xv = 1:nfold
    
    cur_xv_stim1 = stim(xv_inds1{xv},:);
    cur_xv_stimemb1 = makeStimRows(cur_xv_stim1,flen);
    cur_xv_spkbns1 = xv_spkbns1{xv};
    cur_xv_stim2 = stim(xv_inds2{xv},:);
    cur_xv_stimemb2 = makeStimRows(cur_xv_stim2,flen);
    cur_xv_spkbns2 = xv_spkbns2{xv};
    
    cur_tr_stim = stim(tr_inds{xv},:);
    cur_tr_stimemb = makeStimRows(cur_tr_stim,flen);
    cur_tr_spkbns = tr_spkbns{xv};
    
    stimlenxv = size(cur_xv_stimemb1,1);
    stimlen = size(cur_tr_stimemb,1);
    
    cur_tr_stim = stim(tr_inds{xv},:);
    cur_tr_stimemb = makeStimRows(cur_tr_stim,flen);
    cur_tr_spkbns = tr_spkbns{xv};
    
    Robs = convert_to_binned_spks(cur_tr_spkbns,size(cur_tr_stimemb,1))';
    Robsxv1 = convert_to_binned_spks(cur_xv_spkbns1,size(cur_xv_stimemb1,1))';
    Robsxv2 = convert_to_binned_spks(cur_xv_spkbns2,size(cur_xv_stimemb2,1))';    
        
    %FIT NULL MODEL
    avg_rate = mean(Robs);
    xvpred_rate = ones(1,stimlenxv)*avg_rate;
    trpred_rate = ones(1,stimlen)*avg_rate;
    null_LL(xv) = -sum(Robs.*log(trpred_rate) - trpred_rate)/sum(Robs)
    null_xvLL(xv) = -sum(Robsxv2.*log(xvpred_rate) - xvpred_rate)/sum(Robsxv2)
    
    %% FIT ML QUAD MODEL for subspace identification
    clear defmod
    defmod.lambda_d2XT = 200;
    defmod.lambda_L1x = 2;
    
    cur_npos =1;
    cur_nneg = 0;
    cnt = 1;
    kern_types{1} = 'lin';
    init_kerns = 0.1*randn(klen,cur_npos+cur_nneg);
    init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
    quad_mod(xv,cnt) = createGNM(init_kerns,[ones(1,cur_npos) -1*ones(1,cur_nneg)],kern_types,defmod,stim_params);
    quad_mod(xv,cnt) = fitGNM_filters(quad_mod(xv,cnt),cur_tr_stimemb,cur_tr_spkbns,'none',[],1e-4,1e-6);
    quad_xvLL(xv,cnt) = getLL_GNM(quad_mod(xv,cnt),cur_xv_stimemb1,cur_xv_spkbns1,'none');
    
    cur_xvLL = quad_xvLL(xv,cnt);
    xv_dif = -Inf;
    while xv_dif < 0
        cnt = cnt + 1;
        kern_types{cnt} = 'quad';
        init_kerns = 0.1*randn(klen,cur_npos+1+cur_nneg);
        init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
        test_pos = createGNM(init_kerns,[ones(1,cur_npos+1) -1*ones(1,cur_nneg)],kern_types,defmod,stim_params);
        test_pos = fitGNM_filters(test_pos,cur_tr_stimemb,cur_tr_spkbns,'none',[],1e-4,1e-6);
        test_pos_xvLL = getLL_GNM(test_pos,cur_xv_stimemb1,cur_xv_spkbns1,'none');
        
        test_neg = createGNM(init_kerns,[ones(1,cur_npos) -1*ones(1,cur_nneg+1)],kern_types,defmod,stim_params);
        test_neg = fitGNM_filters(test_neg,cur_tr_stimemb,cur_tr_spkbns,'none',[],1e-4,1e-6);
        test_neg_xvLL = getLL_GNM(test_neg,cur_xv_stimemb1,cur_xv_spkbns1,'none');
        
        if test_pos_xvLL <= test_neg_xvLL
            quad_mod(xv,cnt) = test_pos;
            quad_xvLL(xv,cnt) = test_pos_xvLL;
            cur_npos = cur_npos + 1;
        else
            quad_mod(xv,cnt) = test_neg;
            quad_xvLL(xv,cnt) = test_neg_xvLL;
            cur_nneg = cur_nneg + 1;
        end
        fprintf('%d pos %d neg\n',cur_npos,cur_nneg);
        xv_dif = quad_xvLL(xv,cnt)-cur_xvLL;
        cur_xvLL = quad_xvLL(xv,cnt);
    end
    quad_mod_fin(xv) = quad_mod(xv,cnt-1);
    for k = 1:2
        [~, ~, ~, prate, g] = getLL_GNM(quad_mod_fin(xv),cur_tr_stimemb,cur_tr_spkbns,'none');
        quad_mod_fin(xv) = fitGNM_spkNL(quad_mod_fin(xv),g,cur_tr_spkbns,0);
        quad_mod_fin(xv) = fitGNM_filters(quad_mod_fin(xv),cur_tr_stimemb,cur_tr_spkbns,'none',[],1e-4,1e-6);
        quad_fin_xvLL(xv) = getLL_GNM(quad_mod_fin(xv),cur_xv_stimemb2,cur_xv_spkbns2,'none');
    end

    
    %% FIT GNM for subspace identification
    defmod.lambda_d2XT = 200;
    defmod.lambda_L1x = 2;
    defmod.lnl2 = 400;
    defmod.nlx = zeros(23,1);
    defmod.nlmon = 1;
    
    cur_npos = 1;
    cur_nneg = 0;
    cnt = 1;
    kern_types{1} = 'threshlin';
    init_kerns = 0.1*randn(klen,cur_npos+cur_nneg);
    init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
    gnm_mod(xv,cnt) = createGNM(init_kerns,[ones(1,cur_npos) -1*ones(1,cur_nneg)],kern_types,defmod,stim_params);
    gnm_mod(xv,cnt) = fitGNM_filters(gnm_mod(xv,cnt),cur_tr_stimemb,cur_tr_spkbns,'none',[],1e-4,1e-6);
    gnm_mod(xv,cnt) = setGNM_NLBFs(gnm_mod(xv,cnt),cur_tr_stimemb);
    gnm_mod(xv,cnt) = adjust_all_reg(gnm_mod(xv,cnt),'nltype','uncon');
    gnm_mod(xv,cnt) = fitGNM_internal_NLs(gnm_mod(xv,cnt),cur_tr_stimemb,cur_tr_spkbns,1,2);
    gnm_xvLL(xv,cnt) = getLL_GNM(gnm_mod(xv,cnt),cur_xv_stimemb1,cur_xv_spkbns1,'none');
    
    cur_xvLL = gnm_xvLL(xv,cnt);
    xv_dif = -Inf;
    while xv_dif < 0
        cnt = cnt + 1;
        kern_types{cnt} = 'threshlin';
                init_kerns = 0.1*randn(klen,cur_npos+1+cur_nneg);
                init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
                test_pos = createGNM(init_kerns,[ones(1,cur_npos+1) -1*ones(1,cur_nneg)],kern_types,defmod,stim_params);
%         test_pos = add_gnm_filter(gnm_mod(xv,cnt-1),init_kerns,1,'threshlin');
        test_pos = fitGNM_filters(test_pos,cur_tr_stimemb,cur_tr_spkbns,'none',[],1e-4,1e-6);
        test_pos = setGNM_NLBFs(test_pos,cur_tr_stimemb);
        test_pos = adjust_all_reg(test_pos,'nltype','uncon');
        test_pos = fitGNM_internal_NLs(test_pos,cur_tr_stimemb,cur_tr_spkbns,1,2);
        test_pos_xvLL = getLL_GNM(test_pos,cur_xv_stimemb1,cur_xv_spkbns1,'none');
        
        test_neg = createGNM(init_kerns,[ones(1,cur_npos) -1*ones(1,cur_nneg+1)],kern_types,defmod,stim_params);
%         test_neg = add_gnm_filter(gnm_mod(xv,cnt-1),init_kerns,-1,'threshlin');
        test_neg = fitGNM_filters(test_neg,cur_tr_stimemb,cur_tr_spkbns,'none',[],1e-4,1e-6);
        test_neg = setGNM_NLBFs(test_neg,cur_tr_stimemb);
        test_neg = adjust_all_reg(test_neg,'nltype','uncon');
        test_neg = fitGNM_internal_NLs(test_neg,cur_tr_stimemb,cur_tr_spkbns,1,2);
        test_neg_xvLL = getLL_GNM(test_neg,cur_xv_stimemb1,cur_xv_spkbns1,'none');
        
        if test_pos_xvLL <= test_neg_xvLL
            gnm_mod(xv,cnt) = test_pos;
            gnm_xvLL(xv,cnt) = test_pos_xvLL;
            cur_npos = cur_npos + 1;
        else
            gnm_mod(xv,cnt) = test_neg;
            gnm_xvLL(xv,cnt) = test_neg_xvLL;
            cur_nneg = cur_nneg + 1;
        end
        fprintf('%d pos %d neg\n',cur_npos,cur_nneg);
        xv_dif = gnm_xvLL(xv,cnt)-cur_xvLL;
        cur_xvLL = gnm_xvLL(xv,cnt);
    end
    
    gnm_mod_fin(xv) = gnm_mod(xv,cnt-1);
    [~, ~, ~, prate, g] = getLL_GNM(gnm_mod_fin(xv),cur_tr_stimemb,cur_tr_spkbns,'none');
    gnm_mod_fin(xv) = fitGNM_spkNL(gnm_mod_fin(xv),g,cur_tr_spkbns,0);
    gnm_mod_fin(xv) = fitGNM_filters(gnm_mod_fin(xv),cur_tr_stimemb,cur_tr_spkbns,'none',[],1e-4,1e-6);
    gnm_mod_fin(xv) = fitGNM_internal_NLs(gnm_mod_fin(xv),cur_tr_stimemb,cur_tr_spkbns,0,0);
    [~, ~, ~, prate, g] = getLL_GNM(gnm_mod_fin(xv),cur_tr_stimemb,cur_tr_spkbns,'none');
    gnm_mod_fin(xv) = fitGNM_spkNL(gnm_mod_fin(xv),g,cur_tr_spkbns,0);
    gnm_fin_xvLL(xv) = getLL_GNM(gnm_mod_fin(xv),cur_xv_stimemb2,cur_xv_spkbns2,'none');
    
cd ~/James_scripts/GLM/rust_scripts/
save rust_gnm_filt_search_v1 gnm* quad* null*

end

%%

