clear all;
close all;

addpath('~/Timm/rust/SparseFilterSelection/')
addpath('~/Timm/MatlabRepository/')
addpath(genpath('~/James_Scripts'))

pars = load('~/Data/rust/infos/StandardParametersRust.mat');
defmod = pars.defmod;
defmod.h(1:end-1) = [];%restrict PSC term to delta function
flen = pars.flen;
hilen = length(defmod.h);
foff = flen + pars.hilen;

datdir = '~/Data/rust/stcbar/Data/';
cfiles = dir([datdir,'*stc.mat']); ncells = length(cfiles);
fnames = arrayfun(@(x)x.name,cfiles,'UniformOutput',0);
uids   = cellfun(@(x)[x(2:3),'-',x(6:7)],fnames,'UniformOutput',0);
cd ~/Data/blanche/rec_75/matlabdata/
load stdparsRec75.mat

cd /Users/James/James_scripts/stc_sparse_test_figs/stac_allcells
load ./used_stcims.mat


ncells = length(uids);
nn = 5; flen = 14;

cd /Users/James/James_scripts/stc_sparse_test_figs/revised_modfits2/

fprintf('ANALYZING CELL %d OF %d\n\n',nn,ncells);

%% load data
eval(['load ',['~/Data/rust/stcbar/Data/',fnames{nn}]]);
spikebins = convert_to_spikebins(spikes_per_frm(:));

% psth = spikes_per_frm(:);
% rbins = (find(psth>0.5));
% nsp = psth(rbins);
% spikebins =[];
% uvals = unique(psth);
% for i = 1:length(uvals)
%     cur_set = find(psth==uvals(i));
%     spikebins = [spikebins; repmat(cur_set(:),uvals(i),1)];
% end
% spikebins = sort(spikebins);
%
%% create XV data
[stimlen,sdim] = size(stim);
nparts = 60;
partlen = floor(stimlen/nparts);
nfold = 10;
partsperfold = nparts/nfold;
nxvparts = nparts/nfold;

% %boundaries of parts
% pbounds  = [(0:nparts-1)*partlen+1;(1:nparts)*partlen]';
%
% for i = 1:nfold
%     xv_inds{i} = [];
%     xv_spkbns{i} = [];
%     tr_inds{i} = [];
%     tr_spkbns{i} = [];
%
%     cur_perm = randperm(nparts);
%     cur_xv_parts = sort(cur_perm(1:nxvparts));
%     cur_tr_parts = setdiff(1:nparts,cur_xv_parts);
%
%     xv_spks = [];
%     xv_new_inds = nan(stimlen,1);
%     for j = 1:length(cur_xv_parts)
%         cur_start = pbounds(cur_xv_parts(j),1);
%         cur_stop = pbounds(cur_xv_parts(j),2);
%         xv_spks = [xv_spks; spikebins(spikebins >= cur_start & spikebins < cur_stop)];
%
%         cur_inds = (cur_start:cur_stop) - cur_start + length(xv_inds{i}) + 1;
%         xv_new_inds(cur_start:cur_stop) = cur_inds;
%
%         xv_inds{i} = [xv_inds{i} cur_start:cur_stop];
%     end
%     xv_spkbns{i} = xv_new_inds(xv_spks);
%
%     tr_spks = [];
%     tr_new_inds = nan(stimlen,1);
%     for j = 1:length(cur_tr_parts)
%         cur_start = pbounds(cur_tr_parts(j),1);
%         cur_stop = pbounds(cur_tr_parts(j),2);
%         tr_spks = [tr_spks; spikebins(spikebins >= cur_start & spikebins < cur_stop)];
%
%         cur_inds = (cur_start:cur_stop) - cur_start + length(tr_inds{i}) + 1;
%         tr_new_inds(cur_start:cur_stop) = cur_inds;
%
%         tr_inds{i} = [tr_inds{i} cur_start:cur_stop];
%     end
%     tr_spkbns{i} = tr_new_inds(tr_spks);
% end
% [NT,sdim] = size(stim);


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
group_parts = [];
for i = 1:partsperfold
    cur_part_range = (i-1)*nfold + (1:nfold);
    group_assignments = circshift(1:nfold,[0 (i-1)]);
    group_parts = [group_parts; cur_part_range(group_assignments)];
end
group_parts = part_rate_orders(group_parts);
for i = 1:nfold
    group_avg_rate(i) = mean(part_rate(group_parts(:,i)));
end


%create Nfold different sets of XV and TR data
for i = 1:nfold
    xv_inds{i} = [];
    tr_inds{i} = [];
    
    xv_inds{i} = [];
    for j = 1:partsperfold
        xv_inds{i} = [xv_inds{i} part_inds{group_parts(j,i)}];
    end
    tr_set = setdiff(1:nparts,group_parts(:,i));
    
    tr_inds{i} = [];
    for j = 1:length(tr_set)
        tr_inds{i} = [tr_inds{i} part_inds{tr_set(j)}];
    end
    
end

for i = 1:nfold
    tr_spkbns{i} = find(ismember(tr_inds{i},spikebins));
    xv_spkbns{i} = find(ismember(xv_inds{i},spikebins));
    tr_rate(i) = length(tr_spkbns{i})/length(tr_inds{i});
    xv_rate(i) = length(xv_spkbns{i})/length(xv_inds{i});
end

%%
klen = flen*sdim;
stim_params.spatial_dims = 1;
stim_params.sdim = sdim;
stim_params.flen = flen;

for xv = 1:nfold
    
    cur_xv_stim = stim(xv_inds{xv},:);
    cur_xv_stimemb = makeStimRows(cur_xv_stim,flen);
    cur_xv_spkbns = xv_spkbns{xv};
    
    cur_tr_stim = stim(tr_inds{xv},:);
    cur_tr_stimemb = makeStimRows(cur_tr_stim,flen);
    cur_tr_spkbns = tr_spkbns{xv};
    
    stimlen = size(cur_tr_stimemb,1);
    Robs = zeros(1,stimlen);
    ftable = tabulate(cur_tr_spkbns);
    Robs(ftable(:,1)) = ftable(:,2);
    stimlenxv = size(cur_xv_stimemb,1);
    Robsxv = zeros(1,stimlenxv);
    ftable = tabulate(cur_xv_spkbns);
    Robsxv(ftable(:,1)) = ftable(:,2);
    
    
%     %compute STA
%     sta = mean(cur_tr_stimemb(cur_tr_spkbns,:)) - mean(cur_tr_stimemb);
%     
%     %project out STA
%     proj_mat = sta'*inv(sta*sta')*sta;
%     stim_proj = cur_tr_stimemb - cur_tr_stimemb*proj_mat;
%     stvcv = cov(stim_proj(cur_tr_spkbns,:));
%     utvcv = cov(stim_proj);
%     [evecs,evals] = eig(stvcv-utvcv);
%     evs = diag(evals);
%     STCbvs  = evecs;
%     sta = sta';
%     npos = 1; nneg = 5;
%     stcs_compareset  = evecs(:,[1:nneg,length(evs)-npos+1:end]);
%     stcs_compareset  = stcs_compareset(:,end:-1:1);
%     rstcs = fliplr(stcs_compareset); %reversed STC kernels (suppressive first)
%     used_stcs = [sta stcs_compareset(:,1:npos) rstcs(:,1:nneg)];
%     used_stcs = bsxfun(@rdivide,used_stcs,sqrt(sum(used_stcs.^2)));
%     
%     %FIT NULL MODEL
%     avg_rate = mean(Robs);
%     xvpred_rate = ones(1,stimlenxv)*avg_rate;
%     trpred_rate = ones(1,stimlen)*avg_rate;
%     null_LL(xv) = -sum(Robs.*log(trpred_rate) - trpred_rate)/sum(Robs)
%     null_xvLL(xv) = -sum(Robsxv.*log(xvpred_rate) - xvpred_rate)/sum(Robsxv)
%     
%     %FIT STA MODEL
%     kern_types{1} = 'lin';
%     sta_mod = createGNM(sta,1,kern_types,[],stim_params);
%     g = cur_tr_stimemb*sta;
%     sta_mod = fitGNM_spkNL(sta_mod,g,cur_tr_spkbns,0);
%     sta_LL(xv) = getLL_GNM(sta_mod,cur_tr_stimemb,cur_tr_spkbns,'none');
%     sta_xvLL(xv) = getLL_GNM(sta_mod,cur_xv_stimemb,cur_xv_spkbns,'none');
%     
%     %FIT STC MODEL
%     clear kern_types
%     kern_types{1} = 'lin';
%     for i = 2:(size(used_stcs,2)+1)
%         kern_types{i} = 'quad';
%     end
%     stc_mod = createGNM(used_stcs,[1 ones(1,npos) -1*ones(1,nneg)],kern_types,[],stim_params);
%     g_mat = cur_tr_stimemb*get_k_mat(stc_mod);
%     stc_mod = fitGNM_weights(stc_mod,g_mat,cur_tr_spkbns,0);
%     [ll, ~, ~, prate, g] = getLL_GNM(stc_mod,cur_tr_stimemb,cur_tr_spkbns,'none');
%     stc_mod = fitGNM_spkNL(stc_mod,g,cur_tr_spkbns,0);
%     stc_LL(xv) = getLL_GNM(stc_mod,cur_tr_stimemb,cur_tr_spkbns,'none');
%     stc_xvLL(xv) = getLL_GNM(stc_mod,cur_xv_stimemb,cur_xv_spkbns,'none');

%     %FIT GLM
%     fprintf('XV set %d, fitting quadratic model\n',xv);
%     npos = 1; nneg = 0;
%     defmod.lambda_L1x = 10;
%     defmod.lambda_d2XT = 100;
%     clear kern_types
%     kern_types{1} = 'lin';
%     init_kerns = randn(klen,npos+nneg);
%     init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
%     glm(xv) = createGNM(init_kerns,[ones(1,npos) -1*ones(1,nneg)],kern_types,defmod,stim_params);
%     glm(xv) = fitGNM_filters(glm(xv),cur_tr_stimemb,cur_tr_spkbns,'none',[],1e-4,1e-6);    
%     [~, ~, ~, prate, g] = getLL_GNM(glm(xv),cur_tr_stimemb,cur_tr_spkbns,'none');
%     glm(xv) = fitGNM_spkNL(glm(xv),g,cur_tr_spkbns,0);
%     glm_xvLL(xv) = getLL_GNM(glm(xv),cur_xv_stimemb,cur_xv_spkbns,'none');
%     glm_LL(xv) = getLL_GNM(glm(xv),cur_tr_stimemb,cur_tr_spkbns,'none');

%     %FIT ML QUAD MODEL
%     fprintf('XV set %d, fitting quadratic model\n',xv);
%     npos = 3; nneg = 4;
%     defmod.lambda_L1x = 10;
%     defmod.lambda_d2XT = 100;
%     clear kern_types
%     kern_types{1} = 'lin';
%     for i = 2:(size(used_stcs,2)+1)
%         kern_types{i} = 'quad';
%     end
%     init_kerns = randn(klen,npos+nneg);
%     init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
%     quad_mod(xv) = createGNM(init_kerns,[ones(1,npos) -1*ones(1,nneg)],kern_types,defmod,stim_params);
%     for i = (npos+1):(npos+nneg)
%         quad_mod(xv).mods(i).lambda_d2XT = 200; %200
%         quad_mod(xv).mods(i).lambda_L1x = 20; %60
%     end
%     quad_mod(xv) = fitGNM_filters(quad_mod(xv),cur_tr_stimemb,cur_tr_spkbns,'none',[],1e-4,1e-6);
%     g_mat = cur_tr_stimemb*get_k_mat(quad_mod(xv));
%     temp = fitGNM_weights(quad_mod(xv),g_mat,cur_tr_spkbns,0);
%     
%     [~, ~, ~, prate, g] = getLL_GNM(quad_mod(xv),cur_tr_stimemb,cur_tr_spkbns,'none');
%     quad_mod(xv) = fitGNM_spkNL(quad_mod(xv),g,cur_tr_spkbns,0);
%     quad_xvLL(xv) = getLL_GNM(quad_mod(xv),cur_xv_stimemb,cur_xv_spkbns,'none');
%     quad_LL(xv) = getLL_GNM(quad_mod(xv),cur_tr_stimemb,cur_tr_spkbns,'none');
    
    
    %FIT GNM
    fprintf('XV set %d, fitting GNM\n',xv);
    npos = 5; nneg = 7;
    defmod.lambda_L1x = 50;
    defmod.lambda_d2XT = 200;
    clear kern_types
    for i = 1:(npos+nneg)
        kern_types{i} = 'threshlin';
    end
    n_test_gnm = 1;
    for j = 1:n_test_gnm
        init_kerns = randn(klen,npos+nneg);
        init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
        test_gnm(j) = createGNM(init_kerns,[ones(1,npos) -1*ones(1,nneg)],kern_types,defmod,stim_params);
        for i = (npos+1):(npos+nneg)
            %         gnm(xv).mods(i).lambda_d2XT = 450; %200
            %         gnm(xv).mods(i).lambda_L1x = 60; %60
            test_gnm(j).mods(i).lambda_d2XT = 400; %200
            test_gnm(j).mods(i).lambda_L1x = 100; %60
        end
        test_gnm(j) = fitGNM_filters(test_gnm(j),cur_tr_stimemb,cur_tr_spkbns,'none',[],1e-4,1e-6);
        test_xvLL(j) = getLL_GNM(test_gnm(j),cur_xv_stimemb,cur_xv_spkbns,'none');
    end
    [~,best_mod] = min(test_xvLL);
    gnm(xv) = test_gnm(best_mod);
    
    gnm_s(xv) = fitGNM_spkNL_reg(gnm(xv),cur_tr_stimemb,cur_tr_spkbns,30);
    
    %     [~, ~, ~, ~, g] = getLL_GNM(gnm(xv),cur_tr_stimemb,cur_tr_spkbns,'none');
    % %     old_spk_beta = gnm(xv).spk_beta;
    %     gnm(xv) = fitGNM_spkNL(gnm(xv),g,cur_tr_spkbns,0);
    % %     gnm(xv).spk_beta = old_spk_beta;
    %     [gnm(xv),old_scales,new_scales] = fitGNM_filter_scales(gnm(xv),cur_tr_stimemb,cur_tr_spkbns,'iter');
    %     gnm(xv) = rescale_GNM_regparams(gnm(xv),old_scales,new_scales);
    
    gnm_xvLL(xv) = getLL_GNM(gnm_s(xv),cur_xv_stimemb,cur_xv_spkbns,'none');
    gnm_LL(xv) = getLL_GNM(gnm_s(xv),cur_tr_stimemb,cur_tr_spkbns,'none');
    
    gnmr(xv) = gnm_s(xv);
    gnmr(xv) = adjust_all_reg(gnmr(xv),'lnl2',400); %100
    gnmr(xv) = adjust_all_reg(gnmr(xv),'nlmon',1);
    gnmr(xv) = adjust_all_reg(gnmr(xv),'nltype','uncon');
    gnmr(xv) = adjust_all_reg(gnmr(xv),'nlx',zeros(21,1));
    gnmr(xv) = adjust_all_reg(gnmr(xv),'nly',zeros(21,1));
    gnmr(xv) = fitGNM_internal_NLs(gnmr(xv),cur_tr_stimemb,cur_tr_spkbns,0,1);
    gnmr_xvLL(xv) = getLL_GNM(gnmr(xv),cur_xv_stimemb,cur_xv_spkbns,'none');
    gnmr_LL(xv) = getLL_GNM(gnmr(xv),cur_tr_stimemb,cur_tr_spkbns,'none');
    gnmr2(xv) = fitGNM_filters(gnmr(xv),cur_tr_stimemb,cur_tr_spkbns,'none',[],1e-4,1e-6);
    gnmr2_xvLL(xv) = getLL_GNM(gnmr2(xv),cur_xv_stimemb,cur_xv_spkbns,'none');
% 
%     gnmr3(xv) = fitGNM_internal_NLs(gnmr2(xv),cur_tr_stimemb,cur_tr_spkbns,0,1);
%     gnmr3_xvLL(xv) = getLL_GNM(gnmr3(xv),cur_xv_stimemb,cur_xv_spkbns,'none');
% 
%     gnmr4(xv) = fitGNM_filters(gnmr3(xv),cur_tr_stimemb,cur_tr_spkbns,'none',[],1e-4,1e-6);
%     gnmr4_xvLL(xv) = getLL_GNM(gnmr4(xv),cur_xv_stimemb,cur_xv_spkbns,'none');
   
end

%%
sta_rel_xvLL = (null_xvLL - sta_xvLL)/log(2);
stc_rel_xvLL = (null_xvLL - stc_xvLL)/log(2);
glm_rel_xvLL = (null_xvLL - glm_xvLL)/log(2);
quad_rel_xvLL = (null_xvLL - quad_xvLL)/log(2);
gnm_rel_xvLL = (null_xvLL - gnm_xvLL)/log(2);
gnmr_rel_xvLL = (null_xvLL - gnmr_xvLL)/log(2);
gnmr2_rel_xvLL = (null_xvLL - gnmr2_xvLL)/log(2);
gnmr3_rel_xvLL = (null_xvLL - gnmr3_xvLL)/log(2);
gnmr4_rel_xvLL = (null_xvLL - gnmr4_xvLL)/log(2);

% Y = [glm_rel_xvLL(:);sta_rel_xvLL(:);stc_rel_xvLL(:);quad_rel_xvLL(:);gnm_rel_xvLL(:); gnmr2_rel_xvLL(:)];
% G = [ones(nfold,1);2*ones(nfold,1);3*ones(nfold,1);4*ones(nfold,1);5*ones(nfold,1);6*ones(nfold,1)];
% Y = [sta_rel_xvLL(:);stc_rel_xvLL(:);quad_rel_xvLL(:);gnmr2_rel_xvLL(:)];
Y = [glm_rel_xvLL(:);stc_rel_xvLL(:);quad_rel_xvLL(:);gnmr2_rel_xvLL(:)];
G = [ones(nfold,1);2*ones(nfold,1);3*ones(nfold,1);4*ones(nfold,1)];
figure
boxplot(Y,G)

% Y = [sta_rel_xvLL(:) stc_rel_xvLL(:) quad_rel_xvLL(:) gnmr2_rel_xvLL(:)];
Y = [glm_rel_xvLL(:) stc_rel_xvLL(:) quad_rel_xvLL(:) gnmr2_rel_xvLL(:)];
hold on
for i = 1:nfold
    plot(1:4,Y(i,:),'ko-')
end


% stc_rel_xvLL = (sta_xvLL - stc_xvLL)/log(2);
% quad_rel_xvLL = (sta_xvLL - quad_xvLL)/log(2);
% gnm_rel_xvLL = (sta_xvLL - gnm_xvLL)/log(2);
% gnmr_rel_xvLL = (sta_xvLL - gnmr_xvLL)/log(2);
%
% figure
% Y = [stc_rel_xvLL(:);quad_rel_xvLL(:);gnm_rel_xvLL(:); gnmr_rel_xvLL(:)];
% G = [2*ones(nfold,1);3*ones(nfold,1);4*ones(nfold,1);5*ones(nfold,1)];
% boxplot(Y,G)


%%
cd ~/James_scripts/GLM/rust_scripts/
save rust_gnm_examp_v7 glm* quad* gnm* gnmr* null*


%%
use_xv = 4;



