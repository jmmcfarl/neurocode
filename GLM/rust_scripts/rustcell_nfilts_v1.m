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
% cd ~/Data/blanche/rec_75/matlabdata/
% load stdparsRec75.mat

cd ~/James_scripts/stc_sparse_test_figs/stac_allcells
load ./used_stcims.mat


ncells = length(uids);
nn = 5; flen = 14;

cd ~/James_scripts/stc_sparse_test_figs/revised_modfits2/

fprintf('ANALYZING CELL %d OF %d\n\n',nn,ncells);

%% load data
eval(['load ',['~/Data/rust/stcbar/DATA/',fnames{nn}]]);
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
    
    
    %compute STA
    sta = mean(cur_tr_stimemb(cur_tr_spkbns,:)) - mean(cur_tr_stimemb);
    
    %project out STA
    proj_mat = sta'*inv(sta*sta')*sta;
    stim_proj = cur_tr_stimemb - cur_tr_stimemb*proj_mat;
    stvcv = cov(stim_proj(cur_tr_spkbns,:));
    utvcv = cov(stim_proj);
    [evecs,evals] = eig(stvcv-utvcv);
    evs = diag(evals);
    STCbvs  = evecs;
    sta = sta';
    npos = 2; nneg = 6;
    stcs_compareset  = evecs(:,[1:nneg,length(evs)-npos+1:end]);
    stcs_compareset  = stcs_compareset(:,end:-1:1);
    rstcs = fliplr(stcs_compareset); %reversed STC kernels (suppressive first)
    used_stcs = [sta stcs_compareset(:,1:npos) rstcs(:,1:nneg)];
    used_stcs = bsxfun(@rdivide,used_stcs,sqrt(sum(used_stcs.^2)));
    n_sdc_dims = size(used_stcs,2);
    
    %FIT NULL MODEL
    avg_rate = mean(Robs);
    xvpred_rate = ones(1,stimlenxv)*avg_rate;
    trpred_rate = ones(1,stimlen)*avg_rate;
    null_LL(xv) = -sum(Robs.*log(trpred_rate) - trpred_rate)/sum(Robs)
    null_xvLL(xv) = -sum(Robsxv.*log(xvpred_rate) - xvpred_rate)/sum(Robsxv)
        

%FIT ML QUAD MODEL
    fprintf('XV set %d, fitting quadratic model\n',xv);
    for j = 1:6
        npos = j; nneg = j;
        defmod.lambda_L1x = 0.1;
        defmod.lambda_d2XT = 150;
        clear kern_types
        kern_types{1} = 'lin';
        for i = 2:(npos+nneg)
            kern_types{i} = 'quad';
        end
        init_kerns = 0.05*randn(klen,npos+nneg);
        init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
        quad_mod(j,xv) = createGNM(init_kerns,[ones(1,npos) -1*ones(1,nneg)],kern_types,defmod,stim_params);
        for i = (npos+1):(npos+nneg)
            quad_mod(j,xv).mods(i).lambda_d2XT = 300; %200
            quad_mod(j,xv).mods(i).lambda_L1x = 10; %60
        end
        
        n_iter = 2;
        for i = 1:n_iter
            quad_mod(j,xv) = fitGNM_filters(quad_mod(j,xv),cur_tr_stimemb,cur_tr_spkbns,'none',[],1e-4,1e-6);
            [~, ~, ~, prate, g] = getLL_GNM(quad_mod(j,xv),cur_tr_stimemb,cur_tr_spkbns,'none');
            quad_mod(j,xv) = fitGNM_spkNL(quad_mod(j,xv),g,cur_tr_spkbns,0);
            quad_xvLL(j,xv) = getLL_GNM(quad_mod(j,xv),cur_xv_stimemb,cur_xv_spkbns,'none');
            quad_LL(j,xv) = getLL_GNM(quad_mod(j,xv),cur_tr_stimemb,cur_tr_spkbns,'none');
        end
    end

%     stc_klen = 1*n_sdc_dims;
%     stc_stim_params.spatial_dims = 1;
%     stc_stim_params.sdim = n_sdc_dims;
%     stc_stim_params.flen = 1;
%     stc_out = cur_tr_stimemb*used_stcs;
%     xv_stc_out = cur_xv_stimemb*used_stcs;
%     
%     defmod.lambda_L1x = 0;
%     defmod.lambda_d2XT = 0;
%     for j = 1:7
%         for k = 1:7
%             npos = j;
%             nneg = k;
%             clear kern_types
%             for i = 1:(npos+nneg)
%                 kern_types{i} = 'threshlin';
%             end
%             init_kerns = randn(stc_klen,npos+nneg);
%             init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
%             gnm(j,k,xv) = createGNM(init_kerns,[ones(1,npos) -1*ones(1,nneg)],kern_types,defmod,stc_stim_params);
%             gnm(j,k,xv) = fitGNM_filters(gnm(j,k,xv),stc_out,cur_tr_spkbns,'none',[],1e-4,1e-6);
%             gnm_xvLL(j,k,xv) = getLL_GNM(gnm(j,k,xv),xv_stc_out,cur_xv_spkbns,'none');            
%         end
%     end
    
    %%
    mlq_basis = get_k_mat(quad_mod(6,xv));
    n_mlq_dims = size(mlq_basis,2);
    mlq_klen = 1*n_mlq_dims;
    mlq_stim_params.spatial_dims = 1;
    mlq_stim_params.sdim = n_mlq_dims;
    mlq_stim_params.flen = 1;
    mlq_out = cur_tr_stimemb*mlq_basis;
    xv_mlq_out = cur_xv_stimemb*mlq_basis;    
    defmod.lambda_L1x = 0;
    defmod.lambda_d2XT = 0;
    for j = 1:7
        for k = 1:7
            fprintf('%d Pos, %d Neg\n',j,k);
            npos = j;
            nneg = k;
            clear kern_types
            for i = 1:(npos+nneg)
                kern_types{i} = 'threshlin';
            end
            init_kerns = randn(mlq_klen,npos+nneg);
            init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
            gnm_mlq(j,k,xv) = createGNM(init_kerns,[ones(1,npos) -1*ones(1,nneg)],kern_types,defmod,mlq_stim_params);
            gnm_mlq(j,k,xv) = fitGNM_filters(gnm_mlq(j,k,xv),mlq_out,cur_tr_spkbns,'none',[],1e-4,1e-6);
            gnm_mlq_xvLL(j,k,xv) = getLL_GNM(gnm_mlq(j,k,xv),xv_mlq_out,cur_xv_spkbns,'none');
            
            gnmr_mlq(j,k,xv) = adjust_all_reg(gnm_mlq(j,k,xv),'lnl2',100);
            gnmr_mlq(j,k,xv) = setGNM_NLBFs(gnmr_mlq(j,k,xv),mlq_out);
            gnmr_mlq(j,k,xv) = adjust_all_reg(gnmr_mlq(j,k,xv),'nltype','uncon');
            gnmr_mlq(j,k,xv) = adjust_all_reg(gnmr_mlq(j,k,xv),'nlmon',1);
            for i = 1:2
                gnmr_mlq(j,k,xv) = fitGNM_internal_NLs(gnmr_mlq(j,k,xv),mlq_out,cur_tr_spkbns,1,2);
                gnmr_mlq(j,k,xv) = fitGNM_filters(gnmr_mlq(j,k,xv),mlq_out,cur_tr_spkbns,'none',[],1e-4,1e-6);
                [~, ~, ~, prate, g] = getLL_GNM(gnmr_mlq(j,k,xv),mlq_out,cur_tr_spkbns,'none');
                gnmr_mlq(j,k,xv) = fitGNM_spkNL(gnmr_mlq(j,k,xv),g,cur_tr_spkbns,0);
                gnmr_mlq_xvLL(j,k,xv) = getLL_GNM(gnmr_mlq(j,k,xv),xv_mlq_out,cur_xv_spkbns,'none');
                gnmr_mlq_LL(j,k,xv) = getLL_GNM(gnmr_mlq(j,k,xv),mlq_out,cur_tr_spkbns,'none');
            end
            
        end
    end
    
    
    %%
%     %FIT GNM
%     fprintf('XV set %d, fitting GNM\n',xv);
%     npos = 6; nneg = 6;
%     base_mod = gnmr_mlq(6,6,xv);
%     gnmr_ref(xv) = gnmr_mlq(6,6,xv);
%     for i = 1:length(base_mod.mods)
%         gnmr_ref(xv).mods(i).k = mlq_basis*base_mod.mods(i).k;
%     end
%     gnmr_ref(xv).stim_params = quad_mod(1,xv).stim_params;
%     for i = 1:npos
%         gnmr_ref(xv).mods(i).lambda_L1x = 10;
%         gnmr_ref(xv).mods(i).lambda_d2XT = 50;
%     end
%     for i = (npos+1):(npos+nneg)
%         gnmr_ref(xv).mods(i).lambda_d2XT = 50; %200
%         gnmr_ref(xv).mods(i).lambda_L1x = 20; %60
%     end
%     gnmr_ref(xv) = fitGNM_filters(gnmr_ref(xv),cur_tr_stimemb,cur_tr_spkbns,'none',50,1e-4,1e-6,0);
%     gnmr_ref_xvLL(j) = getLL_GNM(gnmr_ref(xv),cur_xv_stimemb,cur_xv_spkbns,'none');
%     
%     gnmr_ref(xv) = fitGNM_internal_NLs(gnmr_ref(xv),cur_tr_stimemb,cur_tr_spkbns,1,2);
%     [~, ~, ~, prate, g] = getLL_GNM(gnmr_ref(xv),cur_tr_stimemb,cur_tr_spkbns,'none');
%     gnmr_ref(xv) = fitGNM_spkNL(gnmr_ref(xv),g,cur_tr_spkbns,0);
%     gnmr_ref_xvLL(j) = getLL_GNM(gnmr_ref(xv),cur_xv_stimemb,cur_xv_spkbns,'none');
% 
    %%
end

%%
gnmr_mlq_xvLL_s = permute(gnmr_mlq_xvLL,[3 1 2]);
gnm_mlq_xvLL_s = permute(gnm_mlq_xvLL,[3 1 2]);
gnm_stc_xvLL_s = permute(gnm_xvLL(:,:,1:4),[3 1 2]);
gnmr_mlq_xvLL_rel = -bsxfun(@minus,gnmr_mlq_xvLL_s,null_xvLL(1:4)');
gnm_mlq_xvLL_rel = -bsxfun(@minus,gnm_mlq_xvLL_s,null_xvLL(1:4)');
gnm_stc_xvLL_rel = -bsxfun(@minus,gnm_stc_xvLL_s,null_xvLL(1:4)');

%%
quad_rel_xvLL = -bsxfun(@minus,quad_xvLL,null_xvLL);

cmap = jet(7);
for i = 1:7
    errorbar(1:7,squeeze(nanmean(gnmr_mlq_xvLL_rel(:,i,:))),squeeze(nanstd(gnmr_mlq_xvLL_rel(:,i,:)))/sqrt(4),'color',cmap(i,:))
    hold on
end

errorbar(1:6,nanmean(quad_rel_xvLLpl(:,1:4),2),nanstd(quad_rel_xvLL(:,1:4),[],2)/sqrt(4),'k','linewidth',2)


%%
cd ~/James_scripts/GLM/rust_scripts/
save rust_gnm_filt_test_subspace gnm* quad* null* stc* xv


