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

%% create XV data
[stimlen,sdim] = size(stim);
nparts = 115;
partlen = floor(stimlen/nparts);
nfold = 5;
partsperfold = nparts/nfold;
nxvparts = nparts/nfold;

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
    tr_set = tr_set(randperm(length(tr_set)));
    
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
tr_durs = round([0.5 1 5 10 20 50]*partlen);

%%
klen = flen*sdim;
stim_params.spatial_dims = 1;
stim_params.sdim = sdim;
stim_params.flen = flen;

for xv = 1:nfold
    cur_xv_stim = stim(xv_inds{xv},:);
    cur_xv_stimemb = makeStimRows(cur_xv_stim,flen);
    cur_xv_spkbns = xv_spkbns{xv};
    
    for tr = 1:length(tr_durs)
        fprintf('Fold %d trdur %d\n',xv,tr);
        cur_tr_inds = tr_inds{xv}(1:tr_durs(tr));
        cur_tr_stim = stim(cur_tr_inds,:);
        cur_tr_stimemb = makeStimRows(cur_tr_stim,flen);
        
        cur_tr_spkbns = tr_spkbns{xv};
        cur_tr_spkbns(cur_tr_spkbns > tr_durs(tr)) = [];
        
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
        
        %don't project out STA
        proj_mat = sta'*inv(sta*sta')*sta;
        %     stim_proj = cur_tr_stimemb - cur_tr_stimemb*proj_mat;
        stim_proj = cur_tr_stimemb;
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
        
        %FIT NULL MODEL
        avg_rate = mean(Robs);
        xvpred_rate = ones(1,stimlenxv)*avg_rate;
        trpred_rate = ones(1,stimlen)*avg_rate;
        null_LL(xv,tr) = -sum(Robs.*log(trpred_rate) - trpred_rate)/sum(Robs)
        null_xvLL(xv,tr) = -sum(Robsxv.*log(xvpred_rate) - xvpred_rate)/sum(Robsxv)
        
        %FIT STA MODEL
        kern_types{1} = 'lin';
        sta_mod = createGNM(sta,1,kern_types,[],stim_params);
        g = cur_tr_stimemb*sta;
        sta_mod = fitGNM_spkNL(sta_mod,g,cur_tr_spkbns,0);
        sta_LL(xv,tr) = getLL_GNM(sta_mod,cur_tr_stimemb,cur_tr_spkbns,'none');
        sta_xvLL(xv,tr) = getLL_GNM(sta_mod,cur_xv_stimemb,cur_xv_spkbns,'none');
        
        %FIT STC MODEL
        clear kern_types
        kern_types{1} = 'lin';
        for i = 2:(size(used_stcs,2)+1)
            kern_types{i} = 'quad';
        end
        stc_mod = createGNM(used_stcs,[1 ones(1,npos) -1*ones(1,nneg)],kern_types,[],stim_params);
        g_mat = cur_tr_stimemb*get_k_mat(stc_mod);
        stc_mod = fitGNM_weights(stc_mod,g_mat,cur_tr_spkbns,0);
        [ll, ~, ~, prate, g] = getLL_GNM(stc_mod,cur_tr_stimemb,cur_tr_spkbns,'none');
        stc_mod = fitGNM_spkNL(stc_mod,g,cur_tr_spkbns,0);
        stc_LL(xv,tr) = getLL_GNM(stc_mod,cur_tr_stimemb,cur_tr_spkbns,'none');
        stc_xvLL(xv,tr) = getLL_GNM(stc_mod,cur_xv_stimemb,cur_xv_spkbns,'none');
        
        %FIT GLM
        fprintf('XV set %d, fitting quadratic model\n',xv);
        npos = 1; nneg = 0;
        defmod.lambda_L1x = 20;
        defmod.lambda_d2XT = 200;
        clear kern_types
        kern_types{1} = 'lin';
        init_kerns = randn(klen,npos+nneg);
        init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
        glm(xv,tr) = createGNM(init_kerns,[ones(1,npos) -1*ones(1,nneg)],kern_types,defmod,stim_params);
        glm(xv,tr) = fitGNM_filters(glm(xv,tr),cur_tr_stimemb,cur_tr_spkbns,'none',[],1e-4,1e-6);
        [~, ~, ~, prate, g] = getLL_GNM(glm(xv,tr),cur_tr_stimemb,cur_tr_spkbns,'none');
        glm(xv,tr) = fitGNM_spkNL(glm(xv,tr),g,cur_tr_spkbns,0);
        glm_xvLL(xv,tr) = getLL_GNM(glm(xv,tr),cur_xv_stimemb,cur_xv_spkbns,'none');
        glm_LL(xv,tr) = getLL_GNM(glm(xv,tr),cur_tr_stimemb,cur_tr_spkbns,'none');
        
        %FIT ML QUAD MODEL
        fprintf('XV set %d, fitting quadratic model\n',xv);
        npos = 3; nneg = 6;
        defmod.lambda_L1x = 1;
        defmod.lambda_d2XT = 200;
        clear kern_types
        kern_types{1} = 'lin';
        for i = 2:(npos+nneg)
            kern_types{i} = 'quad';
        end
        init_kerns = randn(klen,npos+nneg);
        init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
        quad_mod(xv,tr) = createGNM(init_kerns,[ones(1,npos) -1*ones(1,nneg)],kern_types,defmod,stim_params);
        for i = (npos+1):(npos+nneg)
            quad_mod(xv,tr).mods(i).lambda_d2XT = 400; %200
            quad_mod(xv,tr).mods(i).lambda_L1x = 40; %60
        end
        quad_mod(xv,tr) = fitGNM_filters(quad_mod(xv,tr),cur_tr_stimemb,cur_tr_spkbns,'none',[],1e-4,1e-6);
        g_mat = cur_tr_stimemb*get_k_mat(quad_mod(xv,tr));
        temp = fitGNM_weights(quad_mod(xv,tr),g_mat,cur_tr_spkbns,0);
        
        [~, ~, ~, prate, g] = getLL_GNM(quad_mod(xv,tr),cur_tr_stimemb,cur_tr_spkbns,'none');
        quad_mod(xv,tr) = fitGNM_spkNL(quad_mod(xv,tr),g,cur_tr_spkbns,0);
        quad_xvLL(xv,tr) = getLL_GNM(quad_mod(xv,tr),cur_xv_stimemb,cur_xv_spkbns,'none');
        quad_LL(xv,tr) = getLL_GNM(quad_mod(xv,tr),cur_tr_stimemb,cur_tr_spkbns,'none');
        
        %%
        mlq_basis = get_k_mat(quad_mod(xv,tr));
        n_mlq_dims = size(mlq_basis,2);
        mlq_klen = 1*n_mlq_dims;
        mlq_stim_params.spatial_dims = 1;
        mlq_stim_params.sdim = n_mlq_dims;
        mlq_stim_params.flen = 1;
        mlq_out = cur_tr_stimemb*mlq_basis;
        xv_mlq_out = cur_xv_stimemb*mlq_basis;
        defmod.lambda_L1x = 0;
        defmod.lambda_d2XT = 0;
        npos = 6; nneg = 6;
        clear kern_types
        for i = 1:(npos+nneg)
            kern_types{i} = 'threshlin';
        end
        n_iter = 10;
        for ii = 1:n_iter
            fprintf('Test fit %d of %d\n',ii,n_iter);
            init_kerns = randn(mlq_klen,npos+nneg);
            init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
            gnm_test(ii) = createGNM(init_kerns,[ones(1,npos) -1*ones(1,nneg)],kern_types,defmod,mlq_stim_params);
            gnm_test(ii) = fitGNM_filters(gnm_test(ii),mlq_out,cur_tr_spkbns,'none',[],1e-4,1e-6);
            gnm_test_xvLL(ii) = getLL_GNM(gnm_test(ii),xv_mlq_out,cur_xv_spkbns,'none');
        end
        [~,best_mod] = min(gnm_test_xvLL);
        gnm_init_mod = gnm_test(best_mod);
        
        %%
        %FIT GNM
        npos = 6; nneg = 6;
        clear kern_types
        for ii = 1:(npos+nneg)
            kern_types{ii} = 'threshlin';
        end
        init_kerns = mlq_basis*get_k_mat(gnm_init_mod);
        fprintf('XV set %d, fitting GNM\n',xv);
        defmod.lambda_L1x = 40;
        defmod.lambda_d2XT = 400;
        defmod.lnl2 = 1000;
        defmod.nlmon = 1;
        defmod.nlx = zeros(23,1);
        gnm_mod(xv,tr) = createGNM(init_kerns,[ones(1,npos) -1*ones(1,nneg)],kern_types,defmod,stim_params);
        for i = (npos+1):(npos+nneg)
            gnm_mod(xv,tr).mods(i).lambda_d2XT = 800; %200
            gnm_mod(xv,tr).mods(i).lambda_L1x = 80; %60
        end
        gnm_mod(xv,tr) = fitGNM_filters(gnm_mod(xv,tr),cur_tr_stimemb,cur_tr_spkbns,'none',[],1e-4,1e-6);
        gnm_xvLL(xv,tr) = getLL_GNM(gnm_mod(xv,tr),cur_xv_stimemb,cur_xv_spkbns,'none');
        [~, ~, ~, ~, g] = getLL_GNM(gnm_mod(xv,tr),cur_tr_stimemb,cur_tr_spkbns,'none');
        gnm_mod(xv,tr) = fitGNM_spkNL(gnm_mod(xv,tr),g,cur_tr_spkbns,0);
        gnm_xvLL(xv,tr) = getLL_GNM(gnm_mod(xv,tr),cur_xv_stimemb,cur_xv_spkbns,'none');
        
        gnmr_mod(xv,tr) = setGNM_NLBFs(gnm_mod(xv,tr),cur_tr_stimemb);
        gnmr_mod(xv,tr) = adjust_all_reg(gnmr_mod(xv,tr),'nltype','uncon');
        gnmr_mod(xv,tr) = fitGNM_internal_NLs(gnmr_mod(xv,tr),cur_tr_stimemb,cur_tr_spkbns,0,2);
        gnmr_xvLL(xv,tr) = getLL_GNM(gnmr_mod(xv,tr),cur_xv_stimemb,cur_xv_spkbns,'none');
        
        gnmr2_mod(xv,tr) = fitGNM_filters(gnmr_mod(xv,tr),cur_tr_stimemb,cur_tr_spkbns,'none',[],1e-4,1e-6);
        gnmr2_mod(xv,tr) = fitGNM_internal_NLs(gnmr2_mod(xv,tr),cur_tr_stimemb,cur_tr_spkbns,0,0);
        gnmr2_xvLL(xv,tr) = getLL_GNM(gnmr2_mod(xv,tr),cur_xv_stimemb,cur_xv_spkbns,'none');
        %
        %     gnmr3(xv) = fitGNM_internal_NLs(gnmr2(xv),cur_tr_stimemb,cur_tr_spkbns,0,1);
        %     gnmr3_xvLL(xv) = getLL_GNM(gnmr3(xv),cur_xv_stimemb,cur_xv_spkbns,'none');
        %
        %     gnmr4(xv) = fitGNM_filters(gnmr3(xv),cur_tr_stimemb,cur_tr_spkbns,'none',[],1e-4,1e-6);
        %     gnmr4_xvLL(xv) = getLL_GNM(gnmr4(xv),cur_xv_stimemb,cur_xv_spkbns,'none');
        
    end
end

%%
sta_rel_xvLL = (null_xvLL - sta_xvLL)/log(2);
stc_rel_xvLL = (null_xvLL - stc_xvLL)/log(2);
glm_rel_xvLL = (null_xvLL - glm_xvLL)/log(2);
quad_rel_xvLL = (null_xvLL - quad_xvLL)/log(2);
gnm_rel_xvLL = (null_xvLL - gnm_xvLL)/log(2);
gnmr_rel_xvLL = (null_xvLL - gnmr_xvLL)/log(2);
gnmr2_rel_xvLL = (null_xvLL - gnmr2_xvLL)/log(2);

Y = [glm_rel_xvLL(:);stc_rel_xvLL(:);quad_rel_xvLL(:);gnmr2_rel_xvLL(:)];
G = [ones(nfold,1);2*ones(nfold,1);3*ones(nfold,1);4*ones(nfold,1)];
figure
boxplot(Y,G)

Y = [glm_rel_xvLL(:) stc_rel_xvLL(:) quad_rel_xvLL(:) gnmr2_rel_xvLL(:)];
Y = [glm_rel_xvLL(:) stc_rel_xvLL(:) quad_rel_xvLL(:) gnmr2_rel_xvLL(:)];
hold on
for i = 1:nfold
    plot(1:4,Y(i,:),'ko-')
end


stc_rel_xvLL = (glm_xvLL - stc_xvLL)/log(2);shg

quad_rel_xvLL = (glm_xvLL - quad_xvLL)/log(2);
gnmr2_rel_xvLL = (glm_xvLL - gnmr2_xvLL)/log(2);

figure
Y = [stc_rel_xvLL(:);quad_rel_xvLL(:);gnmr2_rel_xvLL(:)];
G = [2*ones(nfold,1);3*ones(nfold,1);4*ones(nfold,1)];
boxplot(Y,G)



