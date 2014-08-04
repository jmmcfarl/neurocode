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
nfold = 10;
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
cd ~/James_scripts/GLM/rust_scripts/
load ./rust_gnm_filt_v4.mat
% for xv = 1:nfold
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
    
        
    %FIT NULL MODEL
    avg_rate = mean(Robs);
    xvpred_rate = ones(1,stimlenxv)*avg_rate;
    trpred_rate = ones(1,stimlen)*avg_rate;
    null_LL(xv) = -sum(Robs.*log(trpred_rate) - trpred_rate)/sum(Robs)
    null_xvLL(xv) = -sum(Robsxv.*log(xvpred_rate) - xvpred_rate)/sum(Robsxv)
    
    %% FIT ML QUAD MODEL for subspace identification
    npos = 4; nneg = 6;
    defmod.lambda_L1x = 0.5;
    defmod.lambda_d2XT = 150;
    clear kern_types
    kern_types{1} = 'lin';
    for i = 2:(npos+nneg)
        kern_types{i} = 'quad';
    end
    init_kerns = 0.1*randn(klen,npos+nneg);
    init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
    quad_mod(xv) = createGNM(init_kerns,[ones(1,npos) -1*ones(1,nneg)],kern_types,defmod,stim_params);
    for i = (npos+1):(npos+nneg)
        quad_mod(xv).mods(i).lambda_d2XT = 400; %200
        quad_mod(xv).mods(i).lambda_L1x = 10; %60
    end
    
    quad_mod(xv) = fitGNM_filters(quad_mod(xv),cur_tr_stimemb,cur_tr_spkbns,'none',[],1e-4,1e-6);
    quad_xvLL(xv) = getLL_GNM(quad_mod(xv),cur_xv_stimemb,cur_xv_spkbns,'none');
    quad_LL(xv) = getLL_GNM(quad_mod(xv),cur_tr_stimemb,cur_tr_spkbns,'none');
    
    for k = 1:2
        [~, ~, ~, prate, g] = getLL_GNM(quad_mod(xv),cur_tr_stimemb,cur_tr_spkbns,'none');
        quad_mod(xv) = fitGNM_spkNL(quad_mod(xv),g,cur_tr_spkbns,0);
        quad_mod(xv) = fitGNM_filters(quad_mod(xv),cur_tr_stimemb,cur_tr_spkbns,'none',[],1e-4,1e-6);
        quad_xvLL(xv) = getLL_GNM(quad_mod(xv),cur_xv_stimemb,cur_xv_spkbns,'none');
    end
%          [~, ~, ~, prate, g,quad_intg] = getLL_GNM(quad_mod(xv),cur_tr_stimemb,cur_tr_spkbns,'none');
   
   %%
    mlq_basis = get_k_mat(quad_mod(xv));
    n_mlq_dims = size(mlq_basis,2);
    mlq_klen = 1*n_mlq_dims;
    mlq_stim_params.spatial_dims = 1;
    mlq_stim_params.sdim = n_mlq_dims;
    mlq_stim_params.flen = 1;
    mlq_out = cur_tr_stimemb*mlq_basis;
    xv_mlq_out = cur_xv_stimemb*mlq_basis;    
    defmod.lambda_L1x = 0;
    defmod.lambda_d2XT = 0;
    for j = 1:8
        for k = 1:8
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

            gnmr_mlq(j,k,xv) = fitGNM_internal_NLs(gnmr_mlq(j,k,xv),mlq_out,cur_tr_spkbns,1,2);
            gnmr_mlq(j,k,xv) = fitGNM_filters(gnmr_mlq(j,k,xv),mlq_out,cur_tr_spkbns,'none',[],1e-4,1e-6);
            [~, ~, ~, prate, g] = getLL_GNM(gnmr_mlq(j,k,xv),mlq_out,cur_tr_spkbns,'none');
            gnmr_mlq(j,k,xv) = fitGNM_spkNL(gnmr_mlq(j,k,xv),g,cur_tr_spkbns,0);
            gnmr_mlq_xvLL(j,k,xv) = getLL_GNM(gnmr_mlq(j,k,xv),xv_mlq_out,cur_xv_spkbns,'none');
            gnmr_mlq_LL(j,k,xv) = getLL_GNM(gnmr_mlq(j,k,xv),mlq_out,cur_tr_spkbns,'none');
            
            gnmr2_mlq(j,k,xv) = fitGNM_internal_NLs(gnmr_mlq(j,k,xv),mlq_out,cur_tr_spkbns,1,2);
            gnmr2_mlq(j,k,xv) = fitGNM_filters(gnmr2_mlq(j,k,xv),mlq_out,cur_tr_spkbns,'none',[],1e-4,1e-6);
            [~, ~, ~, prate, g] = getLL_GNM(gnmr2_mlq(j,k,xv),mlq_out,cur_tr_spkbns,'none');
            gnmr2_mlq(j,k,xv) = fitGNM_spkNL(gnmr2_mlq(j,k,xv),g,cur_tr_spkbns,0);
            gnmr2_mlq_xvLL(j,k,xv) = getLL_GNM(gnmr2_mlq(j,k,xv),xv_mlq_out,cur_xv_spkbns,'none');
            gnmr2_mlq_LL(j,k,xv) = getLL_GNM(gnmr2_mlq(j,k,xv),mlq_out,cur_tr_spkbns,'none');
            
            gnmr3_mlq(j,k,xv) = setGNM_NLBFs(gnmr2_mlq(j,k,xv),mlq_out);
            gnmr3_mlq(j,k,xv) = fitGNM_internal_NLs(gnmr3_mlq(j,k,xv),mlq_out,cur_tr_spkbns,1,0);
            gnmr3_mlq_xvLL(j,k,xv) = getLL_GNM(gnmr3_mlq(j,k,xv),xv_mlq_out,cur_xv_spkbns,'none');
            gnmr3_mlq_LL(j,k,xv) = getLL_GNM(gnmr3_mlq(j,k,xv),mlq_out,cur_tr_spkbns,'none');
            [~, ~, ~, ~, ~,gnm_intg] = getLL_GNM(gnmr3_mlq(j,k,xv),mlq_out,cur_tr_spkbns,'none');
        
        end
    end
    
%%
% mlq_basis = get_k_mat(quad_mod(xv));
% for j = 1:8
%     for k = 1:8
%             fprintf('%d Pos, %d Neg\n',j,k);
%         npos = j; nneg = k;
%         defmod.lambda_L1x = 50;
%         defmod.lambda_d2XT = 50; %200
%         clear kern_types
%         for ii = 1:(npos+nneg)
%             kern_types{ii} = 'threshlin';
%         end
%         init_kerns = mlq_basis*get_k_mat(gnmr3_mlq(j,k,xv));
% %         fgnm(j,k,xv) = createGNM(init_kerns,[ones(1,npos) -1*ones(1,nneg)],kern_types,defmod,stim_params);
% %         fgnm(j,k,xv) = fitGNM_filters(fgnm(j,k,xv),cur_tr_stimemb,cur_tr_spkbns,'none',[],1e-4,1e-6);
% %         [~, ~, ~, prate, g] = getLL_GNM(fgnm(j,k,xv),cur_tr_stimemb,cur_tr_spkbns,'none');
% %         fgnm(j,k,xv) = fitGNM_spkNL(fgnm(j,k,xv),g,cur_tr_spkbns,0);
% %         fgnm_xvLL(j,k,xv) = getLL_GNM(fgnm(j,k,xv),cur_xv_stimemb,cur_xv_spkbns,'none');
% %         
% %         fgnmr(j,k,xv) = adjust_all_reg(fgnm(j,k,xv),'nlx',zeros(25,1));
% %         fgnmr(j,k,xv) = setGNM_NLBFs(fgnmr(j,k,xv),cur_tr_stimemb);
% %         fgnmr(j,k,xv) = adjust_all_reg(fgnmr(j,k,xv),'lnl2',400); %100
% %         fgnmr(j,k,xv) = adjust_all_reg(fgnmr(j,k,xv),'nlmon',1);
% %         fgnmr(j,k,xv) = adjust_all_reg(fgnmr(j,k,xv),'nltype','uncon');
% %         fgnmr(j,k,xv) = fitGNM_internal_NLs(fgnmr(j,k,xv),cur_tr_stimemb,cur_tr_spkbns,1,0);
% %         
% %         [~, ~, ~, prate, g] = getLL_GNM(fgnmr(j,k,xv),cur_tr_stimemb,cur_tr_spkbns,'none');
% %         fgnmr(j,k,xv) = fitGNM_spkNL(fgnmr(j,k,xv),g,cur_tr_spkbns,0);
% %         fgnmr_xvLL(j,k,xv) = getLL_GNM(fgnmr(j,k,xv),cur_xv_stimemb,cur_xv_spkbns,'none');
% %         
% %         fgnmr2(j,k,xv) = fitGNM_filters(fgnmr(j,k,xv),cur_tr_stimemb,cur_tr_spkbns,'none',[],1e-4,1e-6);
% %         fgnmr2(j,k,xv) = fitGNM_internal_NLs(fgnmr2(j,k,xv),cur_tr_stimemb,cur_tr_spkbns,1,0);
% %         [~, ~, ~, prate, g] = getLL_GNM(fgnmr2(j,k,xv),cur_tr_stimemb,cur_tr_spkbns,'none');
% %         fgnmr2(j,k,xv) = fitGNM_spkNL(fgnmr2(j,k,xv),g,cur_tr_spkbns,0);
% %         fgnmr2_xvLL(j,k,xv) = getLL_GNM(fgnmr2(j,k,xv),cur_xv_stimemb,cur_xv_spkbns,'none');
% % 
%         fgnm3(j,k,xv) = createGNM(init_kerns,[ones(1,npos) -1*ones(1,nneg)],kern_types,defmod,stim_params);
%         fgnm3(j,k,xv).spk_alpha = gnmr3_mlq(j,k,xv).spk_alpha;
%         fgnm3(j,k,xv).spk_beta = gnmr3_mlq(j,k,xv).spk_beta;
%         fgnm3(j,k,xv).spk_theta = gnmr3_mlq(j,k,xv).spk_theta;
%         for ii = 1:(npos+nneg)
%            fgnm3(j,k,xv).mods(ii).nlx = gnmr3_mlq(j,k,xv).mods(ii).nlx; 
%            fgnm3(j,k,xv).mods(ii).nly = gnmr3_mlq(j,k,xv).mods(ii).nly; 
%         end
%          fgnm3(j,k,xv) = adjust_all_reg(fgnm3(j,k,xv),'lnl2',400); %100
%         fgnm3(j,k,xv) = adjust_all_reg(fgnm3(j,k,xv),'nlmon',1);
%         fgnm3(j,k,xv) = adjust_all_reg(fgnm3(j,k,xv),'nltype','uncon');
%        getLL_GNM(fgnm3(j,k,xv),cur_xv_stimemb,cur_xv_spkbns,'none')
%         fgnm3(j,k,xv) = fitGNM_filters(fgnm3(j,k,xv),cur_tr_stimemb,cur_tr_spkbns,'none',[],1e-4,1e-6);
%        getLL_GNM(fgnm3(j,k,xv),cur_xv_stimemb,cur_xv_spkbns,'none')
%         fgnm3(j,k,xv) = fitGNM_internal_NLs(fgnm3(j,k,xv),cur_tr_stimemb,cur_tr_spkbns,1,0);
%         [~, ~, ~, prate, g] = getLL_GNM(fgnm3(j,k,xv),cur_tr_stimemb,cur_tr_spkbns,'none');
%         fgnm3(j,k,xv) = fitGNM_spkNL(fgnm3(j,k,xv),g,cur_tr_spkbns,0);
%         fgnm3_xvLL(j,k,xv) = getLL_GNM(fgnm3(j,k,xv),cur_xv_stimemb,cur_xv_spkbns,'none');
%        
% %         gnmr_xvLL(xv) = getLL_GNM(gnmr(xv),cur_xv_stimemb,cur_xv_spkbns,'none');
% %         gnmr_LL(xv) = getLL_GNM(gnmr(xv),cur_tr_stimemb,cur_tr_spkbns,'none');
% %         gnmr2(xv) = fitGNM_filters(gnmr(xv),cur_tr_stimemb,cur_tr_spkbns,'none',[],1e-4,1e-6);
% %         gnmr2_xvLL(xv) = getLL_GNM(gnmr2(xv),cur_xv_stimemb,cur_xv_spkbns,'none');
%     end
% end
cd ~/James_scripts/GLM/rust_scripts/
save rust_gnm_filt_v5_full fgnm* gnm* quad* null*

end

%%

%%
use_xv = 1:10;
xv_imp = -(bsxfun(@minus,permute(gnmr3_mlq_xvLL(:,:,use_xv),[3 1 2]),null_xvLL(use_xv)'))/log(2);
% xv_imp = -(bsxfun(@minus,permute(fgnmr_xvLL(:,:,use_xv),[3 1 2]),null_xvLL(use_xv)'))/log(2);
% xv_imp = -(bsxfun(@minus,permute(fgnm3_xvLL(:,:,use_xv),[3 1 2]),null_xvLL(use_xv)'))/log(2);
% xv_imp2 = -(bsxfun(@minus,permute(gnmr2_mlq_xvLL(:,:,use_xv),[3 1 2]),null_xvLL(use_xv)'))/log(2);

quad_imp = -(quad_xvLL(use_xv)-null_xvLL(use_xv))/log(2);

% xv_imp = -(bsxfun(@minus,permute(gnmr3_mlq_xvLL(:,:,use_xv),[3 1 2]),squeeze(gnmr3_mlq_xvLL(2,1,use_xv))))/log(2);
% quad_imp = -(quad_xvLL(use_xv)-squeeze(gnmr3_mlq_xvLL(2,1,use_xv))')/log(2);

avg_xv_imp = squeeze(mean(xv_imp));
sem_xv_imp = squeeze(std(xv_imp))/sqrt(length(use_xv));
avg_quad_imp = mean(quad_imp);
[XX,YY] = meshgrid(1:8,2:8);
ZZ = avg_xv_imp(2:end,:); 
EE = sem_xv_imp(2:end,:); 

figure
plot3d_errorbars(XX(:),YY(:),ZZ(:),EE(:))
hold on
surf(1:8,2:8,avg_xv_imp(2:end,:))
% alpha(.4)
hold on
% plot3(XX(:),YY(:),ZZ(:),'bo','markersize',12,'linewidth',2)
% surf(repmat(avg_quad_imp,8,8))
plot3(1,8,avg_quad_imp,'ko','markersize',18,'linewidth',2)
%  caxis([0.83 0.89])
%  zlim([0.84 0.89])
xlim([0.8 8.2])
ylim([0.8 8.2])


%%
use_xv = 1:10;
xv_imp = -(bsxfun(@minus,permute(gnmr3_mlq_xvLL(:,:,use_xv),[3 1 2]),null_xvLL(use_xv)'))/log(2);
% xv_imp2 = -(bsxfun(@minus,permute(gnmr2_mlq_xvLL(:,:,use_xv),[3 1 2]),null_xvLL(use_xv)'))/log(2);

quad_imp = -(quad_xvLL(use_xv)-null_xvLL(use_xv))/log(2);

xv_imp = bsxfun(@minus,xv_imp,quad_imp');

avg_xv_imp = squeeze(mean(xv_imp));
sem_xv_imp = squeeze(std(xv_imp))/sqrt(length(use_xv));
[XX,YY] = meshgrid(1:8,2:8);
ZZ = avg_xv_imp(2:end,:); 
EE = sem_xv_imp(2:end,:); 

figure
plot3d_errorbars(XX(:),YY(:),ZZ(:),EE(:))
hold on
surf(1:8,2:8,avg_xv_imp(2:end,:))
hold on
xlim([0.8 8.2])
ylim([0.8 8.2])

avg_xv_imp = squeeze(mean(xv_imp));
sem_xv_imp = squeeze(std(xv_imp))/sqrt(length(use_xv));
[XX,YY] = meshgrid(1:8,1:8);
ZZ = avg_xv_imp(1:end,:); 
EE = sem_xv_imp(1:end,:); 

figure
plot3d_errorbars(XX(:),YY(:),ZZ(:),EE(:))
hold on
surf(1:8,1:8,avg_xv_imp(1:end,:))
hold on
xlim([0.8 8.2])
ylim([0.8 8.2])

%%
close all

% use_xv = 3;
% use_j = 7;exc_set = 1:use_j;
% use_i = 7;inh_set = (use_j+1):(use_j+use_i);
use_xv = 1;
mlq_basis = get_k_mat(quad_mod(use_xv));

use_j = 4;exc_set = 1:use_j;
use_i = 4;inh_set = (use_j+1):(use_j+use_i);

cur_mod = gnmr3_mlq(use_j,use_i,1);


k_mat = get_k_mat(cur_mod);
k_mat = mlq_basis*k_mat;

k_norms = sqrt(sum(k_mat.^2));

[~,exc_ord] = sort(k_norms(exc_set),'descend');
[~,inh_ord] = sort(k_norms(inh_set),'descend');
cur_mod.mods(exc_set) = cur_mod.mods(exc_set(exc_ord));
cur_mod.mods(inh_set) = cur_mod.mods(inh_set(inh_ord));

plotfo1d_nopsc_subspace(cur_mod,3,'centered',mlq_basis,sdim)

% use_j = 4;exc_set = 1:use_j;
% use_i = 4;inh_set = (use_j+1):(use_j+use_i);
use_j = 6;exc_set = 1:use_j;
use_i = 6;inh_set = (use_j+1):(use_j+use_i);

cur_mod = gnmr3_mlq(use_j,use_i,1);


k_mat = get_k_mat(cur_mod);
k_mat = mlq_basis*k_mat;

k_norms = sqrt(sum(k_mat.^2));

[~,exc_ord] = sort(k_norms(exc_set),'descend');
[~,inh_ord] = sort(k_norms(inh_set),'descend');
cur_mod.mods(exc_set) = cur_mod.mods(exc_set(exc_ord));
cur_mod.mods(inh_set) = cur_mod.mods(inh_set(inh_ord));

plotfo1d_nopsc_subspace(cur_mod,3,'centered',mlq_basis,sdim)

%%
for use_xv = 1:10;
    fprintf('XV set %d of %d\n',use_xv,10);
use_gnm = gnmr3_mlq(end,end,use_xv);

cur_xv_stim = stim(xv_inds{use_xv},:);
cur_xv_stimemb = makeStimRows(cur_xv_stim,flen);
cur_xv_spkbns = xv_spkbns{use_xv};

cur_tr_stim = stim(tr_inds{use_xv},:);
cur_tr_stimemb = makeStimRows(cur_tr_stim,flen);
cur_tr_spkbns = tr_spkbns{use_xv};

stimlen = size(cur_tr_stimemb,1);
Robs = zeros(1,stimlen);
ftable = tabulate(cur_tr_spkbns);
Robs(ftable(:,1)) = ftable(:,2);
stimlenxv = size(cur_xv_stimemb,1);
Robsxv = zeros(1,stimlenxv);
ftable = tabulate(cur_xv_spkbns);
Robsxv(ftable(:,1)) = ftable(:,2);

mlq_basis = get_k_mat(quad_mod(use_xv));

cur_k = mlq_basis*get_k_mat(use_gnm);
g_mat = cur_tr_stimemb*cur_k;

xv_mlq_out = cur_xv_stimemb*mlq_basis;

l1_strengths = [0 1 5 10 20 50 100 250 500 1000 5000 10000 50000];
for ll = 1:length(l1_strengths)
    fprintf('L1 %d of %d\n',ll,length(l1_strengths));
    use_gnmf(use_xv,ll) = fitGNM_weights(use_gnm,g_mat,cur_tr_spkbns,0,[],l1_strengths(ll));
    w_set(use_xv,ll,:) = [use_gnmf(use_xv,ll).mods(:).w];
    gnm_l1_xvLL(use_xv,ll) = getLL_GNM(use_gnmf(use_xv,ll),xv_mlq_out,cur_xv_spkbns,'none');
end
end

cur_imp = -bsxfun(@minus,gnm_l1_xvLL,null_xvLL');
figure
errorbar(l1_strengths,mean(cur_imp),std(cur_imp)/sqrt(10));

