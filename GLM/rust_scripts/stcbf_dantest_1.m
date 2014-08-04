clear all;

addpath('~/Timm/rust/SparseFilterSelection/')
addpath('~/Timm/MatlabRepository/')
addpath('~/James_scripts/GLM')

pars = load('~/Data/rust/infos/StandardParametersRust.mat');
defmod = pars.defmod;
defmod.h(1:end-1) = [];%restrict PSC term to delta function
flen = pars.flen;
hilen = length(defmod.h);
foff = flen + pars.hilen;
ooptions = optimset('MaxFunEvals',100000,'MaxIter',10000);

% datdir = '/Users/timm/myrepo/workspace/DataRepository/rust/stcbar/DATA/';
datdir = '~/Data/rust/stcbar/Data/';
cfiles = dir([datdir,'*stc.mat']); ncells = length(cfiles);
fnames = arrayfun(@(x)x.name,cfiles,'UniformOutput',0);
uids   = cellfun(@(x)[x(2:3),'-',x(6:7)],fnames,'UniformOutput',0);

% tuid = '52-18'; hassta=1;  npos=1; nneg=1;
tuid = '44-29'; hassta=0;  npos=12; nneg=12;
% tuid = '43-03'; hassta=0;  npos=6; nneg=6;
% tuid = '43-21'; hassta=0;  npos=6; nneg=6; %was 6,6
% tuid = '33-27';  %% complicated simple
% tuid = '43-03'; %%??
% tuid = '52-18'; %% simple scell
% tuid = '44-29'; %% the complex cell
% tuid = '33-44';  %% ??

%% load STC data
cd ~/Data/rust
dataname = sprintf('stcbf_data-%s.mat',tuid)
eval(['load ' dataname]);
tid =  find(strcmp(tuid,uids));
npos = 12; nneg = 12;

rstcs = fliplr(stcs); %reversed STC kernels (suppressive first)
nposdims = npos; nnegdims = nneg;
posdims = 1:nposdims; negdims = 1:nnegdims;
% STCbvs = [sta' stcs(:,posdims) rstcs(:,negdims)]; %use only expansive subspace
STCbvs = [stcs(:,posdims) rstcs(:,negdims)]; %use only expansive subspace
negdims = negdims + length(posdims);

%% create XV data
[stimlen,sdim] = size(stim);
nparts = 18;
partlen = floor(stimlen/nparts);
nfold = 3;
nxvparts = nparts/nfold;

%boundaries of parts
pbounds  = [(0:nparts-1)*partlen+1;(1:nparts)*partlen]';

for i = 1:nfold
    xv_inds{i} = [];
    xv_spkbns{i} = [];
    tr_inds{i} = [];
    tr_spkbns{i} = [];
    
    cur_perm = randperm(nparts);
    cur_xv_parts = sort(cur_perm(1:nxvparts));
    cur_tr_parts = setdiff(1:nparts,cur_xv_parts);
    
    xv_spks = [];
    xv_new_inds = nan(stimlen,1);
    for j = 1:length(cur_xv_parts)
        cur_start = pbounds(cur_xv_parts(j),1);
        cur_stop = pbounds(cur_xv_parts(j),2);
        xv_spks = [xv_spks; spikebins(spikebins >= cur_start & spikebins < cur_stop)];
        
        cur_inds = (cur_start:cur_stop) - cur_start + length(xv_inds{i}) + 1;
        xv_new_inds(cur_start:cur_stop) = cur_inds;
        
        xv_inds{i} = [xv_inds{i} cur_start:cur_stop];
    end
    xv_spkbns{i} = xv_new_inds(xv_spks);
    
    tr_spks = [];
    tr_new_inds = nan(stimlen,1);
    for j = 1:length(cur_tr_parts)
        cur_start = pbounds(cur_tr_parts(j),1);
        cur_stop = pbounds(cur_tr_parts(j),2);
        tr_spks = [tr_spks; spikebins(spikebins >= cur_start & spikebins < cur_stop)];
        
        cur_inds = (cur_start:cur_stop) - cur_start + length(tr_inds{i}) + 1;
        tr_new_inds(cur_start:cur_stop) = cur_inds;
        
        tr_inds{i} = [tr_inds{i} cur_start:cur_stop];
    end
    tr_spkbns{i} = tr_new_inds(tr_spks);
end

nmods = [12 18 24];
for xv = 1:nfold
    
    cur_tr_stim = stim(tr_inds{xv},:);
    cur_tr_spkbns = tr_spkbns{xv};
    
    cur_xv_stim = stim(xv_inds{xv},:);
    cur_xv_spkbns = xv_spkbns{xv};
    
    %% precompute the stimulus filtered by each STC kernel for training
    %% stim
    [klen,Nstcbf] = size(STCbvs);
    flen = klen/sdim;
    stimlen = size(cur_tr_stim,1);
    Nstcbf = size(STCbvs,2);
    kern_output = zeros(stimlen-flen+1,Nstcbf); %initialize filtered output to be (NT-N_tau+1)xNmods.
    for ikern = 1:Nstcbf;
        %for the given NL module, store the output of the stimulus filtered by
        %the internal kernel
        kern_output(:,ikern)= kfilterInput(cur_tr_stim,STCbvs(:,ikern));
    end
    
    %% for XV data
    [xvstimlen,sdim] = size(cur_xv_stim);
    
    %precompute the stimulus filtered by each STC kernel
    [klen,cur_Nstcbf] = size(STCbvs);
    flen = klen/sdim;
    xvkern_output = zeros(xvstimlen-flen+1,1); %initialize filtered output to be (NT-N_tau+1)xNmods.
    for ikern = 1:Nstcbf;
        %for the given NL module, store the output of the stimulus filtered by
        %the internal kernel
        xvkern_output(:,ikern)= kfilterInput(cur_xv_stim,STCbvs(:,ikern));
    end
    
    
    %% Initialize model
    defmod.SDIM=sdim;
    defmod.locLambda = 50;
    defmod.lh = 0;
    defmod.lh2 = 0;
    defmod.lnl = 0;
    defmod.lnl2 = 0;
    defmod.hcon = 0;defmod.hmon = 0;defmod.nlcon = 0;defmod.nlmon = 0;
    nSTCbvs = size(STCbvs,2);
    
    foff = flen + length(defmod.h);
    tr_cspkbs = cur_tr_spkbns(cur_tr_spkbns>foff & cur_tr_spkbns<stimlen)-foff+2;
    
    dim_signs = nan(size(STCbvs,2),1);
    dim_signs(posdims) = 1;
    dim_signs(negdims) = -1;
    
    
    %%
    for jj = 1:length(nmods)
        
        cur_nmods = nmods(jj);
        mod_signs = nan(cur_nmods,1);
        
        posmods = 1:(cur_nmods/2); negmods = (cur_nmods/2+1):cur_nmods;
        posmod_inds = 1:length(posmods);
        if cur_nmods <= Nstcbf
            negmod_inds = (length(posdims)+1):(length(posdims)+length(negmods));
            unused_stcs = setdiff(1:nSTCbvs,[posmod_inds negmod_inds]);
            mod_signs(posmods) = 1;
            mod_signs(negmods) = -1;
        else
            extra_mods = (Nstcbf+1):cur_nmods;
            mod_signs(posdims) = 1;
            mod_signs(negdims) = -1;
            mod_signs(extra_mods(mod(extra_mods,2)==0)) = 1;
            mod_signs(extra_mods(mod(extra_mods,2)==1)) = -1;
            unused_stcs = [];
        end
        
        %initialize on STC dims
                STCcf_0 = eye(nSTCbvs);
                STCcf_0(:,unused_stcs) = [];
                if cur_nmods > nSTCbvs
                    n_extra = cur_nmods-nSTCbvs;
                    STCcf_0 = [STCcf_0 randn(nSTCbvs,n_extra)];
                end
        %random initialization
%         n_extra = cur_nmods;
%         STCcf_0 = randn(nSTCbvs,n_extra);
        
        %make sure expansive subunits start out in expansive subspace, etc
        STCcf_0(negdims,mod_signs==1) = 0;
        STCcf_0(posdims,mod_signs==-1) = 0;
        
        %double space
        STCcf_0 = [STCcf_0 -STCcf_0];
        mod_signs = [mod_signs; mod_signs];
        
        %normalize
        for i = 1:2*cur_nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
        
        %initialize model
        glm_stcb = createGLM0_stcb(STCbvs,STCcf_0,defmod,mod_signs,dim_signs,'test'); %initialize
        [glm_stcb,norm_vals] = normalizeRFs_STCB(glm_stcb,kern_output);
        
        glm_stcb.lambdaW = 0; %sparseness on model weights
        
        %fit model
        stc_posneg_mod{xv,jj} = fitNLHI_stcb_nonlpsc(glm_stcb,cur_tr_stim,tr_cspkbs,'tots');
        % stc_posneg_mod = fitNLHI_stcb_norefine(glm_stcb,trstim,cspkbs2,'tots');
        
        xvspkbns = cur_xv_spkbns(cur_xv_spkbns>foff & cur_xv_spkbns<xvstimlen)-foff;
        stc_posneg_xvLL(xv,jj) = getLLGLM_STCBF(stc_posneg_mod{xv,jj},xvkern_output,xvspkbns,'none');
        
        
        %%
        f1 = plotfo1d_nopsc(stc_posneg_mod{xv,jj},12);
        %         set(f1,'Position',[1300 1000 500 800]);
        set(f1,'PaperUnits','centimeters');
        set(f1, 'PaperSize', [30 40]);
        set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))])
        pname = ['~\James_scripts\stc_sparse_test_figs\dtest_24d_' tuid sprintf('_%dmods_IT%d',cur_nmods,xv)];
        print('-dpng',pname); close
        
    end
        %% for STC
        cur_nmods = 24;
        mod_signs = nan(cur_nmods,1);
    
        posmods = 1:(cur_nmods/2); negmods = (cur_nmods/2+1):cur_nmods;
        posmod_inds = 1:length(posmods);
        if cur_nmods <= Nstcbf
            negmod_inds = (length(posdims)+1):(length(posdims)+length(negmods));
            unused_stcs = setdiff(1:nSTCbvs,[posmod_inds negmod_inds]);
            mod_signs(posmods) = 1;
            mod_signs(negmods) = -1;
        else
            extra_mods = (Nstcbf+1):cur_nmods;
            mod_signs(posdims) = 1;
            mod_signs(negdims) = -1;
            mod_signs(extra_mods(mod(extra_mods,2)==0)) = 1;
            mod_signs(extra_mods(mod(extra_mods,2)==1)) = -1;
            unused_stcs = [];
        end
    
        %initialize on STC dims
        STCcf_0 = eye(nSTCbvs);
        STCcf_0(:,unused_stcs) = [];
        if cur_nmods > nSTCbvs
            n_extra = cur_nmods-nSTCbvs;
            STCcf_0 = [STCcf_0 randn(nSTCbvs,n_extra)];
        end
        %random initialization
        %     n_extra = cur_nmods;
        %     STCcf_0 = randn(nSTCbvs,n_extra);
    
        %make sure expansive subunits start out in expansive subspace, etc
        STCcf_0(negdims,mod_signs==1) = 0;
        STCcf_0(posdims,mod_signs==-1) = 0;
    
        %normalize
        for i = 1:cur_nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
    
        %initialize model
        glm_stcb = createGLM0_stcb(STCbvs,STCcf_0,defmod,mod_signs,dim_signs,'test'); %initialize
        [glm_stcb,norm_vals] = normalizeRFs_STCB(glm_stcb,kern_output);
    
        glm_stcb.lambdaW = 0; %sparseness on model weights
    
        %fit model
        % stc_posneg_mod = fitNLHI_stcb_nonlpsc(glm_stcb,trstim,cspkbs,'tots');
        stc_posneg_mod_stc{xv} = fitNLHI_stcb_norefine(glm_stcb,cur_tr_stim,tr_cspkbs,'tots');
    
        xvspkbns = cur_xv_spkbns(cur_xv_spkbns>foff & cur_xv_spkbns<xvstimlen)-foff;
        stc_posneg_xvLL_stc(xv) = getLLGLM_STCBF(stc_posneg_mod_stc{xv},xvkern_output,xvspkbns,'none');
    
        f1 = plotfo1d_nopsc(stc_posneg_mod_stc{xv},12);
        %         set(f1,'Position',[1300 1000 500 800]);
        set(f1,'PaperUnits','centimeters');
        set(f1, 'PaperSize', [15 40]);
        set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))])
        pname = ['~\James_scripts\stc_sparse_test_figs\dtest_24d_stc_' tuid sprintf('_IT%d',cur_nmods,xv)];
        print('-dpng',pname); close
    
    
end
%%
for xv = 1:nfold
    for jj = 1:length(nmods)
        cur_nmods = nmods(jj);
        [spatial_profiles, temporal_profiles, weights, mod_type, space_COM, temp_COM] = ...
            compute_mod_stats(stc_posneg_mod{xv,jj});
        [~,space_com_ord] = sort(space_COM);
        stc_posneg_mod{xv,jj}.mods = stc_posneg_mod{xv,jj}.mods(space_com_ord);
        
        f1 = plotfo1d_nopsc(stc_posneg_mod{xv,jj},12);
        %         set(f1,'Position',[1300 1000 500 800]);
        set(f1,'PaperUnits','centimeters');
        set(f1, 'PaperSize', [30 40]);
        set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))])
        pname = ['~\James_scripts\stc_sparse_test_figs\dtest_24d_' tuid sprintf('_%dmods_IT%d',cur_nmods,xv) '_sort'];
        print('-dpng',pname); close
        
        
    end
    
            [spatial_profiles, temporal_profiles, weights, mod_type, space_COM, temp_COM] = ...
                compute_mod_stats(stc_posneg_mod_stc{xv});
            [~,space_com_ord] = sort(space_COM);
            stc_posneg_mod_stc{xv}.mods = stc_posneg_mod_stc{xv}.mods(space_com_ord);
    
                        f1 = plotfo1d_nopsc(stc_posneg_mod_stc{xv},12);
                %         set(f1,'Position',[1300 1000 500 800]);
                set(f1,'PaperUnits','centimeters');
                set(f1, 'PaperSize', [30 40]);
                set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))])
                pname = ['~\James_scripts\stc_sparse_test_figs\dtest_24d_stc_' tuid sprintf('_IT%d',xv) '_sort'];
                print('-dpng',pname); close
    
end

cd ~/James_scripts/stc_sparse_test_figs/
dataname = sprintf('dantest_mods_24stc_%s',tuid);
save(dataname,'nfold','nmods','stc_posneg*') 

%%
% load ./dantest_mods_43-21.mat
load ./dantest_mods_24stc_44-29.mat
% load ./dantest_mods.mat
% stcbf_xv_norm = stc_posneg_xvLL./repmat(stc_posneg_xvLL_stc',1,length(nmods));
stcbf_xv_norm = stc_posneg_xvLL - repmat(stc_posneg_xvLL_stc',1,length(nmods));

figure
subplot(2,1,1)
plot(2*nmods,stc_posneg_xvLL','o-')
hold on

% load ./dantest_mods_rand_43-21
% load ./dantest_mods_randinit.mat
% stcbf_rand_xv_norm = stc_posneg_xvLL./repmat(stc_posneg_xvLL_stc',1,4);
% stcbf_rand_xv_norm = stc_posneg_xvLL - repmat(stc_posneg_xvLL_stc',1,4);
% 
% plot(2*nmods,stc_posneg_xvLL,'o--')
plot(2*nmods,repmat(stc_posneg_xvLL_stc,length(nmods),1),'--')
xlabel('Number of modules','fontsize',14)
ylabel('XV nLL','fontsize',14)
xlim([22 52])

subplot(2,1,2)
plot(2*nmods,stcbf_xv_norm','o-')
hold on
% plot(2*nmods,stcbf_rand_xv_norm,'o--')
xlabel('Number of modules','fontsize',14)
ylabel('Relative XV nLL','fontsize',14)
xlim([22 52])
%%
% load ./dantest_mods_43-21.mat
load ./dantest_mods.mat
sdim = stc_posneg_mod{1,1}.mods(1).SDIM;
xv = 3;
jj = 4;
cur_nmods = nmods(jj);
[spatial_profiles, temporal_profiles, weights, mod_type, space_COM, temp_COM] = ...
    compute_mod_stats(stc_posneg_mod{xv,jj});
% [~,space_com_ord] = sort(space_COM);
% stc_posneg_mod{xv,jj}.mods = stc_posneg_mod{xv,jj}.mods(space_com_ord);
% plotfo1d_nopsc(stc_posneg_mod{xv,jj},8)
spatial_profiles = spatial_profiles./repmat(sum(spatial_profiles,2),1,sdim);
relweights = abs(weights)/sum(abs(weights));
pos_comps = find(weights > 0);
neg_comps = find(weights < 0);
% spatial_profiles = spatial_profiles.*repmat(abs(weights),1,sdim);
subplot(3,1,1)
plot(1:sdim,spatial_profiles','k')
hold on
plot(1:sdim,spatial_profiles(pos_comps,:)','b')
plot(1:sdim,spatial_profiles(neg_comps,:)','r')
plot(1:sdim,relweights'*spatial_profiles,'k','linewidth',2)
plot(1:sdim,mean(spatial_profiles),'k--','linewidth',2)
axis tight
xlabel('Pixels','fontsize',14)
ylabel('Temporal Variance','fontsize',14)
title('RF Fit:STC Init','fontsize',14)

% load ./dantest_mods_rand_43-21
load ./dantest_mods_randinit.mat
% xv = 1;
% jj = 4;
cur_nmods = nmods(jj);
[spatial_profiles_rand, temporal_profiles, weights, mod_type, space_COM_rand, temp_COM] = ...
    compute_mod_stats(stc_posneg_mod{xv,jj});
% [~,space_com_ord] = sort(space_COM_rand);
% stc_posneg_mod{xv,jj}.mods = stc_posneg_mod{xv,jj}.mods(space_com_ord);
spatial_profiles_rand = spatial_profiles_rand./repmat(sum(spatial_profiles_rand,2),1,sdim);
relweights = abs(weights)/sum(abs(weights));
pos_comps = find(weights > 0);
neg_comps = find(weights < 0);
% plotfo1d_nopsc(stc_posneg_mod{xv,jj},8)
subplot(3,1,2)
plot(1:sdim,spatial_profiles_rand','k')
hold on
plot(1:sdim,spatial_profiles_rand(pos_comps,:)','b')
plot(1:sdim,spatial_profiles_rand(neg_comps,:)','r')
hold on
plot(1:sdim,relweights'*spatial_profiles_rand,'k','linewidth',2)
plot(1:sdim,mean(spatial_profiles_rand),'k--','linewidth',2)
axis tight
xlabel('Pixels','fontsize',14)
ylabel('Temporal Variance','fontsize',14)
title('RF Fit:Rand Init','fontsize',14)

[spatial_profiles_stc, temporal_profiles, weights, mod_type, space_COM_stc, temp_COM] = ...
    compute_mod_stats(stc_posneg_mod_stc{xv});
% [~,space_com_ord] = sort(space_COM_stc);
% stc_posneg_mod_stc{xv}.mods = stc_posneg_mod_stc{xv}.mods(space_com_ord);
spatial_profiles_stc = spatial_profiles_stc./repmat(sum(spatial_profiles_stc,2),1,sdim);
relweights = abs(weights)/sum(abs(weights));
pos_comps = find(mod_type > 0);
neg_comps = find(mod_type < 0);

subplot(3,1,3)
plot(1:sdim,spatial_profiles_stc','k')
hold on
plot(1:sdim,spatial_profiles_stc(pos_comps,:)','b')
plot(1:sdim,spatial_profiles_stc(neg_comps,:)','r')
hold on
plot(1:sdim,relweights'*spatial_profiles_stc,'k','linewidth',2)
plot(1:sdim,mean(spatial_profiles_stc),'k--','linewidth',2)
axis tight
xlabel('Pixels','fontsize',14)
ylabel('Temporal Variance','fontsize',14)
title('STC Fit','fontsize',14)

%%
sax = 1:sdim;
figure;
for i = 1:cur_nmods
    subplot(3,1,1);hold on
    plot(sax-space_COM(i),spatial_profiles(i,:),'g')
xlim([-10 10])
ylim([0 0.5])
xlabel('Relative position','fontsize',14)
ylabel('Temporal Power','fontsize',14)

    subplot(3,1,2);hold on
    plot(sax-space_COM_rand(i),spatial_profiles_rand(i,:),'r')
xlim([-10 10])
ylim([0 0.5])
xlabel('Relative position','fontsize',14)
ylabel('Temporal Power','fontsize',14)

end

for i = 1:12
        subplot(3,1,3);hold on
    plot(sax-space_COM_stc(i),spatial_profiles_stc(i,:),'b')
xlim([-10 10])
xlabel('Relative position','fontsize',14)
ylabel('Temporal Power','fontsize',14)

end

%%
cd ~/James_scripts/stc_sparse_test_figs/
load ./dantest_mods.mat
% load ./dantest_mods_randinit.mat
% load ./dantest_mods_rand_43-21.mat
% load ./dantest_mods_43-21.mat
sdim = stc_posneg_mod{1,1}.mods(1).SDIM;
klen = length(stc_posneg_mod{1,1}.mods(1).k);
xv = 2;
jj = 4;
cur_nmods = length(stc_posneg_mod{xv,jj}.mods);
[spatial_profiles, temporal_profiles, weights, mod_type, space_COM, temp_COM] = ...
    compute_mod_stats(stc_posneg_mod{xv,jj});
relweights = abs(weights)/sum(abs(weights));
pos_comps = find(weights > 0);
neg_comps = find(weights < 0);
kmat = nan(klen,cur_nmods);
for imod = 1:cur_nmods
    kmat(:,imod) = stc_posneg_mod{xv,jj}.mods(imod).k;
end

npoints = 50; %define spatially interpolated axis for aligned profiles
[arfsE,profmeansE] = alignRFs(kmat(:,pos_comps),sdim,npoints,0);
[arfsS,profmeansS] = alignRFs(kmat(:,neg_comps),sdim,npoints,0);

[evecsE,pcofsE]    = findPCAs(arfsE,npoints,4,1);
[evecsS,pcofsS]    = findPCAs(arfsS,npoints,4,1);

% [arfs,profmeans] = alignRFs(kmat,sdim,npoints);
% [evecs,pcofs]    = findPCAs(arfs,npoints,4);

ncofsE = normvecs(pcofsE(1:2,:));      ncofsS = normvecs(pcofsS(1:2,:));  %normalized PCA scores (first 2 components) for each filter
nphasE = 180*asin(ncofsE(1,:))/pi;     nphasS = 180*asin(ncofsS(1,:))/pi %compute direction of each filter in 2-PCA space




figure
scatter(space_COM(pos_comps),nphasE,weights(pos_comps)*50)
hold on
figure
scatter(space_COM(neg_comps),nphasS,-weights(neg_comps)*50,'r')

%%
cd ~/James_scripts/stc_sparse_test_figs/
load ./dantest_mods.mat
% load ./dantest_mods_rand_43-21.mat
% load ./dantest_mods_43-21.mat
sdim = stc_posneg_mod{1,1}.mods(1).SDIM;
klen = length(stc_posneg_mod{1,1}.mods(1).k);
xv = 1;
jj = 3;
cur_nmods = length(stc_posneg_mod{xv,jj}.mods);
[spatial_profiles, temporal_profiles, weights, mod_type, space_COM, temp_COM] = ...
    compute_mod_stats(stc_posneg_mod{xv,jj});
relweights = abs(weights)/sum(abs(weights));
pos_comps = find(weights > 0);
neg_comps = find(weights < 0);
kmat = nan(klen,cur_nmods);
for imod = 1:cur_nmods
    kmat(:,imod) = stc_posneg_mod{xv,jj}.mods(imod).k;
end

load ./dantest_mods_randinit.mat
sdim = stc_posneg_mod{1,1}.mods(1).SDIM;
klen = length(stc_posneg_mod{1,1}.mods(1).k);
xv = 1;
jj = 4;
cur_nmods = length(stc_posneg_mod{xv,jj}.mods);
[spatial_profiles, temporal_profiles, weights2, mod_type, space_COM2, temp_COM] = ...
    compute_mod_stats(stc_posneg_mod{xv,jj});
relweights = abs(weights2)/sum(abs(weights2));
pos_comps2 = find(weights2 > 0);
neg_comps2 = find(weights2 < 0);
kmat2 = nan(klen,cur_nmods);
for imod = 1:cur_nmods
    kmat2(:,imod) = stc_posneg_mod{xv,jj}.mods(imod).k;
end

combined_kmat = [kmat kmat2];
space_COM_comb = [space_COM; space_COM2];
weights_comb = [weights; weights2];
pos_comps_comb = find(weights_comb > 0);
neg_comps_comb = find(weights_comb < 0);

npoints = 50; %define spatially interpolated axis for aligned profiles
[arfsE,profmeansE] = alignRFs(combined_kmat(:,pos_comps_comb),sdim,npoints,0);
[arfsS,profmeansS] = alignRFs(combined_kmat(:,neg_comps_comb),sdim,npoints,0);

[evecsE,pcofsE]    = findPCAs(arfsE,npoints,4,0);
[evecsS,pcofsS]    = findPCAs(arfsS,npoints,4,0);

% [arfs,profmeans] = alignRFs(kmat,sdim,npoints);
% [evecs,pcofs]    = findPCAs(arfs,npoints,4);

tempE = pcofsE(1,:) + i*pcofsE(2,:);
tempS = pcofsS(1,:) + i*pcofsS(2,:);
nphasE = angle(tempE);
nphasS = angle(tempS);

ncofsE = normvecs(pcofsE(1:2,:));  
ncofsS = normvecs(pcofsS(1:2,:));  %normalized PCA scores (first 2 components) for each filter
nphasE = 180*asin(ncofsE(1,:))/pi;
nphasS = 180*asin(ncofsS(1,:))/pi %compute direction of each filter in 2-PCA space

% nphasE = mod(nphasE,pi);
% nphasS = mod(nphasS,pi);

figure
scatter(space_COM(pos_comps),nphasE(1:length(pos_comps)),weights(pos_comps)*50)
hold on
scatter(space_COM2(pos_comps2),nphasE((length(pos_comps)+1):end),weights2(pos_comps2)*50,'r')

figure
scatter(space_COM(neg_comps),nphasS(1:length(neg_comps)),-weights(neg_comps)*50)
hold on
scatter(space_COM2(neg_comps2),nphasS((length(neg_comps)+1):end),-weights2(neg_comps2)*50,'r')


% figure
% scatter(space_COM(neg_comps),nphasS,-weights(neg_comps)*50,'r')
