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
% tuid = '44-29'; hassta=0;  npos=6; nneg=6;
% tuid = '43-03'; hassta=0;  npos=6; nneg=6;
% tuid = '43-21'; hassta=0;  npos=6; nneg=6; %was 6,6
% tuid = '33-27';  %% complicated simple
% tuid = '43-03'; %%??
% tuid = '52-18'; %% simple scell
tuid = '44-29'; %% the complex cell
% tuid = '33-44';  %% ??

%% load STC data
cd ~/Data/rust
dataname = sprintf('stcbf_data-%s.mat',tuid)
eval(['load ' dataname]);
tid =  find(strcmp(tuid,uids));

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

%% Fit models for each XV data set
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
    defmod.lh = 0;
    defmod.lh2 = 0;
    defmod.lnl = 0;
    defmod.lnl2 = 0;
    defmod.hcon = 0;defmod.hmon = 0;defmod.nlcon = 0;defmod.nlmon = 1;
    nSTCbvs = size(STCbvs,2);
    
    foff = flen + length(defmod.h);
    %adjust spike bin indices and get rid of spike bins before foff
    tr_cspkbs = cur_tr_spkbns(cur_tr_spkbns>foff & cur_tr_spkbns<stimlen)-foff;
    
    nmods = 10;
    posmods = 1:5; negmods = 6:10;
    posmod_inds = 1:length(posmods);
    negmod_inds = (length(posdims)+1):(length(posdims)+length(negmods));
    unused_stcs = setdiff(1:nSTCbvs,[posmod_inds negmod_inds]);
    mod_signs = nan(nmods,1);
    mod_signs(posmods) = 1;
    mod_signs(negmods) = -1;
    dim_signs = nan(size(STCbvs,2),1);
    dim_signs(posdims) = 1;
    dim_signs(negdims) = -1;
    
%     fprintf('%d modules\n',nmods);
    
    %initialize module filters
    STCcf_0 = eye(nSTCbvs);
    STCcf_0(:,unused_stcs) = [];
    % STCcf_0(:,nmods+1:end) = [];
    if nmods > nSTCbvs
        n_extra = nmods-nSTCbvs;
        STCcf_0 = [STCcf_0 randn(nSTCbvs,n_extra)];
    end
    for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
    
    %initialize model
    glm_stcb = createGLM0_stcb(STCbvs,STCcf_0,defmod,mod_signs,dim_signs,'test'); %initialize
    [glm_stcb,norm_vals] = normalizeRFs_STCB(glm_stcb,kern_output);
    
    glm_stcb.lambdaW = 0; %sparseness on model weights
    
    lambda_vals = [0 1 2 5 10 20 50 100];
    for ii = 1:length(lambda_vals)
        
        for jj = 1:nmods
            glm_stcb.mods(jj).locLambda = lambda_vals(ii);
        end
        
        %fit model
        stc_posneg_mod{xv,ii} = fitNLHI_stcb_nopsc(glm_stcb,cur_tr_stim,tr_cspkbs,'tots');
        % stc_posneg_mod = fitNLHI_stcb_norefine(glm_stcb,trstim,cspkbs2,'tots');
               
        %adjust spike bin indices and get rid of spike bins before foff
        xvspkbns = cur_xv_spkbns(cur_xv_spkbns>foff & cur_xv_spkbns<xvstimlen)-foff;       
        stc_posneg_xvLL(xv,ii) = getLLGLM_STCBF(stc_posneg_mod{xv,ii},xvkern_output,xvspkbns,'none');
        
        %%
        f1 = plotfo1d_nopsc(stc_posneg_mod{xv,ii});
        set(f1,'Position',[1300 1000 500 800]);
        set(f1,'PaperUnits','centimeters');
        set(f1, 'PaperSize', [15 40]);
        set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))])
        pname = ['~\James_scripts\stc_sparse_test_figs\' tuid sprintf('_locL_%d_IT_%d',lambda_vals(ii),xv)];
        print('-dpng',pname); close
        
    end
end

%% Save output
cd ~/James_scripts/stc_sparse_test_figs/
save stcbf_lambda_vals_mods lambda_vals nmods STCbvs stc_posneg_mod stc_posneg_xvLL nfold

%%
% plot(lambda_vals,stc_posneg_xvLL)
% xlabel('Localization Penalty','fontsize',14)
% ylabel('XV nLL','fontsize',14)

norm_stc_posneg_xvLL = stc_posneg_xvLL./repmat(stc_posneg_xvLL(:,1),1,length(lambda_vals));
figure
plot(lambda_vals,norm_stc_posneg_xvLL)
xlabel('Localization Penalty','fontsize',14)
ylabel('XV nLL','fontsize',14)
hold on
plot(lambda_vals,mean(norm_stc_posneg_xvLL),'k','linewidth',2)
set(gca,'xscale','log')

stc_posneg_LL = cellfun(@(x) x.LL,stc_posneg_mod);
% stc_posneg_LP = cellfun(@(x) x.LP,stc_posneg_mod);
% plot(lambda_vals,stc_posneg_LL)
% xlabel('Localization Penalty','fontsize',14)
% ylabel('XV nLL','fontsize',14)
norm_stc_posneg_LL = stc_posneg_LL./repmat(stc_posneg_LL(:,1),1,length(lambda_vals));
figure
plot(lambda_vals,norm_stc_posneg_LL)
xlabel('Localization Penalty','fontsize',14)
ylabel('nLL','fontsize',14)
hold on
plot(lambda_vals,mean(norm_stc_posneg_LL),'k','linewidth',2)
set(gca,'xscale','log')

