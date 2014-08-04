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

%% load XV data
cd ~/Data/rust/celldata/
allp    = load(['XVTRData-',char(tuid),'.mat']);
stim    = allp.atrstims;
sdim = size(stim,2);
fstimlen = size(stim,1);
spikebins = allp.atrsbins{1};
trstim  = stim;

%% precompute the stimulus filtered by each STC kernel
[klen,Nstcbf] = size(STCbvs);
flen = klen/sdim;
stimlen = size(trstim,1);
Nstcbf = size(STCbvs,2);
kern_output = zeros(stimlen-flen+1,Nstcbf); %initialize filtered output to be (NT-N_tau+1)xNmods.
for ikern = 1:Nstcbf;
    %for the given NL module, store the output of the stimulus filtered by
    %the internal kernel
    kern_output(:,ikern)= kfilterInput(trstim,STCbvs(:,ikern));
end

%% XV data
xvstim    = allp.axvstims;
[fstimlen,sdim] = size(xvstim);
xvspikebins = allp.axvsbins{1};

%adjust spike bin indices and get rid of spike bins before foff
xvspkbns = xvspikebins(xvspikebins>foff & xvspikebins<fstimlen)-foff;

%precompute the stimulus filtered by each STC kernel
[klen,cur_Nstcbf] = size(STCbvs);
flen = klen/sdim;
xvkern_output = zeros(fstimlen-flen+1,1); %initialize filtered output to be (NT-N_tau+1)xNmods.
for ikern = 1:cur_Nstcbf;
    %for the given NL module, store the output of the stimulus filtered by
    %the internal kernel
    xvkern_output(:,ikern)= kfilterInput(xvstim,STCbvs(:,ikern));
end

%% Initialize model
defmod.SDIM=sdim;
defmod.locLambda = 10;
defmod.lh = 0;
defmod.lh2 = 0;
defmod.lnl = 0;
defmod.lnl2 = 0;
defmod.hcon = 0;defmod.hmon = 0;defmod.nlcon = 0;defmod.nlmon = 1;
nSTCbvs = size(STCbvs,2);

foff = flen + length(defmod.h);
%adjust spike bin indices and get rid of spike bins before foff
cspkbs2 = spikebins(spikebins>foff & spikebins<stimlen)-foff;

dim_signs = nan(size(STCbvs,2),1);
dim_signs(posdims) = 1;
dim_signs(negdims) = -1;

nmods = [2 4 6 8 10 12 16 20];
n_iter = 3;

%%
for ii = 1:length(nmods)
    for jj = 1:n_iter
        cur_nmods = nmods(ii);
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
        
        if jj == 1 %for first iteration initialize on STC dims
            %initialize module filters
            STCcf_0 = eye(nSTCbvs);
            STCcf_0(:,unused_stcs) = [];
            if cur_nmods > nSTCbvs
                n_extra = cur_nmods-nSTCbvs;
                STCcf_0 = [STCcf_0 randn(nSTCbvs,n_extra)];
            end
        else %for rest use random initialization
            n_extra = cur_nmods;
            STCcf_0 = randn(nSTCbvs,n_extra);
        end
        
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
        stc_posneg_mod{jj,ii} = fitNLHI_stcb_nopsc(glm_stcb,trstim,cspkbs2,'tots');
        % stc_posneg_mod = fitNLHI_stcb_norefine(glm_stcb,trstim,cspkbs2,'tots');
        
        stc_posneg_xvLL(jj,ii) = getLLGLM_STCBF(stc_posneg_mod{jj,ii},xvkern_output,xvspkbns,'none');
        
        %%
        f1 = plotfo1d_nopsc(stc_posneg_mod{jj,ii});
        %         set(f1,'Position',[1300 1000 500 800]);
        set(f1,'PaperUnits','centimeters');
        set(f1, 'PaperSize', [15 4*cur_nmods]);
        set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))])
        pname = ['~\James_scripts\stc_sparse_test_figs\nmods_test\' tuid sprintf('_%dmods_IT%d',cur_nmods,jj)];
        print('-dpng',pname); close 
        
    end
end
