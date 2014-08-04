clear all; clc
addpath('~/Timm/rust/SparseFilterSelection/')
addpath('~/Timm/MatlabRepository/')
addpath('~/James_scripts/GLM')

pars = load('~/Data/rust/infos/StandardParametersRust.mat');
flen = pars.flen;
% foff = flen + pars.hilen;
ooptions = optimset('MaxFunEvals',100000,'MaxIter',10000);
datdir = '~/Data/rust/stcbar/Data/';
cfiles = dir([datdir,'*stc.mat']); ncells = length(cfiles);
fnames = arrayfun(@(x)x.name,cfiles,'UniformOutput',0);
uids   = cellfun(@(x)[x(2:3),'-',x(6:7)],fnames,'UniformOutput',0);

nposdims = 5; nnegdims = 5;
nSTCbvs = nposdims + nnegdims;
nmods = 10;
posmods = 1:5; negmods = 6:10;

defmod = pars.defmod;
%restrict PSC term to delta function
defmod.h(1:end-1) = [];
defmod.locLambda = 10;
defmod.lh = 0;
defmod.lh2 = 0;
defmod.lnl = 10;
defmod.lnl2 = 10;
defmod.hcon = 0;defmod.hmon = 0;defmod.nlcon = 0;defmod.nlmon = 0;


%% cycle over all datasets
for tid = 1:length(fnames)
    %load data
    tuid = char(uids(tid));
    eval(['load ',['~/Data/rust/stcbar/Data/',fnames{tid}]]);
    
    psth = spikes_per_frm(:);
    rbins = (find(psth>0.5));
    nsp = psth(rbins);
    spikebins =[];
    
    %construct vector of spike bins
    for isp =1:length(rbins)
        spikebins= [spikebins; repmat(rbins(isp),nsp(isp),1)];
    end
    
    stcstim = stim;
    stclen = size(stcstim,1);
    sdim = size(stim,2);
    
    %readjust spike binning based on new start and end points
    stcstart =1; stcend = size(stim,1);
    stcspks   = spikebins((spikebins >= stcstart) & (spikebins<stcend))-stcstart +1;
    
    cd ~/Data/rust
    dataname = sprintf('stimspikedat-%s.mat',tuid)
    save(dataname,'stcstim','stcspks','pars','tuid');
%     eval(['load ' dataname]);
    
    %% get STC for a clean part of the sequence
    [sta, stcs, fstim] = getSTCfilters(stim,stcspks,flen,nposdims,nnegdims);
    rstcs = fliplr(stcs); %reversed STC kernels (suppressive first)    
    
    %% load XV data
    cd ~/Data/rust/celldata/
    allp    = load(['XVTRData-',char(tuid),'.mat']);
    stim    = allp.atrstims;
    sdim = size(stim,2);
    defmod.SDIM=sdim;
    fstimlen = size(stim,1);
    spikebins = allp.atrsbins{1};
    
    tslen   = min(2^16,fstimlen); %only up to first 2^16 samples for data-fitting
    trstim  = stim(1:tslen,:);
        
    %% initialize GNM
    posdims = 1:nposdims; negdims = 1:nnegdims;
    STCbvs = [stcs(:,posdims) rstcs(:,negdims)]; %use only expansive subspace
    negdims = negdims + length(posdims);
            
    foff = flen + length(defmod.h);
    %adjust spike bin indices and get rid of spike bins before foff
    cspkbs2 = spikebins(spikebins>foff & spikebins<tslen)-foff;
        
    mod_signs = nan(nmods,1);
    mod_signs(posmods) = 1;
    mod_signs(negmods) = -1;
    dim_signs = nan(size(STCbvs,2),1);
    dim_signs(posdims) = 1;
    dim_signs(negdims) = -1;
    
    fprintf('%d modules\n',nmods);
    
    STCcf_0 = eye(nSTCbvs);
    STCcf_0(:,nmods+1:end) = [];
    if nmods > nSTCbvs
        n_extra = nmods-nSTCbvs;
        STCcf_0 = [STCcf_0 randn(nSTCbvs,n_extra)];
    end
    for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
    
    glm_stcb = createGLM0_stcb(STCbvs,STCcf_0,defmod,mod_signs,dim_signs,'test'); %initialize
    
    %precompute the stimulus filtered by each STC kernel
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
    [glm_stcb,norm_vals] = normalizeRFs_STCB(glm_stcb,kern_output);
    
    stc_posneg_mod(tid) = fitNLHI_stcb(glm_stcb,trstim,cspkbs2,'tots');
    % stc_posneg_mod = fitNLHI_stcb_norefine(glm_stcb,trstim,cspkbs2,'tots');
      
    %now compute XVLL for this model
    xvstim    = allp.axvstims;
    [fstimlen,sdim] = size(xvstim);
    spikebins = allp.axvsbins{1};
    
    %adjust spike bin indices and get rid of spike bins before foff
    xvspkbns = spikebins(spikebins>foff & spikebins<fstimlen)-foff;
    
    %precompute the stimulus filtered by each STC kernel
    [klen,cur_Nstcbf] = size(STCbvs);
    flen = klen/sdim;
    kern_output = zeros(fstimlen-flen+1,1); %initialize filtered output to be (NT-N_tau+1)xNmods.
    for ikern = 1:cur_Nstcbf;
        %for the given NL module, store the output of the stimulus filtered by
        %the internal kernel
        kern_output(:,ikern)= kfilterInput(xvstim,STCbvs(:,ikern));
    end
    
    stc_posneg_xvLL(tid) = getLLGLM_STCBF(stc_posneg_mod(tid),kern_output,xvspkbns,'none');
    
    
    %%
    plotfo1d(stc_posneg_mod(tid))
    set(gcf,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [15 32]);
    set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
    pname = ['~\James_scripts\stc_sparse_test_figs\' tuid '_modfit_10loc_10nl_10nl2_0h_0h2_1h'];
    print('-dpng',pname); close
    
    plot1dfilterbank_v2(STCbvs,sdim);
    set(gcf,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [7 32]);
    set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
    pname = ['~\James_scripts\stc_sparse_test_figs\' tuid '_stc'];
    print('-dpng',pname); close

end