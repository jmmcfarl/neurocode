clear all;

addpath('~/Timm/rust/SparseFilterSelection/')
addpath('~/Timm/MatlabRepository/')

pars = load('~/Data/rust/infos/StandardParametersRust.mat');
flen = pars.flen;
foff = flen + pars.hilen;
ooptions = optimset('MaxFunEvals',100000,'MaxIter',10000);

% datdir = '/Users/timm/myrepo/workspace/DataRepository/rust/stcbar/DATA/';
datdir = '~/Data/rust/stcbar/Data/';
cfiles = dir([datdir,'*stc.mat']); ncells = length(cfiles);
fnames = arrayfun(@(x)x.name,cfiles,'UniformOutput',0);
uids   = cellfun(@(x)[x(2:3),'-',x(6:7)],fnames,'UniformOutput',0);
%stcstart  = 0.8*1e5; ; stcend = 278528;

% tuid = '52-18'; hassta=1;  npos=1; nneg=1;
% tuid = '44-29'; hassta=0;  npos=6; nneg=6;
% tuid = '43-03'; hassta=0;  npos=6; nneg=6;
% tuid = '43-21'; hassta=0;  npos=6; nneg=6; %was 6,6
% tuid = '33-27';  %% complicated simple 
% tuid = '43-03'; %%??
% tuid = '52-18'; %% simple scell
tuid = '44-29'; %% the complex cell
% tuid = '33-44';  %% ??


tid =  find(strcmp(tuid,uids));
eval(['load ',['~/Data/rust/stcbar/Data/',fnames{tid}]]);
% eval(['load ',['~/myrepo/workspace/DataRepository/rust/stcbar/Data/',fnames{tid}]])

psth = spikes_per_frm(:);
rbins = (find(psth>0.5));
nsp = psth(rbins);
spikebins =[];

% for isp =1:length(rbins)
%     spikebins= [spikebins; repmat(rbins(isp),nsp(isp),1)];
% end
% stcstart =1; stcend = size(stim,1);
%
% stcstim   = stim(stcstart:stcend,:);
% stclen = size(stcstim,1);
sdim = size(stim,2);
%
% %readjust spike binning based on new start and end points
% stcspks   = spikebins((spikebins >= stcstart) & (spikebins<stcend))-stcstart +1;
% stcspksK  = stcspks(stcspks>flen)-flen+1; %spike binning adjusted for kernel
% stcspksKH = stcspks(stcspks>foff)-foff+2; %spike binning adjusted for kernel and PSC terms

cd ~/Data/rust
dataname = sprintf('stimspikedat-%s.mat',tuid)
% save(dataname,'stcstim','stcspks','stcspksK','stcspksKH','pars','tuid');
eval(['load ' dataname]);

%% get STC for a clean part of the sequence
[sta, stcs, fstim, evs] = getSTCfilters(stcstim,stcspks,flen,50,50);
rstcs = fliplr(stcs); %reversed STC kernels (suppressive first)

eig_diffs = diff(evs);
npos = find(eig_diffs < 0.005,1,'first');
nneg = find(eig_diffs(end:-1:1) < 0.005,1,'first');

if (hassta) %if want to include STA filter
    effks = [sta',stcs(:,1:npos),rstcs(:,1:nneg)];
else
    effks = [stcs(:,1:npos),rstcs(:,1:nneg)];
end
stcs1 = [stcs(:,1:npos+5) rstcs(:,1:nneg+5)];
plot1dfilterbanks({sta',stcs1,effks},sdim);


%% fit GNM
cd ~/Data/rust/celldata/
allp    = load(['XVTRData-',char(tuid),'.mat']);
stim    = allp.atrstims;
sdim = size(stim,2);
fstimlen = size(stim,1);
spikebins = allp.atrsbins{1};

tslen   = min(2^16,fstimlen); %only up to first 2^16 samples for data-fitting
trstim  = stim(1:tslen,:);

%adjust spike bin indices and get rid of spike bins before foff
cspkbs2 = spikebins(spikebins>foff & spikebins<tslen)-foff;


%%
addpath('~/James_scripts/GLM')
stim    = allp.atrstims;
sdim = size(stim,2);
fstimlen = size(stim,1);
spikebins = allp.atrsbins{1};

nposdims = npos+1; nnegdims = nneg;
posdims = 1:nposdims; negdims = 1:nnegdims;
% STCbvs = [stcs(:,posdims) rstcs(:,negdims)]; %use only expansive subspace
STCbvs = [sta' stcs(:,posdims) rstcs(:,negdims)]; %use only expansive subspace
negdims = negdims + length(posdims);

% STCbvs = [stcs(:,1:4)]; %use only expansive subspace

defmod = pars.defmod; 
%restrict PSC term to delta function
defmod.h(1:end-1) = [];
defmod.SDIM=sdim; 
defmod.locLambda = 10;
% defmod.locLambda = 0;
defmod.lh = 0; 
defmod.lh2 = 0;
defmod.lnl = 0;
defmod.lnl2 = 100; 
defmod.hcon = 0;defmod.hmon = 0;defmod.nlcon = 0;defmod.nlmon = 1;
nSTCbvs = size(STCbvs,2);

foff = flen + length(defmod.h);
%adjust spike bin indices and get rid of spike bins before foff
cspkbs2 = spikebins(spikebins>foff & spikebins<tslen)-foff;


nmods = 9;
posmods = 1:5; negmods = 6:9;
% posmods = []; negmods = [];

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

stc_posneg_mod = fitNLHI_stcb(glm_stcb,trstim,cspkbs2,'tots');
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

stc_posneg_xvLL = getLLGLM_STCBF(stc_posneg_mod,kern_output,xvspkbns,'none');


%%
plotfo1d(stc_posneg_mod)
set(gcf,'PaperUnits','centimeters'); 
set(gcf, 'PaperSize', [15 30]); 
set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
pname = ['~\James_scripts\stc_sparse_test_figs\' tuid]; 
print('-dpng',pname); close
