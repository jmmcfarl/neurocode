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
tuid = '43-21'; hassta=0;  npos=6; nneg=6; %was 6,6


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
[sta, stcs, fstim] = getSTCfilters(stcstim,stcspks,flen,20,20);
rstcs = fliplr(stcs); %reversed STC kernels (suppressive first)

if (hassta) %if want to include STA filter
    effks = [sta',stcs(:,1:npos-1),rstcs(:,1:nneg)];    
else
    effks = [stcs(:,1:npos),rstcs(:,1:nneg)];
end
% plot1dfilterbanks({sta',stcs,effks},sdim);


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

% STCcf_0 = randn(nSTCbvs,nmods);
% for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end; 

%%

% 
% kern_length = size(STCbvs,1);
% kern_h = kern_length/sdim;
% sigma_space = 1;
% sigma_time = 0.1;
% scale = 0.75;
% noise_sigma = 1;
% 
% [Xspace,Xtime] = meshgrid(1:sdim,1:kern_h);
% space_dist = (bsxfun(@minus,Xspace(:),Xspace(:)')).^2;
% time_dist = (bsxfun(@minus,Xtime(:),Xtime(:)')).^2;
% sigma = exp(-scale - 1/2*(space_dist/sigma_space^2 + time_dist/sigma_time^2));
% 
% for i = 1:nmods
% cur_bv = loc_kerns(:,i);
% post_cov = inv(1/noise_sigma^2*eye(kern_length) + inv(sigma));
% post_mean(:,i) = 1/noise_sigma^2*post_cov*cur_bv;
% end
% % figure
% % pcolor(reshape(cur_bv,kern_h,sdim));shading flat; colorbar
% % figure
% % pcolor(reshape(post_mean,kern_h,sdim));shading flat; colorbar
% 
% 
% % opts = optimset('GradObj','on','Algorithm','active-set','Display','iter','MaxIter',2000,'MaxFunEvals',20,'TolFun',1e-4);
% % [k_inferred,fval] = fminunc(@(x) smooth_stcb_posterior(x,cur_bv,noise_sigma,inv(sigma)),cur_bv,opts);


%%
addpath('~/James_scripts/GLM')
STCbvs = stcs(:,1:6); %use only expansive subspace

defmod = pars.defmod; defmod.chmon = -1;
defmod.SDIM=sdim; defmod.locLambda = 10;
defmod.lh = 0; defmod.lnl = 0;
defmod.lnl2 = 10; defmod.lh2 = 10;
nSTCbvs = size(STCbvs,2);
nmods = 4;

% defmod.spatial_sigma = 3; defmod.prior_scale = 1;

% [loc_kerns loc_cfs] = min_entropy_filts_projpurs(STCbvs,nmods,sdim);

% STCcf_0 = randn(nSTCbvs,nmods);
STCcf_0 = eye(nSTCbvs);
STCcf_0(:,nmods+1:end) = [];
if nmods > nSTCbvs
   n_extra = nmods-nSTCbvs;
   STCcf_0 = [STCcf_0 randn(nSTCbvs,n_extra)];
end
for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end; 


% glm_stcb = createGLM0_stcb(STCbvs,loc_cfs,defmod,'test'); %initialize
glm_stcb = createGLM0_stcb(STCbvs,STCcf_0,defmod,'test'); %initialize

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

close all
plotfo1d(glm_stcb)
drawnow;

disp_type = 'tots';
clc;
f_glm_stcb = fitNLHI_stcb(glm_stcb,trstim,cspkbs2,disp_type);
plotfo1d(f_glm_stcb);
drawnow;




