clear all;
addpath('~/Timm/tim1/functions/') %for getMatchProf.m
addpath('~/James_scripts/GLM/')

load /Users/James/Data/rust/infos/sigstcs.mat %load data about numbers of significant E and I STCs
pars = load('~/Data/rust/infos/StandardParametersRust.mat'); %load pars structure
x = load('~/Data/rust/infos/RustSignifSTCs.mat'); %load precomputed STCs (including localized versions)
allcells = x.allcells;
clear x

%load current data
% stim is Nxd stimulus time series
% spikebins is a vector containing the bin of each spike (length Nspikes)
tcellid = 1; %test cell
datadir = '~/Data/rust/celldata/';
cd(datadir);
tuid    = allcells{tcellid}.uid
disp(['*** working on cell ',tuid])
allp    = load(['XVTRData-',char(tuid),'.mat']);
stim    = allp.atrstims; sdim = size(stim,2); fstimlen = size(stim,1);
spikebins = allp.atrsbins{1};

%%
defmod = pars.defmod; %default module parameters
defmod.lnl = 200; %set NL lambda
defmod.lh = 10; %set PSC coefficient lambda
defmod.SDIM = sdim; %set stimulus dimensionality
foff = pars.flen + pars.hilen - 2; %need to leave off this many stim samples.  The minus 2 comes from the fact that the current frame is included in the internal kernel (otherwise it would be -1)
lltolf =1e-8; %LL tolerance for model fit

%adjust spike bin indices and get rid of spike bins before foff or after
%the length of the stim
corrected_spikebins = spikebins(spikebins>foff & spikebins<fstimlen)-foff;

%if you want to visualize the STC coefs being used
%plot1dfilterbank(allcells{tcellid}.effks,sdim)

init_glm      = createGLM0_jmm(allcells{tcellid}.effks,defmod,'sub'); %initialized with UNLOCALIZED STC FILTERS
norm_glm      = normalizeRFs(init_glm,stim); %normalize internal RFs to unit output variance (normalize k's and put in w)

%fit GLM with normalized initialized model, using stim and
%corrected_spikebins.  Use LL tolerance specified by lltolf
fit_glm       = fitNLHI_jmm(norm_glm,stim,corrected_spikebins,lltolf); 

%% T2 method
lltolr = 1e-5; 
lltolf =1e-8; 
tolfun=1e-3;
defmod = pars.defmod; 
defmod.lnl=200; 
defmod.lh=10;

foff = pars.flen + pars.hilen-2; %how much to leave off of stimulus?
flen = pars.flen;

%adjust spike bin indices and get rid of spike bins before foff
cspkbs2 = spikebins(spikebins>foff & spikebins<fstimlen)-foff;

%plot1dfilterbank(allcells{tcellid}.effks,sdim)

% get spiketrain & fit model
defmod.SDIM = sdim;
glmsub      = createGLM0(allcells{tcellid}.effks,defmod,'sub'); %initialized with UNLOCALIZED STC FILTERS
nmod        = normalizeRFs(glmsub,stim); %normalize internal RFs to unit output variance
fmod        = fitNLHI(nmod,stim,cspkbs2,lltolf); 

plotfo1d(fmod)

%% DButts method
