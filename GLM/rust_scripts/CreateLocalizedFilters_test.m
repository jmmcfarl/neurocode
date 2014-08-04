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
tuid = '43-21'; hassta=0;  npos=6; nneg=6; 


tid =  find(strcmp(tuid,uids)); 
eval(['load ',['~/Data/rust/stcbar/Data/',fnames{tid}]]);
% eval(['load ',['~/myrepo/workspace/DataRepository/rust/stcbar/Data/',fnames{tid}]])

psth = spikes_per_frm(:); 
rbins = (find(psth>0.5)); 
nsp = psth(rbins); 
spikebins =[]; 

for isp =1:length(rbins)
    spikebins= [spikebins; repmat(rbins(isp),nsp(isp),1)];
end
stcstart =1; stcend = size(stim,1); 

%% check rate stability
x = cumsum(psth)'; y = linspace(0,sum(psth),length(x)); plot(x-y); hline(0); axis tight;

%% get STC for a clean part of the sequence
stcstim   = stim(stcstart:stcend,:); 
stclen = size(stcstim,1); 
sdim = size(stim,2);

%readjust spike binning based on new start and end points
stcspks   = spikebins((spikebins >= stcstart) & (spikebins<stcend))-stcstart +1;
stcspksK  = stcspks(stcspks>flen)-flen+1; %spike binning adjusted for kernel
stcspksKH = stcspks(stcspks>foff)-foff+2; %spike binning adjusted for kernel and PSC terms

cd ~/Data/rust 
dataname = sprintf('stimspikedat-%s.mat',tuid)
save(dataname,'stcstim','stcspks','stcspksK','stcspksKH','pars','tuid'); 
% 
% x = cumsum(psth(stcstart:stcend)); y = linspace(0,sum(psth(stcstart:stcend)),length(x));
% figure; plot(x'-y); title(tuid); hline(0); axis tight;

[sta, stcs, fstim] = getSTCfilters(stcstim,stcspks,flen,8,8); 
rstcs = fliplr(stcs); %reversed STC kernels (suppressive first)

if (hassta) %if want to include STA filter
	effks = [sta',stcs(:,1:npos-1),rstcs(:,1:nneg)];
else
	effks = [stcs(:,1:npos),rstcs(:,1:nneg)];
end
plot1dfilterbanks({sta',stcs,effks},sdim);


%% get localized filters
%% find and save localized solutions for this cell
npos     = 2*sdim-1; %number of filter center-positions
rfcenters = linspace(1,sdim,npos); %equispaced grid of center positions
ninits = 10; %number of initialization repeats for finding localized filters
[lrfs1,lrfs2] = deal(zeros(sdim*flen,npos)); %???
nbanksExp = 2; fbanksExp  = cell(nbanksExp,1); %initialize a cell array of excitatory filter banks (nbanksExp per grid point)
nbanksSup = 2; fbanksSup  = cell(nbanksSup,1); %same for suppressive filters
tksExp = [stcs(:,1:6)]; %used excitatory STC subspace
tksSup = rstcs(:,1:6); %used suppressive STC subspace

%initialize filter banks to have dimensions (flen*SDIM x npos)
for ibank = 1:nbanksExp; 
    fbanksExp{ibank} = zeros(size(tksExp,1),npos);
end
for ibank = 1:nbanksSup
    fbanksSup{ibank} = zeros(size(tksSup,1),npos);
end

for ipos = 1:npos %for each grid position
	effksExp = tksExp; 	effksSup = tksSup; %set current exc and sup subspaces to be full STC spaces
	tpos  = rfcenters(ipos); disp(sprintf('--- position %i',tpos)); %set current grid position
	for nfilt = 1:nbanksExp %cycle through each filter bank at this position
        disp(sprintf(' -- filter #%i',nfilt)); 
		trf      = getCLOCrf_jmm(tpos,effksExp,sdim,ninits,ooptions); %get most localized filter at this position
		effksExp = normvecs(effksExp - trf*(effksExp'*trf)'); %project this filter out of the STC space
		fbanksExp{nfilt}(:,ipos) = trf; %store the filter
	end	
	for nfilt = 1:nbanksSup
        disp(sprintf(' -- filter #%i',nfilt));		
		trf      = getCLOCrf_jmm(tpos,effksSup,sdim,ninits,ooptions);
		effksSup = normvecs(effksSup - trf*(effksSup'*trf)');
		fbanksSup{nfilt}(:,ipos) = trf;
    end
end

plot1dfilterbanks(fbanksExp,sdim); plot1dfilterbanks(fbanksSup,sdim)
x = getPowerProfs(fbanksExp{2},sdim); figure; plot(x)
x = getPowerProfs(fbanksSup{1},sdim); figure; plot(x)


locmatname = sprintf('localizedRepertoire_jmm-%s.mat',tuid);
cd ~/Data/rust 
save(locmatname,'effks','fbanksExp','fbanksSup','tuid','pars'); 

%% select filters via sparse regression
tuid = '43-21'; nmodsLoc = 10;
cd ~/Data/rust 
locmatname = sprintf('localizedRepertoire_jmm-%s.mat',tuid);
eval(['load ' locmatname]);
addpath('~/James_scripts/glmnet_matlab/')

dataname   = sprintf('stimspikedat-%s.mat',tuid); 
stcdat = load(dataname);
[stclen ,sdim] = size(stcdat.stcstim);
stcstim = stcdat.stcstim;

locks = [fbanksExp{1} fbanksExp{2} fbanksSup{1} fbanksSup{2}]; %compile filter banks
nks    = size(locks,2);  %total number of filters
fouts = zeros(stclen-flen+1,nks); %initialize filter output matrix
%compute output of each filter
for ik = 1:nks
    fouts(:,ik) = kfilterInput(stcstim,locks(:,ik)); 
end
fouts2 = rectify([fouts -fouts]); %create positive and negative rectified versions of filter output

psth   = zeros(stclen-flen+1,1); %intialize vector of spike counts
x = tabulate(stcdat.stcspksK); %compute frequency table from the spike bin vector
psth(x(:,1)) = x(:,2); %psth is the spike count vector in each bin


fit1   = glmnet(fouts2,psth);

% cur_opts.alpha = 0.8;
% opts = glmnetSet(cur_opts);
% fit2   = glmnet(fouts2,psth,'gaussian',opts); 

% family = 'multinomial';
% upoints = 1:500;
% fit2   = glmnet(fouts2(upoints,:),psth(upoints)+1,family); 

glmnetPrint(fit1); 
figure; glmnetPlot(fit1,'lambda')

% %% GLM-fitting method for L1 regression of module weights
% addpath('~/James_scripts/GLM/')
% addpath('~/Dan_Matlab/')
% 
% %initialize parameters
% silent = 0;
% lamrange = [];
% lamrange2 = [];
% Pcon = [];
% Pmono = [];
% hold_const = [];
% NLtype = 0;
% 
% nmods = size(fouts2,2);
%     
% nlamvals = 10;
% lamvals = linspace(min(fit1.lambda),max(fit1.lambda),nlamvals);
% 
% %for L1 regularization of weights
% llist = [1];
% llist = [llist 1:nmods];
% 
% %initialize weight_mat
% weight_mat = nan(length(lamvals),nmods);
% plls = nan(length(lamvals),1);
% 
% used_trange = 1:50e3;
% used_mods = 1:50;
% fouts_used = fouts2(used_trange,used_mods);
% used_spikebins = stcdat.stcspksK;
% used_spikebins(used_spikebins > length(used_trange)) = [];
% for i = 1:length(lamvals)
%     
%     fprintf('Trying %d of %d lambda values\n',i,length(lamvals));
%     W0 = [rand(nmods,1); 0]; %random initialization
%     llist(1) = lamvals(i);
%     [fitp,grad] = GLMsolve_jmm(fouts_used, used_spikebins, W0, silent, lamrange, lamrange2, Pcon, Pmono, llist, hold_const, NLtype );
%     weight_mat(i,:) = fitp.k(1:end-1);
%     plls(i) = fitp.LL;
%     
% end
% 
% 
%% ideally: check via crossvalidation which filters to keep
tlambdas = fit1.lambda(2:end); %vector of lambda values tried 
zrange = [min(locks(:)),max(locks(:))]; %range of filter coefs
ctlambdas = tlambdas(1:20); %select from the 20 largest lambda values tried by glmnet
prevsigs = []; %store indices of sig filts in previous iteration
teffks={}; 
imod=0; 
for ilval = 1:length(ctlambdas); 
	tlambda = ctlambdas(ilval); %current lambda value
	cfs = glmnetCoef(fit1,tlambda); %get module weights for current lambda value
	cfs2 = [cfs(2:(nks+1)),-cfs((nks+2):end)]; %ignore the constant term and flip sign of negative modules
	
	sigfilts = find(sum(abs(cfs2),2)>0); %find filters where the positive and/or negative filters have non-zero coefs 
    nsfilts = length(sigfilts); %number of sig filts
	if not(length(sigfilts)==length(prevsigs)) %if there is a change in the number of sig filters from the previous lambda value
		imod = imod+1; 
        teffks{imod} = locks(:,sigfilts); %create filter bank for sig filters
		prevsigs = sigfilts; 
	elseif not(all(sigfilts==prevsigs)) %or if there is any change in which filters are sig 
		imod = imod+1; 
        teffks{imod} = locks(:,sigfilts);
		prevsigs = sigfilts; 
	end
end; 

nmodsSTC = size(effks,2); 

%% define number of models and add model with many filters from the stationary regularization part 
tlambda = exp(-6);
cfs = glmnetCoef(fit1,tlambda); 
cfs2 = [cfs(2:(nks+1)),-cfs((nks+2):end)]; 
sigfilts = find(sum(abs(cfs2),2)>0); 
bteffks  = locks(:,sigfilts); %extract significant filters at this lambda 
teffks = [teffks(1:(nmodsLoc-1)),bteffks]; 
teffks = teffks(1:nmodsLoc); 
plot1dfilterbanks(teffks,sdim); 


defmod = stcdat.pars.defmod; defmod.chmon = -1; defmod.SDIM=sdim

% comparable sequence of STC models
teffksSTC = {}; for imod =1:nmodsSTC; teffksSTC{imod}=effks(:,1:imod); end; 
plot1dfilterbanks(teffksSTC,sdim)

rmodsLoc = cellfun(@(x)createGLM0(x,defmod,'sub'),teffks,'UniformOutput',0);
  rmodsLoc{end}.lambdaW=10; 
  
  
rmodsSTC = cellfun(@(x)createGLM0(x,defmod,'sub'),teffksSTC,'UniformOutput',0);

modmatname = sprintf('rawSparseMods_jmm-%s.mat',tuid);
save(modmatname,'rmodsSTC','rmodsLoc','tuid')