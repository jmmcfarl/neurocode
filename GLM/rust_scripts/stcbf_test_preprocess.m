clear all;
addpath('~/Timm/rust/SparseFilterSelection/')
addpath('~/Timm/MatlabRepository/')
addpath('~/James_scripts/GLM')

pars = load('~/Data/rust/infos/StandardParametersRust.mat');
flen = pars.flen;
foff = flen + pars.hilen;
ooptions = optimset('MaxFunEvals',100000,'MaxIter',10000);

datdir = '~/Data/rust/stcbar/Data/';
cfiles = dir([datdir,'*stc.mat']); ncells = length(cfiles);
fnames = arrayfun(@(x)x.name,cfiles,'UniformOutput',0);
uids   = cellfun(@(x)[x(2:3),'-',x(6:7)],fnames,'UniformOutput',0);

% tuid = '52-18'; hassta=1;  npos=1; nneg=1;
tuid = '44-29'; hassta=0;  npos=6; nneg=6;
% tuid = '43-03'; hassta=0;  npos=6; nneg=6;
% tuid = '43-21'; hassta=0;  npos=6; nneg=6; %was 6,6
% tuid = '33-27';  %% complicated simple 
% tuid = '43-03'; %%??
% tuid = '52-18'; %% simple scell
% tuid = '44-29'; %% the complex cell
% tuid = '33-44';  %% ??

%% load data
tid =  find(strcmp(tuid,uids));
eval(['load ',['~/Data/rust/stcbar/Data/',fnames{tid}]]);

psth = spikes_per_frm(:);
rbins = (find(psth>0.5));
nsp = psth(rbins);
spikebins =[];

for isp =1:length(rbins)
    spikebins= [spikebins; repmat(rbins(isp),nsp(isp),1)];
end
stcstart =1; stcend = size(stim,1);

stcstim   = stim(stcstart:stcend,:);
stclen = size(stcstim,1);
sdim = size(stim,2);

%readjust spike binning based on new start and end points
stcspks   = spikebins((spikebins >= stcstart) & (spikebins<stcend))-stcstart +1;

%% get STC for a clean part of the sequence
[sta, stcs, fstim, evs] = getSTCfilters(stcstim,stcspks,flen,20,20);
rstcs = fliplr(stcs); %reversed STC kernels (suppressive first)

eig_diffs = diff(evs);
% npos = find(eig_diffs < 0.005,1,'first');
% nneg = find(eig_diffs(end:-1:1) < 0.005,1,'first');
npos = 6;
nneg = 6;

% hassta = 0;
% if (hassta) %if want to include STA filter
%     effks = [sta',stcs(:,1:npos),rstcs(:,1:nneg)];
% else
%     effks = [stcs(:,1:npos),rstcs(:,1:nneg)];
% end
% stcs1 = [stcs(:,1:npos+5) rstcs(:,1:nneg+5)];
% plot1dfilterbanks({sta',stcs1,effks},sdim);

%% save data
cd ~/Data/rust
dataname = sprintf('stcbf_data-%s.mat',tuid)
save(dataname,'npos','nneg','stcs','sta','evs','pars','tuid','stim','spikebins');
