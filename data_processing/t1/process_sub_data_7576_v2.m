clear all; flen =8;

%% natural stimli from recording 75
cd ~/Data/blanche/rec_75
load ./sub_dstimpsRec75
stf75=load('stimfiles75.mat');
cd ~/Data/blanche/rec_76/matlabdata
load ./sub_dstimpsRec76.mat;
stf76=load('stimfileNames76.mat');

allstims = [dstimps75;dstimps76];
cnames   = [stf75.stimfiles,stf76.stimfiles]; cellfun(@disp,cnames)
nconds   = length(cnames);

keys = {'af','wn','pn','pt','ps','ns'}
rids = cellfun(@(tkey)...
    find(strcmp(cellfun(@(x)x(1:2),cnames,'UniformOutput',0),tkey)),...
    keys,'UniformOutput',0);

dt = 19.9920031987*1e-3; %true sweep time.
sweeps_per_stim = 6000;

taf        = rids{1}; %cellfun(@disp,cnames(tpn))
twn        = rids{2};
tpn        = rids{3}; %cellfun(@disp,cnames(tpn))
tpt        = rids{4};
tps        = rids{5}; %cellfun(@disp,cnames(tps))
tns        = rids{6}; %cellfun(@disp,cnames(tns))

is_nat_contrast = ~cellfun(@(x) (x(end)=='5'),cnames);
is_nat_mean = ~cellfun(@(x) (x(18)=='1'),cnames);
stim_contrast = nan(size(is_nat_contrast));
stim_contrast(~is_nat_contrast) = cellfun(@(x) str2num(x(end-2:end)),cnames(~is_nat_contrast));
stim_mean = nan(size(is_nat_mean));
stim_mean(~is_nat_mean) = cellfun(@(x) str2num(x(18:20)),cnames(~is_nat_mean));

% DIFFERENT MEAN GROUPS (POOLED ACROSS NS, PN, AND PS)
% used_conds(1,:) = logical(stim_mean == 125 & ismember(1:nconds,[tns tpn tps]));
% used_conds(2,:) = logical(isnan(stim_mean) & ismember(1:nconds,[tns tpn tps]));
% used_conds(2,:) = logical(isnan(stim_mean) & ~isnan(stim_contrast) & ismember(1:nconds,[tns tpn tps]));
% used_conds = tns(~isnan(stim_mean(tns)) & stim_contrast(tns) ~= 15);
% used_conds = [tpn(2:end) tpt(2:end) tps(2:end)];
% used_conds = [tps(2:end)];
used_conds = [tpn(2:end) tpt(2:end) tps(2:end) tns];
used_conds(isnan(stim_mean(used_conds))) = [];

%%
addpath(genpath('~/James_scripts/'))
selstim   = [];
cur_means = [];
cur_vars = [];
X = [];
for icond=1:length(used_conds)
    tcond = used_conds(icond);
    fprintf('Condition %d of %d\n',icond,length(used_conds));
    cur_stim = allstims{tcond}';
    cur_stim = bsxfun(@minus,cur_stim,mean(cur_stim));
    cur_means = [cur_means; mean(cur_stim(:))];
    cur_vars = [cur_vars; var(cur_stim(:))];
    X = [X; makeStimRows(cur_stim,flen,0)];
end

NT   = size(X,1); SDIM  = size(cur_stim,2); NeK   = flen*SDIM;

%% the spikes
cd ~/Data/blanche/rec_75/matlabdata/
load spksegsRec75.mat;
% load stdparsRec75;
cd ~/Data/blanche/rec_76/matlabdata/
load spksegsRec76.mat;

allspks  = [spksegs75,spksegs76]';  ncells = size(allspks,2);
stimlen = dt*6000;
nspks = cellfun(@length,allspks);

% boxplot(nspks','Labels',cnames,'labelorientation','inline')

aselspks = cell(1,ncells);
tconds = used_conds;
for icell = 1:ncells
    selspks = [];
    for icond=1:length(tconds)
        tcond = tconds(icond);
        selspks=[selspks;((icond-1)*stimlen)+allspks{tcond,icell}];
    end
    aselspks{icell} = selspks;
end

%% save data
cd ~/Data/blanche/matlabdata/
tfilename = 'sub_Xmat2';
save(tfilename,'-v7.3','X','used_conds','aselspks','dt','flen');

% %% WHITENING STIM
% disp(' -- whitening'); drawnow;
% [coefs,scores] = princomp(X);
% scorevars = var(scores);
% 
% tfilename = 'NS_PCScores_ds2';
% save(tfilename,'-v7.3','coefs','scores','scorevars','flen','SDIM');

