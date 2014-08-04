clear all; flen =6;

%% natural stimli from recording 75
cd ~/Data/blanche/rec_75/matlabdata/
load dstimpsRec75.mat;
stf75=load('stimfiles75.mat');
cd ~/Data/blanche/rec_76/matlabdata/
load dstimpsRec76.mat;
stf76=load('stimfileNames76.mat');

allstims = [dstimps75;dstimps76];
cnames   = [stf75.stimfiles,stf76.stimfiles]; cellfun(@disp,cnames)
nconds   = length(cnames);

keys = {'af','wn','pn','pt','ps','ns'}
rids = cellfun(@(tkey)...
    find(strcmp(cellfun(@(x)x(1:2),cnames,'UniformOutput',0),tkey)),...
    keys,'UniformOutput',0);

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
used_conds = tns(~isnan(stim_mean(tns)) & stim_contrast(tns) ~= 15);
% used_conds = [tpn(2:end) tpt(2:end) tps(2:end)];
% used_conds = [tps(2:end)];

%define subpatch of image
[XX,YY] = meshgrid(1:32,1:32);
Xr = XX(:); Yr = YY(:);
used_pixs = find(Xr <= 16 & Yr <= 16); 

selstim   = [];
cur_means = [];
cur_vars = [];
for icond=1:length(used_conds)
    tcond = used_conds(icond)
    cur_stim = allstims{tcond};
    cur_stim = cur_stim(:,used_pixs);
    %         cur_stim = dsamplebstim(cur_stim,2);
    %         cur_stim = zscore(cur_stim); %z-score normalization
    selstim =[selstim;cur_stim];
    cur_means = [cur_means; mean(cur_stim(:))];
    cur_vars = [cur_vars; var(cur_stim(:))];
end

addpath('~/James_scripts/')
NT   = size(selstim,1); SDIM  = size(selstim,2); NeK   = flen*SDIM;
X = makeStimRows(selstim,flen,1);
clear allstims


%% the spikes
cd ~/Data/blanche/rec_75/matlabdata/
load spksegsRec75.mat;
load stdparsRec75;
cd ~/Data/blanche/rec_76/matlabdata/
load spksegsRec76.mat;
allspks  = [spksegs75,spksegs76]';  ncells = size(allspks,2);
stimlen = dt*6000
nspks = cellfun(@length,allspks);

% boxplot(nspks','Labels',cnames,'labelorientation','inline')

aselspks = cell(1,ncells);
tconds = used_conds;
for icell = 1:ncells
    selspks = [];
    for icond=1:length(tconds)
        tcond = tconds(icond)
        selspks=[selspks;((icond-1)*stimlen)+allspks{tcond,icell}];
    end
    aselspks{icell} = selspks;
end

% ecdf(aselspks{1,1}); hold on; for icell =2:10; ecdf(aselspks{1,icell}); end;
% [f,x] = ecdf(aselspks{2,1}); plot(x,f,'r'); for icell =2:10; [f,x] = ecdf(aselspks{2,icell}); plot(x,f,'r'); end; hold off;

% tfilename = 'spks7576-NS_d2.mat'; save(tfilename,'-v7.3','aselspks');

%% save data
cd ~/Data/blanche/matlabdata/
tfilename = 'NS_Xmat_ds2';
save(tfilename,'-v7.3','X','used_conds','aselspks');

%% WHITENING STIM
disp(' -- whitening'); drawnow;
[coefs,scores] = princomp(X);
scorevars = var(scores);

tfilename = 'NS_PCScores_ds2';
save(tfilename,'-v7.3','coefs','scores','scorevars','flen','SDIM');

