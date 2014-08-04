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
used_conds(1,:) = logical(stim_mean == 125 & ismember(1:nconds,[tns tpn tps]));
used_conds(2,:) = logical(isnan(stim_mean) & ismember(1:nconds,[tns tpn tps]));
% used_conds(2,:) = logical(isnan(stim_mean) & ~isnan(stim_contrast) & ismember(1:nconds,[tns tpn tps]));

n_cond_groups = size(used_conds,1);

selstim   = cell(n_cond_groups,1);
cur_means = cell(n_cond_groups,1);
cur_vars = cell(n_cond_groups,1);
for gcond = 1:n_cond_groups
    tconds = find(used_conds(gcond,:));
    for icond=1:length(tconds)
        tcond = tconds(icond)
        cur_stim = allstims{tcond};
        cur_stim = zscore(cur_stim); %z-score normalization
        selstim{gcond} =[selstim{gcond};cur_stim];
        cur_means{gcond} = [cur_means{gcond}; mean(cur_stim(:))];
        cur_vars{gcond} = [cur_vars{gcond}; var(cur_stim(:))];
    end
end

addpath('~/James_scripts/')
for gcond = 1:n_cond_groups;
    NT   = size(selstim{gcond},1); SDIM  = size(selstim{gcond},2); NeK   = flen*SDIM;
    X{gcond} = makeStimRows(selstim{gcond},flen,1);
end
clear allstims


%% the spikes
cd ~/Data/blanche/rec_75/matlabdata/
load spksegsRec75.mat;
load stdparsRec75;
cd ~/Data/blanche/rec_76/matlabdata/
load spksegsRec76.mat;
allspks  = [spksegs75,spksegs76]';  ncells = size(allspks,2)
stimlen = dt*6000
nspks = cellfun(@length,allspks);

% boxplot(nspks','Labels',cnames,'labelorientation','inline')

aselspks = cell(n_cond_groups,ncells);
for gcond = 1:n_cond_groups
    tconds = find(used_conds(gcond,:));
    for icell = 1:ncells
        selspks = [];
        for icond=1:length(tconds)
            tcond = tconds(icond)
            selspks=[selspks;((icond-1)*stimlen)+allspks{tcond,icell}];
        end
        aselspks{gcond,icell} = selspks;
    end
end

% ecdf(aselspks{1,1}); hold on; for icell =2:10; ecdf(aselspks{1,icell}); end;
% [f,x] = ecdf(aselspks{2,1}); plot(x,f,'r'); for icell =2:10; [f,x] = ecdf(aselspks{2,icell}); plot(x,f,'r'); end; hold off;

% tfilename = sprintf('spks7576-%s.mat',ctype); save(tfilename,'-v7.3','aselspks');

%% save data
cd ~/Data/blanche/matlabdata/
tfilename = 'Stimmean_Xmat_z';
save(tfilename,'-v7.3','X','used_conds','aselspks');

%% WHITENING STIM
for gcond = 1:2;
disp(' -- whitening'); drawnow;
[coefs{gcond},scores{gcond}] = princomp(X{gcond});
scorevars{gcond} = var(scores{gcond});

tfilename = 'Stimmean_PCScores_z';
save(tfilename,'-v7.3','coefs','scores','scorevars','flen','SDIM');

end