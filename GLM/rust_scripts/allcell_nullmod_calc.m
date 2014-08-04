clear all;
close all;

addpath('~/Timm/rust/SparseFilterSelection/')
addpath('~/Timm/MatlabRepository/')
addpath('~/James_scripts/GLM/')

pars = load('~/Data/rust/infos/StandardParametersRust.mat');
defmod = pars.defmod;
defmod.h(1:end-1) = [];%restrict PSC term to delta function
flen = pars.flen;
hilen = length(defmod.h);
foff = flen + pars.hilen;

datdir = '~/Data/rust/stcbar/Data/';
cfiles = dir([datdir,'*stc.mat']); ncells = length(cfiles);
fnames = arrayfun(@(x)x.name,cfiles,'UniformOutput',0);
uids   = cellfun(@(x)[x(2:3),'-',x(6:7)],fnames,'UniformOutput',0);
cd ~/Data/blanche/rec_75/matlabdata/
load stdparsRec75.mat

%load precomputed STC filt params
cd /Users/James/James_scripts/stc_sparse_test_figs/stac_allcells
load ./used_stcims.mat

cd ~/James_scripts/stc_sparse_test_figs/
load rust_allcells_modanal

ncells = length(uids)
for nn = 1:length(all_xv_params)
    
    fprintf('ANALYZING CELL %d OF %d\n\n',nn,ncells);
    
    %% load data
    eval(['load ',['~/Data/rust/stcbar/Data/',fnames{nn}]]);
    
    psth = spikes_per_frm(:);
    rbins = (find(psth>0.5));
    nsp = psth(rbins);
    spikebins =[];
    uvals = unique(psth);
    for i = 1:length(uvals)
        cur_set = find(psth==uvals(i));
        spikebins = [spikebins; repmat(cur_set(:),uvals(i),1)];
    end
    spikebins = sort(spikebins);
    
    %% create XV data
    [stimlen,sdim] = size(stim);
    nparts = 18;
    partlen = floor(stimlen/nparts);
    nfold = 3;
    nxvparts = nparts/nfold;
    
    %boundaries of parts
    pbounds  = [(0:nparts-1)*partlen+1;(1:nparts)*partlen]';
    
    for i = 1:nfold
        xv_inds{i} = [];
        xv_spkbns{i} = [];
        tr_inds{i} = [];
        tr_spkbns{i} = [];
        
%         cur_perm = randperm(nparts);
%         cur_xv_parts = sort(cur_perm(1:nxvparts));
%         cur_tr_parts = setdiff(1:nparts,cur_xv_parts);
cur_tr_parts = all_xv_params(nn).cur_tr_parts;
cur_xv_parts = all_xv_params(nn).cur_xv_parts;


        xv_spks = [];
        xv_new_inds = nan(stimlen,1);
        for j = 1:length(cur_xv_parts)
            cur_start = pbounds(cur_xv_parts(j),1);
            cur_stop = pbounds(cur_xv_parts(j),2);
            xv_spks = [xv_spks; spikebins(spikebins >= cur_start & spikebins < cur_stop)];
            
            cur_inds = (cur_start:cur_stop) - cur_start + length(xv_inds{i}) + 1;
            xv_new_inds(cur_start:cur_stop) = cur_inds;
            
            xv_inds{i} = [xv_inds{i} cur_start:cur_stop];
        end
        xv_spkbns{i} = xv_new_inds(xv_spks);
        
        tr_spks = [];
        tr_new_inds = nan(stimlen,1);
        for j = 1:length(cur_tr_parts)
            cur_start = pbounds(cur_tr_parts(j),1);
            cur_stop = pbounds(cur_tr_parts(j),2);
            tr_spks = [tr_spks; spikebins(spikebins >= cur_start & spikebins < cur_stop)];
            
            cur_inds = (cur_start:cur_stop) - cur_start + length(tr_inds{i}) + 1;
            tr_new_inds(cur_start:cur_stop) = cur_inds;
            
            tr_inds{i} = [tr_inds{i} cur_start:cur_stop];
        end
        tr_spkbns{i} = tr_new_inds(tr_spks);
    end
    [NT,sdim] = size(stim);
    
    xv = 1; %use only one XV set
    
    cur_xv_stim = stim(xv_inds{xv},:);
    %     cur_xv_stimemb = makeStimRows(cur_xv_stim,flen);
    cur_xv_spkbns = xv_spkbns{xv};
    
    cur_tr_stim = stim(tr_inds{xv},:);
    %     cur_tr_stimemb = makeStimRows(cur_tr_stim,flen);
    cur_tr_spkbns = tr_spkbns{xv};
    
    %%
    avg_rate(nn) = length(cur_tr_spkbns)/length(tr_inds{xv});
    
    null_pred = ones(size(tr_inds{xv}))*avg_rate(nn);
    
    null_xvLL(nn) = -(sum(log(null_pred(cur_xv_spkbns))) - sum(null_pred))/length(cur_xv_spkbns);
    
end
%%
sta_xv = arrayfun(@(x) x.xvLL,all_sta_mod);
stc_xv = arrayfun(@(x) x.xvLL,all_stc_mod);
rot_xv = arrayfun(@(x) x.xvLL,all_rotmods);
fin_xv = arrayfun(@(x) x.xvLL,all_fin_mod);

sta_imp = (null_xvLL - sta_xv)/log(2);
stc_imp = (null_xvLL - stc_xv)/log(2);
rot_imp = (null_xvLL - rot_xv)/log(2);
fin_imp = (null_xvLL - fin_xv)/log(2);

stc_impsta = (sta_xv - stc_xv)/log(2);
rot_impsta = (sta_xv - rot_xv)/log(2);
fin_impsta = (sta_xv - fin_xv)/log(2);

rot_impstc = (stc_xv - rot_xv)/log(2);
fin_impstc = (stc_xv - fin_xv)/log(2);

