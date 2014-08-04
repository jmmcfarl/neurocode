clear all; pars = load('~/Data/rust/infos/StandardParametersRust.mat'); 
datdir = '~/Data/rust/stcbar/DATA/'; 

datdir = '~/Data/rust/stcbar/Data/';
cfiles = dir([datdir,'*stc.mat']); ncells = length(cfiles);
fnames = arrayfun(@(x)x.name,cfiles,'UniformOutput',0);
uids   = cellfun(@(x)[x(2:3),'-',x(6:7)],fnames,'UniformOutput',0);

cd ~/Data/rust/
load RustSTCsallcells

ncells = length(allcells);
nfilts = 8;

 cd ~/James_scripts/stc_sparse_test_figs/stac_allcells/
for ii = 1:ncells
    sdim = allcells{ii}.sdim;
    flen = length(allcells{ii}.sta)/sdim;
    figure
    subplot(3,nfilts,1)
    imagesc(reshape(allcells{ii}.sta,flen,sdim));
    for i = 1:nfilts
        subplot(3,nfilts,nfilts+i)
        imagesc(reshape(allcells{ii}.stcs(:,i),flen,sdim));
    end
    for i = 1:nfilts
        subplot(3,nfilts,2*nfilts+i)
        imagesc(reshape(allcells{ii}.stcs(:,2*nfilts+1-i),flen,sdim));
        colormap(gray)
    end
    
    print(uids{ii},'-dpng'); close 
    
end

%%
n_used_stcdims = ...
    [7 8;
    6 6;
    4 8;
    4 8;
    2 7;
    8 8;
    2 6;
    7 7;
    3 8;
    6 7;
    8 8;
    3 5;
    1 3;
    8 8;
    6 8;
    6 0;
    8 8;
    6 5;
    6 7;
    6 6;
    6 8;
    7 8;
    3 8;
    8 7;
    7 5;
    2 6;
    8 8;
    8 8;
    3 7;
    5 8;
    5 0;
    1 1;
    5 4;
    5 0;
    1 5;
    4 6;
    4 0;
    7 8;
    8 0;
    4 3;
    7 6;
    2 5;
    3 6;
    6 2;
    3 3;
    1 3;
    6 6;
    4 5;
    4 4;
    5 2;
    3 2;
    4 5;
    4 8;
    3 2;
    3 3;
    2 4;
    3 3;
    8 7;
    3 4];
    
    
    