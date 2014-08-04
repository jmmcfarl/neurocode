%% Load Data

clear all;
addpath('~/James_scripts/GLM')
addpath('~/James_scripts/GLM/t1')
addpath('~/Timm/rust/SparseFilterSelection/')
addpath('~/Timm/MatlabRepository/')
addpath('~/Timm/tim1/functions/')

cd ~/Data/blanche/rec_75/matlabdata/
load stdparsRec75.mat

% load spks7576-ns.mat; load PCScores-ns.mat; stype='ns'
% load spks7576-pn.mat; load PCScores-pn.mat; stype='pn'
% load spks7576-ps.mat; load PCScores-ps.mat; stype='ps'

%load z-scored PC data
cd ~/Data/blanche/matlabdata/
load spks7576-all.mat;
load PCScores_z-all.mat;
stype='all';

stcmods = 4;
ncomps = 6144;
% scorevars(1001:end) = [];
% coefs(:,1001:end) = [];
% scores(:,1001:end) = [];
% datacov = coefs*diag(scorevars)*coefs';
X = scores*coefs';
datacov = cov(X);
clear scores
%% Compute whitened data, and STA/STCs
cd ~/James_scripts/GLM/t1/multirec_stc_compare/
for t = 1:27
%old approach
    tcell = t;
    pids =1:1024;
    tsbs      = 1+floor(aselspks{tcell}/dt);
    
    compids   = 1:ncomps;
    
%     WX        = scores(:,compids)*diag(1./sqrt(scorevars(compids)));
%     spikebins = tsbs(tsbs>flen & tsbs<(size(WX,1)+flen-1))-flen+1 ;
%     WS        = WX(spikebins,:);
%     rsta      = mean(WS) - mean(WX);
%     
%     stvcv = cov(WS);  utvcv = cov(WX);
%     [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
%     
%     npos=6; nneg=6;
%     stcs  = evecs(:,[1:nneg,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
%     
%     rstcs = fliplr(stcs); %reversed STC kernels (suppressive first)    
%     
%     STCbvs = [stcs(:,1:npos) rstcs(:,1:nneg)]; %use only expansive subspace
% 
%     pix_conv_mat = diag(sqrt(scorevars(compids)))*coefs(:,compids)';
%     
%     kimages = [rsta',STCbvs]'*pix_conv_mat;
%     f1 = figure('Name',stype,'visible','off'); plotfilterbank(kimages',32,pids)
%     set(f1,'PaperUnits','centimeters');
%     set(f1, 'PaperSize', [25 50]);
%     set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))])
%     fname = sprintf('Cell%d_OldSTC',t);
%     print(f1,'-dpng',fname);close all
%     f1 = figure('visible','off');
%     set(f1,'PaperUnits','centimeters');
%     set(f1, 'PaperSize', [10 10]);
%     set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))])
%     plot(diag(evals),'.'), axis tight
%     fname = sprintf('Cell%d_OldSTC_evals',t);
%     print('-dpng',fname);close

% %new approach    
%     WX        = scores(:,compids);
%     spikebins = tsbs(tsbs>flen & tsbs<(size(WX,1)+flen-1))-flen+1 ;
%     WS        = WX(spikebins,:);
%     rsta      = mean(WS) - mean(WX);
%     
%     stvcv = cov(WS);  utvcv = cov(WX);
%     [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
%     
%     npos=6; nneg=6;
%     stcs  = evecs(:,[1:nneg,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
%     
%     rstcs = fliplr(stcs); %reversed STC kernels (suppressive first)    
%     
%     STCbvs = [stcs(:,1:npos) rstcs(:,1:nneg)]; %use only expansive subspace
%     kimagesc = [rsta' STCbvs]'*diag(1./sqrt(scorevars(compids)))*coefs(:,compids)';
%     f1 = figure('Name',stype,'visible','off'); plotfilterbank(kimagesc',32,pids)
%     set(f1,'PaperUnits','centimeters');
%     set(f1, 'PaperSize', [25 50]);
%     set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))])
%     fname = sprintf('Cell%d_New2STC',t);
%     print(f1,'-dpng',fname);close all
%     f1 = figure('visible','off');
%     set(f1,'PaperUnits','centimeters');
%     set(f1, 'PaperSize', [10 10]);
%     set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))])
%     fname = sprintf('Cell%d_New2STC_evals',t);
%     plot(diag(evals),'.'),axis tight
%     print('-dpng',fname);close

%new approach    
    spikebins = tsbs(tsbs>flen & tsbs<(size(X,1)+flen-1))-flen+1 ;
    S        = X(spikebins,:);
    rsta      = mean(S);
    
    stvcv = cov(S); 
    [evecs,evals] = eig(stvcv-datacov); evs   = diag(evals);
    
    npos=6; nneg=6;
    stcs  = evecs(:,[1:nneg,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
    
    rstcs = fliplr(stcs); %reversed STC kernels (suppressive first)    
    
    STCbvs = [stcs(:,1:npos) rstcs(:,1:nneg)]; %use only expansive subspace
    kimagesc = [rsta' STCbvs];
    f1 = figure('Name',stype,'visible','off'); plotfilterbank(kimagesc,32,pids)
    set(f1,'PaperUnits','centimeters');
    set(f1, 'PaperSize', [25 50]);
    set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))])
    fname = sprintf('Cell%d_rawSTC',t);
    print(f1,'-dpng',fname);close all
    f1 = figure('visible','off');
    set(f1,'PaperUnits','centimeters');
    set(f1, 'PaperSize', [10 10]);
    set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))])
    fname = sprintf('Cell%d_raw_evals',t);
    plot(diag(evals),'.'),axis tight
    print('-dpng',fname);close
    
    kimagesc = [rsta' STCbvs]'*coefs;
    kimagesc = kimagesc(:,1:600)*diag(1./sqrt(scorevars(1:600)))*coefs(:,1:600)';
    f1 = figure('Name',stype,'visible','off'); plotfilterbank(kimagesc',32,pids)
    set(f1,'PaperUnits','centimeters');
    set(f1, 'PaperSize', [25 50]);
    set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))])
    fname = sprintf('Cell%d_rawSTC_white',t);
    print(f1,'-dpng',fname);close all

        kimagesc = [rsta' STCbvs]'*coefs;
    kimagesc = kimagesc(:,1:600)*coefs(:,1:600)';
    f1 = figure('Name',stype,'visible','off'); plotfilterbank(kimagesc',32,pids)
    set(f1,'PaperUnits','centimeters');
    set(f1, 'PaperSize', [25 50]);
    set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))])
    fname = sprintf('Cell%d_rawSTC_uncor',t);
    print(f1,'-dpng',fname);close all

end

