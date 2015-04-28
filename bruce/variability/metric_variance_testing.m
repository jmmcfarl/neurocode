clear all
close all

nTrials = 200;
nStims = 200;
xx = linspace(0,2*pi,nStims);

alpha = 2;
beta = 0.5;
trialGains = gamrnd(alpha,beta,nTrials,nStims);

rates1 = 3*log(1+3*exp(5+bsxfun(@times,sin(xx),trialGains)));
rates2 = log(1+exp(bsxfun(@times,cos(xx),trialGains)));
% rates1 = bsxfun(@times,(sin(xx)),trialGains);
% rates2 = bsxfun(@times,(cos(xx)),trialGains);
% rates1 = (rates1 +5)*2; rates1(rates1 < 0) = 0;
% rates2 = (rates2 +5)*2; rates2(rates2 < 0) = 0;

yy1 = poissrnd(rates1); yy2 = poissrnd(rates2);

psth1 = mean(yy1); psth2 = mean(yy2);


%%
% tot_xcov = mean((yy1(:)-mean(yy1(:))).*(yy2(:) - mean(yy2(:))));
% psth_xcov = mean((psth1-mean(psth1)).*(psth2-mean(psth2)));

%%
yy1_bar = yy1 - mean(yy1(:));
yy2_bar = yy2 - mean(yy2(:));
% 
% [II,JJ] = meshgrid(1:nTrials);
% upairs = II ~= JJ;
% 
% allX = [];
% allD = [];
% for ii = 1:nStims
%    curD = squareform(pdist(trialGains(:,ii)));
%    yyprod = bsxfun(@times,yy1_bar(:,ii),yy2_bar(:,ii)');
%    allD = cat(1,allD,curD(upairs));
%    allX = cat(1,allX,yyprod(upairs));
% end

% yy1_bar = bsxfun(@minus,yy1,psth1);
% yy2_bar = bsxfun(@minus,yy2,psth2);
[II,JJ] = meshgrid(1:nTrials);
upairs = II > JJ;

allX = [];
allD = [];
for ii = 1:nStims
%    curD = squareform(pdist(yy2_bar(:,ii)));
   curD = squareform(pdist(trialGains(:,ii)));
   yyprod = bsxfun(@times,yy1_bar(:,ii),yy1_bar(:,ii)');
   allD = cat(1,allD,curD(upairs));
   allX = cat(1,allX,yyprod(upairs));
end

%%
n_EP_bins = 100;
EP_bin_edges = prctile(allD,linspace(0,100,n_EP_bins+1));
EP_bin_centers = 0.5*EP_bin_edges(1:end-1) + 0.5*EP_bin_edges(2:end);
EP_xcov = nan(n_EP_bins,1);
for ee = 1:n_EP_bins
    curset = find(allD >= EP_bin_edges(ee) & allD < EP_bin_edges(ee+1));
   EP_xcov(ee) = mean(allX(curset));
end
% poss_Ds = 0:50;
% EP_xcov = nan(length(poss_Ds),1);
% for ee = 1:length(poss_Ds)
%     curset = find(allD == poss_Ds(ee));
%    EP_xcov(ee) = mean(allX(curset));
% end


pair_xcov = mean(allX(allD == 0));
pair_psth_cov = mean(allX);
tot_var = var(yy1(:));
noise_var = tot_var - pair_xcov;


%%
clear all
close all

nTrials = 400;
nStims = 50;
xx = linspace(0,2*pi,nStims);

alpha = 4;
beta = 0.25;
trialGains = gamrnd(alpha,beta,nTrials,nStims);
% trialGains = ones(nTrials,nStims);

rates = 3*log(1+3*exp(5+bsxfun(@times,cos(xx)+0.3,trialGains)));
% rates = bsxfun(@plus,cos(xx),trialGains);

%%
yy = poissrnd(rates); 
% yy = rates;
psth = mean(yy);


%%
% tot_xcov = mean((yy1(:)-mean(yy1(:))).*(yy2(:) - mean(yy2(:))));
% psth_xcov = mean((psth1-mean(psth1)).*(psth2-mean(psth2)));

%%
yy_bar = yy - mean(yy(:));
% yy_bar = bsxfun(@minus,yy,psth);

[II,JJ] = meshgrid(1:nTrials);
upairs = (II ~= JJ);

epsilon = 0.02;
% n_EP_bins = 100;
% EP_bin_edges = prctile(allD,linspace(0,100,n_EP_bins+1));
% EP_bin_centers = 0.5*EP_bin_edges(1:end-1) + 0.5*EP_bin_edges(2:end);

% randpicks = 100;

%%
allD = [];
for ii = 1:nStims
%     curpicks = randi(nTrials,randpicks,1);
    
    curD = squareform(pdist(trialGains(:,ii)));
    curD(~upairs) = nan;
    allD = cat(1,allD,curD(upairs));
    
end

n_splines = 4;
% spline_DS = 0:0.5:4;
spline_DS = prctile(allD,100/n_splines:100/n_splines:(100-100/n_splines));
knot_pts = [0 0 0 0 spline_DS max(allD) max(allD) max(allD) max(allD)];
sp = fastBSpline(knot_pts,ones(length(knot_pts) - 4,1));

ovWeights = sp.getBasis(allD);
ovWeights = mean(ovWeights);

n_splines = length(ovWeights);
%%

allX = [];
allS = [];
allD = [];
allM = [];
allW = [];
for ii = 1:nStims
%     curpicks = randi(nTrials,randpicks,1);
    
    curD = squareform(pdist(trialGains(:,ii)));
    curD(~upairs) = nan;
    
    yyprod = bsxfun(@times,yy_bar(:,ii)',yy_bar(:,ii));
    yyvals = [yy_bar(II,ii) yy_bar(JJ,ii)];

%     weights = nansum(curD < epsilon)/(nTrials-1);
%     WW = repmat(weights,nTrials,1);
    
    ddvals = [trialGains(II,ii) trialGains(JJ,ii)];
    
    WW = nan(nTrials,nTrials,n_splines);
%     ww = nan(nTrials,nTrials,n_splines);
    for jj = 1:nTrials
        proj = sp.getBasis(curD(jj,:));
        corfac = ovWeights./mean(proj);
        weights = bsxfun(@rdivide,proj,corfac);
        WW(jj,:,:) = weights;
%         ww(jj,:,:) = proj;
    end
    WW = reshape(WW,[],n_splines);
%     netweights = sum(weights);
%     WW = repmat(netweights,nTrials,1);
%     allW = cat(1,allW(
    
    allD = cat(1,allD,curD(upairs));
    allX = cat(1,allX,yyprod(upairs));
    allS = cat(1,allS,yyvals(upairs,:));
    allM = cat(1,allM,ddvals(upairs,:));
    allW = cat(1,allW,WW(upairs,:));
end

uu = find(~isnan(allX) & max(allW,[],2) > 0);
rr = lsqlin(allW(uu,:),allX(uu))
% w = allW(uu,:)\allX(uu)
%%
allX = [];
allS = [];
allD = [];
allM = [];
allW = [];
for ii = 1:nStims
%     curpicks = randi(nTrials,randpicks,1);
    
    curD = squareform(pdist(trialGains(:,ii)));
    curD(~upairs) = nan;
    
    yyprod = bsxfun(@times,yy_bar(:,ii)',yy_bar(:,ii));
    yyvals = [yy_bar(II,ii) yy_bar(JJ,ii)];

%     weights = nansum(curD < epsilon)/(nTrials-1);
%     WW = repmat(weights,nTrials,1);
    
    ddvals = [trialGains(II,ii) trialGains(JJ,ii)];
    
    weights = sp.getBasis(curD);
    weights = reshape(weights(:,1),nTrials,nTrials);
    netweights = sum(weights);
    WW = repmat(netweights,nTrials,1);
%     allW = cat(1,allW(
    
    allD = cat(1,allD,curD(upairs));
    allX = cat(1,allX,yyprod(upairs));
    allS = cat(1,allS,yyvals(upairs,:));
    allM = cat(1,allM,ddvals(upairs,:));
    allW = cat(1,allW,WW(upairs));
end



%%
% allB = sp.getBasis(allD);
% rr = (allX./allW)\(allD./allW);


% allW(allW == 0) = 1/nTrials;
eps_set = find(allD < epsilon);

ovW = sp.getBasis(allD); ovW = mean(ovW(eps_set,1)./sum(ovW(eps_set,:),2));
% ovW = length(eps_set)/sum(~isnan(allD));
allWeights = ovW./allW;

%%
avgM = mean(allM,2);
n_xbins = 100;
x_bin_edges = prctile(allM(:,1),linspace(0,100,n_xbins+1));
x_bin_centers = 0.5*x_bin_edges(1:end-1) + 0.5*x_bin_edges(2:end);

eps_set = find(allD < 0.02);

marg_prob = hist(allM(:,1),x_bin_centers);
marg_prob = marg_prob/sum(marg_prob);
cond_prob = histc(reshape(allM(eps_set,1),[],1),x_bin_edges);
cond_prob = cond_prob/sum(cond_prob)*n_xbins;

% [~,binids] = hist(mean(all
%%
n_EP_bins = 200;
EP_bin_edges = prctile(allD,linspace(0,100,n_EP_bins+1));
EP_bin_centers = 0.5*EP_bin_edges(1:end-1) + 0.5*EP_bin_edges(2:end);
EP_xcov = nan(n_EP_bins,1);
EP_xcov2 = nan(n_EP_bins,1);
EP_avg = nan(n_EP_bins,1);
for ee = 1:n_EP_bins
    curset = find(allD >= EP_bin_edges(ee) & allD < EP_bin_edges(ee+1));
%     [cond_prob,binids] = histc(allM(curset,1),x_bin_edges);
%     cond_prob = cond_prob/sum(cond_prob)*n_xbins;
%     EP_xcov(ee) = mean(allX(curset).*cond_prob(binids));
    EP_xcov2(ee) = mean(allX(curset));
    EP_avg(ee) = mean(allS(curset,1));
end

%%
n_xbins = 100;
x_bin_edges = prctile(allM(:),linspace(0,100,n_xbins+1));
x_bin_centers = 0.5*x_bin_edges(1:end-1) + 0.5*x_bin_edges(2:end);

[~,binids] = histc(allM(:,1),x_bin_edges);
for ee = 1:n_xbins
    x_avg(ee) = mean(allS(binids==ee),1);
end

