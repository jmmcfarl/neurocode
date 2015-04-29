% clear all
% close all
% 
% nTrials = 400;
% nStims = 50;
% xx = linspace(0,2*pi,nStims);
% 
% alpha = 2;
% beta = 0.5;
% trialGains = gamrnd(alpha,beta,nTrials,nStims);
% 
% rates1 = 3*log(1+3*exp(5+bsxfun(@times,sin(xx),trialGains)));
% rates2 = log(1+exp(bsxfun(@times,cos(xx),trialGains)));
% % rates1 = bsxfun(@times,(sin(xx)),trialGains);
% % rates2 = bsxfun(@times,(cos(xx)),trialGains);
% % rates1 = (rates1 +5)*2; rates1(rates1 < 0) = 0;
% % rates2 = (rates2 +5)*2; rates2(rates2 < 0) = 0;
% 
% yy1 = poissrnd(rates1); yy2 = poissrnd(rates2);
% 
% psth1 = mean(yy1); psth2 = mean(yy2);
% 
% 
% %%
% % tot_xcov = mean((yy1(:)-mean(yy1(:))).*(yy2(:) - mean(yy2(:))));
% % psth_xcov = mean((psth1-mean(psth1)).*(psth2-mean(psth2)));
% 
% %%
% yy1_bar = yy1 - mean(yy1(:));
% yy2_bar = yy2 - mean(yy2(:));
% % 
% % [II,JJ] = meshgrid(1:nTrials);
% % upairs = II ~= JJ;
% % 
% % allX = [];
% % allD = [];
% % for ii = 1:nStims
% %    curD = squareform(pdist(trialGains(:,ii)));
% %    yyprod = bsxfun(@times,yy1_bar(:,ii),yy2_bar(:,ii)');
% %    allD = cat(1,allD,curD(upairs));
% %    allX = cat(1,allX,yyprod(upairs));
% % end
% 
% % yy1_bar = bsxfun(@minus,yy1,psth1);
% % yy2_bar = bsxfun(@minus,yy2,psth2);
% [II,JJ] = meshgrid(1:nTrials);
% upairs = II > JJ;
% 
% allX = [];
% allD = [];
% for ii = 1:nStims
% %    curD = squareform(pdist(yy2_bar(:,ii)));
%    curD = squareform(pdist(trialGains(:,ii)));
%    yyprod = bsxfun(@times,yy1_bar(:,ii),yy1_bar(:,ii)');
%    allD = cat(1,allD,curD(upairs));
%    allX = cat(1,allX,yyprod(upairs));
% end
% 
% %%
% n_EP_bins = 100;
% EP_bin_edges = prctile(allD,linspace(0,100,n_EP_bins+1));
% EP_bin_centers = 0.5*EP_bin_edges(1:end-1) + 0.5*EP_bin_edges(2:end);
% EP_xcov = nan(n_EP_bins,1);
% for ee = 1:n_EP_bins
%     curset = find(allD >= EP_bin_edges(ee) & allD < EP_bin_edges(ee+1));
%    EP_xcov(ee) = mean(allX(curset));
% end
% % poss_Ds = 0:50;
% % EP_xcov = nan(length(poss_Ds),1);
% % for ee = 1:length(poss_Ds)
% %     curset = find(allD == poss_Ds(ee));
% %    EP_xcov(ee) = mean(allX(curset));
% % end
% 
% 
% pair_xcov = mean(allX(allD == 0));
% pair_psth_cov = mean(allX);
% tot_var = var(yy1(:));
% noise_var = tot_var - pair_xcov;


%%
clear all
close all

nTrials = 200;
nStims = 50;
xx = linspace(0,2*pi,nStims);

alpha = 2;
beta = 1/alpha;
trialGains = gamrnd(alpha,beta,nTrials,nStims);
% trialGains = ones(nTrials,nStims);

rates = 3*log(1+3*exp(5+bsxfun(@times,cos(xx),trialGains)));
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

n_Dbins = 100;
D_bin_edges = prctile(allD,linspace(0,100,n_Dbins+1));
D_bin_centers = 0.5*D_bin_edges(1:end-1) + 0.5*D_bin_edges(2:end);
[margDdist,Dbinids] = histc(allD,D_bin_edges);
Dbinids(Dbinids > n_Dbins) = n_Dbins;
margDdist = margDdist(1:end-1)/length(allD);
%%

allX = [];
allS = [];
allD = [];
allM = [];
allW = [];
all_Coefs = [];
for ii = 1:nStims
%     curpicks = randi(nTrials,randpicks,1);
    
    curD = squareform(pdist(trialGains(:,ii)));
    curD(~upairs) = nan;
    
    yyprod = bsxfun(@times,yy_bar(:,ii)',yy_bar(:,ii));
    yyvals = [yy_bar(II,ii) yy_bar(JJ,ii)];

%     weights = nansum(curD < epsilon)/(nTrials-1);
%     WW = repmat(weights,nTrials,1);
    
    ddvals = [trialGains(II,ii) trialGains(JJ,ii)];
    
    WW = nan(nTrials,nTrials);
    for jj = 1:nTrials
        uset = setdiff(1:nTrials,jj);
%        [curN,curIDs] = histc(curD(jj,:),D_bin_edges);
%        curIDs(curIDs > n_Dbins) = n_Dbins;
%        curN = curN(1:end-1)/(nTrials-1);
%        WW(jj,uset) = curN(curIDs(uset));
WW(jj,uset) = sum(curD(jj,uset) < 0.01)/(nTrials-1);
    end
%     WW = n_Dbins;
%     WW = nan(nTrials,nTrials,n_splines);
% %     ww = nan(nTrials,nTrials,n_splines);
%     for jj = 1:nTrials
%         proj = sp.getBasis(curD(jj,:));
%         corfac = ovWeights./mean(proj);
%         weights = bsxfun(@rdivide,proj,corfac);
%         WW(jj,:,:) = weights;
% %         ww(jj,:,:) = proj;
%     end
%     WW = reshape(WW,[],n_splines);
%     netweights = sum(weights);
%     WW = repmat(netweights,nTrials,1);
%     allW = cat(1,allW(
    
%     cur_coefs = nan(nTrials,n_splines);
%     for jj = 1:nTrials
%         jj
%         otherset = setdiff(1:nTrials,jj);
%     cur_sp = sp.lsqspline(knot_pts,3,curD(jj,otherset),yyprod(jj,otherset));
%     cur_coefs(jj,:) = cur_sp.weights';
%     end
%     cur_sp = sp.lsqspline(knot_pts,3,curD(upairs),yyprod(upairs));
%     cur_coefs = cur_sp.weights';
%     all_Coefs = cat(1,all_Coefs,cur_coefs);
%     
    allD = cat(1,allD,curD(upairs));
    allX = cat(1,allX,yyprod(upairs));
    allS = cat(1,allS,yyvals(upairs,:));
    allM = cat(1,allM,ddvals(upairs,:));
    allW = cat(1,allW,WW(upairs));
end

% uu = find(~isnan(allX) & max(allW,[],2) > 0);
% rr = lsqlin(allW(uu,:),allX(uu))
% w = allW(uu,:)\allX(uu)
%%
% allX = [];
% allS = [];
% allD = [];
% allM = [];
% allW = [];
% for ii = 1:nStims
% %     curpicks = randi(nTrials,randpicks,1);
%     
%     curD = squareform(pdist(trialGains(:,ii)));
%     curD(~upairs) = nan;
%     
%     yyprod = bsxfun(@times,yy_bar(:,ii)',yy_bar(:,ii));
%     yyvals = [yy_bar(II,ii) yy_bar(JJ,ii)];
% 
% %     weights = nansum(curD < epsilon)/(nTrials-1);
% %     WW = repmat(weights,nTrials,1);
%     
%     ddvals = [trialGains(II,ii) trialGains(JJ,ii)];
%     
%     weights = sp.getBasis(curD);
%     weights = reshape(weights(:,1),nTrials,nTrials);
%     netweights = sum(weights);
%     WW = repmat(netweights,nTrials,1);
% %     allW = cat(1,allW(
%     
%     allD = cat(1,allD,curD(upairs));
%     allX = cat(1,allX,yyprod(upairs));
%     allS = cat(1,allS,yyvals(upairs,:));
%     allM = cat(1,allM,ddvals(upairs,:));
%     allW = cat(1,allW,WW(upairs));
% end

%%
n_splines = 5;
% spline_DS = 0:0.5:4;
spline_DS = prctile(allD,100/n_splines:100/n_splines:(100-100/n_splines));
knot_pts = [0 0 0 0 spline_DS max(allD) max(allD) max(allD) max(allD)];
sp = fastBSpline(knot_pts,ones(length(knot_pts) - 4,1));

ovWeights = sp.getBasis(allD);
ovWeights = mean(ovWeights);

n_splines = length(ovWeights);

%% correct for p(detlaX|x) assuming p(x) constant across stimuli
% avgM = mean(allM,2);
n_Dbins = 100;
D_bin_edges = prctile(allD,linspace(0,100,n_Dbins+1));
D_bin_centers = 0.5*D_bin_edges(1:end-1) + 0.5*D_bin_edges(2:end);

% marg_sp = mean(sp.getBasis(allD));

% eps_set = find(allD < 0.02);
% 
[marg_prob,Dbinids] = histc(allD,D_bin_edges);
Dbinids(Dbinids > n_Dbins) = n_Dbins;
marg_prob = marg_prob(1:end-1)/sum(marg_prob);
% cond_prob = histc(reshape(allM(eps_set,1),[],1),x_bin_edges);
% cond_prob = cond_prob/sum(cond_prob)*n_xbins;


n_Xbins = 500;
X_bin_edges = prctile(allM(:,1),linspace(0,100,n_Xbins+1));
X_bin_centers = 0.5*X_bin_edges(1:end-1) + 0.5*X_bin_edges(2:end);

epsilon = prctile(allD,0.5);
eps_set = find(allD <= epsilon);
marg_eprob = length(eps_set)/length(allD);

[~,Xbinids] = histc(allM(:,1),X_bin_edges);
Xbinids(Xbinids == n_Xbins + 1) = n_Xbins;
cond_eprob = nan(n_Xbins,n_Dbins);
for xx = 1:n_Xbins
   curset = find(Xbinids == xx);
   temp = histc(allD(curset),D_bin_edges);
   cond_eprob(xx,:) = temp(1:end-1)/length(curset);
end
% all_weights = marg_eprob./cond_eprob(Xbinids);
% all_weights = marg_prob(1)./cond_prob(Xbinids,1);
% mean(allX(Dbinids == 1).*all_weights(Dbinids==1))
% var(rates(:))
inds = sub2ind([n_Xbins n_Dbins],Xbinids,Dbinids);
pt_cprobs = cond_eprob(inds);
pt_mprobs = marg_prob(Dbinids);
all_weights = pt_mprobs./pt_cprobs;
% all_weights = marg_prob(Dbinids)./cond_eprob(Xbinids,Dbinids);

all_spD = sp.getBasis(allD);
% all_weights = mean(bsxfun(@times,1./cond_spline(Xbinids,:),marg_sp),2);
w = all_spD\(allX.*all_weights)
% 
% eval_xx = linspace(0,2,100);
% sp = fastBSpline(knot_pts,w);
% plot(eval_xx,sp.evalAt(eval_xx))

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

