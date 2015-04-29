Nstims = 50;
rate_range = [0 0.5];
rates = rand(Nstims,1)*range(rate_range) + rate_range(1);

Ntrials = 500;
nrpts = 100;

for rr = 1:nrpts
    %%
    yobs = poissrnd(repmat(rates,1,Ntrials));
    psth = mean(yobs,2);
    psth_var(rr) = var(psth);
    
    %%
    [II,JJ] = meshgrid(1:Ntrials);
    yprods = reshape(bsxfun(@times,yobs,reshape(yobs,Nstims,1,Ntrials)),Nstims,[]);
    upairs = II > JJ;
    
    pair_tvar = mean(mean(yprods(:,upairs)));
    msq = mean(yobs(:)).^2;
    pair_var(rr) = pair_tvar - msq;
end


%%
Ntrials = 25;
Nstims = 100;
xx = rand(Ntrials,Nstims)*2*pi;
phases = rand(1,Nstims)*2*pi;
rates = (sin(bsxfun(@plus,xx,phases)) + 1)*0.5;

nrpts = 50;
nbins = 10;
xbin_edges = linspace(0,2*pi,nbins+1);
    curn = Ntrials/nbins;
%%
psth_var = nan(nrpts,1);
psth_var_cor = nan(nrpts,1);
pair_var = nan(nrpts,1);
for rr = 1:nrpts
    %%
    rr
    yobs = poissrnd(rates);
    
    cur_psth_var = nan(Nstims,1);
    cur_psth_var_cor = nan(Nstims,1);
    cur_noise_var = nan(Nstims,1);
    for ii = 1:Nstims
        [~,binids] = histc(xx(:,ii),xbin_edges);
        psth = nan(nbins,1); atn = nan(nbins,1);
        for bb = 1:nbins
            psth(bb) = nanmean(yobs(binids==bb,ii));
        end
        cur_psth_var(ii) = nanvar(psth);
        cur_noise_var(ii) = nanvar(yobs(:,ii));
        cur_psth_var_cor(ii) = cur_psth_var(ii).*(curn/(curn-1)) - cur_noise_var(ii)'./(curn-1)'; %sahani linden correction for PSTH sampling noise
    end
    psth_var(rr) = nanmean(cur_psth_var);  
%     psth_var_cor(rr) = psth_var(rr).*(curn/(curn-1)) - nanmean(cur_noise_var)'./(curn-1)'; %sahani linden correction for PSTH sampling noise
    psth_var_cor(rr) = nanmean(cur_psth_var_cor);
    
    %%
    eps = 0.2;
    
    knot_pts = [0 0 0 0:0.5:2 2 2 2];
    eval_xx = 0:0.1:2;
    
    [II,JJ] = meshgrid(1:Ntrials);
    upairs = II > JJ;
    
    ybar = yobs - mean(yobs(:));
    
    % cur_pair_var = nan(Nstims,1);
    % cur_spline_var = nan(Nstims,1);
    allD = [];
    allY = [];
    for ii = 1:Nstims
        Dmat = abs(squareform(pdist(xx(:,ii),@(Xi,Xj) circ_dist(Xi,Xj))));
        % uset = Dmat <= eps;
        % uset = uset & upairs;
        
        yprods = bsxfun(@times,ybar(:,ii),ybar(:,ii)');
        allD = cat(1,allD,Dmat(upairs));
        allY = cat(1,allY,yprods(upairs));
    end
    sp = fastBSpline.lsqspline(knot_pts,3,allD,allY);
    spline_var(rr) = sp.evalAt(0);
    
    pair_var(rr) = mean(allY(allD <= eps));
    
    
    % msq = mean(yobs(:)).^2;
    % pair_var(rr) = pair_tvar - msq;
    
    % allD = Dmat(upairs);
    % allY = yprods(upairs);
    % EP_bin_edges = xbin_edges;
    % for bb = 1:length(EP_bin_edges)-1
    %     curset = find(allD >= EP_bin_edges(bb) & allD < EP_bin_edges(bb+1));
    %     yavg(bb) = mean(allY(curset));
    % end
    
end

