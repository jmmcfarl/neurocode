function gmm_fit = fit_GMM(X,n_clusts,varargin)
% gmm_fit = fit_GMM(X,n_clusts,varargin)
% helper function to avoid errors when fitting GMMs
try
    gmm_fit = gmdistribution.fit(X,n_clusts,varargin{:});
catch
    gmm_fit = nan;
end