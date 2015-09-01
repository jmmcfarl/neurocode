function [B] = mvgauss_uds_obslike_robust_duallfp(emiss,hmm,seg_number,min_mllik)

eps = 1e-200;
% min_mllik = -6;
min_mllik = hmm.min_mllik;

%compute likelihood vs time for each state
curT = size(emiss,1);
B = zeros(curT,hmm.K);
for l=1:hmm.K
    if strcmp(hmm.meantype,'variable')
        meanfun = hmm.state(l).meanfun{seg_number};
        meanfun2 = hmm.state(l).meanfun2{seg_number};
        state_meanfun = [meanfun meanfun2];
    else
        state_meanfun = repmat(hmm.state(l).fixedmean,curT,1);
    end
    state_var = hmm.state(l).var;
    d = emiss-state_meanfun;
    B(:,l) = mvnpdf(d,zeros(curT,hmm.p),state_var);
end

B(B < eps) = eps;

max_obs_lik = max(log(B),[],2);

bad_model = find(max_obs_lik < min_mllik);

if ~isempty(bad_model)
    B_r = zeros(curT,hmm.K);
    avg_var = 0;
    for l = 1:hmm.K
        avg_var = avg_var + hmm.state(l).var;
    end
    avg_var = avg_var/hmm.K;
    for l=1:hmm.K
        if strcmp(hmm.meantype,'variable')
            meanfun = hmm.state(l).meanfun{seg_number};
            meanfun2 = hmm.state(l).meanfun2{seg_number};
            state_meanfun = [meanfun meanfun2];
        else
            state_meanfun = repmat(hmm.state(l).fixedmean,curT,1);
        end
        d = emiss-state_meanfun;
        B_r(:,l) = mvnpdf(d,zeros(curT,hmm.p),avg_var);
    end  
    B(bad_model,:) = B_r(bad_model,:);
end