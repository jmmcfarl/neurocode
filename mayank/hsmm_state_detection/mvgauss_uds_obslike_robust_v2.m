function [B,bad_model] = mvgauss_uds_obslike_robust_v2(emiss,hmm,rel_prob,seg_number)

eps = 1e-200;
% min_mllik = -2.5;
% min_mllik = -1.5;
min_mllik = hmm.min_mllik;

%compute likelihood vs time for each state
curT = size(emiss,1);
B = zeros(curT,hmm.K);
for l=1:hmm.K
    if strcmp(hmm.meantype,'variable')
        meanfun = hmm.state(l).meanfun{seg_number};
%         if hmm.p > 1
%             fixedmean = hmm.state(l).fixedmean;
%             state_meanfun = [meanfun repmat(fixedmean,curT,1)];
%         else
            state_meanfun = meanfun;
%         end
    else
        state_meanfun = repmat(hmm.state(l).fixedmean,curT,1);
    end
    state_var = hmm.state(l).var;
    d = emiss-state_meanfun;
    B(:,l) = mvnpdf(d,zeros(curT,hmm.p),state_var);
end

B(B < eps) = eps; %any likelihoods less than epsilon are set to epsilon

max_obs_lik = max(log(B),[],2); %compute the maximum of either states likelihood for the current parameters

bad_model = find(max_obs_lik < min_mllik); %find instances where either states likelihood is far from the model

if ~isempty(bad_model)
    B_r = zeros(curT,hmm.K);
   
    %compute the average covariance matrix across states
    avg_var = 0;
    for l = 1:hmm.K
        avg_var = avg_var + rel_prob(l)*hmm.state(l).var;
    end
    avg_var = avg_var/hmm.K;
    
    %compute the observation likelihood under the model where the state
    %covariance matrices are constrained to be equal 
    for l=1:hmm.K
        if strcmp(hmm.meantype,'variable')
            meanfun = hmm.state(l).meanfun{seg_number};
%             if hmm.p > 1
%                 fixedmean = hmm.state(l).fixedmean;
%                 state_meanfun = [meanfun repmat(fixedmean,curT,1)];
%             else
                state_meanfun = meanfun;
%             end
        else
            state_meanfun = repmat(hmm.state(l).fixedmean,curT,1);
        end
        d = emiss-state_meanfun;
        B_r(:,l) = mvnpdf(d,zeros(curT,hmm.p),avg_var);
    end  
    B(bad_model,:) = B_r(bad_model,:);
end