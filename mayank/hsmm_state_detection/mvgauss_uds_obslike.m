function [B] = mvgauss_uds_obslike(emiss,hmm,seg_number)

%compute likelihood vs time for each state
curT = size(emiss,1);
B = zeros(curT,hmm.K);
for l=1:hmm.K
    if strcmp(hmm.meantype,'variable')
        meanfun = hmm.state(l).meanfun{seg_number};
        if hmm.p > 1
            fixedmean = hmm.state(l).fixedmean;
            state_meanfun = [meanfun repmat(fixedmean,curT,1)];
        else
            state_meanfun = meanfun;
        end
    else
        state_meanfun = repmat(hmm.state(l).fixedmean,curT,1);
    end
    state_var = hmm.state(l).var;
    d = emiss-state_meanfun;
    B(:,l) = mvnpdf(d,zeros(curT,hmm.p),state_var);
end


