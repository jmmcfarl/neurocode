function [B] = mvgauss_obslike_varmean(emiss,hmm)

%compute likelihood vs time for each state
B = zeros(hmm.T,hmm.K);
for l=1:hmm.K
    lf_meanfun = hmm.state(l).lf_meanfun;
    if hmm.p > 1
        hf_mean = hmm.state(l).hf_mean;
        state_meanfun = [lf_meanfun repmat(hf_mean,hmm.T,1)];
    else
        state_meanfun = lf_meanfun;
    end
    state_var = hmm.state(l).var;
    d = emiss-state_meanfun;
    B(:,l) = mvnpdf(d,zeros(hmm.T,hmm.p),state_var);
end


