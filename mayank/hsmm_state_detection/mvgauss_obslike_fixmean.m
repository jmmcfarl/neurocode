function [B] = mvgauss_obslike_fixmean(emiss,hmm)

%compute likelihood vs time for each state
B = zeros(hmm.T,hmm.K);
for l=1:hmm.K
    state_var = hmm.state(l).var;
    d = emiss-repmat(hmm.state(l).mean,hmm.T,1);
    B(:,l) = mvnpdf(d,zeros(hmm.T,hmm.p),state_var);
end


