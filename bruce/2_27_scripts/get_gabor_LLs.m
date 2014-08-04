function LLs = get_gabor_LLs(gabor_model,Robs,X_phase_out,X_energy_out)

[n_phases,NT,K] = size(X_phase_out);

ws = gabor_model.ws(1:n_phases);
wc = gabor_model.ws(n_phases+1);
c = gabor_model.c;

gmat = c*ones(NT,K);
for i = 1:n_phases
    gmat = gmat + squeeze(ws(i)*X_phase_out(i,:,:));
end
gmat = gmat + wc*X_energy_out;

too_large = find(gmat > 100);
expg = exp(gmat);
r = log(1+expg);
r(too_large) = gmat(too_large);

r(r < 1e-20) = 1e-20; %minimum predicted rate

LLs = repmat(Robs,1,K).*log(r)-r;

set1 = find(Robs > 10); %use stirling's approx for log(R!)
set2 = find(Robs <= 10); %use log(R!) directly
LLs(set1,:) = LLs(set1,:) - repmat(Robs(set1).*log(Robs(set1))+Robs(set1),1,K);
LLs(set2,:) = LLs(set2,:) - repmat(log(factorial(Robs(set2))),1,K);

