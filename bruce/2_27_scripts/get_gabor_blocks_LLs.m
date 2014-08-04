function LLs = get_gabor_blocks_LLs(gabor_model,Robs,X_phase_out,X_energy_out,blockids)

[n_phases,NT,K] = size(X_phase_out);
n_blocks = length(unique(blockids));
n_kerns = n_phases+1;

ws = gabor_model.ws(1:2,:);
wc = gabor_model.ws(3,:);
c = gabor_model.c;

gmat = zeros(NT,K);
LLs = zeros(NT,K);
for blockid = 1:n_blocks
    cur_set = find(blockids == blockid);
    
    for i = 1:2
        gmat(cur_set,:) = gmat(cur_set,:) + squeeze(ws(i,blockid)*X_phase_out(i,cur_set,:));
    end
    gmat(cur_set,:) = gmat(cur_set,:) + wc(blockid)*X_energy_out(cur_set,:);
    gmat(cur_set,:) = gmat(cur_set,:) + c(blockid);
end

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
