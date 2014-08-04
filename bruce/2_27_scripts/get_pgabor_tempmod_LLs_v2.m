function LLs = get_pgabor_tempmod_LLs_v2(gabor_params,scales,offsets,Robs,X_gabor_out,blockids)

[~,NT,K] = size(X_gabor_out);

X_energy_out = squeeze(sqrt(X_gabor_out(1,:,:).^2 + X_gabor_out(2,:,:).^2));
gmat = gabor_params(7)*X_energy_out + gabor_params(8)*squeeze(X_gabor_out(1,:,:))...
    + gabor_params(9)*squeeze(X_gabor_out(2,:,:));
n_blocks = length(unique(blockids));

for i = 1:n_blocks
    cur_set = find(blockids == i);
    gmat(cur_set,:) = scales(i)*gmat(cur_set,:)+offsets(i);
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

