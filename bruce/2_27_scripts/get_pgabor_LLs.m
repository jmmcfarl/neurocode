function LLs = get_pgabor_LLs(gabor_params,Robs,X_gabor_out)

[~,NT,K] = size(X_gabor_out);

X_energy_out = squeeze(sqrt(X_gabor_out(1,:,:).^2 + X_gabor_out(2,:,:).^2));

ws = gabor_params(8:9);
wc = gabor_params(7);
c = gabor_params(10);

gmat = c*ones(NT,K);
gmat = gmat + ws(1)*squeeze(X_gabor_out(1,:,:));
gmat = gmat + ws(2)*squeeze(X_gabor_out(2,:,:));
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

