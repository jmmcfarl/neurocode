function LLs = get_pgabor_tempmod_LLs(temp_params,Robs,X_gabor_out,blockids)

[~,NT,K] = size(X_gabor_out);

X_energy_out = squeeze(sqrt(X_gabor_out(1,:,:).^2 + X_gabor_out(2,:,:).^2));

ws = temp_params.ws(2,:);
wc = temp_params.ws(1,:);
c = temp_params.c;

n_blocks = length(unique(blockids));

gmat = zeros(NT,K);
for i = 1:n_blocks
    cur_set = find(blockids == i);
    gmat(cur_set,:) = gmat(cur_set,:) + ws(i)*squeeze(X_gabor_out(1,cur_set,:));
    gmat(cur_set,:) = gmat(cur_set,:) + wc(i)*X_energy_out(cur_set,:);
    gmat(cur_set,:) = gmat(cur_set,:) + c(i);
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

