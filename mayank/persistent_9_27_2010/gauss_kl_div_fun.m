function kl = gauss_kl_div_fun(meandiff,covar1,covar2)

p = size(meandiff,1);
% meandiff = meandiff(:);

dc1 = det(covar1);
dc2 = det(covar2);

% mah_dist = meandiff'*inv(covar2)*meandiff;

mah_dist = zeros(size(meandiff,2),1);
for i = 1:length(mah_dist)
    mah_dist(i) = meandiff(:,i)'*inv(covar2)*meandiff(:,i);
end

kl = 1/2*(log(dc2/dc1) + trace(inv(covar2)*covar1) + mah_dist - p);