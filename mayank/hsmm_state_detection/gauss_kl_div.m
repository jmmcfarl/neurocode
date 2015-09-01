function kl = gauss_kl_div(meandiff,covar1,covar2)

p = size(meandiff,1);
% meandiff = meandiff(:);

dc1 = det(covar1);
dc2 = det(covar2);

% mah_dist = meandiff'*inv(covar2)*meandiff;
if size(meandiff,1) > size(meandiff,2)
    meandiff = meandiff';
end

temp = zeros(size(meandiff,2),1);
for i = 1:length(temp)
    temp(i) = meandiff(:,i)'*inv(covar2)*meandiff(:,i);
end
mah_dist = mean(temp);

kl = 1/2*(log(dc2/dc1) + trace(inv(covar2)*covar1) + mah_dist - p);