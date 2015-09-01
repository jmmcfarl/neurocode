function bd = gauss_bd(meandiff,covar1,covar2)

dc1 = det(covar1);
dc2 = det(covar2);
dc = det((covar1 + covar2)/2);

% mah_dist = meandiff'*inv(covar2)*meandiff;

temp = zeros(size(meandiff,2),1);
for i = 1:length(temp)
    temp(i) = meandiff(:,i)'*inv((covar1+covar2)/2)*meandiff(:,i);
end
mah_dist = mean(temp);

bd = 1/8*mah_dist + 1/2*log(dc/(sqrt(dc1*dc2)));
