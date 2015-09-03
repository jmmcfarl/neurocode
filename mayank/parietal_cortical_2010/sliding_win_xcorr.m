function [mean_xcorr,std_xcorr,lags] = sliding_win_xcorr(sig1,sig2,maxlag,markers1,markers2)

lags = -maxlag:maxlag;

xcMat = zeros(size(markers1,1),length(lags));
for i = 1:size(markers1,1) 
   seg1 = detrend(sig1(markers1(i,1):markers1(i,2)));
   seg2 = detrend(sig2(markers2(i,1):markers2(i,2)));
   xcMat(i,:) = xcov(seg1,seg2,maxlag,'coeff');   
end

mean_xcorr = mean(xcMat);
std_xcorr = std(xcMat);