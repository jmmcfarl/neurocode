function plot_errorbar_quantiles(X,Y,quantiles,color)
% h = plot_errorbar_quantiles(X,Y,quantiles,color)
% make an errorbar plot with center points and errorbounds at specified quantiles

if length(quantiles) ~= 3
    error('must provide three quantiles for plotting');
end
if nargin < 4
    color = 'b';
end


Yquants = prctile(Y,quantiles);
UE = Yquants(3,:) - Yquants(2,:);
LE = Yquants(2,:) - Yquants(1,:);

errorbar(X,Yquants(2,:),LE,UE,color);
hold on
plot(X,Yquants(2,:),'o','color',color);