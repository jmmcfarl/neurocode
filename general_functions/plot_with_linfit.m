function fit_p = plot_with_linfit(xpts,ypts)

fit_p = robustfit(xpts(:),ypts(:),'bisquare',[],'off');
xrange = linspace(min(xpts),max(xpts),100);
% pred_y = fit_p(1) + xrange*fit_p(2);
pred_y =xrange*fit_p(1);
f1 = figure;
plot(xpts,ypts,'k.')
hold on
plot(xrange,pred_y,'r')
