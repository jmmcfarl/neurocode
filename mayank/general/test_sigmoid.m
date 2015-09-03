% test sigmoid fit
% needs my_sigmoid.m
clear all
% test data
t = [-600:0.1:600];
ainit = [30 80 20 4]
xx= ainit(1)*((1.0./(1.0+exp((ainit(2)-t)/ainit(3))))+0.1*(rand(size(t))-0.5))+ainit(4);
init_guess = [2*mean(xx) 10 1.0 min(xx)]

% non-linear fit
[betafit, resid, jacobian] = nlinfit(t,xx,@my_sigmoid,init_guess);
mean_rms_error = mean(sqrt((resid.^2)));

betafit

% predicted 
xx_hat = betafit(1)*(1.0./(1.0 + exp((betafit(2)-t)/betafit(3)))) + betafit(4);

% just some plotting
Fig = figure(1);
clf;
fw = 12;

subplot(2,1,1)
plot(t,xx,'bo','markersize',6)
hold on
plot(t,xx_hat,'rx-','markersize',6,'linewidth',2)
grid on;
axis tight
toh = text(1,1,strcat('mean rms error = ',num2str(mean_rms_error)));
set(toh,'units','normal','position',[0.1 0.8],'fontsize',fw)
lg=legend('measured','estimated');
set(lg,'fontsize',10,'location','SouthEast');
set(gca,'fontsize',fw,'fontweight','bold');
subplot(2,1,2)
bar(t,resid)
axis tight
set(gca,'fontsize',fw,'fontweight','bold');

