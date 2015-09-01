function [crossPoint] = find_GMM_xp(u,sig,t)

x_ax = [-5:.01:5];

g1 = t(1)/sqrt(2*pi*sig(1)^2)*exp(-(x_ax-u(1)).^2./(2*sig(1)^2));
g2 = t(2)/sqrt(2*pi*sig(2)^2)*exp(-(x_ax-u(2)).^2./(2*sig(2)^2));
% plot(x_ax,g1,x_ax,g2,'r');pause;clf
start_point = find(g1>g2,1,'first');
crossPoint = x_ax(start_point+find(g1(start_point:end)<g2(start_point:end),1,'first'));