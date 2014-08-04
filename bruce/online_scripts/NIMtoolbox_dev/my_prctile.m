function percentiles = my_prctile(x,p)
%
% percentiles = my_prctile(x,p)
%
% Calculate pth percentiles of x (p can be vector-valued).

%%
x = x(:);
p = p(:);
xlen = length(x);

x = sort(x);
un_ax = [0 100*(0.5:(xlen-0.5))./xlen 100]'; %uniform spacing from 0 to 100 of length xlen + 2
x = [x(1); x; x(end)]; %make length xlen + 2
percentiles = interp1q(un_ax,x,p); %find pth percentiles
