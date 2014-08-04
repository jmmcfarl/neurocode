function [ll,grad] = const_ll_opt(b,g,gs)

eg = exp(g+b);
r = log(1+eg);

egs = exp(gs+b);
rs = log(1+egs);
rs(rs < 1e-20) = 1e-20;
% Standard point-process likelihood
ll = sum(log(rs)) - sum(r);

wg = eg./(1+eg);
ers = (egs./(1+egs))./rs;
% grad = ers'*gs-wg'*g;             
grad = sum(ers)-sum(wg);  % grad w.r.t. b

ls = length(gs); %normalize by number of spikes
ll=-ll/ls;
grad=-grad'/ls; %'

end
