function beta = fit_dist_transform(y,g,disp)

options.MaxFunEvals = 1e4;
options.MaxIter = 1e4;
options.TolFun = 1e-8;
options.TolX = 1e-11;
un_obs_counts = unique(y);
max_obs_cnt = max(un_obs_counts);
emp_pdf = hist(y,0:max_obs_cnt)/length(y);

LB = [0 -1];
UB = [0.1 2];
init_beta = [.0002 1];
beta = fmincon(@(k) internal_cost_function(k,emp_pdf,g,max_obs_cnt),init_beta,...
    [],[],[],[],LB,UB,[],options);

% A = 0:.00005:.005;
% B = [-0.3:0.05:2.5];
% C = nan(length(A),length(B));
% for i = 1:length(A)
%     i
%     for j = 1:length(B)
%         C(i,j) = internal_cost_function([A(i) B(j)],emp_pdf,g,max_obs_cnt);
%     end
% end

if disp==1
    gt = beta(1)*g + beta(2);
    gt = exp(gt);
    n_reps = 100;
    
    sim_counts = poissrnd(repmat(gt(:),n_reps,1));
    sim_pdf = hist(sim_counts,0:max_obs_cnt)/length(sim_counts);
    figure
    plot(emp_pdf,'o-')
    hold on
    plot(sim_pdf,'ro-')
end
end

function cost = internal_cost_function(beta,emp_pdf,g,max_obs_cnt)

n_reps = 100;

gt = beta(1)*g + beta(2);
gt = exp(gt);

sim_counts = poissrnd(repmat(gt(:),n_reps,1));
% if sum(sim_counts > max_obs_cnt) > length(gt)*0.01
%     cost = inf;
% else
    sim_pdf = hist(sim_counts,0:max_obs_cnt)/length(sim_counts);
    % kls = emp_pdf.*log(emp_pdf./sim_pdf);
    kls = sim_pdf.*log(sim_pdf./emp_pdf);
    kls(sim_pdf == 0) = 0;
    cost = sum(kls);
% end
end