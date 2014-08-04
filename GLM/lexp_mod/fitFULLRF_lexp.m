function [fit1 exitflag] = fitFULLRF_lexp(fit0, X, Robs, max_evals, method, silent,optTol,progTol,sa_params)

options = [];
if silent == 1
    options.Display = 'off';
else
    options.Display = 'iter';
end
options.maxFunEvals = max_evals;
options.maxIter = max_evals;
options.Method = method;
if nargin < 7
    options.optTol = 1e-5;
else
    options.optTol = optTol;
end
if nargin < 8
    options.progTol = 1e-7;
else
    options.progTol = progTol;
end

if nargin < 6
    silent = 1;
end

fsdim = fit0.mods(1).fsdim;
if strcmp(fit0.basis,'white')
    pkern_len = size(fit0.pix_conv_mat,2);
elseif strcmp(fit0.basis,'pix')
    pkern_len = size(X,2);
end
kern_t = pkern_len/fsdim;

stimlen = size(X,1);

nmods = length(fit0.mods);
initial_params = [];
for m = 1:nmods
    cur_pkern = fit0.mods(m).pix';
    initial_params = [initial_params cur_pkern]; %add STCBcoefs to initial param vector
    
    NLx = fit0.mods(m).nlx;
    NL = fit0.mods(m).nly;
    if strcmp(fit0.mods(m).nltype,'uncon')
        %compute derivative of non-linearity
        fpr = zeros(1,length(NLx)-1);
        for n = 1:length(fpr)
            fpr(n) = (NL(n+1)-NL(n))/(NLx(n+1)-NLx(n));
        end
        fprimes{m} = fpr;
    else
        fprimes{m} = [];
    end   
end
initial_params(end+1) = fit0.const; %add constant offset term to params

lambda = [];
for i = 1:nmods
    lambda = [lambda fit0.mods(i).lambda_L1x'.*ones(1,pkern_len)];
end
lambda = [lambda 0]/sum(Robs);

if isempty(sa_params)
    if max(lambda) > 0
        [params LL] = L1General2_PSSas(@(K) FULLBF_lexp_LLinternal(K, Robs, X, fit0,fprimes),initial_params',lambda',options,0);
%         [params LL] = L1General2_PSSgb(@(K) FULLBF_lexp_LLinternal(K, Robs, X, fit0,fprimes),initial_params',lambda',options,0);
        exitflag = 1;
    else
        [params LL exitflag] = minFunc( @(K) FULLBF_lexp_LLinternal(K, Robs, X, fit0,fprimes), initial_params', options );
        % [params LL exitflag] = fminunc( @(K) FULLBF_lexp_LLinternal(K, Robs, X, fit0), initial_params' );
    end
else
    %initialize overall best parameter values
    ov_best_LL = Inf;
    ov_best_params = initial_params;
    
    %compute spike bin representation
    spk_vals = unique(Robs); spk_vals(spk_vals==0) = [];
    spikebins = [];
    for i = 1:length(spk_vals)
        cur_set = find(Robs == spk_vals(i));
        spikebins = [spikebins; repmat(cur_set(:),spk_vals(i),1)];
    end
    
    cooling_schedule = sa_params.T0*sa_params.alpha.^((1:sa_params.n_sa_iter)-1);
    options.maxIter = sa_params.it_per_sa;
    for ii = 1:sa_params.n_sa_iter
        if max(lambda) > 0
            [params LL] = L1General2_PSSas(@(K) FULLBF_lexp_LLinternal(K, Robs, X, fit0),initial_params',lambda',options,0);
            exitflag = 1;
        else
            [params LL exitflag] = minFunc( @(K) FULLBF_lexp_LLinternal(K, Robs, X, fit0), initial_params', options );
        end
        fit_params = params;
        
        k_mat = zeros(nmods,pkern_len);
        for nn = 1:nmods
            k_mat(nn,:) = params((nn-1)*pkern_len+(1:pkern_len));
        end
        
        group_merge = randperm(nmods);
        group_merge = group_merge(1:2);
        net = sum(k_mat(group_merge,:));
        k_mat(group_merge(1),:) = net;
        k_mat(group_merge(2),:) = 0;
        k_mat = k_mat';
        cur_x = k_mat(:);
        cur_x = [cur_x; params(end)]; %add constant back
        
%                 group_flip = ceil(rand*nmods);
%         %fprintf('Flip mod %d\n',group_flip);
%         k_mat(group_flip,:) = -k_mat(group_flip,:);
%         k_mat = k_mat';
%         g_mat = X*k_mat;
% 
%         tempmod = fit0;
%         tempmod.const = params(end);
%         cur_x = k_mat(:);
%         cur_x = [cur_x; tempmod.const]; %add constant back
% %         tempmod.mods(group_flip).w = -tempmod.mods(group_flip).w;
%         for n = 1:nmods
%             cur_kern = cur_x((n-1)*pkern_len+(1:pkern_len));
%             if strcmp(tempmod.basis,'white')
%                 cur_wkern = cur_x'*tempmod.kern_conv_mat;
%                 tempmod.mods(n).k = cur_wkern(:);
%                 tempmod.mods(n).pix = cur_kern(:);
%             elseif strcmp(tempmod.basis,'pix')
% %                 kern_out = X * cur_kern;
%                 tempmod.mods(n).k = cur_kern(:);
%                 tempmod.mods(n).pix = cur_kern(:);
%             end
%         end
%         tempmod = fitWeights_full(tempmod,g_mat,spikebins,1,1:nmods);
        
        test_LL = FULLBF_lexp_LLinternal(cur_x, Robs, X, fit0);
        
        acc_prob = min(1,exp((LL-test_LL)/cooling_schedule(ii)));
        is_acc = rand < acc_prob;
        old_LL = LL;
        fprintf('Current LP: %.3f  Proposal LP: %.3f  Temp: %.3f  Acc Prob: %.3f\n',...
            LL,test_LL,cooling_schedule(ii),acc_prob);
        if is_acc
            params = cur_x;
            LL = test_LL;
            fprintf('Reflection accepted\n');
        else
            fprintf('Reflection rejected\n');
        end
        
        if LL < ov_best_LL
            ov_best_LL = LL;
            ov_best_params = params;
        end
        
                     figure(4)
                     clf
                     subplot(3,1,1)
                     plot(initial_params(1:30))
                     hold on
                     plot(initial_params(31:60),'r')
                     subplot(3,1,2)
                     plot(fit_params(1:30))
                     hold on
                     plot(fit_params(31:60),'r')
                     subplot(3,1,3)
                     plot(params(1:30))
                     hold on
                     plot(params(31:60),'r')
                     fprintf('Cur LL %.3f  Prop LL %.3f  Temp: %.3f  Acc Prob: %.3f\n',old_LL,test_LL,cooling_schedule(ii),acc_prob);
                     input('')
        
        initial_params = params';
        
    end
    params = ov_best_params;
    LL = ov_best_LL;
    
end


%%%%%%%%%%%%%%%%%%%%%%
% if silent == 0
% %     opts = optimset('GradObj','on','Algorithm','active-set','Display','iter','MaxIter',200,'MaxFunEvals',10000,'TolFun',1e-7);
%     opts = optimset('GradObj','on','LargeScale','on','Display','iter','MaxIter',200,'MaxFunEvals',10000,'TolFun',1e-7);
% else
%     opts = optimset('GradObj','on','Algorithm','active-set','MaxIter',200,'Display','off','MaxFunEvals',10000,'TolFun',1e-7);
% end
% [params LL exitflag] = fminunc( @(K) FULLBF2d_LLinternal(K, Robs, X, fit0), initial_params', opts );



fit0.LP_seq = [fit0.LP_seq LL];
fit0.LP_ax = [fit0.LP_ax fit0.LP_ax(end)+1];

% Reassign variables
fit1 = fit0;
fit1.const = params(end);
for n = 1:nmods
    cur_kern = params((n-1)*pkern_len+(1:pkern_len));
    if strcmp(fit0.basis,'white')
        cur_wkern = cur_kern'*fit0.kern_conv_mat;
        fit1.mods(n).k = cur_wkern(:);
        fit1.mods(n).pix = cur_kern(:);
    elseif strcmp(fit0.basis,'pix')
        kern_out = X * cur_kern;
        fit1.mods(n).k = cur_kern(:);
        fit1.mods(n).pix = cur_kern(:);
    end
end
