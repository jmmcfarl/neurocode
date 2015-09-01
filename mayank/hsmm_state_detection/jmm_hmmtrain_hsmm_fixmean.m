function [hmm] = jmm_hmmtrain_hsmm_fixmean(hmm,emiss,num_it)

addpath('C:\Code\maphmmbox\')
%%
K = hmm.K;
p = hmm.p;
T = hmm.T;
Fs = hmm.Fs;
Fs  = hmm.Fs;
total_dur = T/Fs;
t_axis = (1:T)/Fs;
eps = 1e-100;

%%

%initialize forward and backward messages
alpha=zeros(T,K);
beta=zeros(T,K);
gamma=zeros(T,K);

%set model parameters to their initialized values
Pi = hmm.Pi;
A = hmm.A;

lpost=0;
LP = [];

for it = 1:num_it

    %compute likelihood of observations for each model under current
    %parameters
    [B] = mvgauss_obslike_fixmean(emiss,hmm);
    B(B < eps) = eps; %prevent likelihoods from rounding to 0

    scale=zeros(T,1); %initialize rescaling parameters

    %compute rescaled forward messages
    alpha(1,:)=Pi.*B(1,:);
    scale(1)=sum(alpha(1,:));
    alpha(1,:)=alpha(1,:)/scale(1);
    for t=2:T
        alpha(t,:)=(alpha(t-1,:)*A).*B(t,:);
        scale(t)=sum(alpha(t,:)); % P(X_i | X_1 ... X_{i-1})
        alpha(t,:)=alpha(t,:)/scale(t);
    end

    %compute rescaled backward messages
    beta(T,:)=ones(1,K)/scale(T);
    for t=T-1:-1:1
        beta(t,:)=(beta(t+1,:).*B(t+1,:))*(A')/scale(t);
    end

    gamma=(alpha.*beta);
    gamma=rdiv(gamma,rsum(gamma));
    gammasum=sum(gamma);

    %compute chi
    chi=zeros(T-1,K*K);
    for t=1:T-1
        temp=A.*(alpha(t,:)' * (beta(t+1,:).*B(t+1,:)));
        chi(t,:)=temp(:)'/sum(temp(:));
    end

    lscale = log(scale);


    oldlpost=lpost;
    loglik=sum(lscale);
    % Transition Probs and probs of first state node
    LP=[LP; loglik];


    %%%% M STEP
    % transition matrix
    sumchi=sum(chi,1);	% expectation over time
    sumchi=reshape(sumchi,K,K);
    sumchin=sum(sumchi,2);
    A=rdiv(sumchi,sumchin); %update transition matrix

    % initial state
    Pi=zeros(1,K);
    Pi=Pi+gamma(1,:);
    Pi=Pi./sum(Pi);

    var_est = zeros(K,p,p);
    for l=1:K
        mean_est(l,:) = sum(repmat(gamma(:,l),1,p).*emiss)/sum(gamma(:,l));
        mdiff = emiss-repmat(mean_est(l,:),T,1);
        var_est(l,:,:) = ((repmat(gamma(:,l),1,p).*mdiff)'*mdiff)/sum(gamma(:,l));
    end

    %update all parameters in hmm
    hmm.A = A;
    hmm.Pi = Pi;
    for l = 1:K
        new_covar = squeeze(var_est(l,:,:));
        if det(new_covar) > 0
            hmm.state(l).var = squeeze(var_est(l,:,:));
        else
            disp('warning: singular covariance')
        end
        hmm.state(l).mean = mean_est(l,:);
    end

    fprintf('iteration %i log posterior = %f \n',it,loglik);
    %         jmm_print_hmm_params(hmm);
    %     pause

end