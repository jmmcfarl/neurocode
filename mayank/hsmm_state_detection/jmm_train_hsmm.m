function [hmm] = jmm_train_hsmm(hmm,emiss,num_it)

addpath('C:\Code\maphmmbox\')
%%
K = hmm.K;
p = hmm.p;
T = hmm.T;
Fs = hmm.Fs;
A_prior = hmm.A_alpha;
Pi_prior = hmm.Pi_alpha;
Fs  = hmm.Fs;
windowSize = hmm.windowSize;
windowSlide = hmm.windowSlide;
total_dur = T/Fs;
t_axis = (1:T)/Fs;
numWins = floor((total_dur-windowSize)/windowSlide);
win_t = (0:numWins-1)*windowSlide+windowSize/2;
win_t = [0 win_t max(t_axis)];
niqf = Fs/2;
% [filt_b,filt_a] = butter(2,hmm.hcf/niqf,'low');
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
    [B] = mvgauss_obslike_varmean(emiss,hmm);
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
    lprior=dirichlet(Pi,hmm.Pi_alpha,1);
    for l=1:K,
        lprior=[lprior dirichlet(A(l,:),hmm.A_alpha(l,:),1)];
    end
    lpost=(loglik+sum(lprior));		% log posterior
    LP=[LP; loglik lprior];


    %%%% M STEP
    % transition matrix
    sumchi=sum(chi,1);	% expectation over time
    sumchi=reshape(sumchi,K,K);
    sumchi=sumchi+(A_prior-1); %using conjugate dirichlet prior
    sumchin=sum(sumchi,2);
    A=rdiv(sumchi,sumchin); %update transition matrix

    % initial state
    Pi=zeros(1,K);
    Pi=Pi+gamma(1,:);
    Pi=(Pi+Pi_prior-1);
    Pi=Pi./sum(Pi);

    lf_meanfun_est = zeros(T,K);
    hf_mean_est = zeros(K,1);
    var_est = zeros(K,p,p);
    for l=1:K
        temp_meanfun = nan(1,numWins);
        for w = 1:numWins
            begT = round((w-1)*windowSlide*Fs)+1;
            endT = begT + round(windowSize*Fs);
            temp_meanfun(w) = sum(gamma(begT:endT,l).*emiss(begT:endT,1))/sum(gamma(begT:endT,l));
        end
        temp_meanfun = [temp_meanfun(1) temp_meanfun temp_meanfun(end)];
        lf_meanfun_est(:,l) = interp1(win_t,temp_meanfun,t_axis);
        lf_meanfun_est(:,l) = filtfilt(filt_b,filt_a,lf_meanfun_est(:,l));
        hf_mean_est(l) = sum(gamma(:,l).*emiss(:,2))/sum(gamma(:,l));
        state_mean_est = [lf_meanfun_est(:,l) repmat(hf_mean_est(l),T,1)];
        mdiff = emiss-state_mean_est;
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
        hmm.state(l).hf_mean = hf_mean_est(l);
        hmm.state(l).lf_meanfun = lf_meanfun_est(:,l);
    end

    fprintf('iteration %i log posterior = %f \n',it,lpost);
    %         jmm_print_hmm_params(hmm);
    %     pause

end