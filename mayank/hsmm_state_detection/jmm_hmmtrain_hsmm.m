function [hmm, LP] = jmm_hmmtrain_hsmm(hmm,emiss)

addpath('C:\Code\maphmmbox\')
tolerance = 1e-5;
cur_log_post = -Inf;
delta_log_post = Inf;

%%
K = hmm.K;
p = hmm.p;
T = hmm.T;
Fs = hmm.Fs;
Fs  = hmm.Fs;
windowSize = hmm.windowSize;
windowSlide = hmm.windowSlide;
total_dur = T/Fs;
t_axis = (1:T)/Fs;
numWins = floor((total_dur-windowSize)/windowSlide);
win_t = (0:numWins-1)*windowSlide+windowSize/2;
win_t = [0 win_t max(t_axis)];
niqf = Fs/2;
if hmm.hcf > 0
    [filt_b,filt_a] = butter(2,hmm.hcf/niqf,'low');
end
eps = 1e-100;

%%

%initialize forward and backward messages
alpha=zeros(T,K);
beta=zeros(T,K);
gamma=zeros(T,K);

%set model parameters to their initialized values
Pi = hmm.Pi;
A = hmm.A;

LP = [];

% for it = 1:num_it
it = 0;
while delta_log_post > tolerance
    
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
    loglik = sum(lscale); %rabiner eq(103), scale is defined as inverse here though
    LP=[LP; loglik];
    delta_log_post = (loglik - cur_log_post)/abs(loglik);
    cur_log_post = loglik;
   
    
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
    
    lf_meanfun_est = zeros(T,K);
    hf_mean_est = zeros(K,p-1);
    var_est = zeros(K,p,p);
    if strcmp(hmm.meantype,'variable')
        tot_prob = nan(2,numWins);
        p_thresh = 0.05;
        temp_meanfun = nan(2,numWins);
       for l=1:K
            for w = 1:numWins
                begT = round((w-1)*windowSlide*Fs)+1;
                endT = begT + round(windowSize*Fs);
                temp_meanfun(l,w) = sum(gamma(begT:endT,l).*emiss(begT:endT,1))/sum(gamma(begT:endT,l));
                tot_prob(l,w) = sum(gamma(begT:endT,l))/(windowSize*Fs);
            end
            temp_meanfun(l,tot_prob(l,:) < p_thresh) = nan;
        end
        lf_meandiff = temp_meanfun(2,:)-temp_meanfun(1,:);
        mlf_meandiff = nanmean(lf_meandiff);
        temp_meanfun(1,isnan(temp_meanfun(1,:))) = temp_meanfun(2,isnan(temp_meanfun(1,:))) - mlf_meandiff;
        temp_meanfun(2,isnan(temp_meanfun(2,:))) = temp_meanfun(1,isnan(temp_meanfun(2,:))) + mlf_meandiff;
        temp_meanfun = [temp_meanfun(:,1) temp_meanfun temp_meanfun(:,end)]';
        lf_meanfun_est = interp1(win_t,temp_meanfun,t_axis);
        for l = 1:K
            lf_meanfun_est(:,l) = filtfilt(filt_b,filt_a,lf_meanfun_est(:,l));
        end
    else
        for l = 1:K
            lf_mean = sum(gamma(:,l).*emiss(:,1))/sum(gamma(:,l));
            lf_meanfun_est(:,l) = repmat(lf_mean,T,1);
        end
    end
    
    for l = 1:K
        if p > 1
            hf_mean_est(l,:) = sum(repmat(gamma(:,l),1,p-1).*emiss(:,2:end))/sum(gamma(:,l));
            state_mean_est = [lf_meanfun_est(:,l) repmat(hf_mean_est(l,:),T,1)];
        else
            state_mean_est = lf_meanfun_est(:,l);
        end
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
        if p > 1
            hmm.state(l).hf_mean = hf_mean_est(l,:);
        end
        hmm.state(l).lf_meanfun = lf_meanfun_est(:,l);
    end
    it = it + 1;
    fprintf('iteration %i log posterior = %f \n',it,loglik);
    %         jmm_print_hmm_params(hmm);
    %     pause
    
end

hmm.loglik = LP;
