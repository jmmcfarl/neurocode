function [hmm, LP, gamma] = hmm_uds_train(hmm,emiss)

addpath(strcat(hmm.drive_letter,':\Code\maphmmbox\'))

tolerance = 1e-5;
cur_log_post = -Inf;
delta_log_post = Inf;

%%
K = hmm.K;
p = hmm.p;
T = hmm.T;
Fs = hmm.Fs;
if strcmp(hmm.meantype,'variable')
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
end

eps = 1e-100;

%set model parameters to their initialized values
Pi = hmm.Pi;
A = hmm.A;

LP = []; %this will store the log likelihoods of the observations given model parameters across iterations of EM

%%
%initialize forward and backward messages
alpha=zeros(T,K);
beta=zeros(T,K);
gamma=zeros(T,K);

it = 0;
while delta_log_post > tolerance
    
%% compute posterior over hidden states 
    %compute likelihood of observations for each model under current
    %parameters
    [B] = mvgauss_uds_obslike(emiss,hmm);
    B(B < eps) = eps; %prevent likelihoods from rounding to 0
    
    scale=zeros(T,1); %initialize rescaling parameters
    
    %compute rescaled forward messages
    alpha(1,:)=Pi.*B(1,:);
    scale(1)=sum(alpha(1,:));
    alpha(1,:)=alpha(1,:)/scale(1);
    for t=2:T
        alpha(t,:)=(alpha(t-1,:)*A).*B(t,:);
        scale(t)=sum(alpha(t,:)); 
        alpha(t,:)=alpha(t,:)/scale(t);
    end
    
    %compute rescaled backward messages
    beta(T,:)=ones(1,K)/scale(T);
    for t=T-1:-1:1
        beta(t,:)=(beta(t+1,:).*B(t+1,:))*(A')/scale(t);
    end
    
    %compute posteriors over hidden states
    gamma=(alpha.*beta);
    gamma=rdiv(gamma,rsum(gamma));
    
    %compute chi (posterior of two consecutive hidden states)
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
   
    
%% M STEP
    % transition matrix
    expchi=sum(chi,1);
    expchi=reshape(expchi,K,K);
    A=rdiv(expchi,sum(expchi,2)); %updated estimate of transition matrix is expectation over chi
    
    % initial state
    Pi=gamma(1,:)/sum(gamma(1,:));
 
    var_est = zeros(K,p,p);
    if strcmp(hmm.meantype,'fixed')
        fixedmean_est = zeros(K,p);
    elseif p > 1
        fixedmean_est = zeros(K,p-1); %if there is a nonstationary dimension the number of fixed dimensions is p-1
    end
    if strcmp(hmm.meantype,'variable')
        meanfun_est = zeros(T,K);
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
        meandiff = temp_meanfun(2,:)-temp_meanfun(1,:);
        avg_meandiff = nanmean(meandiff);
        temp_meanfun(1,isnan(temp_meanfun(1,:))) = temp_meanfun(2,isnan(temp_meanfun(1,:))) - avg_meandiff;
        temp_meanfun(2,isnan(temp_meanfun(2,:))) = temp_meanfun(1,isnan(temp_meanfun(2,:))) + avg_meandiff;
        temp_meanfun = [temp_meanfun(:,1) temp_meanfun temp_meanfun(:,end)]';
        meanfun_est = interp1(win_t,temp_meanfun,t_axis);
        for l = 1:K
            meanfun_est(:,l) = filtfilt(filt_b,filt_a,meanfun_est(:,l));
        end
        if p > 1
            for l = 1:K
                fixedmean_est(l,:) = sum(repmat(gamma(:,l),1,p-1).*emiss(:,2:end))/sum(gamma(:,l));
            end
        end
    else
        for l = 1:K
            fixedmean_est = sum(repmat(gamma(:,l),1,p).*emiss)/sum(gamma(:,l));
        end
    end
    
    for l = 1:K
        if strcmp(hmm.meantype,'variable') && p > 1
            state_mean_est = [meanfun_est(:,l) repmat(fixedmean_est(l,:),T,1)];
        elseif strcmp(hmm.meantype,'variable')
            state_mean_est = meanfun_est(:,l);
        else
            state_mean_est = repmat(fixedmean_est(l,:),T,1);
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
        if p > 1 || strcmp(hmm.meantype,'fixedmean')
            hmm.state(l).fixedmean = fixedmean_est(l,:);
        end
        if strcmp(hmm.meantype,'variable')
            hmm.state(l).meanfun = meanfun_est(:,l);
        end
    end
    it = it + 1;
    fprintf('iteration %i log posterior = %f \n',it,loglik);
    %         jmm_print_hmm_params(hmm);
    %     pause
    
end

%store the observation likelihood under the ML parameters
hmm.loglik = LP;
