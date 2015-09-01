function [hmm, gamma] = hmm_uds_train_seg(hmm,emiss)

%ADAPTED FROM CODE IN THE HMMBOX, version 3.3, Iead Rezek, Oxford University, February 2001

addpath(strcat(hmm.drive_letter,':\Code\maphmmbox\'))

tolerance = 1e-5; %convergence tolerance for the likelihood (or posterior)
cur_log_post = -Inf; %will store the current iterations likelihood
delta_log_post = Inf; %will store the current change in likelihood from the previous iteration

%%
K = hmm.K; %number of hidden states
p = hmm.p; %dimensionality of observations
T = hmm.T; %number of time samples
Fs = hmm.Fs; %sampling frequency
if strcmp(hmm.meantype,'variable')
    windowSize = hmm.windowSize;
    windowSlide = hmm.windowSlide;
    p_thresh = 0.05; %minimum probability of occupying a given state in order to use direct estimation of that state's mean function within a given time window
end

%these are the indices for the segments of data with UDS
UDS_samples = [];
for i = 1:hmm.Nsegs
    UDS_seg_inds{i} = hmm.UDS_segs(i,1):hmm.UDS_segs(i,2);
    UDS_samples = [UDS_samples hmm.UDS_segs(i,1):hmm.UDS_segs(i,2)];
end

%set model parameters to their initialized values
Pi = hmm.Pi;
A = hmm.A;

LP = []; %this will store the log likelihoods of the observations given model parameters across iterations of EM

%%

it = 0; %iteration number
while delta_log_post > tolerance

    %% compute posterior over hidden states
    %compute likelihood of observations for each model under current
    %parameters
   
    %loop over UDS segments
    for i = 1:hmm.Nsegs
        B = mvgauss_uds_obslike_robust(emiss(UDS_seg_inds{i},:),hmm,i);               
            
        %initialize forward and backward messages
        curT = length(UDS_seg_inds{i});
        alpha{i}=zeros(curT,K);
        beta{i}=zeros(curT,K);
        gamma{i}=zeros(curT,K);

        scale=zeros(curT,1); %initialize rescaling parameters

        %compute rescaled forward messages
        alpha{i}(1,:)=Pi(i,:).*B(1,:);
        scale(1)=sum(alpha{i}(1,:));
        alpha{i}(1,:)=alpha{i}(1,:)/scale(1);
        for t=2:curT
            alpha{i}(t,:)=(alpha{i}(t-1,:)*A).*B(t,:);
            scale(t)=sum(alpha{i}(t,:));
            alpha{i}(t,:)=alpha{i}(t,:)/scale(t);
        end

        %compute rescaled backward messages
        beta{i}(curT,:)=ones(1,K)/scale(curT);
        for t=curT-1:-1:1
            beta{i}(t,:)=(beta{i}(t+1,:).*B(t+1,:))*(A')/scale(t);
        end

        %compute posteriors over hidden states
        gamma{i}=(alpha{i}.*beta{i});
        gamma{i} = gamma{i}./repmat(sum(gamma{i},2),1,K);
%         gamma{i}=rdiv(gamma{i},rsum(gamma{i}));

        %compute chi (posterior of two consecutive hidden states)
        chi{i}=zeros(curT-1,K*K);
        for t=1:curT-1
            temp=A.*(alpha{i}(t,:)' * (beta{i}(t+1,:).*B(t+1,:)));
            chi{i}(t,:)=temp(:)'/sum(temp(:));
        end

        lscale = log(scale);
        loglik(i) = sum(lscale); %rabiner eq(103), scale is defined as inverse here though
    end

    loglik = sum(loglik); %likelihood over all UDS segments
    LP=[LP; loglik];
    delta_log_post = (loglik - cur_log_post)/abs(loglik);
    cur_log_post = loglik;

%% M STEP
    % transition matrix
    expchi = zeros(1,K^2);
    for i = 1:hmm.Nsegs
        expchi = expchi + sum(chi{i},1);
    end
    expchi=reshape(expchi,K,K);
    A = expchi./repmat(sum(expchi,2),1,K);
%     A=rdiv(expchi,sum(expchi,2)); %updated estimate of transition matrix is expectation over chi

    % initial state
    for i = 1:hmm.Nsegs
        Pi(i,:)=gamma{i}(1,:)/sum(gamma{i}(1,:));
    end

    var_est = zeros(K,p,p);
    if strcmp(hmm.meantype,'fixed')
        fixedmean_est = zeros(K,p);
    elseif p > 1
        fixedmean_est = zeros(K,p-1); %if there is a nonstationary dimension the number of fixed dimensions is p-1
    end
    if strcmp(hmm.meantype,'variable')
        meandiff = [];
        for i = 1:hmm.Nsegs
            curT = length(UDS_seg_inds{i});
            total_dur = curT/Fs;
            t_axis{i} = (1:curT)/Fs;
            numWins = floor((total_dur-windowSize)/windowSlide);
            win_t{i} = (0:numWins-1)*windowSlide+windowSize/2;
            win_t{i} = [0 win_t{i} max(t_axis{i})];
            cur_tot_prob = nan(2,numWins);
            temp_meanfun{i} = zeros(hmm.K,numWins);
            for l=1:K
                for w = 1:numWins
                    begT = round((w-1)*windowSlide*Fs)+1;
                    endT = begT + round(windowSize*Fs);
                    %restimate state means within the current data window
                    temp_meanfun{i}(l,w) = sum(gamma{i}(begT:endT,l).*emiss(UDS_seg_inds{i}(begT:endT),1))/sum(gamma{i}(begT:endT,l));
                    cur_tot_prob(l,w) = sum(gamma{i}(begT:endT,l))/(windowSize*Fs); %average probability of occupying state witin window
                end
                temp_meanfun{i}(l,cur_tot_prob(l,:) < p_thresh) = nan; %set any estimates of meanfunction with insufficient posterior probability to nan
            end
            meandiff = [meandiff temp_meanfun{i}(2,cur_tot_prob(l,:) >= p_thresh)-...
                temp_meanfun{i}(1,cur_tot_prob(l,:) >= p_thresh)]; %compute a vector of difference of state meanfunctions
        end 
        avg_meandiff = nanmean(meandiff); %average difference of mean functions
        
        %any instances where there was insufficient posterior probability
        %for one state are given by estimates extrapolated from the other
        %state's mean
        for i = 1:hmm.Nsegs
            temp_meanfun{i}(1,isnan(temp_meanfun{i}(1,:))) = temp_meanfun{i}(2,isnan(temp_meanfun{i}(1,:))) - avg_meandiff;
            temp_meanfun{i}(2,isnan(temp_meanfun{i}(2,:))) = temp_meanfun{i}(1,isnan(temp_meanfun{i}(2,:))) + avg_meandiff;
            temp_meanfun{i} = [temp_meanfun{i}(:,1) temp_meanfun{i} temp_meanfun{i}(:,end)]';
            meanfun_est{i} = interp1(win_t{i},temp_meanfun{i},t_axis{i});
        end

        if p > 1 %for multi-dimensional observations, use fixed means for other dimensions
            for l = 1:K
                o_sum = 0;
                g_sum = 0;
                for i = 1:hmm.Nsegs
                    o_sum = o_sum + sum(repmat(gamma{i}(:,l),1,p-1).*emiss(UDS_seg_inds{i},2:end));
                    g_sum = g_sum + sum(gamma{i}(:,l));
                end
                fixedmean_est(l,:) = o_sum/g_sum;
            end
        end
    else %for fixed state means
        for l = 1:K
            o_sum = 0;
            g_sum = 0;
            for i = 1:hmm.Nsegs
               o_sum = o_sum + sum(repmat(gamma{i}(:,l),1,p).*emiss(UDS_seg_inds{i},:));
               g_sum = g_sum + sum(gamma{i}(:,l));
            end
            fixedmean_est(l) = o_sum/g_sum;
        end
    end

    for l = 1:K %restimate state covariance matrices
        for i = 1:hmm.Nsegs
            curT = length(UDS_seg_inds{i});
            if strcmp(hmm.meantype,'variable') && p > 1
                state_mean_est = [meanfun_est{i}(:,l) repmat(fixedmean_est(l,:),curT,1)];
            elseif strcmp(hmm.meantype,'variable')
                state_mean_est = meanfun_est{i}(:,l);
            else
                state_mean_est = repmat(fixedmean_est(l,:),curT,1);
            end
            mdiff{i} = emiss(UDS_seg_inds{i},:)-state_mean_est;
        end
        o_sum = 0;
        g_sum = 0;
        for i = 1:hmm.Nsegs
            o_sum = o_sum + (repmat(gamma{i}(:,l),1,p).*mdiff{i})'*mdiff{i};
            g_sum = g_sum + sum(gamma{i}(:,l));
        end
        var_est(l,:,:) = o_sum/g_sum;
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
            for i = 1:hmm.Nsegs
                hmm.state(l).meanfun{i} = meanfun_est{i}(:,l);
            end
        end
    end
    it = it + 1;
    fprintf('iteration %i log posterior = %f \n',it,loglik);
    %         jmm_print_hmm_params(hmm);
    %     pause

end

%store the observation likelihood under the ML parameters
hmm.hmm_LP = LP;
