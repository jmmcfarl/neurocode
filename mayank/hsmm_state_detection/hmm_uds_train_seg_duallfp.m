function [hmm, gamma] = hmm_uds_train_seg_duallfp(hmm,emiss)

%modified from hmm_uds_train_seg (originally from MAPHMM Box)

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
    p_thresh = 0.05;
end

%these are the indices for the segments of data with UDS
UDS_samples = [];
for i = 1:hmm.Nsegs
    UDS_seg_inds{i} = hmm.UDS_segs(i,1):hmm.UDS_segs(i,2);
    UDS_samples = [UDS_samples hmm.UDS_segs(i,1):hmm.UDS_segs(i,2)];
end


eps = 1e-200;

%set model parameters to their initialized values
Pi = hmm.Pi;
A = hmm.A;

LP = []; %this will store the log likelihoods of the observations given model parameters across iterations of EM

%%

it = 0;
while delta_log_post > tolerance

    %% compute posterior over hidden states
    %compute likelihood of observations for each model under current
    %parameters
   
    %loop over UDS segments
    for i = 1:hmm.Nsegs
        B = mvgauss_uds_obslike_robust_duallfp(emiss(UDS_seg_inds{i},:),hmm,i);      
%         
%         %find instances where the emissions model fails and correct the
%         %likelihood
%         false_ups = find(emiss(UDS_seg_inds{i},:) < hmm.state(1).meanfun{i} & ...
%             B(:,1) < B(:,2));
%         B(false_ups,1) = 1;
%         false_downs = find(emiss(UDS_seg_inds{i},:) > hmm.state(2).meanfun{i} & ...
%             B(:,1) > B(:,2));
%         B(false_downs,2) = 1;
%         hmm.model_failure{i} = [false_ups; false_downs];
         
            
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
        gamma{i}=rdiv(gamma{i},rsum(gamma{i}));

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
    A=rdiv(expchi,sum(expchi,2)); %updated estimate of transition matrix is expectation over chi

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
        meandiff2 = [];
        for i = 1:hmm.Nsegs
            curT = length(UDS_seg_inds{i});
            total_dur = curT/Fs;
            t_axis{i} = (1:curT)/Fs;
            numWins = floor((total_dur-windowSize)/windowSlide);
            win_t{i} = (0:numWins-1)*windowSlide+windowSize/2;
            win_t{i} = [0 win_t{i} max(t_axis{i})];
            cur_tot_prob = nan(2,numWins);
            temp_meanfun{i} = zeros(hmm.K,numWins);
            temp_meanfun2{i} = zeros(hmm.K,numWins);
            for l=1:K
                for w = 1:numWins
                    begT = round((w-1)*windowSlide*Fs)+1;
                    endT = begT + round(windowSize*Fs);
                    temp_meanfun{i}(l,w) = sum(gamma{i}(begT:endT,l).*emiss(UDS_seg_inds{i}(begT:endT),1))/sum(gamma{i}(begT:endT,l));
                    temp_meanfun2{i}(l,w) = sum(gamma{i}(begT:endT,l).*emiss(UDS_seg_inds{i}(begT:endT),2))/sum(gamma{i}(begT:endT,l));
                    cur_tot_prob(l,w) = sum(gamma{i}(begT:endT,l))/(windowSize*Fs);
                end
                temp_meanfun{i}(l,cur_tot_prob(l,:) < p_thresh) = nan;
                temp_meanfun2{i}(l,cur_tot_prob(l,:) < p_thresh) = nan;
            end
            meandiff = [meandiff temp_meanfun{i}(2,cur_tot_prob(l,:) >= p_thresh)-...
                temp_meanfun{i}(1,cur_tot_prob(l,:) >= p_thresh)];
            meandiff2 = [meandiff2 temp_meanfun2{i}(2,cur_tot_prob(l,:) >= p_thresh)-...
                temp_meanfun2{i}(1,cur_tot_prob(l,:) >= p_thresh)];
        end
        avg_meandiff = nanmean(meandiff);
        avg_meandiff2 = nanmean(meandiff2);
        for i = 1:hmm.Nsegs
            temp_meanfun{i}(1,isnan(temp_meanfun{i}(1,:))) = temp_meanfun{i}(2,isnan(temp_meanfun{i}(1,:))) - avg_meandiff;
            temp_meanfun{i}(2,isnan(temp_meanfun{i}(2,:))) = temp_meanfun{i}(1,isnan(temp_meanfun{i}(2,:))) + avg_meandiff;
            temp_meanfun{i} = [temp_meanfun{i}(:,1) temp_meanfun{i} temp_meanfun{i}(:,end)]';
            meanfun_est{i} = interp1(win_t{i},temp_meanfun{i},t_axis{i});
            temp_meanfun2{i}(1,isnan(temp_meanfun2{i}(1,:))) = temp_meanfun2{i}(2,isnan(temp_meanfun2{i}(1,:))) - avg_meandiff2;
            temp_meanfun2{i}(2,isnan(temp_meanfun2{i}(2,:))) = temp_meanfun2{i}(1,isnan(temp_meanfun2{i}(2,:))) + avg_meandiff2;
            temp_meanfun2{i} = [temp_meanfun2{i}(:,1) temp_meanfun2{i} temp_meanfun2{i}(:,end)]';
            meanfun_est2{i} = interp1(win_t{i},temp_meanfun2{i},t_axis{i});
        end
    else
        for l = 1:K
            o_sum = 0;
            g_sum = 0;
            for i = 1:hmm.Nsegs
               o_sum = o_sum + sum(repmat(gamma{i}(:,l),1,p).*emiss(UDS_seg_inds{i},:));
               g_sum = g_sum + sum(gamma{i}(:,l));
            end
            fixedmean_est = o_sum/g_sum;
        end
    end

    for l = 1:K
        for i = 1:hmm.Nsegs
            curT = length(UDS_seg_inds{i});
            state_mean_est = [meanfun_est{i}(:,l) meanfun_est2{i}(:,l)];
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
        if strcmp(hmm.meantype,'variable')
            for i = 1:hmm.Nsegs
                hmm.state(l).meanfun{i} = meanfun_est{i}(:,l);
                hmm.state(l).meanfun2{i} = meanfun_est2{i}(:,l);
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
