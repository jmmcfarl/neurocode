% Do a sigmoid fit to the up or down transitions in the data and save the numbers:
% rlid = ID in LFP where transition occurred. rltime = synchronized time of transition.
% rlamp = amplitude of transition, rlshift = minimum value of transition.
% rltau = time constant of the transition. rlerror = goodness of
% sigmoid fit.

function [rlid,rltime,rlamp,rlshift,rltau,rlerror,t_10] = get_lfp_wcv_sigmoid_fit_dt(trans_ind,data,data_slope,synct)

dt = mean(diff(synct));

% find the true times of state transitions.
rlid = trans_ind;
rltime = trans_ind;
rlamp = trans_ind;
rltau = trans_ind;
rlshift = trans_ind;
rlerror = trans_ind;
t_10 = trans_ind;

%time axis for visualizing results
t2 = -10*2016:10*2016;

max_t = 2016;

for k=1:length(trans_ind)

    cur_trans = trans_ind(k);

    if cur_trans > max_t & cur_trans < length(data) - max_t

        %find preceeding point where down state ended
        next_up_end = find(data_slope(cur_trans:end) < -0.5,1,'first');

        lookback_t = next_up_end;

        if lookback_t > max_t
            t_axis = -max_t:max_t;
        else
            t_axis = -lookback_t:lookback_t;
        end

        cur_data = data(t_axis+cur_trans);

        min_cur_data = min(cur_data);
        cur_data = cur_data-min_cur_data;

        init_guess = [mean(cur_data) 0 100 mean(cur_data)-std(cur_data)];
        [betafit, resid] = nlinfit(t_axis,cur_data,@my_sigmoid,init_guess);
        
        if max(isnan(betafit)) == 0
            
            if abs(betafit(2)) > max(t_axis)
                m = 0;
            else
                m = round(betafit(2));
            end
            rlid(k)=cur_trans+m;
            rltime(k)=synct(cur_trans+m);
            rlamp(k) = betafit(1);
            rltau(k) = betafit(3);
            rlshift(k) = betafit(4)+min_cur_data;
            rlerror(k) = sqrt(mean(resid.^2));
                        fit_eval = my_sigmoid(betafit,t_axis);
            thresh_10 = min(fit_eval)+range(fit_eval)*0.1;
            temp = find(fit_eval > thresh_10,1,'first')+cur_trans-round(length(fit_eval)/2);
            if ~isempty(temp)
                t_10(k) = temp;
            else
                t_10(k) = nan;
            end

        else
            rlid(k) = nan;
            rltime(k) = nan;
            rlamp(k) = nan;
            rltau(k) = nan;
            rlshift(k) = nan;
            rlerror(k) = nan;
            t_10(k) = nan;
        end

        %     net_betafit(k,:) = betafit;

                %jmm fit vis
%                 if cur_trans > max(t2) & cur_trans < length(data) - max(t2)
%         
%                     cur_taxis_2 = data(t2+cur_trans);
%                     cur_taxis_2 = cur_taxis_2 - min_cur_data;
%         
%                     plot(t2,cur_taxis_2,'linewidth',2)
%                     hold on
%                     plot(t_axis,my_sigmoid(betafit,t_axis),'r','linewidth',2)
%                     pause
%                     clf
%                 end


    else
            rlid(k) = nan;
            rltime(k) = nan;
            rlamp(k) = nan;
            rltau(k) = nan;
            rlshift(k) = nan;
            rlerror(k) = nan;
            t_10(k) = nan;
    end        
end

rltau = rltau*dt;
