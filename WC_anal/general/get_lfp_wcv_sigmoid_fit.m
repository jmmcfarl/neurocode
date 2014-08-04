% Do a sigmoid fit to the up or down transitions in the data and save the numbers:
% rlid = ID in LFP where transition occurred. rltime = synchronized time of transition.
% rlamp = amplitude of transition, rlshift = minimum value of transition.
% rltau = time constant of the transition. rlerror = goodness of
% sigmoid fit.

function [rlid,rltime,rlamp,rlshift,rltau,rlerror] = get_lfp_wcv_sigmoid_fit(trans_ind,data,data_slope,synct,t_axis)

dt = mean(diff(synct));

% find the true times of state transitions.
rlid = trans_ind;
rltime = trans_ind;
rlamp = trans_ind;
rltau = trans_ind;
rlshift = trans_ind;
rlerror = trans_ind;

%time axis for visualizing results
t2 = -10*2016:10*2016;

for k=1:length(trans_ind)

    cur_trans = trans_ind(k);
    
    if cur_trans > max(t_axis) & cur_trans < length(data) - max(t_axis)
 
        cur_taxis = data(t_axis+cur_trans);

        mincur_taxis = min(cur_taxis);
        cur_taxis = cur_taxis-mincur_taxis;

        init_guess = [mean(cur_taxis) 0 100 mean(cur_taxis)-std(cur_taxis)];
        [betafit, resid] = nlinfit(t_axis,cur_taxis,@my_sigmoid,init_guess);
        if abs(betafit(2)) > max(t_axis)
            m = 0;
        else
            m = round(betafit(2));
        end
        rlid(k)=cur_trans+m;
        rltime(k)=synct(cur_trans+m);
        rlamp(k) = betafit(1);
        rltau(k) = betafit(3);
        rlshift(k) = betafit(4)+mincur_taxis;
        rlerror(k) = sqrt(mean(resid.^2));
        %     net_betafit(k,:) = betafit;

        %jmm fit vis
        if cur_trans > max(t2) & cur_trans < length(data) - max(t2)
            
            cur_taxis_2 = data(t2+cur_trans);
            cur_taxis_2 = cur_taxis_2 - mincur_taxis;

            plot(t2,cur_taxis_2,'linewidth',2)
            hold on
            plot(t_axis,my_sigmoid(betafit,t_axis),'r','linewidth',2)
            pause
            clf
        end

    end

end

rltau = rltau*dt;
