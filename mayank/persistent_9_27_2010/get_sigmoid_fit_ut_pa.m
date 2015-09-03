% Do a sigmoid fit to the up or down transitions in the data and save the numbers:
% rlid = ID in LFP where transition occurred. rltime = synchronized time of transition.
% rlamp = amplitude of transition, rlshift = minimum value of transition.
% rltau = time constant of the transition. rlerror = goodness of
% sigmoid fit.

function [rlamp,rlshift,rltau,rsquared,t_50,t_90,t_10,fit_data] = ...
    get_sigmoid_fit_ut_pa(up_trans,down_trans,data,time,Fs)

addpath('F:\Code\WC_anal\general\')

data = data(:);

n_ups = length(up_trans);

rlamp = nan(n_ups,1);
rltau = nan(n_ups,1);
rlshift = nan(n_ups,1);
rsquared = nan(n_ups,1);
t_10 = nan(n_ups,1);
t_90 = nan(n_ups,1);
t_50 = nan(n_ups,1);

target_window = round(0.5*Fs);
state_buffer = round(0.15*Fs);

fit_data = nan(size(data));

for i=1:length(up_trans)
    fprintf('fitting state %d of %d\n',i,length(up_trans));
    cur_trans = up_trans(i);
    
    %if this is not the first up transition, set the beggining edge of the
    %fitting window to be the previous down transition plus a buffer to
    %bypass the transition region
    if i > 1
        prev_down = down_trans(i-1)+state_buffer;
    else
        prev_down = 1;
    end
    
    next_down = down_trans(i);
    %if the current up state and previous down states last for the minimum duration (necessary to fit a
    %sigmoid)
    if (next_down - cur_trans) > state_buffer && (cur_trans - prev_down) > state_buffer
        
        %set the range of the data-fitting window
        start_p = max([prev_down,cur_trans-target_window,1]);
        end_p = min([length(time),cur_trans+target_window,next_down-state_buffer]);

        t_axis = [start_p:end_p]';
        cur_data = data(t_axis);
        cur_mean = mean(cur_data);
        cur_data = cur_data-cur_mean;
        t_axis = (t_axis-cur_trans);

        init_guess = [2*std(cur_data) 0 20 0];
        [betafit, resid] = nlinfit(t_axis,cur_data,@my_sigmoid,init_guess);
        
        %make sure it didn't fit a small segment of the sigmoid
        if abs(betafit(1)) > 100 | abs(betafit(2)) > 100
            betafit(:) = nan;
        end
    else
        betafit = nan;
    end

    if max(isnan(betafit)) == 0

        if abs(betafit(2)) > max(t_axis)
            m = 0;
        else
            m = round(betafit(2));
        end
        rlamp(i) = betafit(1);
        rltau(i) = betafit(3);
        rlshift(i) = betafit(4)+cur_mean;
        t_50(i) = cur_trans+m;
        t_10(i) = round(betafit(2) - betafit(3)*log(9));
        t_90(i) = round(betafit(2) + betafit(3)*log(9));
        fit_eval = my_sigmoid(betafit,t_axis);
        fit_data(start_p:end_p) = fit_eval+cur_mean;
        
        sp = find(t_axis > t_10(i),1,'first');
        ep = find(t_axis > t_90(i),1,'first');
        if isempty(ep)
            ep = length(resid);
        end
        t_10(i) = t_10(i) + cur_trans;
        t_90(i) = t_90(i) + cur_trans;
        ssr = sum(resid(sp:ep).^2);
        sst = sum(cur_data(sp:ep).^2);
        rsquared(i) = 1-ssr/sst;        
    else
        t_50(i) = nan;
        rlamp(i) = nan;
        rltau(i) = nan;
        rlshift(i) = nan;
        rsquared(i) = nan;
        t_90(i) = nan;
        t_10(i) = nan;
    end


% %     jmm fit vis
% %time axis for visualizing results
%     t2 = -3*Fs:3*Fs;
%     if cur_trans > max(t2) && cur_trans < length(data) - max(t2)
% 
%         cur_data_2 = data(t2+cur_trans);
%         cur_data_2 = cur_data_2 - cur_mean;
%         plot(t2,cur_data_2,'linewidth',2), hold on
%         if max(isnan(betafit)) == 0
%             fit_fun = my_sigmoid(betafit,t_axis);
%             plot(t_axis,fit_fun,'r','linewidth',2)
%             plot(t_axis,(cur_data-fit_fun)-2,'k')
%             plot(t_90(i)-cur_trans,0,'ro')
%             plot(t_10(i)-cur_trans,0,'go')
%         end
%         pause
%         clf
%     end




end

%convert to secs
rltau = rltau/Fs;
