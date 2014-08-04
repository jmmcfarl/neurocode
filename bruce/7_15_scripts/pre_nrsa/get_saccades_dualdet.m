function [hf_start_inds,hf_stop_inds,sac_data] = get_saccades_dualdet(rEyeXY,lEyeXY,eye_ts,in_blink)
    
%%
sac_eyespeed = 10;
thresh_eyespeed = 5;
min_intersac_interval = 0.1;

eye_dt = median(diff(eye_ts));
eye_fs = 1/eye_dt;

%%
leye_vel = [0 0; diff(lEyeXY)]/eye_dt;
reye_vel = [0 0; diff(rEyeXY)]/eye_dt;
leye_speed = sqrt(leye_vel(:,1).^2+leye_vel(:,2).^2);
reye_speed = sqrt(reye_vel(:,1).^2+reye_vel(:,2).^2);
% leye_speed(in_blink==1) = 0;
% reye_speed(in_blink==1) = 0;

eye_speed = max(reye_speed,leye_speed);
%%
in_sac = zeros(size(eye_ts));
if max(eye_speed) > sac_eyespeed
    
    proposed_sac_starts = 1+find(eye_speed(1:end-1) < sac_eyespeed & eye_speed(2:end) >= sac_eyespeed);
    sac_start_inds = proposed_sac_starts;
    sac_stop_inds = nan(size(proposed_sac_starts));
    for i = 1:length(proposed_sac_starts)
        temp = find(eye_speed(1:proposed_sac_starts(i)) < thresh_eyespeed,1,'last');
        if ~isempty(temp)
            sac_start_inds(i) = temp;
        end
        temp = find(eye_speed(proposed_sac_starts(i)+1:end) < thresh_eyespeed,1,'first');
        if i < length(proposed_sac_starts)
            if temp > proposed_sac_starts(i+1)-proposed_sac_starts(i)
                temp = nan;
            end
        end
        if ~isempty(temp)
            sac_stop_inds(i) = proposed_sac_starts(i)+temp-1;
        end
    end
    bad_sacs = find(isnan(sac_stop_inds));
    sac_start_inds(bad_sacs) = [];
    sac_stop_inds(bad_sacs) = [];
    
    
    blink_sacs = [];
    for i = 1:length(sac_start_inds)
        cur_inds = sac_start_inds(i):sac_stop_inds(i);
        if any(in_blink(cur_inds)==1)
            blink_sacs = [blink_sacs i];
        end
    end
    sac_start_inds(blink_sacs) = [];
    sac_stop_inds(blink_sacs) = [];
    
    sac_start_times = eye_ts(sac_start_inds);
    sac_stop_times = eye_ts(sac_stop_inds);
    
    isi = diff(sac_start_times);
    too_short_int = find(isi < min_intersac_interval);
    rid_set = unique([too_short_int too_short_int+1]);
    hf_start_inds = sac_start_inds(rid_set);
    hf_stop_inds = sac_stop_inds(rid_set);
    
    sac_start_inds(rid_set) = [];
    sac_stop_inds(rid_set) = [];
    
    sac_start_times = eye_ts(sac_start_inds);
    sac_stop_times = eye_ts(sac_stop_inds);
    
    blink_start_inds = find(in_blink(2:end)==1 & in_blink(1:end-1)==0)+1;
    all_event_start_inds = unique([sac_start_inds; hf_start_inds; blink_start_inds']);
    
    n_sacs = length(sac_start_inds);
    for i = 1:n_sacs
        in_sac(sac_start_inds(i):sac_stop_inds(i)) = 1;
    end
    
    %%
    options.Display = 'off';
    min_rebound_amp = 2.5e-3;
    min_rbounce_amp = 0.25;
    min_bounce_amp = -2.5e-3;
    lb = [-Inf 0 0.05];
    ub = [Inf Inf 100];
    
    exp_fun = @(p,x) p(1) + p(2)*exp(-p(3)*x);
    
    leye_deriv = [0 0; diff(lEyeXY)];
    reye_deriv = [0 0; diff(rEyeXY)];
    
    sm_fac = 3;
    leye_deriv(:,1) = smooth(leye_deriv(:,1),sm_fac);
    leye_deriv(:,2) = smooth(leye_deriv(:,2),sm_fac);
    reye_deriv(:,1) = smooth(reye_deriv(:,1),sm_fac);
    reye_deriv(:,2) = smooth(reye_deriv(:,2),sm_fac);
    
    max_lderiv = max(abs(leye_deriv),[],2);
    max_rderiv = max(abs(reye_deriv),[],2);
    
    left_adj_sac_end = nan(length(sac_start_inds),2);
    right_adj_sac_end = nan(length(sac_start_inds),2);
    
    %%
    warning off
    for i = 1:n_sacs-1
        fprintf('Estimating transients Sac %d of %d\n',i,n_sacs');
        %for horizontal component
        within_sac_inds = sac_start_inds(i):sac_stop_inds(i);
        [pk,pkloc] = max(abs(leye_deriv(within_sac_inds,1)));
        pk_amp = leye_deriv(sac_start_inds(i)-1+pkloc,1);
        next_event = all_event_start_inds(find(all_event_start_inds > sac_start_inds(i),1,'first'));
        until_nextsac_inds = (sac_start_inds(i)+pkloc-1):next_event;
        
        cur_snippet = leye_deriv(until_nextsac_inds,1);
        cur_snippet = -1*cur_snippet*sign(pk_amp);
        
        is_trans_h = 1;
%         snip_range = 1:100;
        [rebound_amp,rebound_peak] = findpeaks(cur_snippet,'minpeakheight',min_rebound_amp,'npeaks',1);
%         rebound_peak = find(cur_snippet(1:end-1) > cur_snippet(2:end),1,'first');
%         if cur_snippet(rebound_peak) > min_rebound_amp
        if ~isempty(rebound_peak)
            decay_signal = cur_snippet(rebound_peak:end);
        else
            is_trans_h = 0;
            decay_signal = -cur_snippet;
            rebound_peak = 1;
        end
        
        if ~isempty(decay_signal)
            
            %check for biphasic 'ringing' response
            is_ring_h = 0;
            first_bounce = find(decay_signal(1:end-1) < decay_signal(2:end),1,'first');
            if ~isempty(first_bounce)
                bounce_amp = decay_signal(first_bounce);
                rel_bounce_amp = abs(decay_signal(first_bounce)/decay_signal(1));
                if rel_bounce_amp > min_rbounce_amp & bounce_amp < min_bounce_amp
                    rebound_peak = first_bounce - 1 + rebound_peak;
                    decay_signal = -decay_signal(first_bounce:end);
                    is_ring_h = 1;
                end
            end
            if length(decay_signal) > 10
                xx = (0:(length(decay_signal)-1))';
                p0 = [median(decay_signal) decay_signal(1) 0.1];
                beta = lsqcurvefit(exp_fun,p0,xx,decay_signal,lb,ub,options);
                if beta(2) > .0025
                    left_adj_sac_end(i,1) = round(sac_start_inds(i) + pkloc + rebound_peak -2 + 2.5/beta(3));
                elseif is_ring_h==1
                    left_adj_sac_end(i,1) = round(sac_start_inds(i) + pkloc + rebound_peak-2);
                end
            end
        end
        
        
        %for vertical component
        [pk,pkloc] = max(abs(leye_deriv(within_sac_inds,2)));
        pk_amp = leye_deriv(sac_start_inds(i)-1+pkloc,2);
        next_event = all_event_start_inds(find(all_event_start_inds > sac_start_inds(i),1,'first'));
        until_nextsac_inds = (sac_start_inds(i)+pkloc-1):next_event;
        
        cur_snippet = leye_deriv(until_nextsac_inds,2);
        cur_snippet = -1*cur_snippet*sign(pk_amp);
        
        is_trans_v = 1;
        [rebound_amp,rebound_peak] = findpeaks(cur_snippet,'minpeakheight',min_rebound_amp,'npeaks',1);
%         rebound_peak = find(cur_snippet(1:end-1) > cur_snippet(2:end),1,'first');
%         if cur_snippet(rebound_peak) > min_rebound_amp
        if ~isempty(rebound_peak)
            decay_signal = cur_snippet(rebound_peak:end);
        else
            is_trans_v = 0;
            decay_signal = -cur_snippet;
            rebound_peak = 1;
        end
        
        if ~isempty(decay_signal)
            
            %check for biphasic 'ringing' response
            is_ring_v = 0;
            first_bounce = find(decay_signal(1:end-1) < decay_signal(2:end),1,'first');
            if ~isempty(first_bounce)
                bounce_amp = decay_signal(first_bounce);
                rel_bounce_amp = abs(decay_signal(first_bounce)/decay_signal(1));
                if rel_bounce_amp > min_rbounce_amp & bounce_amp < min_bounce_amp
                    rebound_peak = first_bounce - 1 + rebound_peak;
                    decay_signal = -decay_signal(first_bounce:end);
                    is_ring_v = 1;
                end
            end
            
            if length(decay_signal) > 10
                xx = (0:(length(decay_signal)-1))';
                p0 = [median(decay_signal) decay_signal(1) 0.1];
                beta = lsqcurvefit(exp_fun,p0,xx,decay_signal,lb,ub,options);
                if beta(2) > .0025
                    left_adj_sac_end(i,2) = round(sac_start_inds(i) + pkloc + rebound_peak -2 + 2.5/beta(3));
                elseif is_ring_v==1
                    left_adj_sac_end(i,2) = round(sac_start_inds(i) + pkloc + rebound_peak-2);
                end
            end
        end
        %%
%         big_snip = (sac_start_inds(i)-20):next_event;
%         big_snip(big_snip < 0) = [];
%         subplot(2,1,1)
%         plot(eye_ts(big_snip),lEyeXY(big_snip,1)-mean(lEyeXY(big_snip,1)),'linewidth',2)
%         hold on
%         plot(eye_ts(big_snip),leye_deriv(big_snip,1)*10,'r')
%         plot(eye_ts(sac_start_inds(i)),lEyeXY(sac_start_inds(i),1)-mean(lEyeXY(big_snip,1)),'go','linewidth',4)
%         plot(eye_ts(sac_stop_inds(i)),lEyeXY(sac_stop_inds(i),1)-mean(lEyeXY(big_snip,1)),'ro','linewidth',4)
%         if ~isnan(left_adj_sac_end(i,1))
%             plot(eye_ts(left_adj_sac_end(i,1)),lEyeXY(left_adj_sac_end(i,1),1)-mean(lEyeXY(big_snip,1)),'ko','linewidth',4)
%         end
%         if is_trans_h==0
%             title('No transient detected')
%         else
%             if is_ring_h == 1
%                 title('Ring detected')
%             end
%         end
%         axis tight
%         subplot(2,1,2)
%         plot(eye_ts(big_snip),lEyeXY(big_snip,2)-mean(lEyeXY(big_snip,2)),'linewidth',2)
%         hold on
%         plot(eye_ts(big_snip),leye_deriv(big_snip,2)*10,'r')
%         plot(eye_ts(sac_start_inds(i)),lEyeXY(sac_start_inds(i),2)-mean(lEyeXY(big_snip,2)),'go','linewidth',4)
%         plot(eye_ts(sac_stop_inds(i)),lEyeXY(sac_stop_inds(i),2)-mean(lEyeXY(big_snip,2)),'ro','linewidth',4)
%         if ~isnan(left_adj_sac_end(i,2))
%             plot(eye_ts(left_adj_sac_end(i,2)),lEyeXY(left_adj_sac_end(i,2),2)-mean(lEyeXY(big_snip,2)),'ko','linewidth',4)
%         end
%         if is_trans_v==0
%             title('No transient detected')
%         else
%             if is_ring_v == 1
%                 title('Ring detected')
%             end
%         end
%         axis tight
%         
%         pause
%         clf





        %for horizontal component
        within_sac_inds = sac_start_inds(i):sac_stop_inds(i);
        [pk,pkloc] = max(abs(reye_deriv(within_sac_inds,1)));
        pk_amp = reye_deriv(sac_start_inds(i)-1+pkloc,1);
        next_event = all_event_start_inds(find(all_event_start_inds > sac_start_inds(i),1,'first'));
        until_nextsac_inds = (sac_start_inds(i)+pkloc-1):next_event;
        
        cur_snippet = reye_deriv(until_nextsac_inds,1);
        cur_snippet = -1*cur_snippet*sign(pk_amp);
        
        is_trans_h = 1;
        [rebound_amp,rebound_peak] = findpeaks(cur_snippet,'minpeakheight',min_rebound_amp,'npeaks',1);
        if ~isempty(rebound_peak)
            decay_signal = cur_snippet(rebound_peak:end);
        else
            is_trans_h = 0;
            decay_signal = -cur_snippet;
            rebound_peak = 1;
        end
        
        if ~isempty(decay_signal)
            
            %check for biphasic 'ringing' response
            is_ring_h = 0;
            first_bounce = find(decay_signal(1:end-1) < decay_signal(2:end),1,'first');
            if ~isempty(first_bounce)
                bounce_amp = decay_signal(first_bounce);
                rel_bounce_amp = abs(decay_signal(first_bounce)/decay_signal(1));
                if rel_bounce_amp > min_rbounce_amp & bounce_amp < min_bounce_amp
                    rebound_peak = first_bounce - 1 + rebound_peak;
                    decay_signal = -decay_signal(first_bounce:end);
                    is_ring_h = 1;
                end
            end
            if length(decay_signal) > 10
                xx = (0:(length(decay_signal)-1))';
                p0 = [median(decay_signal) decay_signal(1) 0.1];
                beta = lsqcurvefit(exp_fun,p0,xx,decay_signal,lb,ub,options);
                if beta(2) > .0025
                    right_adj_sac_end(i,1) = round(sac_start_inds(i) + pkloc + rebound_peak -2 + 2.5/beta(3));
                elseif is_ring_h==1
                    right_adj_sac_end(i,1) = round(sac_start_inds(i) + pkloc + rebound_peak-2);
                end
            end
        end
        
        
        %for vertical component
        [pk,pkloc] = max(abs(reye_deriv(within_sac_inds,2)));
        pk_amp = reye_deriv(sac_start_inds(i)-1+pkloc,2);
        next_event = all_event_start_inds(find(all_event_start_inds > sac_start_inds(i),1,'first'));
        until_nextsac_inds = (sac_start_inds(i)+pkloc-1):next_event;
        
        cur_snippet = reye_deriv(until_nextsac_inds,2);
        cur_snippet = -1*cur_snippet*sign(pk_amp);
        
        is_trans_v = 1;
        [rebound_amp,rebound_peak] = findpeaks(cur_snippet,'minpeakheight',min_rebound_amp,'npeaks',1);
        if ~isempty(rebound_peak)
            decay_signal = cur_snippet(rebound_peak:end);
        else
            is_trans_v = 0;
            decay_signal = -cur_snippet;
            rebound_peak = 1;
        end
        
        if ~isempty(decay_signal)
            
            %check for biphasic 'ringing' response
            is_ring_v = 0;
            first_bounce = find(decay_signal(1:end-1) < decay_signal(2:end),1,'first');
            if ~isempty(first_bounce)
                bounce_amp = decay_signal(first_bounce);
                rel_bounce_amp = abs(decay_signal(first_bounce)/decay_signal(1));
                if rel_bounce_amp > min_rbounce_amp & bounce_amp < min_bounce_amp
                    rebound_peak = first_bounce - 1 + rebound_peak;
                    decay_signal = -decay_signal(first_bounce:end);
                    is_ring_v = 1;
                end
            end
            
            if length(decay_signal) > 10
                xx = (0:(length(decay_signal)-1))';
                p0 = [median(decay_signal) decay_signal(1) 0.1];
                beta = lsqcurvefit(exp_fun,p0,xx,decay_signal,lb,ub,options);
                if beta(2) > .0025
                    right_adj_sac_end(i,2) = round(sac_start_inds(i) + pkloc + rebound_peak -2 + 2.5/beta(3));
                elseif is_ring_v==1
                    right_adj_sac_end(i,2) = round(sac_start_inds(i) + pkloc + rebound_peak-2);
                end
            end
        end
        %%
%         big_snip = (sac_start_inds(i)-20):next_event;
%         big_snip(big_snip < 0) = [];
%         subplot(2,1,1)
%         plot(eye_ts(big_snip),lEyeXY(big_snip,1)-mean(lEyeXY(big_snip,1)),'linewidth',2)
%         hold on
%         plot(eye_ts(big_snip),leye_deriv(big_snip,1)*10,'r')
%         plot(eye_ts(sac_start_inds(i)),lEyeXY(sac_start_inds(i),1)-mean(lEyeXY(big_snip,1)),'go','linewidth',4)
%         plot(eye_ts(sac_stop_inds(i)),lEyeXY(sac_stop_inds(i),1)-mean(lEyeXY(big_snip,1)),'ro','linewidth',4)
%         if ~isnan(left_adj_sac_end(i,1))
%             plot(eye_ts(left_adj_sac_end(i,1)),lEyeXY(left_adj_sac_end(i,1),1)-mean(lEyeXY(big_snip,1)),'ko','linewidth',4)
%         end
%         if is_trans_h==0
%             title('No transient detected')
%         else
%             if is_ring_h == 1
%                 title('Ring detected')
%             end
%         end
%         axis tight
%         subplot(2,1,2)
%         plot(eye_ts(big_snip),lEyeXY(big_snip,2)-mean(lEyeXY(big_snip,2)),'linewidth',2)
%         hold on
%         plot(eye_ts(big_snip),leye_deriv(big_snip,2)*10,'r')
%         plot(eye_ts(sac_start_inds(i)),lEyeXY(sac_start_inds(i),2)-mean(lEyeXY(big_snip,2)),'go','linewidth',4)
%         plot(eye_ts(sac_stop_inds(i)),lEyeXY(sac_stop_inds(i),2)-mean(lEyeXY(big_snip,2)),'ro','linewidth',4)
%         if ~isnan(left_adj_sac_end(i,2))
%             plot(eye_ts(left_adj_sac_end(i,2)),lEyeXY(left_adj_sac_end(i,2),2)-mean(lEyeXY(big_snip,2)),'ko','linewidth',4)
%         end
%         if is_trans_v==0
%             title('No transient detected')
%         else
%             if is_ring_v == 1
%                 title('Ring detected')
%             end
%         end
%         axis tight
%         
%         pause
%         clf

    end
    h_nan = find(isnan(left_adj_sac_end(:,1)));
    left_adj_sac_end(h_nan,1) = sac_stop_inds(h_nan);
    v_nan = find(isnan(left_adj_sac_end(:,2)));
    left_adj_sac_end(v_nan,2) = sac_stop_inds(v_nan);
    
    h_nan = find(isnan(right_adj_sac_end(:,1)));
    right_adj_sac_end(h_nan,1) = sac_stop_inds(h_nan);
    v_nan = find(isnan(right_adj_sac_end(:,2)));
    right_adj_sac_end(v_nan,2) = sac_stop_inds(v_nan);
    
else
    n_sacs = 0;
%     sac_data = struct('start_inds',[],'stop_inds',[],'left_amps',[], 'right_amps',[],'left_camps',[],'right_camps',[]);
end

%compute saccade amps
left_sac_amps = nan(n_sacs,2);
right_sac_amps = nan(n_sacs,2);
left_sac_camps = nan(n_sacs,2);
right_sac_camps = nan(n_sacs,2);
for i = 1:n_sacs
   left_sac_amps(i,:) = lEyeXY(sac_stop_inds(i),:) - lEyeXY(sac_start_inds(i),:);
   right_sac_amps(i,:) = rEyeXY(sac_stop_inds(i),:) - rEyeXY(sac_start_inds(i),:);
   
   left_sac_camps(i,1) = lEyeXY(left_adj_sac_end(i,1),1) - lEyeXY(sac_start_inds(i),1);
   left_sac_camps(i,2) = lEyeXY(left_adj_sac_end(i,2),2) - lEyeXY(sac_start_inds(i),2);
   right_sac_camps(i,1) = rEyeXY(right_adj_sac_end(i,1),1) - rEyeXY(sac_start_inds(i),1);
   right_sac_camps(i,2) = rEyeXY(right_adj_sac_end(i,2),2) - rEyeXY(sac_start_inds(i),2);
end

sac_data = struct('start_inds',mat2cell(sac_start_inds),'stop_inds',mat2cell(sac_stop_inds),...
    'left_amps',mat2cell(left_sac_amps), 'right_amps',mat2cell(right_sac_amps),...
    'left_camps',mat2cell(left_sac_camps),'right_camps',mat2cell(right_sac_camps),...
    'left_stop_inds',mat2cell(left_adj_sac_end),'right_stop_inds',mat2cell(right_adj_sac_end));


%%
% plot(eye_ts-eye_ts(1),in_blink,'k')
% hold on
% plot(eye_ts-eye_ts(1),lEyeXY(:,1))
% plot(eye_ts-eye_ts(1),rEyeXY(:,1),'r')
% % plot(eye_ts(sac_stop_inds)-eye_ts(1),lEyeXY(sac_stop_inds,1),'ko','linewidth',2)
% plot(eye_ts(sac_start_inds)-eye_ts(1),lEyeXY(sac_start_inds,1),'ro','linewidth',2)
% plot(eye_ts(left_adj_sac_end(:,1))-eye_ts(1),lEyeXY(left_adj_sac_end(:,1),1),'go','linewidth',2)
% plot(eye_ts(sac_start_inds)-eye_ts(1),rEyeXY(sac_start_inds,1),'ro','linewidth',2)
% plot(eye_ts(right_adj_sac_end(:,1))-eye_ts(1),rEyeXY(right_adj_sac_end(:,1),1),'go','linewidth',2)