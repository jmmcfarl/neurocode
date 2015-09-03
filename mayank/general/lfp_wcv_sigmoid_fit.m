% Do a sigmoid fit to the up and down transitions save the numbers: 1)
% amplitude of the state. 2) Exact time of the transition. 3) Speed of
% transition. 4) Y-offset of the transition.

load transitions;
load lfpwcv;

% global synct

lf8 = loflf8;
wcv = lofwcv;
%%
% find the duration of LFP up states.
% relevant variables lf8_up_start_time lf8_up_end_time wcv_up_start_time
% wcv_up_end_time

% ensure that the lfp data starts with an UP state in the LFP.
lf8_up_start_id = find(lf8_up_start_id);
lf8_up_end_id = find(lf8_up_end_id);
wcv_up_start_id = find(wcv_up_start_id);
wcv_up_end_id = find(wcv_up_end_id);

while lf8_up_end_time(1) - lf8_up_start_time(1) < 0
    lf8_up_end_time(1)=[];
    lf8_up_end_id(1)=[];
end
while lf8_up_end_time(end)- lf8_up_start_time(end) < 0
    lf8_up_start_time(end)=[];
    lf8_up_start_id(end)=[];
end

while wcv_up_end_time(1)- wcv_up_start_time(1) < 0
    wcv_up_end_time(1)=[];
    wcv_up_end_id(1)=[];
end
while wcv_up_end_time(end)- wcv_up_start_time(end) < 0
    wcv_up_start_time(end)=[];
    wcv_up_start_id(end)=[];
end
while wcv_up_start_time(1)-lf8_up_start_time(1) < 0
    wcv_up_start_time(1) = [];
    wcv_up_end_time(1) = [];
    wcv_up_start_id(1) = [];
    wcv_up_end_id(1) = [];
end

% find the true times of state transitions.
lf8_raw_up_start_id = lf8_up_start_id;
lf8_raw_up_start_time = lf8_up_start_time;
wcv_raw_up_start_id = wcv_up_start_id;
wcv_raw_up_start_time = wcv_up_start_time;
lf8_raw_up_end_id = lf8_up_end_id;
lf8_raw_up_end_time = lf8_up_end_time;
wcv_raw_up_end_id = wcv_up_end_id;
wcv_raw_up_end_time = wcv_up_end_time;

lf8_raw_up_start_tau = lf8_up_start_id;
lf8_raw_up_end_tau = lf8_up_start_id;
lf8_raw_up_start_amp = lf8_up_start_id;
lf8_raw_up_end_amp = lf8_up_start_id;
lf8_raw_up_start_shift = lf8_up_start_id;
lf8_raw_up_end_shift = lf8_up_start_id;
lf8_raw_up_start_error = lf8_up_start_id;
lf8_raw_up_end_error = lf8_up_start_id;

wcv_raw_up_start_tau = wcv_up_start_id;
wcv_raw_up_end_tau = wcv_up_start_id;
wcv_raw_up_start_amp = wcv_up_start_id;
wcv_raw_up_end_amp = wcv_up_start_id;
wcv_raw_up_start_shift = wcv_up_start_id;
wcv_raw_up_end_shift = wcv_up_start_id;
wcv_raw_up_start_error = wcv_up_start_id;
wcv_raw_up_end_error = wcv_up_start_id;

t = -800:800;
t=t';

[lf8_raw_up_start_id,lf8_raw_up_start_time,lf8_raw_up_start_amp,...
    lf8_raw_up_start_shift,lf8_raw_up_start_tau,lf8_raw_up_start_error] ...
    = get_lfp_wcv_sigmoid_fit(lf8_up_start_id,lf8,t);

[lf8_raw_up_end_id,lf8_raw_up_end_time,lf8_raw_up_end_amp,...
    lf8_raw_up_end_shift,lf8_raw_up_end_tau,lf8_raw_up_end_error] ...
    = get_lfp_wcv_sigmoid_fit(lf8_up_end_id,-lf8,t); % Note, pass -1 times lf8

[wcv_raw_up_start_id,wcv_raw_up_start_time,wcv_raw_up_start_amp,...
    wcv_raw_up_start_shift,wcv_raw_up_start_tau,wcv_raw_up_start_error,net_betafit_up] ...
    = get_lfp_wcv_sigmoid_fit(wcv_up_start_id,wcv,t);

[wcv_raw_up_end_id,wcv_raw_up_end_time,wcv_raw_up_end_amp,...
    wcv_raw_up_end_shift,wcv_raw_up_end_tau,wcv_raw_up_end_error,net_betafit_down] ...
    = get_lfp_wcv_sigmoid_fit(wcv_up_end_id,-wcv,t); % Note, pass -1 times wcv
% to get +ve sigmoid fit to down transtions.

break

%     plot(t,xx,'bo',...
%          t,betafit(1)*(1.0./(1.0 + exp((betafit(2)-t)/betafit(3)))) +
%          betafit(4),'r');


%%
% Get rid of transitions that yield very short up or down states


lf8_up_dur = lf8_up_end_time-lf8_up_start_time;
lf8_dn_dur = lf8_up_start_time(2:end)-lf8_up_end_time(1:end-1);
lf8_raw_up_dur = lf8_raw_up_end_time-lf8_raw_up_start_time;
lf8_raw_dn_dur = lf8_raw_up_start_time(2:end)-lf8_raw_up_end_time(1:end-1);

lf8_up_area = lf8_up_dur*0;
lf8_raw_up_area = lf8_up_area;
lf8_dn_area = lf8_dn_dur*0;
lf8_raw_dn_area = lf8_dn_area;

lf8_up_skew = lf8_up_area;
lf8_raw_up_skew = lf8_up_area;
lf8_dn_skew = lf8_dn_area;
lf8_raw_dn_skew = lf8_dn_area;

for k = 1:length(lf8_up_start_id)
    a1 = lf8(lf8_up_start_id(k):lf8_up_end_id(k));
%    a1(a1 < 0)= 0;
    ar1 = lf8(lf8_raw_up_start_id(k):lf8_raw_up_end_id(k));
%    ar1(ar1<0) = 0;
    
    lf8_up_area(k) = sum( a1 );
    lf8_raw_up_area(k) = sum( ar1 );

    lf8_up_mean_amp(k) = mean( a1 );
    lf8_raw_up_mean_amp(k) = mean( ar1 );

    lf8_up_max_amp(k) = max( a1 );
    lf8_raw_up_max_amp(k) = max( ar1 );
    
    [amean,amed,asig,adev,lf8_up_skew(k),askew,akurt,anorkurt] = ...
        allwtstat(1:lf8_up_end_id(k)-lf8_up_start_id(k)+1,a1');
    [amean,amed,asig,adev,lf8_raw_up_skew(k),askew,akurt,anorkurt] = ...
        allwtstat(1:lf8_raw_up_end_id(k)-lf8_raw_up_start_id(k)+1,ar1');

    lf8_up_cor(k) = mycorrcoef(1:lf8_up_end_id(k)-lf8_up_start_id(k)+1,a1');
    lf8_raw_up_cor(k) = mycorrcoef(1:lf8_raw_up_end_id(k)-lf8_raw_up_start_id(k)+1,ar1');

end
for k = 1:length(lf8_up_start_id)-1
    a1 = -lf8(lf8_up_end_id(k):lf8_up_start_id(k+1));
%    a1(a1 < 0) = 0;
    ar1 = -lf8(lf8_raw_up_end_id(k):lf8_raw_up_start_id(k+1));
%    ar1(ar1<0) = 0;
    
    lf8_dn_area(k) = sum( a1 );
    lf8_raw_dn_area(k) = sum( ar1 );

    lf8_dn_mean_amp(k) = mean( a1 );
    lf8_raw_dn_mean_amp(k) = mean( ar1 );

    lf8_dn_max_amp(k) = max( a1 );
    lf8_raw_dn_max_amp(k) = max( ar1 );
    
    [amean,amed,asig,adev,lf8_dn_skew(k),askew,akurt,anorkurt] = ...
        allwtstat(1:lf8_up_start_id(k+1)-lf8_up_end_id(k)+1,a1');

    [amean,amed,asig,adev,lf8_raw_dn_skew(k),askew,akurt,anorkurt] = ...
        allwtstat(1:lf8_raw_up_start_id(k+1)-lf8_raw_up_end_id(k)+1,ar1');    

    lf8_dn_cor(k) = mycorrcoef(1:lf8_up_start_id(k+1)-lf8_up_end_id(k)+1,a1');
    lf8_raw_dn_cor(k) = mycorrcoef(1:lf8_raw_up_start_id(k+1)-lf8_raw_up_end_id(k)+1,ar1');
end

lf8_up_dn_ratio1 = lf8_dn_dur./lf8_up_dur(1:end-1);
lf8_up_dn_ratio2 = lf8_dn_dur./lf8_up_dur(2:end);

wcv_up_dur = wcv_up_end_time-wcv_up_start_time;
wcv_dn_dur = wcv_up_start_time(2:end)-wcv_up_end_time(1:end-1);
wcv_raw_up_dur = wcv_raw_up_end_time-wcv_raw_up_start_time;
wcv_raw_dn_dur = wcv_raw_up_start_time(2:end)-wcv_raw_up_end_time(1:end-1);

wcv_up_area = wcv_up_dur*0;
wcv_raw_up_area = wcv_up_area;
wcv_up_mean_amp = wcv_up_area;
wcv_raw_up_mean_amp = wcv_up_area;
wcv_up_max_amp = wcv_up_area;
wcv_raw_up_max_amp = wcv_up_area;
wcv_up_skew = wcv_up_area;
wcv_raw_up_skew = wcv_up_area;

wcv_dn_area = wcv_dn_dur*0;
wcv_raw_dn_area = wcv_dn_area;
wcv_dn_mean_amp = wcv_dn_area;
wcv_raw_dn_mean_amp = wcv_dn_area;
wcv_dn_max_amp = wcv_dn_area;
wcv_raw_dn_max_amp = wcv_dn_area;
wcv_dn_skew = wcv_dn_area;
wcv_raw_dn_skew = wcv_dn_area;

for k = 1:length(wcv_up_start_id)
    a1 = wcv(wcv_up_start_id(k):wcv_up_end_id(k));
    a1(a1 < 0) = 0;
    ar1 = wcv(wcv_raw_up_start_id(k):wcv_raw_up_end_id(k)); 
    ar1(ar1 < 0) = 0;
    
    wcv_up_area(k) = sum( a1 );
    wcv_raw_up_area(k) = sum( ar1 );
    
    wcv_up_mean_amp(k)= mean(a1);
    wcv_raw_up_mean_amp(k) = mean(ar1);
    
    wcv_up_max_amp(k) = max( a1 );
    wcv_raw_up_max_amp(k) = max( ar1 );
    
    [amean,amed,asig,adev,wcv_up_skew(k),askew,akurt,anorkurt] = ...
        allwtstat(1:wcv_up_end_id(k)-wcv_up_start_id(k)+1,a1');
    [amean,amed,asig,adev,wcv_raw_up_skew(k),askew,akurt,anorkurt] = ...
        allwtstat(1:wcv_raw_up_end_id(k)-wcv_raw_up_start_id(k)+1,ar1');

    wcv_up_cor(k) = mycorrcoef(1:wcv_up_end_id(k)-wcv_up_start_id(k)+1,a1');
    wcv_raw_up_cor(k) = mycorrcoef(1:wcv_raw_up_end_id(k)-wcv_raw_up_start_id(k)+1,ar1');
end
for k = 1:length(wcv_up_start_id)-1
    a1 = -wcv(wcv_up_end_id(k):wcv_up_start_id(k+1));
    ar1 = -wcv(wcv_raw_up_end_id(k):wcv_raw_up_start_id(k+1));
    a1(a1 < 0) = 0;
    ar1(ar1<0) = 0;
    
    wcv_dn_area(k) = sum( a1 );
    wcv_raw_dn_area(k) = sum( ar1 );

    wcv_dn_mean_amp(k) = mean( a1 );
    wcv_raw_dn_mean_amp(k) = mean( ar1 );

    wcv_dn_max_amp(k) = max( a1 );
    wcv_raw_dn_max_amp(k) = max( ar1 );
    
    [amean,amed,asig,adev,wcv_dn_skew(k),askew,akurt,anorkurt] = ...
        allwtstat(1:wcv_up_start_id(k+1)-wcv_up_end_id(k)+1,a1');
    [amean,amed,asig,adev,wcv_raw_up_skew(k),askew,akurt,anorkurt] = ...
        allwtstat(1:wcv_raw_up_start_id(k+1)-wcv_raw_up_end_id(k)+1,ar1');    

    wcv_dn_cor(k) = mycorrcoef(1:wcv_up_start_id(k+1)-wcv_up_end_id(k)+1,a1');
    wcv_raw_dn_cor(k) = mycorrcoef(1:wcv_raw_up_start_id(k+1)-wcv_raw_up_end_id(k)+1,ar1');
end

wcv_up_dn_ratio1 = wcv_dn_dur./wcv_up_dur(1:end-1);
wcv_up_dn_ratio2 = wcv_dn_dur./wcv_up_dur(2:end);

%%
figure(1)
plot(lf8_up_start_time,lf8_up_dur,'b',wcv_up_start_time,wcv_up_dur,'r')
title('Up start versus up duration');
legend('LFP','WCV');

figure(2)
plot(lf8_up_end_time(1:end-1),lf8_dn_dur,'b',wcv_up_end_time(1:end-1),wcv_dn_dur,'r')
title('Up end versus down duration');
legend('LFP','WCV');

wcv_min_lf8_start_time = wcv_up_start_time*0;
wcv_min_lf8_start_id = wcv_min_lf8_start_time;
wcv_min_lf8_start_up_phase = wcv_min_lf8_start_time;
wcv_lf8_up_dur = wcv_min_lf8_start_id;
wcv_min_lf8_start_neg_time = wcv_min_lf8_start_time;
wcv_min_lf8_start_neg_phase = wcv_min_lf8_start_time;
wcv_min_lf8_start_phase = wcv_min_lf8_start_time;
wcv_lf8_dn_dur = wcv_min_lf8_start_time;
wcv_lf8_dur = wcv_min_lf8_start_time;
wcv_lf8_up_area = wcv_min_lf8_start_time;
wcv_lf8_period = wcv_min_lf8_start_time;
wcv_lf8_up_amp = wcv_min_lf8_start_time;
num_lf8_up = length(lf8_up_start_time);
num_lf8_dn = length(lf8_dn_dur);
num_wcv_up = length(wcv_up_start_time);
num_wcv_dn = length(wcv_dn_dur);

for k = 1:length(wcv_up_start_time)
    [tmp,phid] = min(abs(wcv_up_start_time(k)-lf8_up_start_time));
%     while (wcv_up_start_time(k)-lf8_up_start_time(phid) < 0)
%         phid = phid-1;
%     end
    a1 = wcv_up_start_time(k)-lf8_up_start_time(phid);
    wcv_min_lf8_start_time(k) = a1;
    wcv_min_lf8_start_id(k) = phid;
    wcv_lf8_up_area(k) = lf8_up_area(phid);
    wcv_min_lf8_start_up_phase(k) =  a1/lf8_up_dur(phid);
    wcv_lf8_up_dur(k)=lf8_up_dur(phid);
    wcv_lf8_up_amp(k)=wcv_raw_up_max_amp(phid);
    wcv_lf8_period(k) = lf8_up_dur(phid)+lf8_dn_dur(min(phid,num_lf8_dn));
    if phid < length(lf8_up_start_time)
        wcv_min_lf8_start_neg_time(k) = lf8_up_start_time(phid+1)-wcv_up_start_time(k);
        wcv_min_lf8_start_neg_phase(k) = wcv_min_lf8_start_neg_time(k)/(lf8_up_dur(phid)+lf8_dn_dur(phid));
        wcv_min_lf8_start_phase(k) =  a1/(lf8_up_dur(phid)+lf8_dn_dur(phid));
        wcv_lf8_dn_dur(k)=lf8_dn_dur(phid);
        wcv_lf8_dur(k) = lf8_up_dur(phid)+lf8_dn_dur(phid);
    else
        wcv_min_lf8_start_neg_time(k) = 0;
        wcv_min_lf8_start_phase(k) = 0;
        wcv_min_lf8_start_neg_phase(k) = 0;
        wcv_lf8_dn_dur(k)=0;
        wcv_lf8_dur(k) = 0;
    end
end

figure(3)
plot(wcv_up_start_time,wcv_min_lf8_start_neg_time,'r-o')
title('WCV up start time with respect to LF8 up start time');

'Raw Lf8 duration up-down or down-up'
[mycorrcoef(lf8_raw_up_dur(1:end-1),lf8_raw_dn_dur(1:end)) ...
    mycorrcoef(lf8_raw_up_dur(2:end),lf8_raw_dn_dur(1:end))]

'Raw lf8 area  up-down or down-up'
[mycorrcoef(lf8_raw_up_area(1:end-1),lf8_raw_dn_area(1:end)) ...
    mycorrcoef(lf8_raw_up_area(2:end),lf8_raw_dn_area(1:end))]

'Raw lf8 mean_amp  up-down or down-up'
[mycorrcoef(lf8_raw_up_mean_amp(1:end-1),lf8_raw_dn_mean_amp(1:end)) ...
    mycorrcoef(lf8_raw_up_mean_amp(2:end),lf8_raw_dn_mean_amp(1:end))]

'Raw lf8 max_amp  up-down or down-up'
[mycorrcoef(lf8_raw_up_max_amp(1:end-1),lf8_raw_dn_max_amp(1:end)) ...
    mycorrcoef(lf8_raw_up_max_amp(2:end),lf8_raw_dn_max_amp(1:end))]

'lf8 area  up-down or down-up'
[mycorrcoef(lf8_up_area(1:end-1),lf8_dn_area(1:end)) ...
    mycorrcoef(lf8_up_area(2:end),lf8_dn_area(1:end))]

'Raw wcv duration up-down or down-up'
[mycorrcoef(wcv_raw_up_dur(1:end-1),wcv_raw_dn_dur(1:end)) ...
    mycorrcoef(wcv_raw_up_dur(2:end),wcv_raw_dn_dur(1:end))]

'Raw wcv area  up-down or down-up'
[mycorrcoef(wcv_raw_up_area(1:end-1),wcv_raw_dn_area(1:end)) ...
    mycorrcoef(wcv_raw_up_area(2:end),wcv_raw_dn_area(1:end))]

'Raw wcv mean_amp  up-down or down-up'
[mycorrcoef(wcv_raw_up_mean_amp(1:end-1),wcv_raw_dn_mean_amp(1:end)) ...
    mycorrcoef(wcv_raw_up_mean_amp(2:end),wcv_raw_dn_mean_amp(1:end))]

'Raw wcv max_amp  up-down or down-up'
[mycorrcoef(wcv_raw_up_max_amp(1:end-1),wcv_raw_dn_max_amp(1:end)) ...
    mycorrcoef(wcv_raw_up_max_amp(2:end),wcv_raw_dn_max_amp(1:end))]

[lf8_raw_up_dn_dur_xcor,alag]=xcov(lf8_raw_up_dur(1:end-1),lf8_raw_dn_dur,300,'coeff');
[lf8_raw_up_dn_area_xcor,alag]=xcov(lf8_raw_up_area(1:end-1),lf8_raw_dn_area,300,'coeff');
[lf8_raw_up_dn_meanamp_xcor,alag]=xcov(lf8_raw_up_mean_amp(1:end-1),lf8_raw_dn_mean_amp,300,'coeff');
[lf8_raw_up_dn_maxamp_xcor,alag]=xcov(lf8_raw_up_max_amp(1:end-1),lf8_raw_dn_max_amp,300,'coeff');
figure(4)
plot(alag,lf8_raw_up_dn_dur_xcor,'r',...
     alag,lf8_raw_up_dn_area_xcor,'b',...
     alag,lf8_raw_up_dn_meanamp_xcor,'g',...
     alag,lf8_raw_up_dn_maxamp_xcor,'k');
grid;
legend('duration','area','mean amp','max amp');
xlabel('#LF8 state lag');
ylabel('Coefficient');

[wcv_raw_up_dn_dur_xcor,alag]=xcov(wcv_raw_up_dur(1:end-1),wcv_raw_dn_dur,300,'coeff');
[wcv_raw_up_dn_area_xcor,alag]=xcov(wcv_raw_up_area(1:end-1),wcv_raw_dn_area,300,'coeff');
[wcv_raw_up_dn_meanamp_xcor,alag]=xcov(wcv_raw_up_mean_amp(1:end-1),wcv_raw_dn_mean_amp,300,'coeff');
[wcv_raw_up_dn_maxamp_xcor,alag]=xcov(wcv_raw_up_max_amp(1:end-1),wcv_raw_dn_max_amp,300,'coeff');
figure(5)
plot(alag,wcv_raw_up_dn_dur_xcor,'r',...
     alag,wcv_raw_up_dn_area_xcor,'b',...
     alag,wcv_raw_up_dn_meanamp_xcor,'g',...
     alag,wcv_raw_up_dn_maxamp_xcor,'k');
grid;
legend('duration','area','mean amp','max amp');
xlabel('#WCV state lag');
ylabel('Coefficient');

[ldn1,x1]=hist(lf8_raw_dn_cor,-1:.01:1);
[lun1,x1]=hist(lf8_raw_up_cor,-1:.01:1);
[wdn1,x1]=hist(wcv_raw_dn_cor,-1:.01:1);
[wun1,x1]=hist(wcv_raw_up_cor,-1:.01:1);
figure(6)
plot(x1,fgsmooth(wun1,2),'r',x1,fgsmooth(lun1,2),'b');legend('wcv raw up','lf8 raw up');shg
figure(7)
plot(x1,fgsmooth(wdn1,2),'r',x1,fgsmooth(ldn1,2),'b');legend('wcv raw dn','lf8 raw dn');shg

