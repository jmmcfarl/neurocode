 up = 10:28;
 backlag = 5;
 forlag = 20;
 
 cc = 102;
 cur_Robs = Robs_mat(:,cc);
 %%
temp = find(all(X(1:end,up) == 0,2));
temp = find(ismember(used_inds,temp));
[aa,la] = get_event_trig_avg(cur_Robs,temp,backlag,forlag);

 
 temp = find(all(X(1:end-1,up) == 0,2) & all(X(2:end,up)==0,2));
temp = find(ismember(used_inds,temp));
[aa2,la] = get_event_trig_avg(cur_Robs,temp,backlag,forlag);

temp = find(all(X(1:end-2,up) == 0,2) & all(X(2:end-1,up)==0,2) & all(X(3:end,up)==0,2));
temp = find(ismember(used_inds,temp));
[aa3,la] = get_event_trig_avg(cur_Robs,temp,backlag,forlag);


temp = find(all(X(1:end-3,up) == 0,2) & all(X(2:end-2,up)==0,2) & all(X(3:end-1,up)==0,2) ...
    & all(X(4:end,up)==0,2));
temp = find(ismember(used_inds,temp));
[aa4,la] = get_event_trig_avg(cur_Robs,temp,backlag,forlag);

temp = find(all(X(1:end-4,up) == 0,2) & all(X(2:end-3,up)==0,2) & all(X(3:end-2,up)==0,2) ...
    & all(X(4:end-1,up)==0,2) & all(X(5:end,up) == 0,2));
temp = find(ismember(used_inds,temp));
[aa5,la] = get_event_trig_avg(cur_Robs,temp,backlag,forlag);

%%
figure;
plot(la,aa,la,aa2,'r',la,aa3,'g',la,aa4,'k',la,aa5,'m')
