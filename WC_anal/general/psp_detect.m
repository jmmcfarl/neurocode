function [psp_times,psp_sign,psp_det] = psp_detect(signal,threshold)

%% create template function
dsf = 5;
dt = 0.2*dsf;
Fsd = 1/dt*1000;
template_dur = 80; %in ms

temp_taxis = linspace(0,template_dur,round(template_dur/dt)-1);
temp_length = length(temp_taxis);

tau_R = 20;
tau_D = 40;

e = (1-exp(-temp_taxis/tau_R)).*exp(-temp_taxis/tau_D);

%% initialize calculation
dsignal = downsample(signal,dsf);

pdet = zeros(length(dsignal)-temp_length,1);

temp_sum = sum(e);
temp_sq_sum = sum(e.^2);
norm_fact = temp_sq_sum-temp_sum^2/temp_length;

cur_dsignal = dsignal(1:temp_length);
cur_sig_sum = sum(cur_dsignal);
cur_sig_sq_sum = sum(cur_dsignal.^2);

%% compute psp detection
for i = 1:length(pdet)
    cur_dsignal = dsignal(i:i+temp_length-1);
    %cur_sig_sum = sum(cur_dsignal);
    %cur_sig_sq_sum = sum(cur_dsignal.^2);
    curdot = dot(e,cur_dsignal);
   scale = (curdot - temp_sum*cur_sig_sum/temp_length)/norm_fact;
   offset = (cur_sig_sum-scale*temp_sum)/temp_length;
   fitfun = scale*e+offset;
   %SSE = sum((fitfun'-cur_dsignal).^2);
   SSE = cur_sig_sq_sum+scale^2*temp_sq_sum+temp_length*offset^2 ...
       - 2*(scale*curdot+offset*cur_sig_sum-scale*offset*temp_sum);
   SSE = sqrt(SSE/(temp_length-1));
   pdet(i) = scale/SSE;
   cur_sig_sum = cur_sig_sum - cur_dsignal(1)+dsignal(i+temp_length);
   cur_sig_sq_sum = cur_sig_sq_sum-cur_dsignal(1)^2+dsignal(i+temp_length)^2;
   
end

psp_det = interp1(1:length(pdet),pdet',linspace(1,length(pdet),length(signal)-temp_length*dsf));
if max(pdet) > threshold
    [dummy,epsp_times] = findpeaks(pdet,'minpeakheight',threshold);
    diff_etimes = [0 diff(epsp_times)];
    epsp_times(find(diff_etimes < temp_length)) =  [];
else
    epsp_times = [];
end
if min(pdet) < -threshold
    [dummy,ipsp_times] = findpeaks(-pdet,'minpeakheight',threshold);
    diff_itimes = [0 diff(ipsp_times)];
    ipsp_times(find(diff_itimes < temp_length)) = [];
else
    ipsp_times = [];
end

psp_sign = [ones(1,length(epsp_times)) zeros(1,length(ipsp_times))];
psp_times = [epsp_times ipsp_times];
[psp_times,order] = sort(psp_times);
psp_sign = psp_sign(order);
psp_sign = logical(psp_sign);
psp_times = psp_times*dsf;


