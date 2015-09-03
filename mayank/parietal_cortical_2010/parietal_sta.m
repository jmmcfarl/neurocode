clear all
% close all

load G:\WC_Germany\parietal_cortical_2010\parietal_cortical_2010

Fs = 2016;
niqf = 2016/2;
lcf = 10/niqf;
[b,a] = butter(2,lcf,'high');
maxlag = round(0.25*Fs);
lags = -maxlag:maxlag;

%get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];

frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');

% thomas_el = find_struct_field_vals(sess_data,'thom_elec',1);
% sess_data = sess_data(thomas_el);
% desynch_times = desynch_times(thomas_el);

n = length(sess_data);
set = [14 16 17 20];
for d = set

    to_dir = sess_data(d).directory;
    to_dir(1) = 'G';
    disp(sprintf('session %d',d))
    cd(to_dir);

    load used_data wcv 
    load spike_time_jmm

%     lf8_f = filtfilt(b,a,lf8)/sess_data(d).gains(8);
%     lf5_f = filtfilt(b,a,lf5)/sess_data(d).gains(5);
    wcv_f = filtfilt(b,a,wcv)/sess_data(d).gains(5);

    if sess_data(d).gains(6) ~= 0
%         lf7_f = filtfilt(b,a,lf7)/sess_data(d).gains(7);
%         lf6_f = filtfilt(b,a,lf6)/sess_data(d).gains(6);
%         csd = (2*lf7_f-lf6_f-lf8_f);
    end

    spkid(spkid < maxlag | spkid > length(wcv_f) - maxlag) = [];

    n_spikes = length(spkid);
%     stmat_lf8 = zeros(n_spikes,length(lags));
%     stmat_lf5 = stmat_lf8;
    stmat_wcv = zeros(n_spikes,length(lags));
%     if sess_data(d).gains(6) ~= 0
%         stmat_csd = zeros(n_spikes,length(lags));
%     end

    for i = 1:n_spikes
%         stmat_lf8(i,:) = lf8_f(spkid(i)-maxlag:spkid(i)+maxlag);
%         stmat_lf5(i,:) = lf5_f(spkid(i)-maxlag:spkid(i)+maxlag);
        stmat_wcv(i,:) = wcv_f(spkid(i)-maxlag:spkid(i)+maxlag);
%         if sess_data(d).gains(6) ~= 0
%             stmat_csd(i,:) = csd(spkid(i)-maxlag:spkid(i)+maxlag);
% 
%         end
    end

    figure
    plot(lags/Fs,mean(stmat_wcv)), hold on
%     plot(lags/Fs,mean(stmat_lf5),'r')
%     if sess_data(d).gains(6) ~= 0
%         plot(lags/Fs,mean(stmat_csd),'k')
%     end
    xlim([-maxlag/Fs maxlag/Fs])
    xlabel('Time (s)','fontsize',16)
    ylabel('Amplitude (mV)','fontsize',16)
%     legend('LF8','LF5','CSD')
    pause
    close all

end
%%
