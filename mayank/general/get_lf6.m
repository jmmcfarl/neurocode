% Compute the shift between the LFP and MP
% clear all

if exist('all_eeg_data.mat')
    load all_eeg_data
elseif exist('part1_eeg_data.mat')
    load part1_eeg_data
elseif exist('part2_eeg_data.mat')
    load part2_eeg_data
else
    disp('ERROR DATA DOES NOT EXIST')
    return
end

clear CSC4* 
clear *NumberValidSamples *ChannelNumbers

set(0,'DefaultAxesFontSize',15);
set(0,'DefaultTextFontSize',15);

load sync_times
if exist('spike_time.mat')
    load spike_time
elseif exist('spike_time_br.mat')
    load spike_time_br
else
    disp('ERROR NO SPIKE TIME DATA')
end

global synct

[len,widwcv] = size(CSC1_Samples);
[widlfp] = length(CSC2_Samples(1,:));

% Get the data.

clear CSC1_Samples;

lf6 = reshape(CSC6_Samples,len*widlfp,1);
lf6 = lf6(synct2id);

% clear *Samples *TimeStamps *SampleFrequencies
% clear CSC2_Samples CSC5_Samples CSC8_Samples;
clear *_Samples

lf6 = detrend(lf6);


save used_data_lf6 lf6
% dataViewer(lf8,wcv_minus_spike,synct,dt)


% spktim = zeros(size(wcv));
% spktim(spkid) = 1;
% 
% % Get low passband filtered data.
% f=CSC8_SampleFrequencies(1);
% lof = .1; hif = 1; nyqf = f/2;
% lof = lof/nyqf; hif = hif/nyqf;
% [b,a] = butter(2, [lof hif]);
% 
% loflf8 = filtfilt(b,a,lf8);
% lofwcv = filtfilt(b,a,wcv_minus_spike);
% 
% lof = .1; hif = 2.5; nyqf = f/2;
% lof = lof/nyqf; hif = hif/nyqf;
% [b,a] = butter(2, [lof hif]);
% 
% loflf8_1 = filtfilt(b,a,lf8);
% lofwcv_1 = filtfilt(b,a,wcv_minus_spike);
% 
% difloflf8_1 = [0; diff(loflf8_1)];
% diflofwcv_1 = [0; diff(lofwcv_1)];
% 
% % loflf2 = filtfilt(b,a,lf2);
% difloflf8 = [0; diff(loflf8)];
% diflofwcv = [0; diff(lofwcv)];

%****************MAYANK UDS TRANSITION CATCHER**********************
% poslf8std = 0;
% neglf8std = 0;
%
% tmp = loflf8;
% tmp(tmp < poslf8std) = 0;
% tmp(tmp > 0) = 1;
% upon = tmp;
% tmp = [0; diff(tmp)];
% lf8_up_start_id = (tmp > 0);
% lf8_up_start_time = synct(lf8_up_start_id);
% lf8_up_end_id = (tmp < 0);
% lf8_up_end_time = synct(lf8_up_end_id);
% upspk = spktim & upon;
% dnspk = spktim & ~upon;
%
% tmp = lofwcv;
% tmp(tmp < poslf8std) = 0;
% tmp(tmp > 0) = 1;
% upon = tmp;
% tmp = [0; diff(tmp)];
% wcv_up_start_id = (tmp > 0);
% wcv_up_start_time = synct(wcv_up_start_id);
% wcv_up_end_id = (tmp < 0);
% wcv_up_end_time = synct(wcv_up_end_id);
%clear tmp
%*********************************************


%*********************JAMES UDS TRANS EXTRACTION**************
%set these thresholds depending on amp of up modulation.

% poslf8std = 0.4;
% neglf8std = -0.4;
% 
% tmp = loflf8;
% loflf8z = tmp;
% loflf8z = loflf8z - mean(loflf8z);
% loflf8z = loflf8z/std(loflf8z);
% tmp(loflf8z > poslf8std) = 1;
% tmp(loflf8z < neglf8std) = -2;
% tmp(loflf8z < poslf8std & loflf8z > neglf8std) = 0;
% tmp = [0; diff(tmp)];
% lf8_up_start_id = (tmp == 1);
% lf8_up_end_id = (tmp == -2);
% lf8_up_start_points = find(tmp == 1);
% lf8_up_end_points = find(tmp == -2);
% 
% %eliminate instances where up and down transitions aren't alternating
% %properly
% 
% %make sure we start on an up transition and end on a down state
% while min(lf8_up_end_points) < min(lf8_up_start_points)
%     lf8_up_end_id(lf8_up_end_points(1)) = 0;
%     lf8_up_end_points(1) = [];
% end
% while max(lf8_up_start_points) > max(lf8_up_end_points)
%     lf8_up_start_id(lf8_up_start_points(end)) = 0;
%     lf8_up_start_points(end) = [];
% end
% %make sure there is only one up transition before the first down
% firstUps = find(lf8_up_start_points < lf8_up_end_points(1));
% if length(firstUps) > 1
%     lf8_up_start_id(lf8_up_start_points(firstUps(2:end))) = 0;
%     lf8_up_start_points(firstUps(2:end)) = [];
% end
% %make sure only one down transition after last up
% lastDowns = find(lf8_up_end_points > lf8_up_start_points(end));
% if length(lastDowns) > 1
%     lf8_up_end_id(lf8_up_end_points(lastDowns(2:end))) = 0;
%     lf8_up_end_points(lastDowns(2:end)) = [];
% end
% clear firstUps lastDowns
% 
% %cur state is 1 for up and 0 for down
% curID = 1;
% curPoint = lf8_up_start_points(1);
% trashStartPoints = [];
% trashEndPoints = [];
% lastEndPoint = lf8_up_end_points(end);
% 
% % subP = synct(1);
% 
% while curPoint < lastEndPoint
%     
%     %if current state is up
%     if curID
%         %find next down transition
% 
%         nextDown = lf8_up_end_points(find(lf8_up_end_points > curPoint,1,'first'));
%         %mark as throw away any up transitions detected in between these
%         %points
%         trashStartPoints = [trashStartPoints; ...
%             lf8_up_start_points(lf8_up_start_points>curPoint&lf8_up_start_points<nextDown)];
%         %set state to down and fix current point
%         curID = 0;
%         curPoint = nextDown;
% % synct(curPoint) - subP
%         %else if current state is down
%     else
%         %find next up transition
%         nextUp = lf8_up_start_points(find(lf8_up_start_points > curPoint,1,'first'));
%         %mark as throw away any down transitions that occur between points
%         trashEndPoints = [trashEndPoints; ...
%             lf8_up_end_points(lf8_up_end_points>curPoint&lf8_up_end_points<nextUp)];
%         %set state to up and fix current point
%         curID = 1;
%         curPoint = nextUp;
% % synct(curPoint) - subP
%     end
% 
% end
% 
% %get rid of trash points
% lf8_up_start_id(trashStartPoints) = 0;
% lf8_up_end_id(trashEndPoints) = 0;
% for i = 1:length(trashStartPoints)
%     trashStartID(i) = find(lf8_up_start_points == trashStartPoints(i));
% end
% for i = 1:length(trashEndPoints)
%     trashEndID(i) = find(lf8_up_end_points == trashEndPoints(i));
% end
% lf8_up_start_time = synct(lf8_up_start_id);
% lf8_up_end_time = synct(lf8_up_end_id);
% lf8_up_start_points(trashStartID) = [];
% lf8_up_end_points(trashEndID) = [];
% 
% 
% %now go through and find spikes during up and down states
% upon = zeros(size(lf8));%initialize 
% for i = 1:length(lf8_up_start_points)
%    
%     upon(lf8_up_start_points(i):lf8_up_end_points(i)) = 1; 
%     
% end
% upspk = spktim & upon;
% dnspk = spktim & ~upon;
% 
% %********************now for wcv
% 
% poswcvstd = 0.4;
% negwcvstd = -0.4;
% 
% tmp = lofwcv;
% lofwcvz = tmp;
% lofwcvz = lofwcvz - mean(lofwcvz);
% lofwcvz = lofwcvz/std(lofwcvz);
% tmp(lofwcvz > poswcvstd) = 1;
% tmp(lofwcvz < negwcvstd) = -2;
% tmp(lofwcvz < poswcvstd & lofwcvz > negwcvstd) = 0;
% tmp = [0; diff(tmp)];
% wcv_up_start_id = (tmp == 1);
% wcv_up_end_id = (tmp == -2);
% 
% wcv_up_start_points = find(tmp == 1);
% wcv_up_end_points = find(tmp == -2);
% 
% 
% %eliminate instances where up and down transitions aren't alternating
% %properly
% 
% %make sure we start on an up transition and end on a down transition
% while min(wcv_up_end_points) < min(wcv_up_start_points)
%     wcv_up_end_id(wcv_up_end_points(1)) = 0;
%     wcv_up_end_points(1) = [];
% end
% while max(wcv_up_start_points) > max(wcv_up_end_points)
%     wcv_up_start_id(wcv_up_start_points(end)) = 0;
%     wcv_up_start_points(end) = [];
% end
% %make sure there is only one up transition before the first down
% firstUps = find(wcv_up_start_points < wcv_up_end_points(1));
% if length(firstUps) > 1
%     wcv_up_start_id(wcv_up_start_points(firstUps(2:end))) = 0;
%     wcv_up_start_points(firstUps(2:end)) = [];
% end
% %make sure only one down transition after last up
% lastDowns = find(wcv_up_end_points > wcv_up_start_points(end));
% if length(lastDowns) > 1
%     wcv_up_end_id(wcv_up_end_points(lastDowns(2:end))) = 0;
%     wcv_up_end_points(lastDowns(2:end)) = [];
% end
% 
% %cur state is 1 for up and 0 for down
% curID = 1;
% curPoint = wcv_up_start_points(1);
% upon = zeros(size(wcv));
% trashStartPoints = [];
% trashEndPoints = [];
% lastEndPoint = wcv_up_end_points(end);
% 
% % subP = synct(1);
% 
% while curPoint < lastEndPoint
%     
%     %if current state is up
%     if curID
%         %find next down transition
% 
%         nextDown = wcv_up_end_points(find(wcv_up_end_points > curPoint,1,'first'));
%         if isempty(nextDown)
%             break
%         end
%         %mark as throw away any up transitions detected in between these
%         %points
%         trashStartPoints = [trashStartPoints; ...
%             wcv_up_start_points(wcv_up_start_points>curPoint&wcv_up_start_points<nextDown)];
%         %set state to down and fix current point
%         curID = 0;
%         curPoint = nextDown;
% % synct(curPoint) - subP
%         %else if current state is down
%     else
%         %find next up transition
%         nextUp = wcv_up_start_points(find(wcv_up_start_points > curPoint,1,'first'));
%         if isempty(nextUp)
%             break
%         end
%         %mark as throw away any down transitions that occur between points
%         trashEndPoints = [trashEndPoints; ...
%             wcv_up_end_points(wcv_up_end_points>curPoint&wcv_up_end_points<nextUp)];
%         %set state to up and fix current point
%         curID = 1;
%         curPoint = nextUp;
% % synct(curPoint) - subP
%     end
% 
% end
% 
% wcv_up_start_id(trashStartPoints) = 0;
% wcv_up_end_id(trashEndPoints) = 0;
% clear trashStartID trashEndID
% for i = 1:length(trashStartPoints)
%     trashStartID(i) = find(wcv_up_start_points == trashStartPoints(i));
% end
% for i = 1:length(trashEndPoints)
%     trashEndID(i) = find(wcv_up_end_points == trashEndPoints(i));
% end
% wcv_up_start_time = synct(wcv_up_start_id);
% wcv_up_end_time = synct(wcv_up_end_id);
% wcv_up_start_points(trashStartID) = [];
% wcv_up_end_points(trashEndID) = [];
% 
% clear tmp tmp2
% %**************************************************************
% 
% 
% 
% % save transitions upspk dnspk lf8_up_start_id lf8_up_start_time ...
% %     lf8_up_end_id lf8_up_end_time wcv_up_start_id wcv_up_start_time ...
% %     wcv_up_end_id wcv_up_end_time
% 
% 
% 
% % filter the data to remove all signals above 30 Hz, specially the 60
% % cycles noise.
% 
% f=CSC8_SampleFrequencies(1);
% lof = .1; hif = 10; nyqf = f/2;
% lof = lof/nyqf; hif = hif/nyqf;
% [b,a] = butter(2, [lof hif]);
% loflf8 = filtfilt(b,a,lf8);
% lofwcv = filtfilt(b,a,wcv_minus_spike);
% loflf2 = filtfilt(b,a,lf2);

% save lfpwcv lofwcv loflf8 synct dt








% % do ripple analysis
% if 0
%     lof = 80; hif = 200; nyqf = f/2;
%     lof = lof/nyqf; hif = hif/nyqf;
%     [b,a] = butter(2, [lof hif]);
%     hifwcv = filtfilt(b,a,wcv);
%     rip = filtfilt(b,a,lf3);
%     rip = sqrt(rip.^2);
%     ripthresh = mean(rip)+2*std(rip);
%     rip(rip < ripthresh) = 0;
% end

% compute transition-time-triggerred averaged LFP and MP

% translen = 3000;
% lf8_up_start_id = find(lf8_up_start_id);
% lf8_up_end_id = find(lf8_up_end_id);
% % up_trans_lf8 = zeros([length(lf8_up_start_id)-9,2*translen+1]);
% % up_trans_wcv = size(up_trans_lf8);
% % dn_trans_lf8 = zeros([length(lf8_up_end_id)-9,2*translen+1]);
% % dn_trans_wcv = size(dn_trans_lf8);
%
%
% % normalize lf8 and wcv in short segments to have less drift and z-scored
% % values.
%
% % wcv = wcv-mean(wcv);
% % wcv = wcv/std(wcv);
% %
% % lf8 = lf8-mean(lf8);
% % lf8 = lf8/mean(lf8);
% %
% % lf2 = lf2-mean(lf2);
% % lf2 = lf2/std(lf2);
% %
% % lofwcv = lofwcv-mean(lofwcv);
% % lofwcv = lofwcv/std(lofwcv);
% %
% % loflf8 = loflf8-mean(loflf8);
% % loflf8 = loflf8/mean(loflf8);
% %
% % loflf2 = loflf2-mean(loflf2);
% % loflf2 = loflf2/std(loflf2);
%
% nseg = 20; % number of segments in which the data should be split
% len = length(wcv);
% seglen = floor(len/(nseg));
% histax = -3.5:.01:3.5;
%
% m = 0;
%
% lf8minlf2 = zeros([1,nseg]);
% abslf8minabslf2 = lf8minlf2;
% loflf8minloflf2 = lf8minlf2;
% bimodlofwcv = lf8minlf2;
% bimodloflf8 = lf8minlf2;
% bimodloflf2 = lf8minlf2;
%
%
% for k = 1:nseg
%     a = round(seglen*(k-1)+(1:seglen));
%     m = m+1;
%
%     awcv = wcv(a);
%     awcv = detrend(awcv);
%     awcv = awcv-mean(awcv);
%     awcv = awcv/std(awcv);
%     wcv(a) = awcv;
%
%     alf8 = lf8(a);
%     alf8 = detrend(alf8);
%     alf8 = alf8-mean(alf8);
%     alf8 = alf8/std(alf8);
%     lf8(a) = alf8;
%
%     alf2 = lf2(a);
%     alf2 = detrend(alf2);
%     alf2 = alf2-mean(alf2);
%     alf2 = alf2/std(alf2);
%     lf2(a) = alf2;
%
%     alofwcv = lofwcv(a);
%     alofwcv = detrend(alofwcv);
%     alofwcv = alofwcv-mean(alofwcv);
%     alofwcv = alofwcv/std(alofwcv);
%     lofwcv(a) = alofwcv;
%
%     aloflf8 = loflf8(a);
%     aloflf8 = detrend(aloflf8);
%     aloflf8 = aloflf8-mean(aloflf8);
%     aloflf8 = aloflf8/std(aloflf8);
%     loflf8(a) = aloflf8;
%
%     aloflf2 = loflf2(a);
%     aloflf2 = detrend(aloflf2);
%     aloflf2 = aloflf2-mean(aloflf2);
%     aloflf2 = aloflf2/std(aloflf2);
%     loflf2(a) = aloflf2;
%
%     anormwcvhist(k,:) = hist(awcv,histax)/length(awcv);
%     anormlf8hist(k,:) = hist(alf8,histax)/length(alf8);
%     anormlf2hist(k,:) = hist(alf2,histax)/length(alf2);
%     anormlofwcvhist(k,:) = hist(alofwcv,histax)/length(alofwcv);
%     anormloflf8hist(k,:) = hist(aloflf8,histax)/length(aloflf8);
%     anormloflf2hist(k,:) = hist(aloflf2,histax)/length(aloflf2);
%
%     lf8minlf2(k) = mean(alf8-alf2);
%     abslf8minabslf2(k) = mean(abs(alf8)-abs(alf2));
%     loflf8minloflf2(k) = mean(aloflf8-aloflf2);
%
%     alofwcv = alofwcv(abs(alofwcv) < 2.5);
%     aloflf8 = aloflf8(abs(aloflf8) < 2.5);
%     aloflf2 = aloflf2(abs(aloflf2) < 2.5);
%
%     bimodlofwcv(k) = (1+skewness(alofwcv)^2)/(kurtosis(alofwcv)+3);
%     bimodloflf8(k) = (1+skewness(aloflf8)^2)/(kurtosis(aloflf8)+3);
%     bimodloflf2(k) = (1+skewness(aloflf2)^2)/(kurtosis(aloflf2)+3);
%
% end
%
% fig = 0;
% mnormlofwcvhist = mean(anormlofwcvhist);
% snormlofwcvhist = std(anormlofwcvhist);
%
% mnormloflf8hist = mean(anormloflf8hist);
% snormloflf8hist = std(anormloflf8hist)/sqrt(nseg);
%
% mnormloflf2hist = mean(anormloflf2hist);
% snormloflf2hist = std(anormloflf2hist)/sqrt(nseg);
%
% % plot the histograms of the normalized membrane potential and the lfp.
%         if 1
% fig = fig+1; figure(fig);
% clf;
% hold on;
% h=plot(histax,fgsmooth(mnormloflf8hist,0),'b', ...
%        histax,fgsmooth(mnormlofwcvhist,0),'r');
% set(h,'LineWidth',2);
% hlen = (length(histax)-1)/2;
% errorbar(histax(hlen),mnormlofwcvhist(hlen),snormlofwcvhist(hlen),'r');
% errorbar(histax(hlen),mnormloflf8hist(hlen),snormloflf8hist(hlen),'b');
% axis([-3.1 3.1 0 Inf]);
% xlabel('Amplitude (z)'); ylabel('#'); legend('LFP','MP');
% box on;
% drawnow; shg;
% print -djpeg100 -r300 norm_hist_wcv_lf8;
%         end
%
% k2=0;
% for k = 5:length(lf8_up_start_id)-5 % skip the first and last few cycles
%     k2 = k2+1;
%     k1 = lf8_up_start_id(k);
%     lf = lf8(k1-translen:k1+translen)';
%     lf = lf - mean(lf);
%     lf = lf/std(lf);
%     wc = wcv(k1-translen:k1+translen)';
%     wc = wc - mean(wc);
%     wc = wc/std(wc);
%     sp = spktim(k1-translen:k1+translen)';
%     up_trans_lf8(k2,:) = lf;
%     up_trans_wcv(k2,:) = wc;
%     up_trans_spk(k2,:) = sp;
% end
%
% mup_trans_lf8 = mean(up_trans_lf8);
% mup_trans_wcv = mean(up_trans_wcv);
% mup_trans_spk = mean(up_trans_spk);
% sup_trans_lf8 = std(up_trans_lf8);
% sup_trans_wcv = std(up_trans_wcv);
% sup_trans_spk = std(up_trans_spk);
%
% k2 = 0;
% for k = 5:length(lf8_up_end_id)-5 % skip the
%     k2 = k2+1;
%     k1 = lf8_up_end_id(k);
%     lf = lf8(k1-translen:k1+translen)';
%     lf = lf - mean(lf);
%     lf = lf/std(lf);
%     wc = wcv(k1-translen:k1+translen)';
%     wc = wc - mean(wc);
%     wc = wc/std(wc);
%     sp = spktim(k1-translen:k1+translen)';
%     dn_trans_lf8(k2,:) = lf;
%     dn_trans_wcv(k2,:) = wc;
%     dn_trans_spk(k2,:) = sp;
% end
%
% mdn_trans_lf8 = mean(dn_trans_lf8);
% mdn_trans_wcv = mean(dn_trans_wcv);
% mdn_trans_spk = mean(dn_trans_spk);
% sdn_trans_lf8 = std(dn_trans_lf8);
% sdn_trans_wcv = std(dn_trans_wcv);
% sdn_trans_spk = std(dn_trans_spk);
%
% atax = -translen:translen;
% atax = atax*dt;
%
%         if 1
%
% fig = fig+1; figure(fig);
% h=plot(atax,fgsmooth(mup_trans_wcv,4),'b',atax,fgsmooth(mup_trans_lf8,4),'r');
% set(h,'LineWidth',2);
% legend('MP','LFP');
% %a1 = translen + 500;
% %errorbar(atax(a1),mup_trans_wcv(a1),sup_trans_wcv(a1),'b');
% %errorbar(atax(a1),mup_trans_lf8(a1),sup_trans_lf8(a1),'r');
% xlabel('Latency (s)');
% ylabel('Mean amplitude (z)');
% shg
% print -djpeg100 -r300 lfp_up_trig_lfp_wcv.jpg
%
% fig = fig+1; figure(fig);
% h=plot(atax,fgsmooth(mdn_trans_wcv,4),'b',atax,fgsmooth(mdn_trans_lf8,4),'r');
% set(h,'LineWidth',2);
% legend('MP','LFP');
% xlabel('Latency (s)');
% ylabel('Mean amplitude (z)')
% drawnow; shg;
% print -djpeg100 -r300 lfp_dn_trig_lfp_wcv.jpg
%
% fig = fig+1; figure(fig);
% h=plot(atax,fgsmooth(mup_trans_spk,4),'b',atax,0.01*fgsmooth(mup_trans_lf8,4),'r');
% set(h,'LineWidth',2);
% legend('AP','LFP');
% %a1 = translen + 500;
% %errorbar(atax(a1),mup_trans_wcv(a1),sup_trans_wcv(a1),'b');
% %errorbar(atax(a1),mup_trans_lf8(a1),sup_trans_lf8(a1),'r');
% xlabel('Latency (s)');
% drawnow; shg;
% ylabel('Mean amplitude (z)')
% print -djpeg100 -r300 lfp_up_trig_lfp_spk.jpg
%
% fig = fig+1; figure(fig);
% h=plot(atax,fgsmooth(mdn_trans_spk,4),'b',atax,0.01*fgsmooth(mdn_trans_lf8,4),'r');
% set(h,'LineWidth',2);
% legend('MP','LFP');
% xlabel('Latency (s)');
% drawnow; shg;
% ylabel('Mean amplitude (z)')
% print -djpeg100 -r300 lfp_dn_trig_lfp_spk.jpg
%
%         end
% % compute all sorts of correlations between LFP and MP
%
%
% len = length(wcv);
% delay = 10000;
% m = 0;
%
% for k = 1:nseg
%     a = round(seglen*(k-1)+(1:seglen));
%     m = m+1;
%     [wclf8cor(:,m)] = xcov(wcv(a),lf8(a),delay,'coeff');
%     if (sum(spktim(a)) > 0)
%         [spkloflf8cor(:,m)] = xcov(spktim(a),loflf8(a),delay,'coeff');
%         [spkacor(:,m)] = xcov(spktim(a),delay,'coeff');
%     else
%         spkloflf8cor(:,m) = zeros(size(wclf8cor(:,1)));
%         spkacor(:,m) = spkloflf8cor(:,m);
%     end
%     [lf2lf8cor(:,m)] = xcov(lf2(a),lf8(a),delay,'coeff');
%     [wclf2cor(:,m)] = xcov(wcv(a),lf2(a),delay,'coeff');
%     [diflofwcvdifloflf8cor(:,m)] = xcov(diflofwcv(a),difloflf8(a),delay,'coeff');
%     [lofwcloflf8cor(:,m)] = xcov(lofwcv(a),loflf8(a),delay,'coeff');
%     [lofwcloflf2cor(:,m)] = xcov(lofwcv(a),loflf2(a),delay,'coeff');
%     [loflf2loflf8cor(:,m),lags] = xcov(loflf2(a),loflf8(a),delay,'coeff');
%     [wcacor(:,m)] = xcov(wcv(a),delay,'coeff');
%     [lf8acor(:,m)] = xcov(lf8(a),delay,'coeff');
%     [lofwcacor(:,m)] = xcov(lofwcv(a),delay,'coeff');
%     [loflf8acor(:,m)] = xcov(loflf8(a),delay,'coeff');
%
% end
%
% [h,siglofwcloflf8cor] = ttest(lofwcloflf8cor(delay,:));
%
% spkacor(isnan(spkacor)) = 0;
% spkloflf8cor(isnan(spkloflf8cor)) = 0;
%
% dt = median(diff(synct));
% tax = lags*dt;
%
% mlf2lf8cor = mean(lf2lf8cor,2);
% mlofwcloflf8cor = mean(lofwcloflf8cor,2);
% mlofwcloflf2cor = mean(lofwcloflf2cor,2);
% mloflf2loflf8cor = mean(loflf2loflf8cor,2);
%
%     if 1
% fig = fig+1; figure(fig);
% clf;
% hold on;
% h=plot(tax,mlofwcloflf8cor,'r');
% set(h,'LineWidth',2);
% errorbar(tax(delay),mean(lofwcloflf8cor(delay,:)),...
%     std(lofwcloflf8cor(delay,:))/sqrt(nseg),'m-');
%     xlabel('Lag (s)'); ylabel('Covariance');
% axis([-2.5 2.5 1.1*min(mlofwcloflf8cor) 1.1*max(mlofwcloflf8cor)]);
% hold off;
% box on;
% shg
% print -djpeg100 -r300 lofwc_loflf8_cor1;
%
% fig = fig+1; figure(fig);
% clf;
% hold on;
% h=plot(tax,mean(spkloflf8cor,2),'r');
% set(h,'LineWidth',2);
% errorbar(tax(delay),mean(spkloflf8cor(delay,:)),...
%     std(spkloflf8cor(delay,:))/sqrt(nseg),'m-');
%     xlabel('Lag (s)'); ylabel('Covariance');
% axis([-2.5 2.5 -Inf Inf]);
% hold off;
% box on;
% shg
% print -djpeg100 -r300 spk_loflf8_cor;
%
%      end
%
% mwclf8cor = mean(wclf8cor,2);
% [maxwclf8cor,i] = min(mwclf8cor);
% shif_wclf8cor = tax(i);
%
% mwclf2cor = mean(wclf2cor,2);
% [maxwclf2cor,i] = min(mwclf2cor);
% shif_wclf2cor = tax(i);
%
% [maxlofwcloflf8cor,i] = max(mlofwcloflf8cor);
% max_shif_lofwcloflf8cor = tax(i);
%
% [minlofwcloflf8cor,i] = min(mlofwcloflf8cor);
% min_shif_lofwcloflf8cor = tax(i);
%
% [maxlofwcloflf2cor,i] = max(mlofwcloflf2cor);
% max_shif_lofwcloflf2cor = tax(i);
%
% [minlofwcloflf2cor,i] = min(mlofwcloflf2cor);
% min_shif_lofwcloflf2cor = tax(i);
%
% mdiflofwcvdifloflf8cor = mean(diflofwcvdifloflf8cor,2);
% sdiflofwcvdifloflf8cor = std(diflofwcvdifloflf8cor,0,2);
% [maxdiflofwcvdifloflf8cor,i] = max(mdiflofwcvdifloflf8cor);
% max_shif_diflofwcvdifloflf8cor = tax(i);
% [mindiflofwcvdifloflf8cor,i] = min(mdiflofwcvdifloflf8cor);
% min_shif_diflofwcvdifloflf8cor = tax(i);
%
% [maxlofwcloflf8cor  minlofwcloflf8cor]
% [max_shif_lofwcloflf8cor min_shif_lofwcloflf8cor]
%
% nspk = sum(spktim);
% spkrate = nspk/(max(synct)-min(synct));
% siglofwcloflf8cor
%
% save lfp_wcv_summary_new mwclf8cor mlofwcloflf8cor mdiflofwcvdifloflf8cor sdiflofwcvdifloflf8cor ...
%     maxwclf8cor  maxlofwcloflf8cor  minlofwcloflf8cor maxdiflofwcvdifloflf8cor siglofwcloflf8cor ...
%     shif_wclf8cor  max_shif_lofwcloflf8cor min_shif_lofwcloflf8cor mindiflofwcvdifloflf8cor ...
%     max_shif_diflofwcvdifloflf8cor mup_trans_lf8 mup_trans_wcv mdn_trans_lf8 mdn_trans_wcv ...
%     mnormlofwcvhist mnormloflf8hist lf8minlf2 loflf8minloflf2 bimodlofwcv bimodloflf8 ...
%     bimodloflf2 mnormloflf2hist snormlofwcvhist snormloflf8hist snormloflf2hist abslf8minabslf2 ...
%     mlofwcloflf2cor mwclf2cor mloflf2loflf8cor mlf2lf8cor min_shif_diflofwcvdifloflf8cor ...
%     lf8acor wcacor loflf8acor lofwcacor spkloflf8cor spkacor mup_trans_spk mdn_trans_spk nspk spkrate
%
%
