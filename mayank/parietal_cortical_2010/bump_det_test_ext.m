function [theta_freqs,theta_phases] = bump_det_test_ext(lfp,wcv,lfp_state_seq,wcv_state_seq,hmm,hmm2,Fs,winsize,isplot)

t_axis = (1:length(lfp))/Fs;
niqf = Fs/2;
maxlag = round(Fs*0.5);

%extract the index values of LF8 up transitions
down_inds = [];
UDS_segs = (hmm.UDS_segs-1)*5+1;
for i = 1:length(lfp_state_seq)
    up_trans = UDS_segs(i,1)+find(lfp_state_seq{i}(1:end-1) == 1 & lfp_state_seq{i}(2:end) == 2);
    down_trans = UDS_segs(i,1)+find(lfp_state_seq{i}(1:end-1) == 2 & lfp_state_seq{i}(2:end) == 1);
    up_trans(up_trans < down_trans(1)) = [];
    down_trans(down_trans > up_trans(end)) = [];
    cdown_inds(1,:) = down_trans;
    cdown_inds(2,:) = up_trans;
    down_inds = [down_inds cdown_inds];
    clear cdown*
end

%create a vector containing the MP state sequence over the whole recording
UDS_segs = (hmm2.UDS_segs-1)*5+1;
cwcv_state_seq = ones(UDS_segs(1,1)-1,1);
for i = 1:length(wcv_state_seq)
    cwcv_state_seq= [cwcv_state_seq; wcv_state_seq{i}];
    if length(wcv_state_seq) > i
        cwcv_state_seq= [cwcv_state_seq; ones(UDS_segs(i+1,1)-UDS_segs(i,2),1)];
    end
end

down_durs = (down_inds(2,:)-down_inds(1,:))/Fs;

F = linspace(2,20,60);
minf = find(F > 2,1,'first');
% min_pow = 1e-3;
min_dur = 1.;
min_cor = 0.15;
min_lag = round(1/5*Fs);

down_inds(:,down_durs < min_dur) = [];

extra_dur = round(Fs*0.2);
% rem_dur = round(Fs*0.15);
rem_dur = 0;

theta_freqs = [];
theta_phases = [];
% all_theta_phases = [];

winsize = round(Fs*winsize);
num_downs = size(down_inds,2);
for i = 1:num_downs
%     fprintf('%d of %d\n',i,num_downs)
    cur_down = (down_inds(1,i)+rem_dur):(down_inds(2,i) - rem_dur);
    cur_sig = wcv(cur_down);
%     [Pxx,F] = pwelch(cur_sig,length(cur_sig),[],F,Fs);
%     [maxpow,maxloc] = findpeaks(Pxx);
%     maxpow(maxloc < minf) = [];
%     maxloc(maxloc < minf) = [];
%     [dummym,biggest] = max(maxpow);
%     maxpow = maxpow(biggest);
%     maxloc = maxloc(biggest);
%     maxfreq = F(maxloc);
%     Pxx = Pxx/sum(Pxx);
%     Pent = sum(-Pxx.*log(Pxx));

[Pcor,lags] = xcov(cur_sig,maxlag,'coeff');
Pcor(lags < 0) = [];
lags(lags < 0) = [];
[maxcorr,maxloc] = findpeaks(Pcor(min_lag:end));
[dummy,biggest] = max(maxcorr);
if ~isempty(maxloc)
maxcorr = maxcorr(biggest(1));
maxloc = maxloc(biggest(1));
maxloc = maxloc + min_lag-1;
maxfreq = (1/(maxloc/Fs));
else
    maxfreq = [];
end
%     if maxpow > min_pow & length(cur_down)/Fs > min_dur & ~isempty(maxfreq)
    if maxcorr > min_cor & length(cur_down)/Fs > min_dur & ~isempty(maxfreq)
        [b,a] = butter(2,[(maxfreq-1)/niqf (maxfreq+1)/niqf]);
%         [b,a] = butter(2,[1.5/niqf 4/niqf]);
        cur_sig_f = filtfilt(b,a,cur_sig);
        cur_l = length(cur_sig_f) + extra_dur;
%         [dpkvals,dpklocs] = findpeaks(-cur_sig_f);
        [uppkvals,uppklocs] = findpeaks(cur_sig_f);
%         uppklocs(uppklocs < dpklocs(1)) = [];
%         pklocs = sort([dpklocs uppklocs]);
      pklocs = uppklocs;
        hf_pklocs = [];
        cur_hf_sig = -wcv(cur_down);
        for j = 1:length(pklocs)-1           
            [dummy,cur_pkloc] = max(cur_hf_sig(pklocs(j):pklocs(j+1)));
            hf_pklocs = [hf_pklocs (cur_pkloc + pklocs(j))];                 
        end
            [dummy,cur_pkloc] = max(cur_hf_sig(pklocs(end):end));
            hf_pklocs = [hf_pklocs (cur_pkloc + pklocs(end))];
            hf_pklocs(hf_pklocs > length(cur_down)) = [];
            [dummy,cur_pkloc] = max(cur_hf_sig(1:pklocs(1)));
            hf_pklocs = [cur_pkloc hf_pklocs];
            
% %             hf_pklocs(end) = [];
%             mdiff = mean(diff(hf_pklocs));
%             hf_pklocs = round([hf_pklocs hf_pklocs(end) + mdiff]);
%             hf_pklocs = round([hf_pklocs hf_pklocs(end) + mdiff]);
            
            hf_pklocs(hf_pklocs > length(cur_down) + extra_dur) = [];
            
        epeakphase = (cur_l-pklocs(end))*maxfreq/Fs*2*pi;
     
        pis = pi*((1:length(pklocs))-1);
        pis = [pis pis(end) + epeakphase];
        cur_taxis = 1:cur_l;
        lin_phase = interp1q([pklocs(:); cur_l],pis(:),cur_taxis(:));
        lin_phase = mod(lin_phase,2*pi)-pi;
        
        theta_freqs = [theta_freqs maxfreq];
        ncur_down = (down_inds(1,i)+rem_dur):(down_inds(2,i)+extra_dur - rem_dur);
        if max(cur_down) < length(cwcv_state_seq)
            last_mp_up = find(cwcv_state_seq(ncur_down(1:end-1)) == 1 & cwcv_state_seq(ncur_down(2:end)) == 2,1,'last');
            if ~isempty(last_mp_up)
                theta_phases = [theta_phases lin_phase(last_mp_up)];
            end
        end
    end

    if isplot
        plot_sig = (down_inds(1,i) - winsize):(down_inds(2,i)+winsize);
        plot_sig(plot_sig < 1 | plot_sig > length(lfp)) = [];
        subplot(3,1,[1 2])
        plot(t_axis(plot_sig),lfp(plot_sig),'k','linewidth',2), hold on
        plot(t_axis(plot_sig),wcv(plot_sig),'b')
        
%         if maxpow > min_pow & length(cur_down)/Fs > min_dur & ~isempty(maxfreq)
    if maxcorr > min_cor & length(cur_down)/Fs > min_dur & ~isempty(maxfreq)
            if ~isempty(hf_pklocs)
            plot(t_axis(ncur_down(hf_pklocs)),wcv(ncur_down(hf_pklocs)),'r.','markersize',16)
%             plot(t_axis(ncur_down(hf_pklocs(end-1))),wcv(ncur_down(hf_pklocs(end-1))),'g.','markersize',16)
            plot(t_axis(ncur_down(hf_pklocs(end))),wcv(ncur_down(hf_pklocs(end))),'g.','markersize',16)
            end
        end
            
        plot(t_axis(plot_sig),cwcv_state_seq(plot_sig),'r')
%         if maxpow > min_pow & length(cur_down)/Fs > min_dur & ~isempty(maxfreq)
    if maxcorr > min_cor & length(cur_down)/Fs > min_dur & ~isempty(maxfreq)
%             plot(t_axis(ncur_down),lin_phase/2,'g')
            if ~isempty(last_mp_up)
                plot(t_axis(ncur_down(last_mp_up)),lin_phase(last_mp_up)/2,'ro')
            end
            plot(t_axis(cur_down),cur_sig_f*3,'c')
%             plot(t_axis(cur_down(pklocs)),cur_sig_f(pklocs),'r.')
        end
        plot(t_axis(cur_down),lfp(cur_down),'r')
        axis tight
        subplot(3,1,3)
%         plot(F,Pxx)
%         title(num2str(Pent))
plot(lags/Fs,Pcor)
        hold on
%         plot(F(maxloc),Pxx(maxloc),'ro')
plot(lags(maxloc)/Fs,Pcor(maxloc),'ro')
        axis tight
        pause
        clf
    end
end

