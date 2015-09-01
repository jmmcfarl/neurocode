function [theta_freqs,theta_phases] = bump_det_test(lfp,wcv,lfp_state_seq,wcv_state_seq,hmm,hmm2,Fs,winsize,isplot)

t_axis = (1:length(lfp))/Fs;
niqf = Fs/2;

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

F = linspace(2,4.5,30);
minf = find(F > 2,1,'first');
min_pow = 1e-3;
min_dur = 1.;

extra_dur = round(Fs*0.1);
% rem_dur = 0;

theta_freqs = [];
theta_phases = [];
% all_theta_phases = [];

winsize = round(Fs*winsize);
num_downs = size(down_inds,2);
for i = 1:num_downs
%     fprintf('%d of %d\n',i,num_downs)

    cur_down = down_inds(1,i):(down_inds(2,i));
    cur_sig = detrend(lfp(cur_down));
    [Pxx,F] = pwelch(cur_sig,length(cur_sig),[],F,Fs);
    [maxpow,maxloc] = findpeaks(Pxx);
    maxpow(maxloc < minf) = [];
    maxloc(maxloc < minf) = [];
    [dummym,biggest] = max(maxpow);
    maxpow = maxpow(biggest);
    maxloc = maxloc(biggest);
    maxfreq = F(maxloc);
    if maxpow > min_pow & length(cur_down)/Fs > min_dur & ~isempty(maxfreq)
        [b,a] = butter(2,[(maxfreq-1)/niqf (maxfreq+1)/niqf]);
        cur_sig_f = filtfilt(b,a,cur_sig);
        cur_l = length(cur_sig_f) + extra_dur;
        [pkvals,pklocs] = findpeaks(-cur_sig_f);
        epeakphase = (cur_l-pklocs(end))*maxfreq/Fs*2*pi;
        [ppkvals,ppklocs] = findpeaks(cur_sig_f);
        ppklocs(ppklocs < pklocs(1)) = [];
        pklocs = sort([pklocs ppklocs]);
%         lin_phase = zeros(length(cur_sig_f),1);
%         dp = 2*pi*maxfreq/Fs;
%         for j = 1:length(pklocs)-1
%             cur_l = pklocs(j+1)-pklocs(j)+1;
%             lin_phase(pklocs(j):pklocs(j+1)) = linspace(-pi,pi,cur_l);
%         end
%         cur_l = length(lin_phase)-pklocs(end);
%         lin_phase(pklocs(end)+1:end) = -pi:dp:(-pi+dp*cur_l-dp);
        
        pis = pi*((1:length(pklocs))-1);
        pis = [pis pis(end) + epeakphase];
        cur_taxis = 1:cur_l;
        lin_phase = interp1q([pklocs(:); cur_l],pis(:),cur_taxis(:));
        lin_phase = mod(lin_phase,2*pi)-pi;
        
        
%         lfp_analytic = hilbert(cur_sig_f);
%         lfp_phase = angle(lfp_analytic);
%         extra_phase = (lfp_phase(end)+pi):dp:(lfp_phase(end)+pi+extra_dur*dp);
%         extra_phase = mod(extra_phase,2*pi)-pi;
%         lfp_phase = [lfp_phase; extra_phase(:)];
        theta_freqs = [theta_freqs F(maxloc)];
%         all_theta_phases = [all_theta_phases; lfp_phase(1:end-extra_dur)];
        ncur_down = down_inds(1,i):(down_inds(2,i)+extra_dur);
%         lfp_phase(length(ncur_down)+1:end) = [];
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
        plot(t_axis(plot_sig),cwcv_state_seq(plot_sig),'r')
        if maxpow > min_pow & length(cur_down)/Fs > min_dur & ~isempty(maxfreq)
            plot(t_axis(ncur_down),lin_phase/2,'g')
            if ~isempty(last_mp_up)
                plot(t_axis(ncur_down(last_mp_up)),lin_phase(last_mp_up)/2,'ro')
            end
            plot(t_axis(cur_down),cur_sig_f*3,'c')
        end
        plot(t_axis(cur_down),lfp(cur_down),'r')
        axis tight
        subplot(3,1,3)
        plot(F,Pxx)
        hold on
        plot(F(maxloc),Pxx(maxloc),'ro')
        axis tight
        pause
        clf
    end
end

