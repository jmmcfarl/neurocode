function [theta_freqs,ltheta_phases,theta_phases] = bump_det_test_ext_f(lfp,wcv,lfp_state_seq,wcv_state_seq,hmm8,hmm,Fs,winsize,isplot)

t_axis = (1:length(lfp))/Fs;
niqf = Fs/2;

[b,a] = butter(2,[2/niqf 40/niqf]);
[b2,a2] = butter(2,[2/niqf 5/niqf]);

%extract the index values of LF8 up transitions
down_inds = [];
UDS_segs = (hmm8.UDS_segs-1)*5+1;
clfp_state_seq = ones(UDS_segs(1,1)-1,1);
for i = 1:length(lfp_state_seq)
    up_trans = UDS_segs(i,1)+find(lfp_state_seq{i}(1:end-1) == 1 & lfp_state_seq{i}(2:end) == 2);
    down_trans = UDS_segs(i,1)+find(lfp_state_seq{i}(1:end-1) == 2 & lfp_state_seq{i}(2:end) == 1);
    up_trans(up_trans < down_trans(1)) = [];
    down_trans(down_trans > up_trans(end)) = [];
    cdown_inds(1,:) = down_trans;
    cdown_inds(2,:) = up_trans;
    down_inds = [down_inds cdown_inds];
    clear cdown*
    clfp_state_seq= [clfp_state_seq; lfp_state_seq{i}];
    if length(lfp_state_seq) > i
        clfp_state_seq= [clfp_state_seq; ones(UDS_segs(i+1,1)-UDS_segs(i,2),1)];
    end

end

%create a vector containing the MP state sequence over the whole recording
UDS_segs = (hmm.UDS_segs-1)*5+1;
cwcv_state_seq = ones(UDS_segs(1,1)-1,1);
for i = 1:length(wcv_state_seq)
    cwcv_state_seq= [cwcv_state_seq; wcv_state_seq{i}];
    if length(wcv_state_seq) > i
        cwcv_state_seq= [cwcv_state_seq; ones(UDS_segs(i+1,1)-UDS_segs(i,2),1)];
    end
end


down_durs = (down_inds(2,:)-down_inds(1,:))/Fs;

min_dur = 1.0;
min_cor = 0.15;
min_pow = 0.05;
down_inds(:,down_durs < min_dur) = [];

rem_dur = round(Fs*0);
extra_dur = round(Fs*1);

F = linspace(2,4.,30);
maxlag = round(0.5*Fs);
min_lag = round(0.2*Fs);
theta_freqs = [];
theta_phases = [];

winsize = round(Fs*winsize);
num_downs = size(down_inds,2);

theta_phases = [];
ltheta_phases = [];
theta_freqs = [];

for i = 1:num_downs

    cur_down = down_inds(1,i):(down_inds(2,i));
    cur_sig = lfp(cur_down);
    cur_sig_f1 = filtfilt(b,a,cur_sig);
    %     cur_down2 = down_inds(1,i):(down_inds(2,i) - rem_dur);
    %     cur_sig2 = lfp(cur_down2);
    [Pxx,F] = pwelch(cur_sig_f1,length(cur_sig_f1),[],F,Fs);
    Pxx = Pxx/sum(Pxx);

    [maxpow,maxloc] = findpeaks(Pxx);
    [dummym,biggest] = max(maxpow);
    maxpow = maxpow(biggest);
    maxloc = maxloc(biggest);
    maxfreq = F(maxloc);
    %     [Pcor,lags] = xcov(cur_sig_f1,maxlag,'coeff');
    %     Pcor(lags < 0) = [];
    %     lags(lags < 0) = [];
    %     [maxcorr,maxloc] = findpeaks(Pcor(min_lag:end));
    %     [dummy,biggest] = max(maxcorr);
    %     if ~isempty(maxloc)
    %         maxcorr = maxcorr(biggest(1));
    %         maxloc = maxloc(biggest(1));
    %         maxloc = maxloc + min_lag-1;
    %         maxfreq = (1/(maxloc/Fs));
    %     else
    %         maxfreq = [];
    %     end

    if ~isempty(maxfreq) & maxpow >= min_pow
        [b2,a2] = butter(2,[(maxfreq-1)/niqf (maxfreq+1)/niqf]);
        cur_sig_f2 = filtfilt(b2,a2,cur_sig);
        dp = maxfreq/Fs*2*pi;

        [dummy,pklocs] = findpeaks(cur_sig_f2);
        %         pis = 2*pi*((1:length(pklocs))-1);
        %         cur_l = length(cur_sig_f2)+extra_dur;
        %         epeakphase = (cur_l-pklocs(end))*maxfreq/Fs*2*pi;
        %         pis = [pis pis(end) + epeakphase];
        %         cur_taxis = 1:cur_l;
        %         lin_phase = interp1q([pklocs(:); cur_l],pis(:),cur_taxis(:));
        %         lin_phase = mod(lin_phase,2*pi)-pi;

        hf_pklocs = [];
        for j = 1:length(pklocs)-1
            [dummy,cur_pkloc] = max(-cur_sig_f1(pklocs(j):pklocs(j+1)));
            hf_pklocs = [hf_pklocs (cur_pkloc + pklocs(j))];
        end
                if ~isempty(hf_pklocs)

        [dummy,cur_pkloc] = max(-cur_sig_f1(pklocs(end):end));
        hf_pklocs = [hf_pklocs (cur_pkloc + pklocs(end))];
        %             hf_pklocs(hf_pklocs > length(cur_down2)) = [];
        [dummy,cur_pkloc] = max(-cur_sig_f1(1:pklocs(1)));
        hf_pklocs = [cur_pkloc hf_pklocs];

        hf_pklocs = unique(hf_pklocs);
        epeak = hf_pklocs(end);
        hf_pklocs(end) = [];

        mdiff = round(mean(diff(hf_pklocs)));
        hf_pklocs = [hf_pklocs (hf_pklocs(end) + (1:2)*mdiff)];
        hf_pklocs = [1 hf_pklocs];
        %             gpeak = hf_pklocs(end) + mdiff;

        phase_down = down_inds(1,i):(down_inds(2,i)+extra_dur);
        phase_down(phase_down > length(clfp_state_seq) | phase_down > length(cwcv_state_seq)) = [];
%         hf_pklocs(hf_pklocs > length(phase_down)) = [];

        if ~isempty(hf_pklocs)
            pis = 2*pi*(1:length(hf_pklocs));
            cur_l = length(phase_down);
            cur_taxis = 1:cur_l;
            %             cur_taxis(cur_taxis > hf_pklocs(end)) = [];
            lin_phase = interp1q(hf_pklocs(:),pis(:),cur_taxis(:));
            lin_phase = mod(lin_phase,2*pi)-pi;


            %         cur_down3 = down_inds(1,i):(down_inds(2,i) - rem_dur+extra_dur);
            %        cur_down3(cur_down3 > length(cwcv_state_seq)) = [];
            last_mp_up = find(cwcv_state_seq(phase_down(1:end-1)) == 1 & cwcv_state_seq(phase_down(2:end)) == 2,1,'last');
            last_mp_up(last_mp_up > length(phase_down)) = [];
            if ~isempty(last_mp_up)
                if isnan(lin_phase(last_mp_up))
                    disp(i)
                end
                theta_phases = [theta_phases lin_phase(last_mp_up)];
                theta_freqs = [theta_freqs maxfreq];
            end
            last_lfp_up = find(clfp_state_seq(phase_down(1:end-1)) == 1 & clfp_state_seq(phase_down(2:end)) == 2,1,'last');
            last_lfp_up(last_lfp_up > length(phase_down)) = [];
            if ~isempty(last_lfp_up)
                ltheta_phases = [ltheta_phases lin_phase(last_lfp_up)];
            end

            if isplot
                plot_sig = (down_inds(1,i) - winsize):(down_inds(2,i)+winsize);
                plot_sig(plot_sig < 1 | plot_sig > length(lfp)) = [];
                cur_sig = lfp(phase_down);
                cur_sig_f1 = filtfilt(b,a,cur_sig);

                %                         subplot(3,1,[1 2])
                %             [bt,at] = butter(2,10/niqf,'low');
                %             chunk = -filtfilt(b,a,lfp(plot_sig));
                %             [dummy,lmindiff] = findpeaks(chunk(1:winsize+down_inds(2,i)-down_inds(1,i)));
                %             temploc = max(lmindiff);
                plot(t_axis(plot_sig),lfp(plot_sig),'k','linewidth',2), hold on
                %             title(num2str(maxcorr))
                plot(t_axis(plot_sig),wcv(plot_sig))
                %             plot(t_axis(down_inds(2,i)),lfp(down_inds(2,i)),'ro','markersize',12,'linewidth',2)
                %             plot(t_axis(plot_sig(temploc)),lfp(plot_sig(temploc)),'go','markersize',12,'linewidth',2)
                plot(t_axis(phase_down(last_mp_up)),lin_phase(last_mp_up)/10,'ro','markersize',12,'linewidth',2)
                plot(t_axis(phase_down(last_lfp_up)),lin_phase(last_lfp_up)/10,'go','markersize',12,'linewidth',2)
                plot(t_axis(phase_down),cur_sig_f1*2+0.1,'r')
                %             plot(t_axis(cur_down2),cur_sig_f2*2,'g')
                plot(t_axis(phase_down(hf_pklocs)),cur_sig_f1(hf_pklocs)*2+0.1,'k.')
                plot(t_axis(phase_down(hf_pklocs(end-1:end))),cur_sig_f1(hf_pklocs(end-1:end))*2+0.1,'g.')
                %             plot(t_axis(cur_down(epeak)),cur_sig_f1(epeak)*2+0.1,'b.')
                %             plot(t_axis(cur_down(gpeak)),cur_sig_f1(gpeak)*2+0.1,'r.')
                plot(t_axis(phase_down),lin_phase/20,'g')
                xlim([t_axis(plot_sig(1)) t_axis(plot_sig(end))])
                %                         subplot(3,1,3)
                % %                         plot(lags/Fs,Pcor), hold on
                % %                                      subplot(4,1,4)
                %            plot(F,Pxx,'r')
                %             %             hold on
                %             %             plot(lags(maxloc)/Fs,Pcor(maxloc),'ro')
                %             % %             plot(F,Pxx), hold on
                % %             plot(F(maxloc),Pxx(maxloc),'ro')
                pause
                clf
            end
        end
                end
    end
end

theta_phases(isnan(theta_phases)) = [];
ltheta_phases(isnan(ltheta_phases)) = [];