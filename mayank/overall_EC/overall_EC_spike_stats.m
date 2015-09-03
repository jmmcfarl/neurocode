
clear all
close all

load F:\WC_Germany\overall_EC\overall_EC_dir.mat
drive_letter = 'F';
dsf = 8;
Fsd = 2016/dsf;

isi_bins = logspace(log10(2),log10(1e3),50);

for d = 1:length(sess_data)
    fprintf('session %d\n',d)
    cdir = sess_data(d).directory;
    cdir(1) = drive_letter;
    cd(cdir)
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);

    load spike_time_jmm 
    spkid = round(spkid/dsf);
    load ec_hmm_state_seq
    mp_state_seq = hmm_bbstate_seq;
    UDS_segs_d = round(hmm.UDS_segs*5);
    seg_dur = (UDS_segs_d(:,2) - UDS_segs_d(:,1))/Fsd;
    T = sum(seg_dur);
    
    used_spikes = [];
    up_spikes = [];
    up_t = 0;
    for i = 1:hmm.Nsegs
        cur_spikes = spkid(spkid >= UDS_segs_d(i,1) & spkid <= UDS_segs_d(i,2));
       used_spikes = [used_spikes; cur_spikes(:)]; 
       cur_spikes_rel = cur_spikes - UDS_segs_d(i,1)+1; 
       cur_up_spikes = find(mp_state_seq{i}(cur_spikes_rel) == 2);
       up_spikes = [up_spikes; cur_spikes(cur_up_spikes(:))];
       up_t = up_t + length(find(mp_state_seq{i}==2))/Fsd;
    end
    
    overall_rate(d) = length(used_spikes)/T;
    up_rate(d) = length(up_spikes)/up_t;
    
    isis = diff(spkid/Fsd)*1e3;
    isi_h = histc(isis,isi_bins);
    
%     figure('visible','off')
%     stairs(isi_bins,isi_h)
%     set(gca,'xscale','log'), xlim([2 1000])
%     t_names = ['F:\WC_Germany\overall_EC\isi_dist\' s_name];
%     print(t_names,'-dpng')
%     close
    
end

cd F:\WC_Germany\overall_EC\
save overall_EC_spike_data overall_rate up_rate