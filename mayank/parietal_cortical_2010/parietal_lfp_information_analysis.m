clear all
close all
cd G:\WC_Germany\parietal_cortical_2010

load F:\WC_Germany\parietal_cortical_2010\parietal_cortical_2010
load F:\WC_Germany\parietal_cortical_2010\desynch_times_mp_lf8
addpath('F:\WC_Germany\parietal_cortical_2010')

raw_Fs = 2016;
dsf = 16;
Fsd = raw_Fs/dsf;

frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');

sess_data = sess_data(parietal);
desynch_start_times = desynch_start_times(parietal);
desynch_stop_times = desynch_stop_times(parietal);
%get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];
desynch_start_times(interneurons) = [];
desynch_stop_times(interneurons) = [];

d=1;
% save lfp_information_data_unsm
while  d <= length(sess_data)
    
    disp(sprintf('session %d',d))
    cd(sess_data(d).directory);

    load used_data wcv_minus_spike lf8
    load hsmm_state_seq_lf_pert15

    [lf8_w,new_t,Fsd] = get_lf_features(lf8,raw_Fs,Fsd,[0.05 40]);
    [lf8_lf,new_t,Fsd] = get_lf_features(lf8,raw_Fs,Fsd,[0.05 5]);
    [hf_features,t_axis,Fs] = get_hf_features(lf8,raw_Fs,Fsd,[20 80]);
    [scal_features,wfreqs,new_t,Fsd] = get_full_scalogram_features(lf8,raw_Fs,Fsd);

    mp_state_seq =  hsmm_bbstate_seq;
    cdsf = round(252/Fsd);
    mp_state_seq = downsample(mp_state_seq,cdsf);
    %     old_t = (1:length(hsmm_bbstate_seq))/Fs_bb;
    %     mp_state_seq = round(interp1(old_t,mp_state_seq,new_t));

    bad_samps = [];
    for i = 1:length(desynch_start_times{d})
        bad_start = find(new_t >= desynch_start_times{d}(i),1,'first');
        bad_stop = find(new_t >= desynch_stop_times{d}(i),1,'first');
        bad_samps = [bad_samps bad_start:bad_stop];
    end
    bad_samps = [bad_samps find(isnan(mp_state_seq))];
    
    lf8_w(bad_samps) = [];
    lf8_lf(bad_samps) = [];
    hf_features(bad_samps) = [];
    scal_features(bad_samps,:) = [];    
    mp_state_seq(bad_samps) = [];
    lscal_features = log(scal_features + 1e-5);
    
%%
    mp_state_seq = mp_state_seq(:);
    [lf_response, nt] = buildr(mp_state_seq, lf8_lf);
    [w_response, nt] = buildr(mp_state_seq, lf8_w);
    [hf_response, nt] = buildr(mp_state_seq, hf_features);
%     opts.nt = nt;
%     opts.method = 'gs';
%     opts.bias = 'gsb';
%     info_lf(d) = information(lf_response,opts,'I');
%     info_w(d) = information(w_response,opts,'I');
%     info_hf(d) = information(hf_response,opts,'I');

opts.method = 'dr';
    opts.bias = 'pt';
    opts.nt = nt;
    
    lscal_features = [lscal_features lf8_lf];
    scal_response = zeros(length(wfreqs)+1,max(nt),2);
    scal_response(:,1:nt(1),1) = lscal_features(mp_state_seq==1,:)';
    scal_response(:,1:nt(2),2) = lscal_features(mp_state_seq==2,:)';   

    joint_info_scal{d} = zeros(length(wfreqs),length(wfreqs));
    joint_info_syn{d} = zeros(length(wfreqs),length(wfreqs));
    info_scal(d,:) = zeros(length(wfreqs),1);
    joint_info_lf_scal(d,:) = zeros(length(wfreqs),1);
    joint_info_lf_syn(d,:) = zeros(length(wfreqs),1);
    for i = 1:length(wfreqs)
        for j = i+1:length(wfreqs)
            cf = [i j];
            nb = 50;
            scal_response_d = binr(scal_response(cf,:,:),nt,nb,'eqpop');
            [joint_info_scal{d}(i,j),joint_info_syn{d}(i,j)] = ...
                information(scal_response_d,opts,'I','Syn');
            joint_info_scal{d}(j,i) = joint_info_scal{d}(i,j);
            joint_info_syn{d}(j,i) = joint_info_syn{d}(i,j);
        end
        nb = 50;
        scal_response_d = binr(scal_response(i,:,:),nt,nb,'eqpop');
        info_scal(d,i) = information(scal_response_d,opts,'I');
        nb = 50;
        cf = [i length(wfreqs)+1];
        scal_response_d = binr(scal_response(cf,:,:),nt,nb,'eqpop');
        [joint_info_lf_scal(d,i),joint_info_lf_syn(d,i)] = ...
            information(scal_response_d,opts,'I','Syn');
    end

    %%
    nb = 50;
    opts.method = 'dr';
    opts.bias = 'pt';
    opts.nt = nt;
    lf_response = binr(lf_response, nt, nb, 'eqpop');
    info_lf_bin(d) = information(lf_response,opts,'I');
    w_response = binr(w_response,nt,nb,'eqpop');
    info_w_bin(d) = information(w_response,opts,'I');
    hf_response = binr(hf_response,nt,nb,'eqpop');
    info_hf_bin(d) = information(hf_response,opts,'I');
    
%     info_scal_bin(d,:) = zeros(length(wfreqs),1);
%     for i = 1:length(wfreqs)
%         scal_response_d = binr(scal_response(i,:,:),nt,nb,'eqpop');
%         info_scal_bin(d,i) = information(scal_response_d,opts,'I');
%     end
    
%     lscal_features = downsample(lscal_features',3)';
%     dwfreqs = downsample(wfreqs,3);
%     scal_response = zeros(length(dwfreqs),max(nt),2);
%     scal_response(:,1:nt(1),1) = lscal_features(mp_state_seq==1,:)';
%     scal_response(:,1:nt(2),2) = lscal_features(mp_state_seq==2,:)';
%     opts.method = 'gs';
%     opts.bias = 'gsb';
%     total_scal_info(d) = information(scal_response,opts,'I');
%     opts.method = 'dr';
%     opts.bias = 'pt';
%     nb = 8;
%     scal_response = binr(scal_response,nt,nb,'eqpop');
%     total_scal_info_bin(d) = information(scal_response,opts,'I');

    d = d+1;
 cd G:\WC_Germany\parietal_cortical_2010\
save lfp_information_data_unsm
clear all   
load lfp_information_data_unsm
end


%%
mean_ji_scal = zeros(length(wfreqs),length(wfreqs));
mean_ji_syn = zeros(length(wfreqs),length(wfreqs));
for i = 1:length(sess_data)
    mean_ji_scal = mean_ji_scal + joint_info_scal{i};
    mean_ji_syn = mean_ji_syn + joint_info_syn{i};
end
mean_ji_scal = mean_ji_scal/length(sess_data);
mean_ji_syn = mean_ji_syn/length(sess_data);

%%
figure
pcolor(wfreqs,wfreqs,mean_ji_scal);shading flat
colorbar
shg
caxis([0.4 0.65])
figure
pcolor(wfreqs,wfreqs,mean_ji_syn);shading flat
colorbar
shg
figure
errorbar(wfreqs,mean(info_scal),std(info_scal)/sqrt(21))
shg
figure
errorbar(wfreqs,mean(joint_info_lf_scal),std(joint_info_lf_scal)/sqrt(21))
figure
errorbar(wfreqs,mean(joint_info_lf_syn),std(joint_info_lf_syn)/sqrt(21))