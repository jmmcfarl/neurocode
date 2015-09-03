clear all
close all
%%
load G:\WC_Germany\overall_EC\overall_allcells_dir
addpath('G:\Code\WC_anal\general\')
addpath('G:\WC_Germany\Overall_EC\')
addpath('G:\Code\Chronux\spectral_analysis\continuous\')
addpath('G:\WC_Germany\hsmm_state_detection\')

drive_letter = 'G';

%%
dsf = 32;
Fsd = 2016/dsf;
minSegLength = 60;
maxLag = 5*Fsd;
niqf = 2016/2;
lcf = .05/niqf;
hcf = 2/niqf;
[b,a] = butter(2,[lcf hcf]);

for d = 1:200
    
    cdir = sess_data(d).directory;
    cdir(1) = 'G';
    disp(sprintf('session %d',d))
    cd(cdir);
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);
    load used_data lf8 lf3 
    
    lf8_d = filtfilt(b,a,lf8);
    lf8_d = downsample(lf8_d,dsf);
    lf3_d = filtfilt(b,a,lf3);
    lf3_d = downsample(lf3_d,dsf);
    
    t_axis = (1:length(lf3_d))/Fsd;
    load ec_hmm_state_seq8
    sMarkers = resample_uds_seg_inds(hmm8.UDS_segs,hmm8.Fs,Fsd,length(lf8_d));
    seg_durs = diff(sMarkers')'/Fsd;

    cnt = 0;
    for i = 1:size(sMarkers,1)
        cnt = cnt+1;
        [lf8_acorr(cnt,:),lags] = xcov(lf8_d(sMarkers(i,1):sMarkers(i,2)),maxLag,'coeff');
        [lf3_acorr(cnt,:),lags] = xcov(lf3_d(sMarkers(i,1):sMarkers(i,2)),maxLag,'coeff');
        [x83(cnt,:),lags] = xcov(lf8_d(sMarkers(i,1):sMarkers(i,2)),lf3_d(sMarkers(i,1):sMarkers(i,2)),maxLag,'coeff');
    end
    
    weighting = seg_durs/sum(seg_durs);
    tot_lf8_acorr(d,:) = sum(lf8_acorr.*repmat(weighting(:),1,length(lags)),1);
    tot_lf3_acorr(d,:) = sum(lf3_acorr.*repmat(weighting(:),1,length(lags)),1);
    tot_83_x(d,:) = sum(x83.*repmat(weighting(:),1,length(lags)),1);

    clear lf8* lf3* x83

end

%%
cd G:\WC_Germany\overall_EC\
save overall_allcell_timedomain tot_* lags Fsd
%%
nolf3nans = find(~any(isnan(tot_83_x),2));
[COEFF, SCORE, LATENT] = princomp(tot_83_x(nolf3nans,:));
obj = gmdistribution.fit(SCORE(:,1:3),2);
idx = cluster(obj,SCORE(:,1:3));

[temp,temploc] = max(abs(tot_83_x),[],2);
for i = 1:length(temp)
   maxval(i) = tot_83_x(i,temploc(i)); 
end
good_lf3 = find(maxval > 0);
bad_lf3 = setdiff(1:200,good_lf3);

ca1 = find_struct_field_vals(sess_data,'region','ca1');
ca3 = find_struct_field_vals(sess_data,'region','ca3');
dg = find_struct_field_vals(sess_data,'region','dg');
pyr = find_struct_field_vals(sess_data,'cell_type','pyr');
int = find_struct_field_vals(sess_data,'cell_type','int');
ca1pyr = intersect(ca1,pyr);
ca1int = intersect(ca1,int);

good_ca1pyr = intersect(good_lf3,ca1pyr);
good_ca1int = intersect(good_lf3,ca1int);
good_ca3 = intersect(good_lf3,ca3);
good_dg = intersect(good_lf3,dg);

