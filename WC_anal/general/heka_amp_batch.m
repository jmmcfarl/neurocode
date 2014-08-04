
clear all
cd 'F:\EC_MPasci'
load E:/WC_Germany/jmM_Analysis/used_sessions_heka.mat
dsf = 10;
Fs = 5000;
% data_down = detrend(data_down);

Fsd = Fs/dsf;
seg_length = 10*Fsd;
niqf = Fsd/2;
hif = 20/niqf;

[b,a] = butter(2,hif,'low');

for d = 1:length(used_sess)
    load(used_sess{d},[used_sess{d} '_MP'])
    eval(['data = ',used_sess{d},'_MP;']);
    if d < 16
        f_names{d} = ['stell_' used_sess{d}];
    else
        f_names{d} = ['pyr_' used_sess{d}];
    end
    clear A*
    data_down = downsample(data,dsf);

    data_filt = filtfilt(b,a,data_down);


    num_segs = floor(length(data_down)/seg_length);

    for i = 1:num_segs

        begPt = (i-1)*seg_length+1;
        endPt = begPt+seg_length-1;

        dataSeg = data_filt(begPt:endPt);

        [u(i,:),sig(i,:),t(i,:),iter] = fit_mix_gaussian(dataSeg,2);

    end

    mu = mean(u);

    detu = detrend(u);

    detu(:,1) = detu(:,1)+mu(1);
    detu(:,2) = detu(:,2)+mu(2);

    amp_Down{d} = detu(:,1);
    amp_Up{d} = detu(:,2);
    sig_Down{d} = sig(:,1);
    sig_Up{d} = sig(:,2);
    
        Fig = figure(1)
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    errorbar([1:num_segs],detu(:,1),sig(:,1),'b')
    hold on
    errorbar([1:num_segs],detu(:,2),sig(:,2),'r')
    legend('Down State Amp','Up State Amp')
    xlabel('Time (10 s)','FontSize',14)
    ylabel('Amplitude (.1 V)','FontSize',14)
    tname = ['E:\WC_Germany\JMM_Analysis\heka_amp\test' f_names{d}];
    print('-dpng',tname);

    clear data* u sig t mu num_segs detu
    
%     save E:\WC_Germany\JMM_Analysis\amp_data amp* sig*
    
    
    
    
end