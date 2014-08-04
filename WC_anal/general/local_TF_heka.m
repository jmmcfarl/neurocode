load('C:\WC_Germany\JMM_analysis_ste\stellate_heka_dir.mat')
% load('C:\WC_Germany\JMM_analysis_pyr\pyr_heka_dir.mat')

d = 5;
load(f_loc{d})
dat_name = f_loc{d}(25:end);

eval(['data = ' dat_name '_MP;'])
eval(['time = ' dat_name '_sampletimes;'])
Fs = 5000;
niqf = Fs/2;
[b,a] = butter(2,12/niqf,'low');
[b2,a2] = butter(2,[7/niqf 20/niqf]);

sweep_t = find(diff(time) < 0);
low = -0.8;
high = -0.4;

%initialize
deltaT_f1 = [];
deltaA_f1 = [];
deltaM_f1 = [];
deltaT_f2 = [];
deltaA_f2 = [];
deltaM_f2 = [];

% for i = 1:length(sweep_t)-1
for i = 1:50
    dataseg = data(sweep_t(i)+1:sweep_t(i+1));
    filtseg = filtfilt(b,a,dataseg);
    filtseg2 = filtfilt(b2,a2,dataseg);

%% optional plotting
    plot(time(sweep_t(i)+1:sweep_t(i+1)),data(sweep_t(i)+1:sweep_t(i+1)))
    hold on
    plot(time(sweep_t(i)+1:sweep_t(i+1)),filtseg,'r','linewidth',2)
%     plot(time(sweep_t(i)+1:sweep_t(i+1)),filtseg2*3+mean(dataseg),'g','linewidth',2)
    ylim([low high])
    grid
    pause
    clf
%%

    %for filtered signal 1
    [u_hts_f1,u_lcs_f1] = findpeaks(filtseg);
    [d_hts_f1,d_lcs_f1] = findpeaks(-filtseg);
    d_hts_f1 = -d_hts_f1;
    npeaks = length(u_hts_f1)+length(d_hts_f1);
    cur_pk_hts = zeros(npeaks,1);
    cur_pk_lcs = zeros(npeaks,1);

    %if first peak is upwards
    if u_lcs_f1(1) < d_lcs_f1(1)
        odds = linspace(1,2*length(u_hts_f1)-1,length(u_hts_f1));
        evens = linspace(2,2*length(d_hts_f1),length(d_hts_f1));
        cur_pk_hts(odds) = u_hts_f1;
        cur_pk_hts(evens) = d_hts_f1;
        cur_pk_lcs(odds) = u_lcs_f1;
        cur_pk_lcs(evens) = d_lcs_f1;
    else
        %if first peak is downwards
        odds = linspace(1,2*length(d_hts_f1)-1,length(d_hts_f1));
        evens = linspace(2,2*length(u_hts_f1),length(u_hts_f1));
        cur_pk_hts(evens) = u_hts_f1;
        cur_pk_hts(odds) = d_hts_f1;
        cur_pk_lcs(evens) = u_lcs_f1;
        cur_pk_lcs(odds) = d_lcs_f1;
    end

    cur_deltaT = diff(cur_pk_lcs)/5000;
    cur_deltaA = diff(cur_pk_hts);

    %now find average signal inbetween each peak
    cur_deltaM = zeros(length(cur_pk_lcs)-1,1);
    for p = 1:length(cur_pk_lcs)-1
        cur_deltaM(p) = mean(dataseg(cur_pk_lcs(p):cur_pk_lcs(p+1)));
    end

deltaT_f1 = [deltaT_f1;cur_deltaT];
deltaA_f1 = [deltaA_f1;cur_deltaA];
deltaM_f1 = [deltaM_f1;cur_deltaM];


    %for filtered signal 2
    [u_hts_f2,u_lcs_f2] = findpeaks(filtseg2);
    [d_hts_f2,d_lcs_f2] = findpeaks(-filtseg2);
    d_hts_f2 = -d_hts_f2;
    npeaks = length(u_hts_f2)+length(d_hts_f2);
    cur_pk_hts = zeros(npeaks,1);
    cur_pk_lcs = zeros(npeaks,1);

    %if first peak is upwards
    if u_lcs_f2(1) < d_lcs_f2(1)
        odds = linspace(1,2*length(u_hts_f2)-1,length(u_hts_f2));
        evens = linspace(2,2*length(d_hts_f2),length(d_hts_f2));
        cur_pk_hts(odds) = u_hts_f2;
        cur_pk_hts(evens) = d_hts_f2;
        cur_pk_lcs(odds) = u_lcs_f2;
        cur_pk_lcs(evens) = d_lcs_f2;
    else
        %if first peak is downwards
        odds = linspace(1,2*length(d_hts_f2)-1,length(d_hts_f2));
        evens = linspace(2,2*length(u_hts_f2),length(u_hts_f2));
        cur_pk_hts(evens) = u_hts_f2;
        cur_pk_hts(odds) = d_hts_f2;
        cur_pk_lcs(evens) = u_lcs_f2;
        cur_pk_lcs(odds) = d_lcs_f2;
    end

    cur_deltaT = diff(cur_pk_lcs)/5000;
    cur_deltaA = diff(cur_pk_hts);

    %now find average signal inbetween each peak
    cur_deltaM = zeros(length(cur_pk_lcs)-1,1);
    for p = 1:length(cur_pk_lcs)-1
        cur_deltaM(p) = mean(dataseg(cur_pk_lcs(p):cur_pk_lcs(p+1)));
    end

deltaT_f2 = [deltaT_f2;cur_deltaT];
deltaA_f2 = [deltaA_f2;cur_deltaA];
deltaM_f2 = [deltaM_f2;cur_deltaM];

i
    
end





