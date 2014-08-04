%% cortical pyramidal neurons
cd C:\WC_Germany\Cortical_analysis
clear all


%% L2/3 Pyrmidal Parietal  (N=11)
cur_tot = 0;
n_23_pyr_par = 11;

sess_data(cur_tot+1).directory = 'D:\wc_data\2004-11-04_CWC_LFP';
sess_data(cur_tot+2).directory = 'D:\wc_data\2005-06\2005-12-05_CWC_LFP_C';
sess_data(cur_tot+3).directory = 'D:\wc_data\2005-06\2005-12-06_CWC_LFP_A';
sess_data(cur_tot+4).directory = 'D:\wc_data\2005-06\2005-12-06_CWC_LFP_B';
sess_data(cur_tot+5).directory = 'D:\wc_data\2005-06\2005-12-13_CWC_LFP_A';
sess_data(cur_tot+6).directory = 'D:\wc_data\2005-06\2005-12-13_CWC_LFP_D';
sess_data(cur_tot+7).directory = 'D:\wc_data\2005-06\2005-12-14_CWC_LFP_C';
sess_data(cur_tot+8).directory = 'D:\wc_data\2005-06\2005-12-14_CWC_LFP_D';
sess_data(cur_tot+9).directory = 'D:\wc_data\2005-06\Cortex_data\2007-04-16_CWC_LFP_B\2007-4-16_14-8-43';
sess_data(cur_tot+10).directory = 'D:\wc_data\2005-06\Cortex_data\2007-04-17_CWC_LFP_B\2007-4-17_14-49-45';
sess_data(cur_tot+11).directory = 'D:\wc_data\2005-06\Cortex_data\2007-04-18_CWC_LFP_B\2007-4-18_14-28-16';

sess_data(cur_tot+1).name = '2004-11-04';
sess_data(cur_tot+2).name = '2005-12-05_C';
sess_data(cur_tot+3).name = '2005-12-06_A';
sess_data(cur_tot+4).name = '2005-12-06_B';
sess_data(cur_tot+5).name = '2005-12-13_A';
sess_data(cur_tot+6).name = '2005-12-13_D';
sess_data(cur_tot+7).name = '2005-12-14_C';
sess_data(cur_tot+8).name = '2005-12-14_D';
sess_data(cur_tot+9).name = '2007-04-16_B';
sess_data(cur_tot+10).name = '2007-04-17_B';
sess_data(cur_tot+11).name = '2007-04-18_B';

for i = 1:n_23_pyr_par
    sess_data(cur_tot+i).layer = '23';
    sess_data(cur_tot+i).cell_type = 'pyramidal';
    sess_data(cur_tot+i).region = 'parietal';
end

sess_data(cur_tot+1).depth = 240;
sess_data(cur_tot+1).ant_post = -1900;
sess_data(cur_tot+1).lateral = 1800;

sess_data(cur_tot+2).depth = 160;
sess_data(cur_tot+2).ant_post = -2300;
sess_data(cur_tot+2).lateral = 1500;

sess_data(cur_tot+3).depth = 260;
sess_data(cur_tot+3).ant_post = -2100;
sess_data(cur_tot+3).lateral = 700;

sess_data(cur_tot+4).depth = 240;
sess_data(cur_tot+4).ant_post = -2700;
sess_data(cur_tot+4).lateral = 1150;

sess_data(cur_tot+5).depth = 180;
sess_data(cur_tot+5).ant_post = -2100;
sess_data(cur_tot+5).lateral = 1200;

sess_data(cur_tot+6).depth = 200;
sess_data(cur_tot+6).ant_post = -2400;
sess_data(cur_tot+6).lateral = 1530;

sess_data(cur_tot+7).depth = 180;
sess_data(cur_tot+7).ant_post = -2700;
sess_data(cur_tot+7).lateral = 1500;

sess_data(cur_tot+8).depth = 300;
sess_data(cur_tot+8).ant_post = -2200;
sess_data(cur_tot+8).lateral = 1000;

sess_data(cur_tot+9).depth = 290;
sess_data(cur_tot+9).ant_post = -1800;
sess_data(cur_tot+9).lateral = 640;

sess_data(cur_tot+10).depth = 190;
sess_data(cur_tot+10).ant_post = -1700;
sess_data(cur_tot+10).lateral = 300;

sess_data(cur_tot+11).depth = 160;
sess_data(cur_tot+11).ant_post = -1800;
sess_data(cur_tot+11).lateral = 510;

%% Parietal Interneuron (N=3)
n_int_par = 3;
cur_tot = cur_tot + n_23_pyr_par;

sess_data(cur_tot+1).directory = 'D:\wc_data\2005-06\2005-06-13_CWC_LFP\2005-6-13_17-30-12'; %L2/3
sess_data(cur_tot+2).directory = 'D:\wc_data\2005-06\2005-12-13_CWC_LFP_B'; %L2/3
sess_data(cur_tot+3).directory = 'D:\wc_data\2005-06\2005-12-14_CWC_LFP_B'; %maybe pyramidal
% sess_data(cur_tot+4).directory = 'D:\wc_data\2005-06\Cortex_data\2007-04-05_CWC_LFP_B\2007-4-5_17-55-51'; %L5

sess_data(cur_tot+1).name = '2005-06-13';
sess_data(cur_tot+2).name = '2005-12-13_B';
sess_data(cur_tot+3).name = '2005-12-14_B';
% sess_data(cur_tot+4).name = '2007-04-05_B';

sess_data(cur_tot+1).layer = '23';
sess_data(cur_tot+2).layer = '23';
sess_data(cur_tot+3).layer = 'U';
% sess_data(cur_tot+4).layer = '5'; %wierd issue with the data, dropped
% packetU

for i = 1:n_int_par
    sess_data(cur_tot+i).cell_type = 'interneuron';
    sess_data(cur_tot+i).region = 'parietal';
end


sess_data(cur_tot+1).depth = 290;
sess_data(cur_tot+1).ant_post = -1900;
sess_data(cur_tot+1).lateral = 610;

sess_data(cur_tot+2).depth = 220;
sess_data(cur_tot+2).ant_post = -2000;
sess_data(cur_tot+2).lateral = 880;

sess_data(cur_tot+3).depth = 350;
sess_data(cur_tot+3).ant_post = -2600;
sess_data(cur_tot+3).lateral = 800;


%% L5/6 Pyramidal Parietal (N=10)
n_56_pyr_par = 10;
cur_tot = cur_tot + n_int_par;

sess_data(cur_tot+1).directory = 'D:\wc_data\2005-06\2005-12-14_CWC_LFP_A';
sess_data(cur_tot+2).directory = 'D:\wc_data\2005-06\Cortex_data\2007-04-04_CWC_LFP_A\2007-4-4_15-6-52';%maybe L6
sess_data(cur_tot+3).directory = 'D:\wc_data\2005-06\Cortex_data\2007-04-10_CWC_LFP_A\2007-4-10_18-43-7';%border L2/3
sess_data(cur_tot+4).directory = 'D:\wc_data\2005-06\Cortex_data\2007-04-11_CWC_LFP_A\2007-4-11_16-43-19';
sess_data(cur_tot+5).directory = 'D:\wc_data\2005-06\Cortex_data\2007-04-12_CWC_LFP_B\2007-4-12_16-55-57';
sess_data(cur_tot+6).directory = 'D:\wc_data\2005-06\Cortex_data\2007-04-13_CWC_LFP_C\2007-4-13_16-2-15';
sess_data(cur_tot+7).directory = 'D:\wc_data\2005-06\Cortex_data\2007-04-16_CWC_LFP_A\2007-4-16_13-13-21';
sess_data(cur_tot+8).directory = 'D:\wc_data\2005-06\Cortex_data\2007-04-17_CWC_LFP_A\2007-4-17_13-36-49';%maybe L6
sess_data(cur_tot+9).directory = 'D:\wc_data\2005-06\Cortex_data\2007-04-18_CWC_LFP_A\2007-4-18_13-40-23';%maybe L6
sess_data(cur_tot+10).directory = 'D:\wc_data\2005-06\Cortex_data\2007-01-12_CWC_LFP\2007-1-12_16-56-12';%L6
% sess_data(cur_tot+11).directory =
% 'D:\wc_data\2005-06\Cortex_data\2007-01-11_CWC_LFP\2007-1-11_17-44-23';
% %lower sampling rate for wcv

sess_data(cur_tot+1).name = '2005-12-14_A';
sess_data(cur_tot+2).name = '2007-04-04_A';
sess_data(cur_tot+3).name = '2007-04-10_A';
sess_data(cur_tot+4).name = '2007-04-11_A';
sess_data(cur_tot+5).name = '2007-04-12_B';
sess_data(cur_tot+6).name = '2007-04-13_C';
sess_data(cur_tot+7).name = '2007-04-16_A';
sess_data(cur_tot+8).name = '2007-04-17_A';
sess_data(cur_tot+9).name = '2007-04-18_A';
sess_data(cur_tot+10).name = '2007-01-12';
% sess_data(cur_tot+11).name = '2007-01-11';

for i = 1:n_56_pyr_par
    sess_data(cur_tot+i).layer = '56';
    sess_data(cur_tot+i).cell_type = 'pyramidal';
    sess_data(cur_tot+i).region = 'parietal';
end

sess_data(cur_tot+1).depth = 360;
sess_data(cur_tot+1).ant_post = -2400;
sess_data(cur_tot+1).lateral = 1350;

sess_data(cur_tot+2).depth = 670;
sess_data(cur_tot+2).ant_post = -2000;
sess_data(cur_tot+2).lateral = 1050;

sess_data(cur_tot+3).depth = 400;
sess_data(cur_tot+3).ant_post = -2100;
sess_data(cur_tot+3).lateral = 700;

sess_data(cur_tot+4).depth = 490;
sess_data(cur_tot+4).ant_post = -1900;
sess_data(cur_tot+4).lateral = 580;

sess_data(cur_tot+5).depth = 550;
sess_data(cur_tot+5).ant_post = -2000;
sess_data(cur_tot+5).lateral = 640;

sess_data(cur_tot+6).depth = 500;
sess_data(cur_tot+6).ant_post = -2000;
sess_data(cur_tot+6).lateral = 670;

sess_data(cur_tot+7).depth = 450;
sess_data(cur_tot+7).ant_post = -1800;
sess_data(cur_tot+7).lateral = 1050;

sess_data(cur_tot+8).depth = 630;
sess_data(cur_tot+8).ant_post = -1600;
sess_data(cur_tot+8).lateral = 900;

sess_data(cur_tot+9).depth = 680;
sess_data(cur_tot+9).ant_post = -1700;
sess_data(cur_tot+9).lateral = 1100;

sess_data(cur_tot+10).depth = 680;
sess_data(cur_tot+10).ant_post = -1900;
sess_data(cur_tot+10).lateral = 1200;


%% L2/3 Pyramidal Frontal (N=4)
n_23_pyr_fro = 4;
cur_tot = cur_tot + n_56_pyr_par;

sess_data(cur_tot+1).directory = 'D:\wc_data\2005-06\2005-12-05_CWC_LFP_D';
sess_data(cur_tot+2).directory = 'D:\wc_data\2005-06\2005-12-06_CWC_LFP_C';
sess_data(cur_tot+3).directory = 'D:\wc_data\2005-06\2005-12-13_CWC_LFP_E';
sess_data(cur_tot+4).directory = 'D:\wc_data\2005-06\Cortex_data\2007-04-16_CWC_LFP_E\2007-4-16_18-34-42';

sess_data(cur_tot+1).name = '2005-12-05_D';
sess_data(cur_tot+2).name = '2005-12-06_C';
sess_data(cur_tot+3).name = '2005-12-13_E';
sess_data(cur_tot+4).name = '2007-04-16_E';

for i = 1:n_23_pyr_fro
    sess_data(cur_tot+i).layer = '23';
    sess_data(cur_tot+i).cell_type = 'pyramidal';
    sess_data(cur_tot+i).region = 'frontal';
end

sess_data(cur_tot+1).depth = 250;
sess_data(cur_tot+1).ant_post = 1300;
sess_data(cur_tot+1).lateral = 1375;

sess_data(cur_tot+2).depth = 330;
sess_data(cur_tot+2).ant_post = 1200;
sess_data(cur_tot+2).lateral = 930;

sess_data(cur_tot+3).depth = 360;
sess_data(cur_tot+3).ant_post = 1300;
sess_data(cur_tot+3).lateral = 1200;

sess_data(cur_tot+4).depth = 350;
sess_data(cur_tot+4).ant_post = 1400;
sess_data(cur_tot+4).lateral = 800;


%% L5 Pyramidal Frontal (N=3)
n_5_pyr_fro = 3;
cur_tot = cur_tot+ n_23_pyr_fro;

sess_data(cur_tot+1).directory = 'D:\wc_data\2005-06\Cortex_data\2007-04-11_CWC_LFP_B\2007-4-11_17-59-41';
sess_data(cur_tot+2).directory = 'D:\wc_data\2005-06\Cortex_data\2007-04-12_CWC_LFP_C\2007-4-12_18-6-33';
sess_data(cur_tot+3).directory = 'D:\wc_data\2005-06\Cortex_data\2007-04-13_CWC_LFP_D\2007-4-13_17-41-15';

sess_data(cur_tot+1).name = '2007-04-11_B';
sess_data(cur_tot+2).name = '2007-04-12_C';
sess_data(cur_tot+3).name = '2007-04-13_D';

for i = 1:n_5_pyr_fro
    sess_data(cur_tot+i).layer = '5';
    sess_data(cur_tot+i).cell_type = 'pyramidal';
    sess_data(cur_tot+i).region = 'frontal';
end

sess_data(cur_tot+1).depth = 550;
sess_data(cur_tot+1).ant_post = 800;
sess_data(cur_tot+1).lateral = 850;

sess_data(cur_tot+2).depth = 650;
sess_data(cur_tot+2).ant_post = 900;
sess_data(cur_tot+2).lateral = 860;

sess_data(cur_tot+3).depth = 670;
sess_data(cur_tot+3).ant_post = 1700;
sess_data(cur_tot+3).lateral = 900;

%% L2/3 Pyramidal Prefrontal (N=5)
n_23_pyr_pre = 5;
cur_tot = cur_tot + n_5_pyr_fro;

sess_data(cur_tot+1).directory = 'D:\wc_data\2005-06\2005-12-05_CWC_LFP_A';
sess_data(cur_tot+2).directory = 'D:\wc_data\2005-06\2005-12-13_CWC_LFP_C'; %border L5
sess_data(cur_tot+3).directory = 'D:\wc_data\2005-06\Cortex_data\2007-04-05_CWC_LFP_C\2007-4-5_19-31-0';
sess_data(cur_tot+4).directory = 'D:\wc_data\2005-06\Cortex_data\2007-04-10_CWC_LFP_B\2007-4-10_19-35-5';%maybe L5
sess_data(cur_tot+5).directory = 'D:\wc_data\2005-06\Cortex_data\2007-04-17_CWC_LFP_C\2007-4-17_15-41-51';

sess_data(cur_tot+1).name = '2005-12-05_A';
sess_data(cur_tot+2).name = '2005-12-13_C';
sess_data(cur_tot+3).name = '2007-04-05_C';
sess_data(cur_tot+4).name = '2007-04-10_B';
sess_data(cur_tot+5).name = '2007-04-17_C';

for i = 1:n_23_pyr_pre
    sess_data(cur_tot+i).layer = '23';
    sess_data(cur_tot+i).cell_type = 'pyramidal';
    sess_data(cur_tot+i).region = 'prefrontal';
end

sess_data(cur_tot+1).depth = 290;
sess_data(cur_tot+1).ant_post = 2800;
sess_data(cur_tot+1).lateral = 1200;

sess_data(cur_tot+2).depth = 430;
sess_data(cur_tot+2).ant_post = 2500;
sess_data(cur_tot+2).lateral = 1050;

sess_data(cur_tot+3).depth = 300;
sess_data(cur_tot+3).ant_post = 2900;
sess_data(cur_tot+3).lateral = 1300;

sess_data(cur_tot+4).depth = 450;
sess_data(cur_tot+4).ant_post = 2600;
sess_data(cur_tot+4).lateral = 1000;

sess_data(cur_tot+5).depth = 250;
sess_data(cur_tot+5).ant_post = 2500;
sess_data(cur_tot+5).lateral = 1150;

%% L5 Pyramidal Prefrontal (N=8)
n_5_pyr_pre = 8;
cur_tot = cur_tot + n_23_pyr_pre;

sess_data(cur_tot+1).directory = 'D:\wc_data\2005-06\Cortex_data\2007-04-04_CWC_LFP_B\2007-4-4_16-0-35';
sess_data(cur_tot+2).directory = 'D:\wc_data\2005-06\Cortex_data\2007-04-13_CWC_LFP_B\2007-4-13_13-45-34';
sess_data(cur_tot+3).directory = 'D:\wc_data\2005-06\Cortex_data\2007-04-16_CWC_LFP_C\2007-4-16_15-18-23';
sess_data(cur_tot+4).directory = 'D:\wc_data\2005-06\Cortex_data\2007-07-27_CWC_LFP_B'; %slender !no LF8!
sess_data(cur_tot+5).directory = 'D:\wc_data\2005-06\Cortex_data\2007-08-14_CWC_LFP'; %thick
sess_data(cur_tot+6).directory = 'D:\wc_data\2005-06\Cortex_data\2007-08-15_CWC_LFP';  %!no LF8!
sess_data(cur_tot+7).directory = 'D:\wc_data\2005-06\Cortex_data\2007-08-17_CWC_LFP_A\2007-8-17_16-33-25'; %!no LF8!
% sess_data(cur_tot+8).directory = 'D:\wc_data\2005-06\Cortex_data\2007-08-17-CWC_LFP_B'; %!no LF8!
sess_data(cur_tot+8).directory = 'D:\wc_data\2005-06\Cortex_data\2007-08-21_CWC_LFP\2007-8-21_14-4-42'; %!no LF8!

sess_data(cur_tot+1).name = '2007-04-04_B';
sess_data(cur_tot+2).name = '2007-04-13_B';
sess_data(cur_tot+3).name = '2007-04-16_C';
sess_data(cur_tot+4).name = '2007-07-27_B';
sess_data(cur_tot+5).name = '2007-08-14'; %best spikelet example
sess_data(cur_tot+6).name = '2007-08-15'; %2nd half of recording all theta, no UDS
sess_data(cur_tot+7).name = '2007-08-17_A';
% sess_data(cur_tot+8).name = '2007-08-17_B';
sess_data(cur_tot+8).name = '2007-08-21';

for i = 1:n_5_pyr_pre
    sess_data(cur_tot+i).layer = '5';
    sess_data(cur_tot+i).cell_type = 'pyramidal';
    sess_data(cur_tot+i).region = 'prefrontal';
end


sess_data(cur_tot+1).depth = 550;
sess_data(cur_tot+1).ant_post = 2800;
sess_data(cur_tot+1).lateral = 800;

sess_data(cur_tot+2).depth = 440;
sess_data(cur_tot+2).ant_post = 2800;
sess_data(cur_tot+2).lateral = 800;

sess_data(cur_tot+3).depth = 500;
sess_data(cur_tot+3).ant_post = 2600;
sess_data(cur_tot+3).lateral = 850;

sess_data(cur_tot+4).depth = 1000;
sess_data(cur_tot+4).ant_post = 2700;
sess_data(cur_tot+4).lateral = 500;

sess_data(cur_tot+5).depth = 1250;
sess_data(cur_tot+5).ant_post = 2200;
sess_data(cur_tot+5).lateral = 450;

sess_data(cur_tot+6).depth = 1600;
sess_data(cur_tot+6).ant_post = 2100;
sess_data(cur_tot+6).lateral = 850;

sess_data(cur_tot+7).depth = 1500;
sess_data(cur_tot+7).ant_post = 2500;
sess_data(cur_tot+7).lateral = 600;

sess_data(cur_tot+8).depth = 1500;
sess_data(cur_tot+8).ant_post = 2100;
sess_data(cur_tot+8).lateral = 540;


%% Thalamic (N = 5)
n_thal = 5;
cur_tot = cur_tot + n_5_pyr_pre;

sess_data(cur_tot+1).directory = 'D:\wc_data\2005-06\TH_LFP\2006-02-28_TH_LFP';
sess_data(cur_tot+2).directory = 'D:\wc_data\2005-06\TH_LFP\2006-03-10_TH_LFP';
sess_data(cur_tot+3).directory = 'D:\wc_data\2005-06\TH_LFP\2007-01-15_TH_LFP\2007-1-15_19-2-10';
sess_data(cur_tot+4).directory = 'D:\wc_data\2005-06\TH_LFP\2007-01-16_TH_LFP\2007-1-16_16-07-47';
sess_data(cur_tot+5).directory = 'D:\wc_data\2005-06\TH_LFP\2007-01-19_TH_LFP\2007-1-19_19-42-24';

sess_data(cur_tot+1).name = '2006-02-28-TH';
sess_data(cur_tot+2).name = '2006-03-10-TH';
sess_data(cur_tot+3).name = '2007-01-15-TH';
sess_data(cur_tot+4).name = '2007-01-16-TH';
sess_data(cur_tot+5).name = '2007-01-19-TH';

for i = 1:n_thal
    sess_data(cur_tot+i).layer = 'U';
    sess_data(cur_tot+i).cell_type = 'U';
    sess_data(cur_tot+i).region = 'thalamus';
end

sess_data(cur_tot+1).depth = nan;
sess_data(cur_tot+1).ant_post = nan;
sess_data(cur_tot+1).lateral = nan;

sess_data(cur_tot+2).depth = nan;
sess_data(cur_tot+2).ant_post = nan;
sess_data(cur_tot+2).lateral = nan;

sess_data(cur_tot+3).depth = nan;
sess_data(cur_tot+3).ant_post = nan;
sess_data(cur_tot+3).lateral = nan;

sess_data(cur_tot+4).depth = nan;
sess_data(cur_tot+4).ant_post = nan;
sess_data(cur_tot+4).lateral = nan;

sess_data(cur_tot+5).depth = nan;
sess_data(cur_tot+5).ant_post = nan;
sess_data(cur_tot+5).lateral = nan;


%% Barrel Cortex (N = 4)
n_barrel = 4;
cur_tot = cur_tot + n_thal;

sess_data(cur_tot+1).directory = 'D:\wc_data\2005-06\2006-12-17_CWC_LFP_B\2006-12-17_18-47-4';
sess_data(cur_tot+2).directory = 'D:\wc_data\2005-06\2006-12-19_CWC_LFP_A\2006-12-19_18-36-3';
sess_data(cur_tot+3).directory = 'D:\wc_data\2005-06\2006-12-19_CWC_LFP_B\2006-12-19_20-41-15';
sess_data(cur_tot+4).directory = 'D:\wc_data\2005-06\2006-12-20_CWC_LFP\2006-12-20_18-14-19';

sess_data(cur_tot+1).name = '2006-12-17_B';
sess_data(cur_tot+2).name = '2006-12-19_A';
sess_data(cur_tot+3).name = '2006-12-19_B';
sess_data(cur_tot+4).name = '2006-12-20';

sess_data(cur_tot+1).cell_type = 'pyramidal';
sess_data(cur_tot+2).cell_type = 'pyramidal';
sess_data(cur_tot+3).cell_type = 'spiny stellate';
sess_data(cur_tot+4).cell_type = 'pyramidal';


for i = 1:n_barrel
    sess_data(cur_tot+i).layer = '34';
    sess_data(cur_tot+i).region = 'barrel';
end

sess_data(cur_tot+1).depth = nan;
sess_data(cur_tot+1).ant_post = nan;
sess_data(cur_tot+1).lateral = nan;

sess_data(cur_tot+2).depth = nan;
sess_data(cur_tot+2).ant_post = nan;
sess_data(cur_tot+2).lateral = nan;

sess_data(cur_tot+3).depth = nan;
sess_data(cur_tot+3).ant_post = nan;
sess_data(cur_tot+3).lateral = nan;

sess_data(cur_tot+4).depth = nan;
sess_data(cur_tot+4).ant_post = nan;
sess_data(cur_tot+4).lateral = nan;


%% V1 (N = 1)
n_v1 = 1;
cur_tot = cur_tot+n_barrel;

sess_data(cur_tot+1).directory = 'D:\wc_data\2005-06\Cortex_data\2007-05-30_CWC_LFP_B';

sess_data(cur_tot+1).name = '2007-05-30_B';

for i = 1:n_v1
    sess_data(cur_tot+i).cell_type = 'pyramidal';
    sess_data(cur_tot+i).layer = '23';
    sess_data(cur_tot+i).region = 'V1';
end

sess_data(cur_tot+1).depth = 230;
sess_data(cur_tot+1).ant_post = -4800;
sess_data(cur_tot+1).lateral = 2650;

%% Basal Ganglia (N=2)
n_bg = 2;
cur_tot = cur_tot + n_v1;

sess_data(cur_tot+1).directory = 'D:\wc_data\2005-06\BG_WC_LFP\2006-05-05_BG_LFP_A';
sess_data(cur_tot+2).directory = 'D:\wc_data\2005-06\BG_WC_LFP\2006-05-05_BG_LFP_B';

sess_data(cur_tot+1).name = '2006-05-05_A';
sess_data(cur_tot+2).name = '2006-05-05_B';

for i = 1:n_bg
    sess_data(cur_tot+i).cell_type = 'U';
    sess_data(cur_tot+i).layer = 'U';
    sess_data(cur_tot+i).region = 'Basal_Ganglia';
end

sess_data(cur_tot+1).depth = nan;
sess_data(cur_tot+1).ant_post = nan;
sess_data(cur_tot+1).lateral = nan;

sess_data(cur_tot+2).depth = nan;
sess_data(cur_tot+2).ant_post = nan;
sess_data(cur_tot+2).lateral = nan;
