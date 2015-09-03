%% cortical pyramidal neurons
cd G:\WC_Germany\parietal_cortical_2010\
clear all


%% L2/3 Pyrmidal Parietal  (N=11)
cur_tot = 0;
n_23_pyr_par = 11;

sess_data(cur_tot+1).directory = 'G:\wc_data\2004-11-04_CWC_LFP';
sess_data(cur_tot+2).directory = 'G:\wc_data\2005-12-05_CWC_LFP_C';
sess_data(cur_tot+3).directory = 'G:\wc_data\2005-12-06_CWC_LFP_A';
sess_data(cur_tot+4).directory = 'G:\wc_data\2005-12-06_CWC_LFP_B';
sess_data(cur_tot+5).directory = 'G:\wc_data\2005-12-13_CWC_LFP_A';
sess_data(cur_tot+6).directory = 'G:\wc_data\2005-12-13_CWC_LFP_D';
sess_data(cur_tot+7).directory = 'G:\wc_data\2005-12-14_CWC_LFP_C';
sess_data(cur_tot+8).directory = 'G:\wc_data\2005-12-14_CWC_LFP_D';
sess_data(cur_tot+9).directory = 'G:\wc_data\2007-04-16_CWC_LFP_B\2007-4-16_14-8-43';
sess_data(cur_tot+10).directory = 'G:\wc_data\2007-04-17_CWC_LFP_B\2007-4-17_14-49-45';
sess_data(cur_tot+11).directory = 'G:\wc_data\2007-04-18_CWC_LFP_B\2007-4-18_14-28-16';

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

sess_data(cur_tot+1).an_con =.14/14;
sess_data(cur_tot+1).depth = 240;
sess_data(cur_tot+1).ant_post = -1900;
sess_data(cur_tot+1).lateral = 1800;
sess_data(cur_tot+1).gains = [60;2e3;5e3;2e3;2e3;2e3;2e3;2e3];
sess_data(cur_tot+1).thom_elec = 0;

sess_data(cur_tot+2).an_con =.13/15;
sess_data(cur_tot+2).depth = 160;
sess_data(cur_tot+2).ant_post = -2300;
sess_data(cur_tot+2).lateral = 1500;
sess_data(cur_tot+2).gains = [80;1e3;2e3;1e3;1e3;0;2e3;1e3];
sess_data(cur_tot+2).thom_elec = 0;

sess_data(cur_tot+3).an_con =.13/15.2;
sess_data(cur_tot+3).depth = 260;
sess_data(cur_tot+3).ant_post = -2100;
sess_data(cur_tot+3).lateral = 700;
sess_data(cur_tot+3).gains = [60;1e3;2e3;1e3;1e3;0;1e3;1e3];
sess_data(cur_tot+3).thom_elec = 0;

sess_data(cur_tot+4).an_con =.13/15.6;
sess_data(cur_tot+4).depth = 240;
sess_data(cur_tot+4).ant_post = -2700;
sess_data(cur_tot+4).lateral = 1150;
sess_data(cur_tot+4).gains = [70;1e3;2e3;1e3;1e3;0;1e3;1e3];
sess_data(cur_tot+4).thom_elec = 0;

sess_data(cur_tot+5).an_con =.14/16.2;
sess_data(cur_tot+5).depth = 180;
sess_data(cur_tot+5).ant_post = -2100;
sess_data(cur_tot+5).lateral = 1200;
sess_data(cur_tot+5).gains = [50;1e3;2e3;1e3;1e3;0;2e3;2e3];
sess_data(cur_tot+5).thom_elec = 0;

sess_data(cur_tot+6).an_con =.15/17;
sess_data(cur_tot+6).depth = 200;
sess_data(cur_tot+6).ant_post = -2400;
sess_data(cur_tot+6).lateral = 1530;
sess_data(cur_tot+6).gains = [70;1e3;2e3;1e3;1e3;0;2e3;2e3];
sess_data(cur_tot+6).thom_elec = 0;

sess_data(cur_tot+7).an_con =.15/17.3;
sess_data(cur_tot+7).depth = 180;
sess_data(cur_tot+7).ant_post = -2700;
sess_data(cur_tot+7).lateral = 1500;
sess_data(cur_tot+7).gains = [60;1e3;2e3;1e3;1e3;0;2e3;2e3];
sess_data(cur_tot+7).thom_elec = 0;

sess_data(cur_tot+8).an_con =.15/17.3;
sess_data(cur_tot+8).depth = 300;
sess_data(cur_tot+8).ant_post = -2200;
sess_data(cur_tot+8).lateral = 1000;
sess_data(cur_tot+8).gains = [60;1e3;2e3;1e3;1e3;0;2e3;2e3];
sess_data(cur_tot+8).thom_elec = 0;

sess_data(cur_tot+9).an_con =.15/17.5;
sess_data(cur_tot+9).depth = 290;
sess_data(cur_tot+9).ant_post = -1800;
sess_data(cur_tot+9).lateral = 640;
sess_data(cur_tot+9).gains = [70;2e3;2e3;2e3;3e3;3e3;3e3;3e3];
sess_data(cur_tot+9).thom_elec = 1;

sess_data(cur_tot+10).an_con =.16/17.4;
sess_data(cur_tot+10).depth = 190;
sess_data(cur_tot+10).ant_post = -1700;
sess_data(cur_tot+10).lateral = 300;
sess_data(cur_tot+10).gains = [100;2e3;2e3;2e3;3e3;3e3;3e3;3e3];
sess_data(cur_tot+10).thom_elec = 1;

sess_data(cur_tot+11).an_con =.15/16.2;
sess_data(cur_tot+11).depth = 160;
sess_data(cur_tot+11).ant_post = -1800;
sess_data(cur_tot+11).lateral = 510;
sess_data(cur_tot+11).gains = [100;2e3;2e3;2e3;3e3;3e3;3e3;3e3];
sess_data(cur_tot+11).thom_elec = 1;

%% Parietal Interneuron (N=3)
n_int_par = 3;
cur_tot = cur_tot + n_23_pyr_par;

sess_data(cur_tot+1).directory = 'G:\wc_data\2005-06-13_CWC_LFP\2005-6-13_17-30-12'; %L2/3
sess_data(cur_tot+2).directory = 'G:\wc_data\2005-12-13_CWC_LFP_B'; %L2/3
sess_data(cur_tot+3).directory = 'G:\wc_data\2005-12-14_CWC_LFP_B'; %maybe pyramidal
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
sess_data(cur_tot+1).gains = [50;3e3;5e3;3e3;3e3;3e3;3e3;3e3];
sess_data(cur_tot+1).thom_elec = 0;

sess_data(cur_tot+2).depth = 220;
sess_data(cur_tot+2).ant_post = -2000;
sess_data(cur_tot+2).lateral = 880;
sess_data(cur_tot+2).gains = [70;1e3;2e3;1e3;1e3;0;2e3;2e3];
sess_data(cur_tot+2).thom_elec = 0;

sess_data(cur_tot+3).depth = 350;
sess_data(cur_tot+3).ant_post = -2600;
sess_data(cur_tot+3).lateral = 800;
sess_data(cur_tot+3).gains = [50;1e3;2e3;1e3;1e3;0;2e3;2e3];
sess_data(cur_tot+3).thom_elec = 0;


%% L5/6 Pyramidal Parietal (N=10)
n_56_pyr_par = 10;
cur_tot = cur_tot + n_int_par;

% TO BE ADDED! SOME TRACHEOSTOMY
% data_dir{1} = 'G:\wc_data\2005-08-30_CWC_LFP_B\2005-8-30_18-20-39';
% data_dir{2} = 'G:\wc_data\2005-08-31_CWC_LFP\2005-8-31_16-5-40';
% data_dir{3} = 'G:\wc_data\2005-10-06_CWC_LFP\2005-10-6_17-35-50';
% data_dir{4} = 'G:\wc_data\2005-10-10_CWC_LFP';
% data_dir{5} = 'G:\wc_data\2007-01-11_CWC_LFP\2007-1-11_17-44-23';
% data_dir{7} = 'G:\wc_data\2007-04-16_CWC_LFP_D\2007-4-16_17-34-6';
% data_dir{8} = 'G:\wc_data\2007-04-05_CWC_LFP_A\2007-4-5_14-49-1';

sess_data(cur_tot+1).directory = 'G:\wc_data\2005-12-14_CWC_LFP_A';
sess_data(cur_tot+2).directory = 'G:\wc_data\2007-04-04_CWC_LFP_A\2007-4-4_15-6-52';%maybe L6
sess_data(cur_tot+3).directory = 'G:\wc_data\2007-04-10_CWC_LFP_A\2007-4-10_18-43-7';%border L2/3
sess_data(cur_tot+4).directory = 'G:\wc_data\2007-04-11_CWC_LFP_A\2007-4-11_16-43-19';
sess_data(cur_tot+5).directory = 'G:\wc_data\2007-04-12_CWC_LFP_B\2007-4-12_16-55-57';
sess_data(cur_tot+6).directory = 'G:\wc_data\2007-04-13_CWC_LFP_C\2007-4-13_16-2-15';
sess_data(cur_tot+7).directory = 'G:\wc_data\2007-04-16_CWC_LFP_A\2007-4-16_13-13-21';
sess_data(cur_tot+8).directory = 'G:\wc_data\2007-04-17_CWC_LFP_A\2007-4-17_13-36-49';%maybe L6
sess_data(cur_tot+9).directory = 'G:\wc_data\2007-04-18_CWC_LFP_A\2007-4-18_13-40-23';%maybe L6
sess_data(cur_tot+10).directory = 'G:\wc_data\2007-01-12_CWC_LFP\2007-1-12_16-56-12';%L6
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

sess_data(cur_tot+1).an_con =.16/18.8;
sess_data(cur_tot+1).depth = 360;
sess_data(cur_tot+1).ant_post = -2400;
sess_data(cur_tot+1).lateral = 1350;
sess_data(cur_tot+1).gains = [50;1e3;2e3;1e3;1e3;0;2e3;2e3];
sess_data(cur_tot+1).thom_elec = 0;

sess_data(cur_tot+2).an_con =.16/17.3;
sess_data(cur_tot+2).depth = 670;
sess_data(cur_tot+2).ant_post = -2000;
sess_data(cur_tot+2).lateral = 1050;
sess_data(cur_tot+2).gains = [50;2e3;2e3;2e3;2e3;2e3;2e3;2e3];
sess_data(cur_tot+2).thom_elec = 0;

sess_data(cur_tot+3).an_con =.18/18;
sess_data(cur_tot+3).depth = 400;
sess_data(cur_tot+3).ant_post = -2100;
sess_data(cur_tot+3).lateral = 700;
sess_data(cur_tot+3).gains = [100;2e3;2e3;2e3;3e3;3e3;3e3;3e3];
sess_data(cur_tot+3).thom_elec = 1;

sess_data(cur_tot+4).an_con =.2/23;
sess_data(cur_tot+4).depth = 490;
sess_data(cur_tot+4).ant_post = -1900;
sess_data(cur_tot+4).lateral = 580;
sess_data(cur_tot+4).gains = [80;2e3;2e3;2e3;3e3;3e3;3e3;3e3];
sess_data(cur_tot+4).thom_elec = 1;

sess_data(cur_tot+5).an_con =.14/15.8;
sess_data(cur_tot+5).depth = 550;
sess_data(cur_tot+5).ant_post = -2000;
sess_data(cur_tot+5).lateral = 640;
sess_data(cur_tot+5).gains = [80;2e3;2e3;2e3;3e3;3e3;3e3;3e3];
sess_data(cur_tot+5).thom_elec = 1;

sess_data(cur_tot+6).an_con =.15/16.5;
sess_data(cur_tot+6).depth = 500;
sess_data(cur_tot+6).ant_post = -2000;
sess_data(cur_tot+6).lateral = 670;
sess_data(cur_tot+6).gains = [100;2e3;2e3;2e3;3e3;3e3;3e3;3e3];
sess_data(cur_tot+6).thom_elec = 1;

sess_data(cur_tot+7).an_con =.15/17.5;
sess_data(cur_tot+7).depth = 450;
sess_data(cur_tot+7).ant_post = -1800;
sess_data(cur_tot+7).lateral = 1050;
sess_data(cur_tot+7).gains = [80;2e3;2e3;2e3;3e3;3e3;3e3;3e3];
sess_data(cur_tot+7).thom_elec = 1;

sess_data(cur_tot+8).an_con =.16/17.4;
sess_data(cur_tot+8).depth = 630;
sess_data(cur_tot+8).ant_post = -1600;
sess_data(cur_tot+8).lateral = 900;
sess_data(cur_tot+8).gains = [100;2e3;2e3;2e3;3e3;3e3;3e3;3e3];
sess_data(cur_tot+8).thom_elec = 1;

sess_data(cur_tot+9).an_con =.15/16.2;
sess_data(cur_tot+9).depth = 680;
sess_data(cur_tot+9).ant_post = -1700;
sess_data(cur_tot+9).lateral = 1100;
sess_data(cur_tot+9).gains = [100;2e3;2e3;2e3;3e3;3e3;3e3;3e3];
sess_data(cur_tot+9).thom_elec = 1;

sess_data(cur_tot+10).an_con =.19/21.2;
sess_data(cur_tot+10).depth = 680;
sess_data(cur_tot+10).ant_post = -1900;
sess_data(cur_tot+10).lateral = 1200;
sess_data(cur_tot+10).gains = [80;2e3;2e3;2e3;3e3;3e3;3e3;3e3];
sess_data(cur_tot+10).thom_elec = 0;


%% L2/3 Pyramidal Frontal (N=4)
n_23_pyr_fro = 4;
cur_tot = cur_tot + n_56_pyr_par;

sess_data(cur_tot+1).directory = 'G:\wc_data\2005-12-05_CWC_LFP_D';
sess_data(cur_tot+2).directory = 'G:\wc_data\2005-12-06_CWC_LFP_C';
sess_data(cur_tot+3).directory = 'G:\wc_data\2005-12-13_CWC_LFP_E';
sess_data(cur_tot+4).directory = 'G:\wc_data\2007-04-16_CWC_LFP_E\2007-4-16_18-34-42';

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
sess_data(cur_tot+1).gains = [80;1e3;2e3;1e3;1e3;0;2e3;1e3];
sess_data(cur_tot+1).thom_elec = 0;

sess_data(cur_tot+2).depth = 330;
sess_data(cur_tot+2).ant_post = 1200;
sess_data(cur_tot+2).lateral = 930;
sess_data(cur_tot+2).gains = [70;1e3;2e3;1e3;1e3;0;1e3;1e3];
sess_data(cur_tot+2).thom_elec = 0;

sess_data(cur_tot+3).depth = 360;
sess_data(cur_tot+3).ant_post = 1300;
sess_data(cur_tot+3).lateral = 1200;
sess_data(cur_tot+3).gains = [70;1e3;2e3;1e3;1e3;0;2e3;2e3];
sess_data(cur_tot+3).thom_elec = 0;

sess_data(cur_tot+4).depth = 350;
sess_data(cur_tot+4).ant_post = 1400;
sess_data(cur_tot+4).lateral = 800;
sess_data(cur_tot+4).gains = [100;2e3;2e3;2e3;3e3;3e3;3e3;3e3];
sess_data(cur_tot+4).thom_elec = 1;


%% L5 Pyramidal Frontal (N=3)
n_5_pyr_fro = 3;
cur_tot = cur_tot+ n_23_pyr_fro;

sess_data(cur_tot+1).directory = 'G:\wc_data\2007-04-11_CWC_LFP_B\2007-4-11_17-59-41';
sess_data(cur_tot+2).directory = 'G:\wc_data\2007-04-12_CWC_LFP_C\2007-4-12_18-6-33';
sess_data(cur_tot+3).directory = 'G:\wc_data\2007-04-13_CWC_LFP_D\2007-4-13_17-41-15';

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
sess_data(cur_tot+1).gains = [90;2e3;2e3;2e3;3e3;3e3;3e3;3e3];
sess_data(cur_tot+1).thom_elec = 1;

sess_data(cur_tot+2).depth = 650;
sess_data(cur_tot+2).ant_post = 900;
sess_data(cur_tot+2).lateral = 860;
sess_data(cur_tot+2).gains = [80;2e3;2e3;2e3;3e3;3e3;3e3;3e3];
sess_data(cur_tot+2).thom_elec = 1;

sess_data(cur_tot+3).depth = 670;
sess_data(cur_tot+3).ant_post = 1700;
sess_data(cur_tot+3).lateral = 900;
sess_data(cur_tot+3).gains = [100;2e3;2e3;2e3;3e3;3e3;3e3;3e3];
sess_data(cur_tot+3).thom_elec = 1;

%% L2/3 Pyramidal Prefrontal (N=5)
n_23_pyr_pre = 5;
cur_tot = cur_tot + n_5_pyr_fro;

sess_data(cur_tot+1).directory = 'G:\wc_data\2005-12-05_CWC_LFP_A';
sess_data(cur_tot+2).directory = 'G:\wc_data\2005-12-13_CWC_LFP_C'; %border L5
sess_data(cur_tot+3).directory = 'G:\wc_data\2007-04-05_CWC_LFP_C\2007-4-5_19-31-0';
sess_data(cur_tot+4).directory = 'G:\wc_data\2007-04-10_CWC_LFP_B\2007-4-10_19-35-5';%maybe L5
sess_data(cur_tot+5).directory = 'G:\wc_data\2007-04-17_CWC_LFP_C\2007-4-17_15-41-51';

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
sess_data(cur_tot+1).gains = [70;1e3;2e3;1e3;1e3;0;2e3;1e3];
sess_data(cur_tot+1).thom_elec = 0;

sess_data(cur_tot+2).depth = 430;
sess_data(cur_tot+2).ant_post = 2500;
sess_data(cur_tot+2).lateral = 1050;
sess_data(cur_tot+2).gains = [80;1e3;2e3;1e3;1e3;0;2e3;2e3];
sess_data(cur_tot+2).thom_elec = 0;

sess_data(cur_tot+3).depth = 300;
sess_data(cur_tot+3).ant_post = 2900;
sess_data(cur_tot+3).lateral = 1300;
sess_data(cur_tot+3).gains = [150;2e3;2e3;1e3;3e3;3e3;3e3;3e3];
sess_data(cur_tot+3).thom_elec = 1;

sess_data(cur_tot+4).depth = 450;
sess_data(cur_tot+4).ant_post = 2600;
sess_data(cur_tot+4).lateral = 1000;
sess_data(cur_tot+4).gains = [100;2e3;2e3;2e3;3e3;3e3;3e3;3e3];
sess_data(cur_tot+4).thom_elec = 1;

sess_data(cur_tot+5).depth = 250;
sess_data(cur_tot+5).ant_post = 2500;
sess_data(cur_tot+5).lateral = 1150;
sess_data(cur_tot+5).gains = [80;2e3;2e3;2e3;3e3;3e3;3e3;3e3];
sess_data(cur_tot+5).thom_elec = 1;

%% L5 Pyramidal Prefrontal (N=8)
n_5_pyr_pre = 8;
cur_tot = cur_tot + n_23_pyr_pre;

sess_data(cur_tot+1).directory = 'G:\wc_data\2007-04-04_CWC_LFP_B\2007-4-4_16-0-35';
sess_data(cur_tot+2).directory = 'G:\wc_data\2007-04-13_CWC_LFP_B\2007-4-13_13-45-34';
sess_data(cur_tot+3).directory = 'G:\wc_data\2007-04-16_CWC_LFP_C\2007-4-16_15-18-23';
sess_data(cur_tot+4).directory = 'G:\wc_data\2007-07-27_CWC_LFP_B'; %slender !no LF8!
sess_data(cur_tot+5).directory = 'G:\wc_data\2007-08-14_CWC_LFP'; %thick
sess_data(cur_tot+6).directory = 'G:\wc_data\2007-08-15_CWC_LFP';  %!no LF8!
sess_data(cur_tot+7).directory = 'G:\wc_data\2007-08-17_CWC_LFP_A\2007-8-17_16-33-25'; %!no LF8!
% sess_data(cur_tot+8).directory = 'D:\wc_data\2005-06\Cortex_data\2007-08-17-CWC_LFP_B'; %!no LF8!
sess_data(cur_tot+8).directory = 'G:\wc_data\2007-08-21_CWC_LFP\2007-8-21_14-4-42'; %!no LF8!

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
sess_data(cur_tot+1).gains = [50;3e3;3e3;3e3;3e3;3e3;3e3;3e3];
sess_data(cur_tot+1).thom_elec = 0;

sess_data(cur_tot+2).depth = 440;
sess_data(cur_tot+2).ant_post = 2800;
sess_data(cur_tot+2).lateral = 800;
sess_data(cur_tot+2).gains = [100;2e3;2e3;2e3;3e3;3e3;3e3;3e3];
sess_data(cur_tot+2).thom_elec = 1;

sess_data(cur_tot+3).depth = 500;
sess_data(cur_tot+3).ant_post = 2600;
sess_data(cur_tot+3).lateral = 850;
sess_data(cur_tot+3).gains = [90;2e3;2e3;2e3;3e3;3e3;3e3;3e3];
sess_data(cur_tot+3).thom_elec = 1;

sess_data(cur_tot+4).depth = 1000;
sess_data(cur_tot+4).ant_post = 2700;
sess_data(cur_tot+4).lateral = 500;
sess_data(cur_tot+4).gains = [100;2e3;2e3;0;4e3;4e3;4e3;4e3];
sess_data(cur_tot+4).thom_elec = 0;

sess_data(cur_tot+5).depth = 1250;
sess_data(cur_tot+5).ant_post = 2200;
sess_data(cur_tot+5).lateral = 450;
sess_data(cur_tot+5).gains = [100;1e3;2e3;0;3e3;3e3;3e3;4e3];
sess_data(cur_tot+5).thom_elec = 0;

sess_data(cur_tot+6).depth = 1600;
sess_data(cur_tot+6).ant_post = 2100;
sess_data(cur_tot+6).lateral = 850;
sess_data(cur_tot+6).gains = [100;1e3;2e3;0;3e3;3e3;3e3;4e3];
sess_data(cur_tot+6).thom_elec = 0;

sess_data(cur_tot+7).depth = 1500;
sess_data(cur_tot+7).ant_post = 2500;
sess_data(cur_tot+7).lateral = 600;
sess_data(cur_tot+7).gains = [100;2e3;2e3;0;3e3;3e3;3e3;4e3];
sess_data(cur_tot+7).thom_elec = 0;

sess_data(cur_tot+8).depth = 1500;
sess_data(cur_tot+8).ant_post = 2100;
sess_data(cur_tot+8).lateral = 540;
sess_data(cur_tot+8).gains = [100;2e3;3e3;0;4e3;4e3;4e3;5e3];
sess_data(cur_tot+8).thom_elec = 0;

%% NEW ADDITIONS
cur_tot = cur_tot + n_5_pyr_pre;

% sess_data(cur_tot+1).directory = 'G:\wc_data\2007-04-04_CWC_LFP_B\2007-4-4_16-0-35';
% sess_data(cur_tot+1).name = '2007-04-04_B';
% sess_data(cur_tot+1).layer = '23';
% sess_data(cur_tot+1).cell_type = 'pyramidal';
% sess_data(cur_tot+1).region = 'parietal';
% sess_data(cur_tot+1).depth = 350;
% sess_data(cur_tot+1).ant_post = -2600;
% sess_data(cur_tot+1).lateral = 800;
% sess_data(cur_tot+1).gains = [50;1e3;2e3;1e3;1e3;0;2e3;2e3];
% sess_data(cur_tot+1).thom_elec = 0;

sess_data(cur_tot+1).directory = 'G:\wc_data\2007-01-11_CWC_LFP\2007-1-11_17-44-23';
sess_data(cur_tot+1).name = '2007-01-11';
sess_data(cur_tot+1).layer = '5';
sess_data(cur_tot+1).cell_type = 'pyramidal';
sess_data(cur_tot+1).region = 'parietal';
sess_data(cur_tot+1).depth = 500;
sess_data(cur_tot+1).ant_post = -1600;
sess_data(cur_tot+1).lateral = 1300;
sess_data(cur_tot+1).gains = [50;2e3;1e3;2e3;2e3;0;2e3;3e3];
sess_data(cur_tot+1).thom_elec = 0;

sess_data(cur_tot+2).directory = 'G:\wc_data\2007-04-05_CWC_LFP_A\2007-4-5_14-49-1';
sess_data(cur_tot+2).name = '2007-04-05_A';
sess_data(cur_tot+2).layer = '6';
sess_data(cur_tot+2).cell_type = 'pyramidal';
sess_data(cur_tot+2).region = 'parietal';
sess_data(cur_tot+2).depth = 720;
sess_data(cur_tot+2).ant_post = -1900;
sess_data(cur_tot+2).lateral = 850;
sess_data(cur_tot+2).gains = [100;2e3;2e3;1e3;3e3;3e3;3e3;3e3];
sess_data(cur_tot+2).thom_elec = 1;

sess_data(cur_tot+3).directory = 'G:\wc_data\2007-04-16_CWC_LFP_D\2007-4-16_17-34-6';
sess_data(cur_tot+3).name = '2007-04-16_D';
sess_data(cur_tot+3).layer = '5';
sess_data(cur_tot+3).cell_type = 'pyramidal';
sess_data(cur_tot+3).region = 'parietal';
sess_data(cur_tot+3).depth = 640;
sess_data(cur_tot+3).ant_post = -1700;
sess_data(cur_tot+3).lateral = 990;
sess_data(cur_tot+3).gains = [90;2e3;2e3;2e3;3e3;3e3;3e3;3e3];
sess_data(cur_tot+3).thom_elec = 1;

% sess_data(cur_tot+4).directory = 'G:\wc_data\2007-08-17-CWC_LFP_B';
% sess_data(cur_tot+4).name = '2007-08-17_B';
% sess_data(cur_tot+4).layer = '5';
% sess_data(cur_tot+4).cell_type = 'pyramidal';
% sess_data(cur_tot+4).region = 'prefrontal';
% sess_data(cur_tot+4).depth = 1500;
% sess_data(cur_tot+4).ant_post = 1900;
% sess_data(cur_tot+4).lateral = 650;
% sess_data(cur_tot+4).gains = [90;2e3;2e3;2e3;3e3;3e3;3e3;3e3];
% sess_data(cur_tot+4).thom_elec = 1;
