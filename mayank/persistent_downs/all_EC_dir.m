%% overall EC database
clear all
drive_letter = 'C';
root_dir = strcat(drive_letter,':\WC_Germany\overall_EC');
cd(root_dir)
addpath('C:\WC_Germany\parietal_cortical_2010\')

%% L3 MEC pyramidal
cur_tot = 0;
n_mec_l3_pyr = 23;

%neuralynx data directory
sess_data(cur_tot+1).directory = strcat(drive_letter,':\wc_data\2006-04-07_CWC_LFP_C\2006-4-7_18-7-23');
sess_data(cur_tot+2).directory = strcat(drive_letter,':\wc_data\2006-04-07_CWC_LFP_D\2006-4-7_20-49-36');
sess_data(cur_tot+3).directory = strcat(drive_letter,':\wc_data\2006-04-08_CWC_LFP_A\2006-4-8_18-51-3');
sess_data(cur_tot+4).directory = strcat(drive_letter,':\wc_data\2007-05-24_CWC_LFP_B\2007-5-24_16-21-57');
sess_data(cur_tot+5).directory = strcat(drive_letter,':\wc_data\2007-05-30_CWC_LFP_C\2007-5-30_19-0-24');
sess_data(cur_tot+6).directory = strcat(drive_letter,':\wc_data\2007-05-31_CWC_LFP\2007-5-31_17-59-31');
sess_data(cur_tot+7).directory = strcat(drive_letter,':\wc_data\2007-06-03_CWC_LFP_B\2007-6-3_20-29-2');
sess_data(cur_tot+8).directory = strcat(drive_letter,':\wc_data\2007-06-04_CWC_LFP_B\2007-6-4_16-24-12');
sess_data(cur_tot+9).directory = strcat(drive_letter,':\wc_data\2007-06-26_CWC_LFP_C\2007-06-26_CWC_LFP_C\2007-6-26_19-52-31');
sess_data(cur_tot+10).directory = strcat(drive_letter,':\wc_data\2007-06-29_CWC_LFP\2007-06-29_CWC_LFP\2007-6-29_19-9-41');
sess_data(cur_tot+11).directory = strcat(drive_letter,':\wc_data\2007-07-03_CWC_LFP_A\2007-7-3_17-32-10');
sess_data(cur_tot+12).directory = strcat(drive_letter,':\wc_data\2007-07-03_CWC_LFP_B\2007-07-03_CWC_LFP_B\2007-7-3_18-37-6');
sess_data(cur_tot+13).directory = strcat(drive_letter,':\wc_data\2007-07-04_CWC_LFP_B\2007-07-04_CWC_LFP_B\2007-7-4_14-44-34');
sess_data(cur_tot+14).directory = strcat(drive_letter,':\wc_data\2007-07-04_CWC_LFP_C\2007-07-04_CWC_LFP_C\2007-7-4_16-57-57');
sess_data(cur_tot+15).directory = strcat(drive_letter,':\wc_data\2007-07-04_CWC_LFP_D\2007-07-04_CWC_LFP_D\2007-7-4_18-8-36');
sess_data(cur_tot+16).directory = strcat(drive_letter,':\wc_data\2007-07-05_CWC_LFP\2007-07-05_CWC_LFP\2007-7-5_17-41-58');
sess_data(cur_tot+17).directory = strcat(drive_letter,':\wc_data\2007-08-29_CWC_LFP_B\2007-08-29_CWC_LFP_B\2007-8-29_18-8-39');
sess_data(cur_tot+18).directory = strcat(drive_letter,':\wc_data\2009-04-07\2009-4-7-19');
sess_data(cur_tot+19).directory = strcat(drive_letter,':\wc_data\2009-04-13_A\2009-4-13-18');
sess_data(cur_tot+20).directory = strcat(drive_letter,':\wc_data\2009-04-13_B');
sess_data(cur_tot+21).directory = strcat(drive_letter,':\wc_data\2007-05-23_CWC_LFP_B\2007-5-23_19-4-42');
sess_data(cur_tot+22).directory = strcat(drive_letter,':\wc_data\2007-05-28_CWC_LFP_B\2007-5-28_19-38-35');
sess_data(cur_tot+23).directory = strcat(drive_letter,':\wc_data\2007-10-10_CWC_LFP_B\2007-10-10_CWC_LFP_B\2007-10-10_16-9-21');

%heka data directory
sess_data(cur_tot+1).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2006_04_07_CWC_LFP_C');
sess_data(cur_tot+2).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2006_04_07_CWC_LFP_D');
sess_data(cur_tot+3).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2006_04_08_CWC_LFP_A');
sess_data(cur_tot+4).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_05_24_CWC_LFP_B');
sess_data(cur_tot+5).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_05_30_CWC_LFP_C');
sess_data(cur_tot+6).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_05_31_CWC_LFP');
sess_data(cur_tot+7).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_06_03_CWC_LFP_B');
sess_data(cur_tot+8).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_06_04_CWC_LFP_B');
sess_data(cur_tot+9).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_06_26_CWC_LFP_C');
sess_data(cur_tot+10).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_06_29_CWC_LFP');
sess_data(cur_tot+11).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_07_03_CWC_LFP_A');
sess_data(cur_tot+12).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_07_03_CWC_LFP_B');
sess_data(cur_tot+13).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_07_04_CWC_LFP_B');
sess_data(cur_tot+14).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_07_04_CWC_LFP_C');
sess_data(cur_tot+15).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_07_04_CWC_LFP_D');
sess_data(cur_tot+16).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_07_05_CWC_LFP');
sess_data(cur_tot+17).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_08_29_CWC_LFP_B');
sess_data(cur_tot+18).heka_dir = strcat(drive_letter,':\wc_data\MPascii\2009-04-07_CWC_LFP');
sess_data(cur_tot+19).heka_dir = strcat(drive_letter,':\wc_data\MPascii\2009-04-13_CWC_LFP_A');
sess_data(cur_tot+20).heka_dir = strcat(drive_letter,':\wc_data\MPascii\2009-04-13_CWC_LFP_B');
sess_data(cur_tot+21).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007-05-23_CWC_LFP_B');
sess_data(cur_tot+22).heka_dir = '';
sess_data(cur_tot+23).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_10_10_CWC_LFP_B');

% cell locations
for i = 1:n_mec_l3_pyr
    sess_data(cur_tot+i).layer = '3';
    sess_data(cur_tot+i).cell_type = 'pyramidal';
    sess_data(cur_tot+i).region = 'MEC';
end

%data short name
sess_data(cur_tot+1).name = '2006-04-07_C';
sess_data(cur_tot+1).depth = 230;
sess_data(cur_tot+1).ant_post = nan;
sess_data(cur_tot+1).lateral = nan;
sess_data(cur_tot+1).dors_vent = nan;
sess_data(cur_tot+1).dist_prh = nan;
sess_data(cur_tot+1).gains = [50; 2e3; 2e3; 2e3; 2e3; nan; 3e3; 3e3];
sess_data(cur_tot+1).heka_type = 'sweep';
sess_data(cur_tot+1).class_certainty = 1;

sess_data(cur_tot+2).name = '2006-04-07_D';
sess_data(cur_tot+2).depth = 260;
sess_data(cur_tot+2).ant_post = nan;
sess_data(cur_tot+2).lateral = nan;
sess_data(cur_tot+2).dors_vent = nan;
sess_data(cur_tot+2).dist_prh = nan;
sess_data(cur_tot+2).gains = [50; 2e3; 2e3; 2e3; 2e3; nan; 3e3; 3e3];
sess_data(cur_tot+2).heka_type = 'sweep';
sess_data(cur_tot+2).class_certainty = 1;

sess_data(cur_tot+3).name = '2006-04-08_A';
sess_data(cur_tot+3).depth = 260;
sess_data(cur_tot+3).ant_post = nan;
sess_data(cur_tot+3).lateral = nan;
sess_data(cur_tot+3).dors_vent = nan;
sess_data(cur_tot+3).dist_prh = nan;
sess_data(cur_tot+3).gains = [40; 2e3; 2e3; 2e3; 2e3; nan; 3e3; 3e3];
sess_data(cur_tot+3).heka_type = 'sweep';
sess_data(cur_tot+3).class_certainty = 1;

sess_data(cur_tot+4).name = '2007-05-24_B';
sess_data(cur_tot+4).depth = 390;
sess_data(cur_tot+4).ant_post = -5;
sess_data(cur_tot+4).lateral = 3.7;
sess_data(cur_tot+4).dors_vent = 3;
sess_data(cur_tot+4).dist_prh = 500;
sess_data(cur_tot+4).gains = [100; 3e3; 3e3; nan; 4e3; 4e3; 4e3; 4e3];
sess_data(cur_tot+4).heka_type = 'sweep';
sess_data(cur_tot+4).class_certainty = 1;

sess_data(cur_tot+5).name = '2007-05-30_C';
sess_data(cur_tot+5).depth = 220;
sess_data(cur_tot+5).ant_post = -5;
sess_data(cur_tot+5).lateral = 3.85;
sess_data(cur_tot+5).dors_vent = 2.5;
sess_data(cur_tot+5).dist_prh = 200;
sess_data(cur_tot+5).gains = [100; 3e3; 3e3; nan; 4e3; 4e3; 4e3; 4e3];
sess_data(cur_tot+5).heka_type = 'sweep';
sess_data(cur_tot+5).class_certainty = 1;

sess_data(cur_tot+6).name = '2007-05-31';
sess_data(cur_tot+6).depth = 210;
sess_data(cur_tot+6).ant_post = -5;
sess_data(cur_tot+6).lateral = 3.55;
sess_data(cur_tot+6).dors_vent = 2.5;
sess_data(cur_tot+6).dist_prh = 500;
sess_data(cur_tot+6).gains = [100; 3e3; 3e3; nan; 4e3; 4e3; 4e3; 4e3];
sess_data(cur_tot+6).heka_type = 'sweep';
sess_data(cur_tot+6).class_certainty = 1;

sess_data(cur_tot+7).name = '2007-06-03_B';
sess_data(cur_tot+7).depth = 240;
sess_data(cur_tot+7).ant_post = -5;
sess_data(cur_tot+7).lateral = 3.4;
sess_data(cur_tot+7).dors_vent = 3;
sess_data(cur_tot+7).dist_prh = 600;
sess_data(cur_tot+7).gains = [90; 3e3; 3e3; nan; 4e3; 4e3; 4e3; 4e3];
sess_data(cur_tot+7).heka_type = 'sweep';
sess_data(cur_tot+7).class_certainty = 1;

sess_data(cur_tot+8).name = '2007-06-04_B';
sess_data(cur_tot+8).depth = 310;
sess_data(cur_tot+8).ant_post = -5;
sess_data(cur_tot+8).lateral = 3.4;
sess_data(cur_tot+8).dors_vent = 2.5;
sess_data(cur_tot+8).dist_prh = 500;
sess_data(cur_tot+8).gains = [120; 3e3; 3e3; nan; 4e3; 4e3; 4e3; 4e3];
sess_data(cur_tot+8).heka_type = 'sweep';
sess_data(cur_tot+8).class_certainty = 1;

sess_data(cur_tot+9).name = '2007-06-26_C';
sess_data(cur_tot+9).depth = 280;
sess_data(cur_tot+9).ant_post = -5;
sess_data(cur_tot+9).lateral = 3.4;
sess_data(cur_tot+9).dors_vent = 3;
sess_data(cur_tot+9).dist_prh = 600;
sess_data(cur_tot+9).gains = [100; 2e3; 2e3; nan; 5e3; 5e3; 5e3; 5e3];
sess_data(cur_tot+9).heka_type = 'sweep';
sess_data(cur_tot+9).class_certainty = 1;

sess_data(cur_tot+10).name = '2007-06-29';
sess_data(cur_tot+10).depth = 360;
sess_data(cur_tot+10).ant_post = -5;
sess_data(cur_tot+10).lateral = 4;
sess_data(cur_tot+10).dors_vent = 2.5;
sess_data(cur_tot+10).dist_prh = 300;
sess_data(cur_tot+10).gains = [80; 2e3; 2e3; nan; 5e3; 5e3; 5e3; 5e3];
sess_data(cur_tot+10).heka_type = 'sweep';
sess_data(cur_tot+10).class_certainty = 1;

sess_data(cur_tot+11).name = '2007-07-03_A';
sess_data(cur_tot+11).depth = 270;
sess_data(cur_tot+11).ant_post = -5;
sess_data(cur_tot+11).lateral = 3.85;
sess_data(cur_tot+11).dors_vent = 2.5;
sess_data(cur_tot+11).dist_prh = 300;
sess_data(cur_tot+11).gains = [100; 3e3; 3e3; nan; 4e3; 4e3; 4e3; 4e3];
sess_data(cur_tot+11).heka_type = 'sweep';
sess_data(cur_tot+11).class_certainty = 1;

sess_data(cur_tot+12).name = '2007-07-03_B';
sess_data(cur_tot+12).depth = 370;
sess_data(cur_tot+12).ant_post = -5;
sess_data(cur_tot+12).lateral = 3.4;
sess_data(cur_tot+12).dors_vent = 3;
sess_data(cur_tot+12).dist_prh = 500;
sess_data(cur_tot+12).gains = [70; 3e3; 3e3; nan; 4e3; 4e3; 4e3; 4e3];
sess_data(cur_tot+12).heka_type = 'sweep';
sess_data(cur_tot+12).class_certainty = 1;

sess_data(cur_tot+13).name = '2007-07-04_B';
sess_data(cur_tot+13).depth = 420;
sess_data(cur_tot+13).ant_post = -5;
sess_data(cur_tot+13).lateral = 3.7;
sess_data(cur_tot+13).dors_vent = 2.5;
sess_data(cur_tot+13).dist_prh = 50;
sess_data(cur_tot+13).gains = [80; 2e3; 2e3; nan; 3e3; 3e3; 3e3; 3e3];
sess_data(cur_tot+13).heka_type = 'sweep';
sess_data(cur_tot+13).class_certainty = 1;

sess_data(cur_tot+14).name = '2007-07-04_C';
sess_data(cur_tot+14).depth = 320;
sess_data(cur_tot+14).ant_post = -5;
sess_data(cur_tot+14).lateral = 3.85;
sess_data(cur_tot+14).dors_vent = 2.5;
sess_data(cur_tot+14).dist_prh = 300;
sess_data(cur_tot+14).gains = [100; 2e3; 2e3; nan; 3e3; 3e3; 3e3; 3e3];
sess_data(cur_tot+14).heka_type = 'sweep';
sess_data(cur_tot+14).class_certainty = 1;

sess_data(cur_tot+15).name = '2007-07-04_D';
sess_data(cur_tot+15).depth = 330;
sess_data(cur_tot+15).ant_post = -5;
sess_data(cur_tot+15).lateral = 3.4;
sess_data(cur_tot+15).dors_vent = 2.5;
sess_data(cur_tot+15).dist_prh = 300;
sess_data(cur_tot+15).gains = [80; 2e3; 2e3; nan; 3e3; 3e3; 3e3; 3e3];
sess_data(cur_tot+15).heka_type = 'sweep';
sess_data(cur_tot+15).class_certainty = 1;

sess_data(cur_tot+16).name = '2007-07-05';
sess_data(cur_tot+16).depth = 430;
sess_data(cur_tot+16).ant_post = -5;
sess_data(cur_tot+16).lateral = 3.7;
sess_data(cur_tot+16).dors_vent = 3;
sess_data(cur_tot+16).dist_prh = 600;
sess_data(cur_tot+16).gains = [100; 2e3; 2e3; nan; 4e3; 4e3; 4e3; 4e3];
sess_data(cur_tot+16).heka_type = 'sweep';
sess_data(cur_tot+16).class_certainty = 1;

sess_data(cur_tot+17).name = '2007-08-29_B';
sess_data(cur_tot+17).depth = 440;
sess_data(cur_tot+17).ant_post = -5;
sess_data(cur_tot+17).lateral = 4;
sess_data(cur_tot+17).dors_vent = 3;
sess_data(cur_tot+17).dist_prh = 500;
sess_data(cur_tot+17).gains = [90; 1e3; 2e3; nan; 4e3; 4e3; 4e3; 5e3];
sess_data(cur_tot+17).heka_type = 'sweep';
sess_data(cur_tot+17).class_certainty = 1;

sess_data(cur_tot+18).name = '2009_04-07';
sess_data(cur_tot+18).depth = 200;
sess_data(cur_tot+18).ant_post = nan;
sess_data(cur_tot+18).lateral = nan;
sess_data(cur_tot+18).dors_vent = nan;
sess_data(cur_tot+18).dist_prh = nan;
sess_data(cur_tot+18).gains = [50; 1e3; 2e3; 2e3; 2e3; 2e3; 2e3; 3e3];
sess_data(cur_tot+18).heka_type = 'cont';
sess_data(cur_tot+18).class_certainty = 1;

sess_data(cur_tot+19).name = '2009_04-13_A';
sess_data(cur_tot+19).depth = 250;
sess_data(cur_tot+19).ant_post = nan;
sess_data(cur_tot+19).lateral = nan;
sess_data(cur_tot+19).dors_vent = nan;
sess_data(cur_tot+19).dist_prh = nan;
sess_data(cur_tot+19).gains = [50; 1e3; 2e3; 2e3; 2e3; 2e3; 2e3; 3e3];
sess_data(cur_tot+19).heka_type = 'cont';
sess_data(cur_tot+19).class_certainty = 1;

sess_data(cur_tot+20).name = '2009_04-13_B';
sess_data(cur_tot+20).depth = 350;
sess_data(cur_tot+20).ant_post = nan;
sess_data(cur_tot+20).lateral = nan;
sess_data(cur_tot+20).dors_vent = nan;
sess_data(cur_tot+20).dist_prh = nan;
sess_data(cur_tot+20).gains = [50; 1e3; 2e3; 2e3; 2e3; 2e3; 2e3; 3e3];
sess_data(cur_tot+20).heka_type = 'cont';
sess_data(cur_tot+20).class_certainty = 1;

sess_data(cur_tot+21).name = '2007-05-23_B';
sess_data(cur_tot+21).depth = 210;
sess_data(cur_tot+21).ant_post = nan;
sess_data(cur_tot+21).lateral = nan;
sess_data(cur_tot+21).dors_vent = nan;
sess_data(cur_tot+21).dist_prh = nan;
sess_data(cur_tot+21).gains = [90; 3e3; 3e3; nan; 4e3; 4e3; 4e3; 4e3];
sess_data(cur_tot+21).heka_type = 'sweep';
sess_data(cur_tot+21).class_certainty = 1;
sess_data(cur_tot+21).cell_type = 'multipolar';

sess_data(cur_tot+22).name = '2007-05-28_B';
sess_data(cur_tot+22).depth = 210;
sess_data(cur_tot+22).ant_post = nan;
sess_data(cur_tot+22).lateral = nan;
sess_data(cur_tot+22).dors_vent = nan;
sess_data(cur_tot+22).dist_prh = nan;
sess_data(cur_tot+22).gains = [90; 3e3; 3e3; nan; 4e3; 4e3; 4e3; 4e3];
sess_data(cur_tot+22).heka_type = 'nan';
sess_data(cur_tot+22).class_certainty = 1;
sess_data(cur_tot+22).cell_type = 'multipolar';

sess_data(cur_tot+23).name = '2007-10-10_B';
sess_data(cur_tot+23).depth = 250;
sess_data(cur_tot+23).ant_post = nan;
sess_data(cur_tot+23).lateral = nan;
sess_data(cur_tot+23).dors_vent = nan;
sess_data(cur_tot+23).dist_prh = nan;
sess_data(cur_tot+23).gains = [100; nan; 1e3; nan; 4e3; 4e3; 4e3; 5e3];
sess_data(cur_tot+23).heka_type = 'sweep';
sess_data(cur_tot+23).class_certainty = 1;
sess_data(cur_tot+23).cell_type = 'multipolar';


%% MEC Layer 2 stellate 
n_mec_2_ste = 16;
cur_tot = cur_tot + n_mec_l3_pyr;

sess_data(cur_tot+1).directory = strcat(drive_letter,':\wc_data\2006-04-05_CWC_LFP_A\2006-4-5_19-10-27');
sess_data(cur_tot+2).directory = strcat(drive_letter,':\wc_data\2006-08-29_CWC_LFP\2006-8-29_20-46-14');
sess_data(cur_tot+3).directory = strcat(drive_letter,':\wc_data\2006-09-03_CWC_LFP\2006-9-3_20-2-20');
sess_data(cur_tot+4).directory = strcat(drive_letter,':\wc_data\2007-05-07_CWC_LFP\2007-5-7_17-52-33');
sess_data(cur_tot+5).directory = strcat(drive_letter,':\wc_data\2007-05-08_CWC_LFP_A\2007-5-8_13-21-1');
sess_data(cur_tot+6).directory = strcat(drive_letter,':\wc_data\2007-05-23_CWC_LFP_A\2007-5-23_18-10-12');
sess_data(cur_tot+7).directory = strcat(drive_letter,':\wc_data\2007-05-24_CWC_LFP_A\2007-5-24_12-48-28');
sess_data(cur_tot+8).directory = strcat(drive_letter,':\wc_data\2007-05-28_CWC_LFP_A\2007-5-28_18-27-1');
sess_data(cur_tot+9).directory = strcat(drive_letter,':\wc_data\2007-05-30_CWC_LFP_A\2007-5-30_13-48-19');
sess_data(cur_tot+10).directory = strcat(drive_letter,':\wc_data\2007-06-01_CWC_LFP_B\2007-6-1_15-56-41');
sess_data(cur_tot+11).directory = strcat(drive_letter,':\wc_data\2007-06-01_CWC_LFP_C\2007-6-1_19-7-5');
sess_data(cur_tot+12).directory = strcat(drive_letter,':\wc_data\2007-06-03_CWC_LFP_A\2007-6-3_19-1-5');
sess_data(cur_tot+13).directory = strcat(drive_letter,':\wc_data\2007-06-04_CWC_LFP_A\2007-6-4_14-46-17');
sess_data(cur_tot+14).directory = strcat(drive_letter,':\wc_data\2007-06-26_CWC_LFP_A\2007-6-26_14-23-29');
sess_data(cur_tot+15).directory = strcat(drive_letter,':\wc_data\2007-06-26_CWC_LFP_B\2007-06-26_CWC_LFP_B\2007-6-26_18-47-13');
sess_data(cur_tot+16).directory = strcat(drive_letter,':\wc_data\2007-07-04_CWC_LFP_A\2007-07-04_CWC_LFP_A\2007-7-4_13-21-34');

sess_data(cur_tot+1).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2006_04_05_CWC_LFP_A');
sess_data(cur_tot+2).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2006_08_29_CWC_LFP');
sess_data(cur_tot+3).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2006_09_03_CWC_LFP');
sess_data(cur_tot+4).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_05_07_CWC_LFP');
sess_data(cur_tot+5).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_05_08_CWC_LFP_A');
sess_data(cur_tot+6).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_05_23_CWC_LFP_A');
sess_data(cur_tot+7).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_05_24_CWC_LFP_A');
sess_data(cur_tot+8).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_05_28_CWC_LFP_A');
sess_data(cur_tot+9).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_05_30_CWC_LFP_A');
sess_data(cur_tot+10).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_06_01_CWC_LFP_B');
sess_data(cur_tot+11).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_06_01_CWC_LFP_C');
sess_data(cur_tot+12).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_06_03_CWC_LFP_A');
sess_data(cur_tot+13).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_06_04_CWC_LFP_A');
sess_data(cur_tot+14).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_06_26_CWC_LFP_A');
sess_data(cur_tot+15).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_06_26_CWC_LFP_B');
sess_data(cur_tot+16).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_07_04_CWC_LFP_A');

% cell locations
for i = 1:n_mec_2_ste
    sess_data(cur_tot+i).layer = '2';
    sess_data(cur_tot+i).cell_type = 'stellate';
    sess_data(cur_tot+i).region = 'MEC';
end

%data short name
sess_data(cur_tot+1).name = '2006-04-05_A';
sess_data(cur_tot+1).depth = 160;
sess_data(cur_tot+1).ant_post = nan;
sess_data(cur_tot+1).lateral = nan;
sess_data(cur_tot+1).dors_vent = nan;
sess_data(cur_tot+1).dist_prh = nan;
sess_data(cur_tot+1).gains = [50; 1e3; 2e3; 2e3; 2e3; nan; 3e3; 3e3];
sess_data(cur_tot+1).heka_type = 'sweep';
sess_data(cur_tot+1).class_certainty = 1;

sess_data(cur_tot+2).name = '2006-08-29';
sess_data(cur_tot+2).depth = 160;
sess_data(cur_tot+2).ant_post = nan;
sess_data(cur_tot+2).lateral = nan;
sess_data(cur_tot+2).dors_vent = nan;
sess_data(cur_tot+2).dist_prh = nan;
sess_data(cur_tot+2).gains = [100; 1e3; 2e3; 2e3; 2e3; nan; 3e3; 3e3];
sess_data(cur_tot+2).heka_type = 'sweep';
sess_data(cur_tot+2).class_certainty = 0;

sess_data(cur_tot+3).name = '2006-09-03';
sess_data(cur_tot+3).depth = 230;
sess_data(cur_tot+3).ant_post = nan;
sess_data(cur_tot+3).lateral = nan;
sess_data(cur_tot+3).dors_vent = nan;
sess_data(cur_tot+3).dist_prh = nan;
sess_data(cur_tot+3).gains = [100; 2e3; 2e3; 2e3; 2e3; nan; 3e3; 3e3];
sess_data(cur_tot+3).heka_type = 'sweep';
sess_data(cur_tot+3).class_certainty = 1;
sess_data(cur_tot+3).layer = '3';

sess_data(cur_tot+4).name = '2007-05-07';
sess_data(cur_tot+4).depth = 110;
sess_data(cur_tot+4).ant_post = -5;
sess_data(cur_tot+4).lateral = 4;
sess_data(cur_tot+4).dors_vent = 2.5;
sess_data(cur_tot+4).dist_prh = 100;
sess_data(cur_tot+4).gains = [100; 3e3; 3e3; nan; 4e3; 4e3; 4e3; 4e3];
sess_data(cur_tot+4).heka_type = 'sweep';
sess_data(cur_tot+4).class_certainty = 1;

sess_data(cur_tot+5).name = '2007-05-08_A';
sess_data(cur_tot+5).depth = 170;
sess_data(cur_tot+5).ant_post = -5;
sess_data(cur_tot+5).lateral = 3.85;
sess_data(cur_tot+5).dors_vent = 2.5;
sess_data(cur_tot+5).dist_prh = 200;
sess_data(cur_tot+5).gains = [100; 3e3; 3e3; nan; 4e3; 4e3; 4e3; 4e3];
sess_data(cur_tot+5).heka_type = 'sweep';
sess_data(cur_tot+5).class_certainty = 1;

sess_data(cur_tot+6).name = '2007-05-23_A';
sess_data(cur_tot+6).depth = 220;
sess_data(cur_tot+6).ant_post = -5;
sess_data(cur_tot+6).lateral = 4;
sess_data(cur_tot+6).dors_vent = 2.5;
sess_data(cur_tot+6).dist_prh = 200;
sess_data(cur_tot+6).gains = [150; 3e3; 3e3; nan; 4e3; 4e3; 4e3; 4e3];
sess_data(cur_tot+6).heka_type = 'sweep';
sess_data(cur_tot+6).class_certainty = 0;

sess_data(cur_tot+7).name = '2007-05-24_A';
sess_data(cur_tot+7).depth = 120;
sess_data(cur_tot+7).ant_post = -5;
sess_data(cur_tot+7).lateral = 3.85;
sess_data(cur_tot+7).dors_vent = 2.5;
sess_data(cur_tot+7).dist_prh = 200;
sess_data(cur_tot+7).gains = [100; 3e3; 3e3; nan; 4e3; 4e3; 4e3; 4e3];
sess_data(cur_tot+7).heka_type = 'sweep';
sess_data(cur_tot+7).class_certainty = 1;

sess_data(cur_tot+8).name = '2007-05-28_A';
sess_data(cur_tot+8).depth = 140;
sess_data(cur_tot+8).ant_post = -5;
sess_data(cur_tot+8).lateral = 4;
sess_data(cur_tot+8).dors_vent = 2.5;
sess_data(cur_tot+8).dist_prh = 100;
sess_data(cur_tot+8).gains = [150; 3e3; 3e3; nan; 4e3; 4e3; 4e3; 4e3];
sess_data(cur_tot+8).heka_type = 'sweep';
sess_data(cur_tot+8).class_certainty = 1;

sess_data(cur_tot+9).name = '2007-05-30_A';
sess_data(cur_tot+9).depth = 120;
sess_data(cur_tot+9).ant_post = -5;
sess_data(cur_tot+9).lateral = 3.85;
sess_data(cur_tot+9).dors_vent = 2.5;
sess_data(cur_tot+9).dist_prh = 400;
sess_data(cur_tot+9).gains = [120; 3e3; 3e3; nan; 4e3; 4e3; 4e3; 4e3];
sess_data(cur_tot+9).heka_type = 'sweep';
sess_data(cur_tot+9).class_certainty = 1;

sess_data(cur_tot+10).name = '2007-06-01_B';
sess_data(cur_tot+10).depth = 170;
sess_data(cur_tot+10).ant_post = -5;
sess_data(cur_tot+10).lateral = 3.55;
sess_data(cur_tot+10).dors_vent = 2.5;
sess_data(cur_tot+10).dist_prh = 500;
sess_data(cur_tot+10).gains = [150; 3e3; 3e3; nan; 5e3; 5e3; 5e3; 5e3];
sess_data(cur_tot+10).heka_type = 'sweep';
sess_data(cur_tot+10).class_certainty = 1;

sess_data(cur_tot+11).name = '2007-06-01_C';
sess_data(cur_tot+11).depth = 180;
sess_data(cur_tot+11).ant_post = -5;
sess_data(cur_tot+11).lateral = 3.55;
sess_data(cur_tot+11).dors_vent = 2.5;
sess_data(cur_tot+11).dist_prh = 100;
sess_data(cur_tot+11).gains = [120; 3e3; 3e3; nan; 5e3; 5e3; 5e3; 5e3];
sess_data(cur_tot+11).heka_type = 'sweep';
sess_data(cur_tot+11).class_certainty = 1;

sess_data(cur_tot+12).name = '2007-06-03_A';
sess_data(cur_tot+12).depth = 180;
sess_data(cur_tot+12).ant_post = -5;
sess_data(cur_tot+12).lateral = 4;
sess_data(cur_tot+12).dors_vent = 2.5;
sess_data(cur_tot+12).dist_prh = 100;
sess_data(cur_tot+12).gains = [100; 3e3; 3e3; nan; 4e3; 4e3; 4e3; 4e3];
sess_data(cur_tot+12).heka_type = 'sweep';
sess_data(cur_tot+12).class_certainty = 1;

sess_data(cur_tot+13).name = '2007-06-04_A';
sess_data(cur_tot+13).depth = 170;
sess_data(cur_tot+13).ant_post = -5;
sess_data(cur_tot+13).lateral = 4;
sess_data(cur_tot+13).dors_vent = 2.5;
sess_data(cur_tot+13).dist_prh = 300;
sess_data(cur_tot+13).gains = [100; 3e3; 3e3; nan; 4e3; 4e3; 4e3; 4e3];
sess_data(cur_tot+13).heka_type = 'sweep';
sess_data(cur_tot+13).class_certainty = 1;

sess_data(cur_tot+14).name = '2007-06-26_A';
sess_data(cur_tot+14).depth = 200;
sess_data(cur_tot+14).ant_post = -5;
sess_data(cur_tot+14).lateral = 4;
sess_data(cur_tot+14).dors_vent = 2.5;
sess_data(cur_tot+14).dist_prh = 200;
sess_data(cur_tot+14).gains = [80; 3e3; 3e3; nan; 4e3; 4e3; 4e3; 4e3];
sess_data(cur_tot+14).heka_type = 'sweep';
sess_data(cur_tot+14).class_certainty = 0;
sess_data(cur_tot+14).layer = '3';
sess_data(cur_tot+14).cell_type = 'pyramidal';

sess_data(cur_tot+15).name = '2007-06-26_B';
sess_data(cur_tot+15).depth = 180;
sess_data(cur_tot+15).ant_post = -5;
sess_data(cur_tot+15).lateral = 3.85;
sess_data(cur_tot+15).dors_vent = 2.5;
sess_data(cur_tot+15).dist_prh = 300;
sess_data(cur_tot+15).gains = [100; 2e3; 2e3; nan; 5e3; 5e3; 5e3; 5e3];
sess_data(cur_tot+15).heka_type = 'sweep';
sess_data(cur_tot+15).class_certainty = 1;

sess_data(cur_tot+16).name = '2007-07-04_A';
sess_data(cur_tot+16).depth = 200;
sess_data(cur_tot+16).ant_post = -5;
sess_data(cur_tot+16).lateral = 4.15;
sess_data(cur_tot+16).dors_vent = 2.5;
sess_data(cur_tot+16).dist_prh = 100;
sess_data(cur_tot+16).gains = [110; 2e3; 2e3; nan; 3e3; 3e3; 3e3; 3e3];
sess_data(cur_tot+16).heka_type = 'sweep';
sess_data(cur_tot+16).class_certainty = 1;


%% Layer 5 MEC
n_mec_5_pyr = 19;
cur_tot = cur_tot + n_mec_2_ste;
c = 1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2007-08-27_CWC_LFP\2007-8-27_17-8-10');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_08_27_CWC_LFP');
sess_data(cur_tot+c).name = '2007-08-27';
sess_data(cur_tot+c).depth = nan;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [100; 1e3; 1e3; nan; 4e3; 4e3; 4e3; 5e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '3/5';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).class_certainty = 0;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2007-08-28_CWC_LFP_A\2007-8-28_16-5-16');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_08_28_CWC_LFP_A');
sess_data(cur_tot+c).name = '2007-08-28_A';
sess_data(cur_tot+c).depth = 520;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [100; 2e3; 2e3; nan; 4e3; 4e3; 4e3; 5e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '5/6';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2007-08-28_CWC_LFP_B\2007-8-28_18-37-38');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_08_28_CWC_LFP_B');
sess_data(cur_tot+c).name = '2007-08-28_B';
sess_data(cur_tot+c).depth = 540;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [100; 2e3; 2e3; nan; 4e3; 4e3; 4e3; 5e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '3/5';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2007-09-25_CWC_LFP_A\2007-9-25_14-53-25');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_09_25_CWC_LFP_A');
sess_data(cur_tot+c).name = '2007-09-25_A';
sess_data(cur_tot+c).depth = 750;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [100; 1e3; 2e3; nan; 4e3; 4e3; 4e3; 5e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '5/6';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2007-09-26_CWC_LFP_B\2007-9-26_17-16-35');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_09_26_CWC_LFP_B');
sess_data(cur_tot+c).name = '2007-09-26_B';
sess_data(cur_tot+c).depth = 700;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [100; 1e3; 2e3; nan; 4e3; 4e3; 4e3; 5e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '5/6';
sess_data(cur_tot+c).cell_type = 'unclear';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).class_certainty = 0;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2007-09-27_CWC_LFP\2007-9-27_18-6-44');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_09_27_CWC_LFP');
sess_data(cur_tot+c).name = '2007-09-27';
sess_data(cur_tot+c).depth = nan;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [100; 1e3; 2e3; nan; 4e3; 4e3; 4e3; 5e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '5/6';
sess_data(cur_tot+c).cell_type = 'unclear';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).class_certainty = 0;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2007-09-28_CWC_LFP_A\2007-9-28_15-11-31');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_09_28_CWC_LFP_A');
sess_data(cur_tot+c).name = '2007-09-28_A';
sess_data(cur_tot+c).depth = 880;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [100; nan; nan; nan; 4e3; 4e3; 4e3; 5e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '5/6';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2007-10-02_CWC_LFP\2007-10-2_17-13-4');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_10_02_CWC_LFP');
sess_data(cur_tot+c).name = '2007-10-02';
sess_data(cur_tot+c).depth = 590;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [100; nan; 1e3; nan; 4e3; 4e3; 4e3; 5e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '3/5';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2007-10-09_CWC_LFP_A\2007-10-9_14-27-18');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_10_09_CWC_LFP_A');
sess_data(cur_tot+c).name = '2007-10-09_A';
sess_data(cur_tot+c).depth = 680;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [90; nan; 1e3; nan; 4e3; 4e3; 4e3; 5e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '3/5';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2007-10-09_CWC_LFP_B\2007-10-9_17-14-16');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_10_09_CWC_LFP_B');
sess_data(cur_tot+c).name = '2007-10-09_B';
sess_data(cur_tot+c).depth = 500;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [100; nan; 1e3; nan; 4e3; 4e3; 4e3; 5e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '3/5';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2007-10-10_CWC_LFP_A\2007-10-10_15-1-25');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_10_10_CWC_LFP_A');
sess_data(cur_tot+c).name = '2007-10-10_A';
sess_data(cur_tot+c).depth = nan;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [100; nan; 1e3; nan; 4e3; 4e3; 4e3; 5e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '5/6';
sess_data(cur_tot+c).cell_type = 'unclear';
sess_data(cur_tot+c).region = 'isocortex';
sess_data(cur_tot+c).class_certainty = 0;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2007-10-12_CWC_LFP\2007-10-12_18-52-11');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_10_12_CWC_LFP');
sess_data(cur_tot+c).name = '2007-10-12';
sess_data(cur_tot+c).depth = 580;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [100; nan; 1e3; nan; 4e3; 4e3; 4e3; 5e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '5/6';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2007-10-23_CWC_LFP_A\2007-10-23_12-56-23');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_10_23_CWC_LFP_A');
sess_data(cur_tot+c).name = '2007-10-23_A';
sess_data(cur_tot+c).depth = 810;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [80; 1e3; 1e3; nan; 4e3; 4e3; 4e3; 5e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '5/6';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2007-10-23_CWC_LFP_B\2007-10-23_16-4-41');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_10_23_CWC_LFP_B');
sess_data(cur_tot+c).name = '2007-10-23_B';
sess_data(cur_tot+c).depth = 630;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [100; 1e3; 1e3; nan; 4e3; 4e3; 4e3; 5e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '3/5';
sess_data(cur_tot+c).cell_type = 'horizontal';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2007-10-25_CWC_LFP_A\2007-10-25_13-26-46');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_10_25_CWC_LFP_A');
sess_data(cur_tot+c).name = '2007-10-25_A';
sess_data(cur_tot+c).depth = 690;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [100; 1e3; 1e3; nan; 4e3; 4e3; 4e3; 5e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '5/6';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).class_certainty = 0;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2007-10-25_CWC_LFP_B\2007-10-25_14-13-6');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_10_25_CWC_LFP_B');
sess_data(cur_tot+c).name = '2007-10-25_B';
sess_data(cur_tot+c).depth = nan;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [100; 1e3; 1e3; nan; 4e3; 4e3; 4e3; 5e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = 'unclear';
sess_data(cur_tot+c).cell_type = 'unclear';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).class_certainty = 0;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2007-10-25_CWC_LFP_C\2007-10-25_16-16-54');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_10_25_CWC_LFP_C');
sess_data(cur_tot+c).name = '2007-10-25_C';
sess_data(cur_tot+c).depth = 780;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [100; 1e3; 1e3; nan; 4e3; 4e3; 4e3; 5e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '5/6';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2007-10-26_CWC_LFP_A\2007-10-26_14-29-11');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_10_26_CWC_LFP_A');
sess_data(cur_tot+c).name = '2007-10-26_A';
sess_data(cur_tot+c).depth = 710;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [90; 1e3; 1e3; nan; 4e3; 4e3; 4e3; 5e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '5/6';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2007-10-26_CWC_LFP_B\2007-10-26_17-20-40');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007_10_26_CWC_LFP_B');
sess_data(cur_tot+c).name = '2007-10-26_B';
sess_data(cur_tot+c).depth = 600;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [90; 1e3; 1e3; nan; 4e3; 4e3; 4e3; 5e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '5/6';
sess_data(cur_tot+c).cell_type = 'unclear';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).class_certainty = 0;
c = c+1;

%% LEC
n_lec = 21;
cur_tot = cur_tot + n_mec_5_pyr;
c = 1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2005-01-31_HC_LFP\2005-1-31_17-55-1');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2005-01-31_HC_LFP');
sess_data(cur_tot+c).name = '2005_01-31';
sess_data(cur_tot+c).depth = 275;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [90; 2e3; 2e3; 2e3; 2e3; 2e3; 2e3; 2e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '3';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'LEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2005-02-03_HC_LFP\2005-2-3_17-11-9');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2005-02-03_HC_LFP');
sess_data(cur_tot+c).name = '2005_02-03';
sess_data(cur_tot+c).depth = 100;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [80; 2e3; 2e3; 2e3; 2e3; 2e3; 2e3; 2e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '2';
sess_data(cur_tot+c).cell_type = 'fan';
sess_data(cur_tot+c).region = 'LEC';
sess_data(cur_tot+c).class_certainty = 0;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2005-11-29_CWC_LFP_C\2005-11-29_18-12-45');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2005-11-29_CWC_LFP_C');
sess_data(cur_tot+c).name = '2005_11-29_C';
sess_data(cur_tot+c).depth = 470;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [80; 1e3; 2e3; 1e3; 1e3; nan; 1e3; 1e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '3/5';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'LEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2005-12-09_CWC_LFP_C\2005-12-9_20-43-55');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2005-12-09_CWC_LFP_C');
sess_data(cur_tot+c).name = '2005_12-09_C';
sess_data(cur_tot+c).depth = nan;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [70; 1e3; 2e3; 1e3; 1e3; nan; 2e3; 2e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = 'nan';
sess_data(cur_tot+c).cell_type = 'nan';
sess_data(cur_tot+c).region = 'LEC';
sess_data(cur_tot+c).class_certainty = 0;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2005-12-09_CWC_LFP_D\2005-12-9_22-4-2');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2005-12-09_CWC_LFP_D');
sess_data(cur_tot+c).name = '2005_12-09_D';
sess_data(cur_tot+c).depth = nan;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [70; 1e3; 2e3; 1e3; 1e3; nan; 2e3; 2e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = 'nan';
sess_data(cur_tot+c).cell_type = 'nan';
sess_data(cur_tot+c).region = 'LEC';
sess_data(cur_tot+c).class_certainty = 0;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2005-12-12_CWC_LFP_C\2005-12-12_18-53-18');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2005-12-12_CWC_LFP_C');
sess_data(cur_tot+c).name = '2005_12-12_C';
sess_data(cur_tot+c).depth = 280;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [50; 1e3; 2e3; 1e3; 1e3; nan; 2e3; 2e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '3';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'LEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2005-11-25_CWC_LFP_A\2005-11-25_19-58-1');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2005-11-25_CWC_LFP_A');
sess_data(cur_tot+c).name = '2005_11-25_A';
sess_data(cur_tot+c).depth = 240;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [80; 2e3; 2e3; 2e3; 2e3; nan; 2e3; 2e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '3';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'LEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2005-12-09_CWC_LFP_A\2005-12-9_12-43-2');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2005-12-09_CWC_LFP_A');
sess_data(cur_tot+c).name = '2005_12-09_A';
sess_data(cur_tot+c).depth = 260;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [100; 1e3; 2e3; 1e3; 1e3; nan; 1e3; 1e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '3';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'LEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2005-12-12_CWC_LFP_A\2005-12-12_13-46-24');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2005-12-12_CWC_LFP_A');
sess_data(cur_tot+c).name = '2005_12-12_A';
sess_data(cur_tot+c).depth = 230;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [80; 1e3; 2e3; 1e3; 1e3; nan; 2e3; 2e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '3';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'LEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2005-12-09_CWC_LFP_B\2005-12-9_14-8-51');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2005-12-09_CWC_LFP_B');
sess_data(cur_tot+c).name = '2005_12-09_B';
sess_data(cur_tot+c).depth = 190;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [60; 1e3; 2e3; 1e3; 1e3; nan; 1e3; 1e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '2';
sess_data(cur_tot+c).cell_type = 'multipolar';
sess_data(cur_tot+c).region = 'LEC';
sess_data(cur_tot+c).class_certainty = 0;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2005-12-12_CWC_LFP_B\2005-12-12_17-41-40');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2005-12-12_CWC_LFP_B');
sess_data(cur_tot+c).name = '2005_12-12_B';
sess_data(cur_tot+c).depth = 260;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [50; 1e3; 2e3; 1e3; 1e3; nan; 2e3; 2e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '3';
sess_data(cur_tot+c).cell_type = 'multipolar';
sess_data(cur_tot+c).region = 'LEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2005-11-28_CWC_LFP\2005-11-28_17-11-19');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2005-11-28_CWC_LFP');
sess_data(cur_tot+c).name = '2005_11-28';
sess_data(cur_tot+c).depth = 170;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [80; 1e3; 2e3; 1e3; 1e3; nan; 1e3; 1e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '2';
sess_data(cur_tot+c).cell_type = 'fan';
sess_data(cur_tot+c).region = 'LEC';
sess_data(cur_tot+c).class_certainty = 0;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2005-11-29_CWC_LFP_A\2005-11-29_14-27-37');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2005-11-29_CWC_LFP_A');
sess_data(cur_tot+c).name = '2005_11-29_A';
sess_data(cur_tot+c).depth = 160;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [80; 1e3; 2e3; 1e3; 1e3; nan; 1e3; 1e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '2';
sess_data(cur_tot+c).cell_type = 'multipolar';
sess_data(cur_tot+c).region = 'LEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2005-11-29_CWC_LFP_B\2005-11-29_17-47-36');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2005-11-29_CWC_LFP_B');
sess_data(cur_tot+c).name = '2005_11-29_B';
sess_data(cur_tot+c).depth = 190;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [80; 1e3; 2e3; 1e3; 1e3; nan; 1e3; 1e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '2';
sess_data(cur_tot+c).cell_type = 'multipolar';
sess_data(cur_tot+c).region = 'LEC';
sess_data(cur_tot+c).class_certainty = 0;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2007-01-17_CWC_LFP\2007-1-17_17-42-58');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007-01-17_CWC_LFP');
sess_data(cur_tot+c).name = '2007_01-17';
sess_data(cur_tot+c).depth = 110;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [100; 2e3; 2e3; 2e3; 2e3; 3e3; 3e3; 3e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '2';
sess_data(cur_tot+c).cell_type = 'fan';
sess_data(cur_tot+c).region = 'LEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2009-04-21_CWC_LFP\2009-4-21_22-20-11');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\2009-04-21_CWC_LFP');
sess_data(cur_tot+c).name = '2009_04-21';
sess_data(cur_tot+c).depth = 260;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [50; 1e3; 2e3; 2e3; 2e3; 2e3; 2e3; 3e3];
sess_data(cur_tot+c).heka_type = 'cont';
sess_data(cur_tot+c).layer = '3';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'LEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2009-05-10_CWC_LFP_A\2009-5-10_16-56-59');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\2009-05-10_CWC_LFP_A');
sess_data(cur_tot+c).name = '2009_05-10_A';
sess_data(cur_tot+c).depth = 300;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [50; 1e3; 2e3; 2e3; 2e3; 2e3; 2e3; 3e3];
sess_data(cur_tot+c).heka_type = 'cont';
sess_data(cur_tot+c).layer = '3';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'LEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2009-05-10_CWC_LFP_B\2009-5-10_18-6-22');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\2009-05-10_CWC_LFP_B');
sess_data(cur_tot+c).name = '2009_05-10_B';
sess_data(cur_tot+c).depth = 440;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [50; 1e3; 2e3; 2e3; 2e3; 2e3; 2e3; 3e3];
sess_data(cur_tot+c).heka_type = 'cont';
sess_data(cur_tot+c).layer = '3';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'LEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2009-05-16_CWC_LFP\2009-5-16_16-43-50');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\2009-05-16_CWC_LFP');
sess_data(cur_tot+c).name = '2009_05-16';
sess_data(cur_tot+c).depth = 340;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [50; 1e3; 2e3; 2e3; 2e3; 2e3; 2e3; 3e3];
sess_data(cur_tot+c).heka_type = 'cont';
sess_data(cur_tot+c).layer = '3';
sess_data(cur_tot+c).cell_type = 'multipolar';
sess_data(cur_tot+c).region = 'LEC';
sess_data(cur_tot+c).class_certainty = 0;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2009-04-24_CWC_LFP\2009-4-24_22-16-59');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\2009-04-24_CWC_LFP');
sess_data(cur_tot+c).name = '2009_04-24';
sess_data(cur_tot+c).depth = 180;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [50; 1e3; 2e3; 2e3; 2e3; 2e3; 2e3; 3e3];
sess_data(cur_tot+c).heka_type = 'cont';
sess_data(cur_tot+c).layer = '2';
sess_data(cur_tot+c).cell_type = 'multipolar';
sess_data(cur_tot+c).region = 'LEC';
sess_data(cur_tot+c).class_certainty = 0;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2009-04-29_CWC_LFP_A\2009-4-29_15-57-12');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\2009-04-29_CWC_LFP_A');
sess_data(cur_tot+c).name = '2009_04-29_A';
sess_data(cur_tot+c).depth = 190;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [50; 1e3; 2e3; 2e3; 2e3; 2e3; 2e3; 3e3];
sess_data(cur_tot+c).heka_type = 'cont';
sess_data(cur_tot+c).layer = '2';
sess_data(cur_tot+c).cell_type = 'multipolar';
sess_data(cur_tot+c).region = 'LEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2009-04-29_CWC_LFP_B\2009-4-29_16-56-26');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\2009-04-29_CWC_LFP_B');
sess_data(cur_tot+c).name = '2009_04-29_B';
sess_data(cur_tot+c).depth = 150;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [50; 1e3; 2e3; 2e3; 2e3; 2e3; 2e3; 3e3];
sess_data(cur_tot+c).heka_type = 'cont';
sess_data(cur_tot+c).layer = '2';
sess_data(cur_tot+c).cell_type = 'fan';
sess_data(cur_tot+c).region = 'LEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2009-12-01_CWC_LFP');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\2009-12-01_CWC_LFP');
sess_data(cur_tot+c).name = '2009_12-01';
sess_data(cur_tot+c).depth = 490;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [50; 1e3; 2e3; 2e3; 2e3; 2e3; 2e3; 3e3];
sess_data(cur_tot+c).heka_type = 'cont';
sess_data(cur_tot+c).layer = '3/5';
sess_data(cur_tot+c).cell_type = 'multipolar';
sess_data(cur_tot+c).region = 'LEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2009-12-02_CWC_LFP');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\2009-12-02_CWC_LFP');
sess_data(cur_tot+c).name = '2009_12-02';
sess_data(cur_tot+c).depth = 300;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [50; 1e3; 2e3; 2e3; 2e3; 2e3; 2e3; 3e3];
sess_data(cur_tot+c).heka_type = 'cont';
sess_data(cur_tot+c).layer = '3';
sess_data(cur_tot+c).cell_type = 'multipolar';
sess_data(cur_tot+c).region = 'LEC';
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2010-05-16_CWC_LFP');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\2010-05-16_CWC_LFP');
sess_data(cur_tot+c).name = '2010-05-16';
sess_data(cur_tot+c).depth = 270;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [50; 1e3; 2e3; 2e3; 2e3; 2e3; 2e3; 3e3];
sess_data(cur_tot+c).heka_type = 'cont';
sess_data(cur_tot+c).layer = '2/3';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'isocortex';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2010-05-26_CWC_LFP');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\2010-05-26_CWC_LFP');
sess_data(cur_tot+c).name = '2010-05-26';
sess_data(cur_tot+c).depth = 100;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [50; 1e3; 2e3; 2e3; 2e3; 2e3; 2e3; 3e3];
sess_data(cur_tot+c).heka_type = 'cont';
sess_data(cur_tot+c).layer = '2';
sess_data(cur_tot+c).cell_type = 'fan';
sess_data(cur_tot+c).region = 'LEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2010-05-27_CWC_LFP_A');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\2010-05-27_CWC_LFP_A');
sess_data(cur_tot+c).name = '2010-05-27_A';
sess_data(cur_tot+c).depth = 250;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [50; 1e3; 2e3; 2e3; 2e3; 2e3; 2e3; 3e3];
sess_data(cur_tot+c).heka_type = 'cont';
sess_data(cur_tot+c).layer = '3';
sess_data(cur_tot+c).cell_type = 'unclear';
sess_data(cur_tot+c).region = 'LEC';
sess_data(cur_tot+c).class_certainty = 0;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2010-05-27_CWC_LFP_B');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\2010-05-27_CWC_LFP_B');
sess_data(cur_tot+c).name = '2010-05-27_B';
sess_data(cur_tot+c).depth = 190;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [50; 1e3; 2e3; 2e3; 2e3; 2e3; 2e3; 3e3];
sess_data(cur_tot+c).heka_type = 'cont';
sess_data(cur_tot+c).layer = '3';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'LEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2010-05-29_CWC_LFP');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\2010-05-29_CWC_LFP');
sess_data(cur_tot+c).name = '2010-05-29';
sess_data(cur_tot+c).depth = 230;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [50; 2e3; 2e3; 2e3; 2e3; 2e3; 2e3; 3e3];
sess_data(cur_tot+c).heka_type = 'cont';
sess_data(cur_tot+c).layer = '3';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'LEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2010-07-31_CWC_LFP_A\2010-7-31_16-49-36');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\2010-07-31_CWC_LFP_A');
sess_data(cur_tot+c).name = '2010-07-31_A';
sess_data(cur_tot+c).depth = 80;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [50; 2e3; 5e3; 2e3; 2e3; 2e3; 2e3; 3e3];
sess_data(cur_tot+c).heka_type = 'cont';
sess_data(cur_tot+c).layer = '2';
sess_data(cur_tot+c).cell_type = 'unclear';
sess_data(cur_tot+c).region = 'LEC';
sess_data(cur_tot+c).class_certainty = 0;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2010-07-31_CWC_LFP_B\2010-7-31_17-36-25');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\2010-07-31_CWC_LFP_B');
sess_data(cur_tot+c).name = '2010-07-31_B';
sess_data(cur_tot+c).depth = 130;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [50; 2e3; 5e3; 2e3; 2e3; 2e3; 2e3; 3e3];
sess_data(cur_tot+c).heka_type = 'cont';
sess_data(cur_tot+c).layer = '2';
sess_data(cur_tot+c).cell_type = 'fan';
sess_data(cur_tot+c).region = 'LEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2010-07-31_CWC_LFP_C\2010-7-31_20-55-55');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\2010-07-31_CWC_LFP_C');
sess_data(cur_tot+c).name = '2010-07-31_C';
sess_data(cur_tot+c).depth = 220;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [50; 2e3; 2e3; 2e3; 2e3; 2e3; 2e3; 3e3];
sess_data(cur_tot+c).heka_type = 'cont';
sess_data(cur_tot+c).layer = '3';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'LEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2010-08-01_CWC_LFP_A\2010-8-1_13-56-56');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\2010-08-01_CWC_LFP_A');
sess_data(cur_tot+c).name = '2010-08-01_A';
sess_data(cur_tot+c).depth = 280;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [50; 2e3; 5e3; 2e3; 2e3; 2e3; 2e3; 3e3];
sess_data(cur_tot+c).heka_type = 'cont';
sess_data(cur_tot+c).layer = '3';
sess_data(cur_tot+c).cell_type = 'unclear';
sess_data(cur_tot+c).region = 'LEC';
sess_data(cur_tot+c).class_certainty = 0;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2010-08-01_CWC_LFP_B\2010-8-1_16-28-28');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\2010-08-01_CWC_LFP_B');
sess_data(cur_tot+c).name = '2010-08-01_B';
sess_data(cur_tot+c).depth = 180;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [50; 2e3; 5e3; 2e3; 2e3; 2e3; 2e3; 3e3];
sess_data(cur_tot+c).heka_type = 'cont';
sess_data(cur_tot+c).layer = '2';
sess_data(cur_tot+c).cell_type = 'multipolar';
sess_data(cur_tot+c).region = 'LEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2010-08-05_CWC_LFP_A\2010-8-5_17-21-0');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\2010-08-05_CWC_LFP_A');
sess_data(cur_tot+c).name = '2010-08-05_A';
sess_data(cur_tot+c).depth = 180;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [50; 2e3; 5e3; 2e3; 2e3; 2e3; 2e3; 3e3];
sess_data(cur_tot+c).heka_type = 'cont';
sess_data(cur_tot+c).layer = '2';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'LEC';
sess_data(cur_tot+c).class_certainty = 0;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2010-08-05_CWC_LFP_B\2010-8-5_18-30-4');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\2010-08-05_CWC_LFP_B');
sess_data(cur_tot+c).name = '2010-08-05_B';
sess_data(cur_tot+c).depth = 410;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [50; 2e3; 5e3; 2e3; 2e3; 2e3; 2e3; 3e3];
sess_data(cur_tot+c).heka_type = 'cont';
sess_data(cur_tot+c).layer = '3';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'LEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2010-08-12_CWC_LFP_A\2010-8-12_14-8-8');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\2010-08-12_CWC_LFP_A');
sess_data(cur_tot+c).name = '2010-08-12_A';
sess_data(cur_tot+c).depth = 250;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [50; 2e3; 5e3; 2e3; 2e3; 2e3; 2e3; 3e3];
sess_data(cur_tot+c).heka_type = 'cont';
sess_data(cur_tot+c).layer = '3';
sess_data(cur_tot+c).cell_type = 'unclear';
sess_data(cur_tot+c).region = 'LEC';
sess_data(cur_tot+c).class_certainty = 0;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2010-08-12_CWC_LFP_B\2010-8-12_14-55-7');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\2010-08-12_CWC_LFP_B');
sess_data(cur_tot+c).name = '2010-08-12_B';
sess_data(cur_tot+c).depth = 300;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [50; 2e3; 5e3; 2e3; 2e3; 2e3; 2e3; 3e3];
sess_data(cur_tot+c).heka_type = 'cont';
sess_data(cur_tot+c).layer = '3';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'LEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2010-08-17_CWC_LFP_A\2010-8-17_17-58-5');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\2010-08-17_CWC_LFP_A');
sess_data(cur_tot+c).name = '2010-08-17_A';
sess_data(cur_tot+c).depth = 140;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [50; 2e3; 5e3; 2e3; 2e3; 2e3; 2e3; 3e3];
sess_data(cur_tot+c).heka_type = 'cont';
sess_data(cur_tot+c).layer = '2';
sess_data(cur_tot+c).cell_type = 'fan';
sess_data(cur_tot+c).region = 'LEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2010-08-17_CWC_LFP_B\2010-8-17_19-7-58');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\2010-08-17_CWC_LFP_B');
sess_data(cur_tot+c).name = '2010-08-17_B';
sess_data(cur_tot+c).depth = 230;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [50; 2e3; 5e3; 2e3; 2e3; 2e3; 2e3; 3e3];
sess_data(cur_tot+c).heka_type = 'cont';
sess_data(cur_tot+c).layer = '3';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'LEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

%% ADJACENT CORTEX/ UNCLEAR
n_adj_unc = 25;
cur_tot = cur_tot + c -1;
c = 1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2005-01-28_HC_LFP\2005-1-28_18-24-8');
sess_data(cur_tot+c).heka_dir = '';
sess_data(cur_tot+c).name = '2005_01-28';
sess_data(cur_tot+c).depth = 350;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [100; 2e3; 2e3; 2e3; 2e3; 2e3; 2e3; 2e3];
sess_data(cur_tot+c).heka_type = 'nan';
sess_data(cur_tot+c).layer = '2/3';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'ectorhinal';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2005-11-25_CWC_LFP_B\2005-11-25_21-1-52');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2005-11-25_CWC_LFP_B');
sess_data(cur_tot+c).name = '2005_11-25_B';
sess_data(cur_tot+c).depth = 180;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [80; 1e3; 2e3; 1e3; 1e3; nan; 1e3; 1e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '2';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'isocortex/LEC';
sess_data(cur_tot+c).class_certainty = 0;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2006-03-24_CWC_LFP');
sess_data(cur_tot+c).heka_dir = '';
sess_data(cur_tot+c).name = '2006_03-24';
sess_data(cur_tot+c).depth = nan;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [nan; nan; nan; nan; nan; nan; nan; nan];
sess_data(cur_tot+c).heka_type = 'nan';
sess_data(cur_tot+c).layer = '2/3';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'perirhinal';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2006-03-26_CWC_LFP_A\2006-3-26_14-13-58');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2006-03-26_CWC_LFP_A');
sess_data(cur_tot+c).name = '2006_03-26_A';
sess_data(cur_tot+c).depth = nan;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [50; 1e3; 2e3; 2e3; 2e3; nan; 3e3; 3e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '2';
sess_data(cur_tot+c).cell_type = 'stellate';
sess_data(cur_tot+c).region = 'isocortex/MEC';
sess_data(cur_tot+c).class_certainty = 0;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2006-03-26_CWC_LFP_B\2006-3-26_21-13-56');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2006-03-26_CWC_LFP_B');
sess_data(cur_tot+c).name = '2006_03-26_B';
sess_data(cur_tot+c).depth = nan;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [50; 1e3; 2e3; 2e3; 2e3; nan; 3e3; 3e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '2';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).class_certainty = 0;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2006-03-29_CWC_LFP\2006-3-29_23-1-1');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2006-03-29_CWC_LFP');
sess_data(cur_tot+c).name = '2006_03-29';
sess_data(cur_tot+c).depth = nan;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [50; 1e3; 2e3; 2e3; 2e3; nan; 2e3; 2e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '2';
sess_data(cur_tot+c).cell_type = 'unclear';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).class_certainty = 0;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2006-03-31_CWC_LFP\2006-3-31_15-12-51');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2006-03-31_CWC_LFP');
sess_data(cur_tot+c).name = '2006_03-31';
sess_data(cur_tot+c).depth = nan;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [100; 1e3; 2e3; 2e3; 2e3; nan; 2e3; 2e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '3';
sess_data(cur_tot+c).cell_type = 'stellate';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).class_certainty = 0;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2006-04-05_CWC_LFP_B\2006-4-5_21-23-32');
sess_data(cur_tot+c).heka_dir = '';
sess_data(cur_tot+c).name = '2006_04-05_B';
sess_data(cur_tot+c).depth = nan;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [100; 1e3; 2e3; 2e3; 2e3; nan; 3e3; 3e3];
sess_data(cur_tot+c).heka_type = 'nan';
sess_data(cur_tot+c).layer = '2';
sess_data(cur_tot+c).cell_type = 'unclear';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).class_certainty = 0;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2006-04-06_CWC_LFP_A\2006-4-6_15-24-23');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2006-04-06_CWC_LFP_A');
sess_data(cur_tot+c).name = '2006_04-06_A';
sess_data(cur_tot+c).depth = nan;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [100; 2e3; 2e3; 2e3; 2e3; nan; 3e3; 3e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '3';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2006-04-06_CWC_LFP_B\2006-4-6_17-30-56');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2006-04-06_CWC_LFP_B');
sess_data(cur_tot+c).name = '2006_04-06_B';
sess_data(cur_tot+c).depth = nan;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [40; 2e3; 2e3; 2e3; 2e3; nan; 3e3; 3e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '3';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'isocortex/MEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2006-04-07_CWC_LFP_A\2006-4-7_14-39-42');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2006-04-07_CWC_LFP_A');
sess_data(cur_tot+c).name = '2006_04-07_A';
sess_data(cur_tot+c).depth = nan;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [100; 2e3; 2e3; 2e3; 2e3; nan; 3e3; 3e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '2';
sess_data(cur_tot+c).cell_type = 'interneuron';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2006-04-07_CWC_LFP_B\2006-4-7_15-24-5');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2006-04-07_CWC_LFP_B');
sess_data(cur_tot+c).name = '2006_04-07_B';
sess_data(cur_tot+c).depth = nan;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [80; 2e3; 2e3; 2e3; 2e3; nan; 3e3; 3e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '2';
sess_data(cur_tot+c).cell_type = 'stellate';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).class_certainty = 0;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2006-04-08_CWC_LFP_B\2006-4-8_21-11-53');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2006-04-08_CWC_LFP_B');
sess_data(cur_tot+c).name = '2006_04-08_B';
sess_data(cur_tot+c).depth = nan;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [80; 2e3; 2e3; 2e3; 2e3; nan; 3e3; 3e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '3';
sess_data(cur_tot+c).cell_type = 'interneuron';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2006-04-08_CWC_LFP_C\2006-4-8_22-16-53');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2006-04-08_CWC_LFP_C');
sess_data(cur_tot+c).name = '2006_04-08_C';
sess_data(cur_tot+c).depth = nan;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [50; 2e3; 2e3; 2e3; 2e3; nan; 3e3; 3e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '3';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'isocortex/MEC';
sess_data(cur_tot+c).class_certainty = 0;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2006-08-11_CWC_LFP\2006-8-11_22-30-42');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2006-08-11_CWC_LFP');
sess_data(cur_tot+c).name = '2006_08-11';
sess_data(cur_tot+c).depth = nan;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [100; 2e3; 2e3; 2e3; 2e3; nan; 3e3; 3e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '2/3';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'ectorhinal';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2006-08-30_CWC_LFP\2006-8-30_20-37-38');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2006-08-30_CWC_LFP');
sess_data(cur_tot+c).name = '2006_08-30';
sess_data(cur_tot+c).depth = nan;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [80; 2e3; 2e3; 2e3; 2e3; nan; 3e3; 3e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '2/3';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'isocortex';
sess_data(cur_tot+c).class_certainty = 0;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2006-09-01_CWC_LFP\2006-9-1_18-48-41');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2006-09-01_CWC_LFP');
sess_data(cur_tot+c).name = '2006_09-01';
sess_data(cur_tot+c).depth = nan;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [80; 2e3; 2e3; 2e3; 2e3; nan; 3e3; 3e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '2';
sess_data(cur_tot+c).cell_type = 'stellate';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).class_certainty = 0;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2006-09-04_CWC_LFP\2006-9-4_21-37-20');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2006-09-04_CWC_LFP');
sess_data(cur_tot+c).name = '2006_09-04';
sess_data(cur_tot+c).depth = nan;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [70; 2e3; 2e3; 2e3; 2e3; nan; 3e3; 3e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = 'nan';
sess_data(cur_tot+c).cell_type = 'nan';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).class_certainty = 0;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2006-09-05_CWC_LFP\2006-9-5_20-14-38');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2006-09-05_CWC_LFP');
sess_data(cur_tot+c).name = '2006_09-05';
sess_data(cur_tot+c).depth = nan;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [150; 2e3; 2e3; 2e3; 2e3; nan; 3e3; 3e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = 'nan';
sess_data(cur_tot+c).cell_type = 'nan';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).class_certainty = 0;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2007-05-08_CWC_LFP_B\2007-5-8_16-6-47');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007-05-08_CWC_LFP_B');
sess_data(cur_tot+c).name = '2007_05-08_B';
sess_data(cur_tot+c).depth = 100;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [100; 3e3; 3e3; nan; 4e3; 4e3; 4e3; 4e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '2';
sess_data(cur_tot+c).cell_type = 'unclear';
sess_data(cur_tot+c).region = 'isocortex/MEC';
sess_data(cur_tot+c).class_certainty = 0;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2007-05-25_CWC_LFP_A\2007-5-25_19-20-12');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007-05-25_CWC_LFP_A');
sess_data(cur_tot+c).name = '2007_05-25_A';
sess_data(cur_tot+c).depth = 350;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [100; 3e3; 3e3; nan; 4e3; 4e3; 4e3; 4e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '3';
sess_data(cur_tot+c).cell_type = 'unclear';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).class_certainty = 0;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2007-05-25_CWC_LFP_B\2007-5-25_20-24-13');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007-05-25_CWC_LFP_B');
sess_data(cur_tot+c).name = '2007_05-25_B';
sess_data(cur_tot+c).depth = 340;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [80; 3e3; 3e3; nan; 4e3; 4e3; 4e3; 4e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '3';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2007-05-30_CWC_LFP_B');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007-05-30_CWC_LFP_B');
sess_data(cur_tot+c).name = '2007_05-30_B';
sess_data(cur_tot+c).depth = 250;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [120; 3e3; 3e3; nan; 4e3; 4e3; 4e3; 4e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = '2/3';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'perirhinal';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2007-06-01_CWC_LFP_A\2007-6-1_15-0-4');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\A2007-06-01_CWC_LFP_A');
sess_data(cur_tot+c).name = '2007_06-01_A';
sess_data(cur_tot+c).depth = nan;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [150; 3e3; 3e3; nan; 5e3; 5e3; 5e3; 5e3];
sess_data(cur_tot+c).heka_type = 'sweep';
sess_data(cur_tot+c).layer = 'nan';
sess_data(cur_tot+c).cell_type = 'nan';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).class_certainty = 0;
c = c+1;


sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2007-09-25_CWC_LFP_B');
sess_data(cur_tot+c).heka_dir = '';
sess_data(cur_tot+c).name = '2007_09-25_B';
sess_data(cur_tot+c).depth = nan;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [100; 1e3; 2e3; nan; 4e3; 4e3; 4e3; 5e3];
sess_data(cur_tot+c).heka_type = 'nan';
sess_data(cur_tot+c).layer = 'nan';
sess_data(cur_tot+c).cell_type = 'nan';
sess_data(cur_tot+c).region = 'perirhinal';
sess_data(cur_tot+c).class_certainty = 0;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2009-04-05_CWC_LFP\2009-4-5_18-26-31');
sess_data(cur_tot+c).heka_dir = strcat(drive_letter,':\wc_data\MPascii\2009-04-05_CWC_LFP');
sess_data(cur_tot+c).name = '2009_04-05';
sess_data(cur_tot+c).depth = 520;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [100; 1e3; 2e3; 2e3; 2e3; 2e3; 2e3; 2e3];
sess_data(cur_tot+c).heka_type = 'cont';
sess_data(cur_tot+c).layer = '3/5';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'isocortex/MEC';
sess_data(cur_tot+c).class_certainty = 1;
c = c+1;

% ATROPINE EXPERIMENTS
cur_tot = cur_tot + c -1;
c = 1;
sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2007-11-22_CWC_LFP_A\2007-11-22_15-14-24');
sess_data(cur_tot+c).heka_dir = '';
sess_data(cur_tot+c).name = '2007_11-22_A';
sess_data(cur_tot+c).depth = nan;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [];
sess_data(cur_tot+c).heka_type = 'nan';
sess_data(cur_tot+c).layer = '3';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).atropine = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2007-11-22_CWC_LFP_B\2007-11-22_17-57-8');
sess_data(cur_tot+c).heka_dir = '';
sess_data(cur_tot+c).name = '2007_11-22_B';
sess_data(cur_tot+c).depth = nan;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [];
sess_data(cur_tot+c).heka_type = 'nan';
sess_data(cur_tot+c).layer = '3';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).atropine = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2007-11-26_CWC_LFP_B\2007-11-26_15-21-49');
sess_data(cur_tot+c).heka_dir = '';
sess_data(cur_tot+c).name = '2007_11-26_B';
sess_data(cur_tot+c).depth = nan;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [];
sess_data(cur_tot+c).heka_type = 'nan';
sess_data(cur_tot+c).layer = '2';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).atropine = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2007-11-27_CWC_LFP_A\2007-11-27_13-56-32');
sess_data(cur_tot+c).heka_dir = '';
sess_data(cur_tot+c).name = '2007_11-27_A';
sess_data(cur_tot+c).depth = nan;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [];
sess_data(cur_tot+c).heka_type = 'nan';
sess_data(cur_tot+c).layer = '2';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).atropine = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2007-11-27_CWC_LFP_C\2007-11-27_20-39-32');
sess_data(cur_tot+c).heka_dir = '';
sess_data(cur_tot+c).name = '2007_11-27_C';
sess_data(cur_tot+c).depth = nan;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [];
sess_data(cur_tot+c).heka_type = 'nan';
sess_data(cur_tot+c).layer = '3';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).atropine = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2007-11-27_CWC_LFP_D\2007-11-27_21-21-18');
sess_data(cur_tot+c).heka_dir = '';
sess_data(cur_tot+c).name = '2007_11-27_D';
sess_data(cur_tot+c).depth = nan;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [];
sess_data(cur_tot+c).heka_type = 'nan';
sess_data(cur_tot+c).layer = '3';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).atropine = 1;
c = c+1;

sess_data(cur_tot+c).directory = strcat(drive_letter,':\wc_data\2007-11-28_CWC_LFP_A\2007-11-28_13-25-26');
sess_data(cur_tot+c).heka_dir = '';
sess_data(cur_tot+c).name = '2007_11-28_A';
sess_data(cur_tot+c).depth = nan;
sess_data(cur_tot+c).ant_post = nan;
sess_data(cur_tot+c).lateral = nan;
sess_data(cur_tot+c).dors_vent = nan;
sess_data(cur_tot+c).dist_prh = nan;
sess_data(cur_tot+c).gains = [];
sess_data(cur_tot+c).heka_type = 'nan';
sess_data(cur_tot+c).layer = '3';
sess_data(cur_tot+c).cell_type = 'pyramidal';
sess_data(cur_tot+c).region = 'MEC';
sess_data(cur_tot+c).atropine = 1;
c = c+1;

%% L2/3 Pyrmidal Parietal  (N=11)
n_23_pyr_par = 11;

sess_data(cur_tot+1).directory = 'C:\wc_data\2004-11-04_CWC_LFP';
sess_data(cur_tot+2).directory = 'C:\wc_data\2005-12-05_CWC_LFP_C';
sess_data(cur_tot+3).directory = 'C:\wc_data\2005-12-06_CWC_LFP_A';
sess_data(cur_tot+4).directory = 'C:\wc_data\2005-12-06_CWC_LFP_B';
sess_data(cur_tot+5).directory = 'C:\wc_data\2005-12-13_CWC_LFP_A';
sess_data(cur_tot+6).directory = 'C:\wc_data\2005-12-13_CWC_LFP_D';
sess_data(cur_tot+7).directory = 'C:\wc_data\2005-12-14_CWC_LFP_C';
sess_data(cur_tot+8).directory = 'C:\wc_data\2005-12-14_CWC_LFP_D';
sess_data(cur_tot+9).directory = 'C:\wc_data\2007-04-16_CWC_LFP_B\2007-4-16_14-8-43';
sess_data(cur_tot+10).directory = 'C:\wc_data\2007-04-17_CWC_LFP_B\2007-4-17_14-49-45';
sess_data(cur_tot+11).directory = 'C:\wc_data\2007-04-18_CWC_LFP_B\2007-4-18_14-28-16';

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

sess_data(cur_tot+1).directory = 'C:\wc_data\2005-06-13_CWC_LFP\2005-6-13_17-30-12'; %L2/3
sess_data(cur_tot+2).directory = 'C:\wc_data\2005-12-13_CWC_LFP_B'; %L2/3
sess_data(cur_tot+3).directory = 'C:\wc_data\2005-12-14_CWC_LFP_B'; %maybe pyramidal
% sess_data(cur_tot+4).directory = 'C:\wc_data\Cortex_data\2007-04-05_CWC_LFP_B\2007-4-5_17-55-51'; %L5

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

sess_data(cur_tot+1).directory = 'C:\wc_data\2005-12-14_CWC_LFP_A';
sess_data(cur_tot+2).directory = 'C:\wc_data\2007-04-04_CWC_LFP_A\2007-4-4_15-6-52';%maybe L6
sess_data(cur_tot+3).directory = 'C:\wc_data\2007-04-10_CWC_LFP_A\2007-4-10_18-43-7';%border L2/3
sess_data(cur_tot+4).directory = 'C:\wc_data\2007-04-11_CWC_LFP_A\2007-4-11_16-43-19';
sess_data(cur_tot+5).directory = 'C:\wc_data\2007-04-12_CWC_LFP_B\2007-4-12_16-55-57';
sess_data(cur_tot+6).directory = 'C:\wc_data\2007-04-13_CWC_LFP_C\2007-4-13_16-2-15';
sess_data(cur_tot+7).directory = 'C:\wc_data\2007-04-16_CWC_LFP_A\2007-4-16_13-13-21';
sess_data(cur_tot+8).directory = 'C:\wc_data\2007-04-17_CWC_LFP_A\2007-4-17_13-36-49';%maybe L6
sess_data(cur_tot+9).directory = 'C:\wc_data\2007-04-18_CWC_LFP_A\2007-4-18_13-40-23';%maybe L6
sess_data(cur_tot+10).directory = 'C:\wc_data\2007-01-12_CWC_LFP\2007-1-12_16-56-12';%L6
% sess_data(cur_tot+11).directory =
% 'C:\wc_data\2005-06\Cortex_data\2007-01-11_CWC_LFP\2007-1-11_17-44-23';
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

sess_data(cur_tot+1).directory = 'C:\wc_data\2005-12-05_CWC_LFP_D';
sess_data(cur_tot+2).directory = 'C:\wc_data\2005-12-06_CWC_LFP_C';
sess_data(cur_tot+3).directory = 'C:\wc_data\2005-12-13_CWC_LFP_E';
sess_data(cur_tot+4).directory = 'C:\wc_data\2007-04-16_CWC_LFP_E\2007-4-16_18-34-42';

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

sess_data(cur_tot+1).directory = 'C:\wc_data\2007-04-11_CWC_LFP_B\2007-4-11_17-59-41';
sess_data(cur_tot+2).directory = 'C:\wc_data\2007-04-12_CWC_LFP_C\2007-4-12_18-6-33';
sess_data(cur_tot+3).directory = 'C:\wc_data\2007-04-13_CWC_LFP_D\2007-4-13_17-41-15';

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

sess_data(cur_tot+1).directory = 'C:\wc_data\2005-12-05_CWC_LFP_A';
sess_data(cur_tot+2).directory = 'C:\wc_data\2005-12-13_CWC_LFP_C'; %border L5
sess_data(cur_tot+3).directory = 'C:\wc_data\2007-04-05_CWC_LFP_C\2007-4-5_19-31-0';
sess_data(cur_tot+4).directory = 'C:\wc_data\2007-04-10_CWC_LFP_B\2007-4-10_19-35-5';%maybe L5
sess_data(cur_tot+5).directory = 'C:\wc_data\2007-04-17_CWC_LFP_C\2007-4-17_15-41-51';

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

sess_data(cur_tot+1).directory = 'C:\wc_data\2007-04-04_CWC_LFP_B\2007-4-4_16-0-35';
sess_data(cur_tot+2).directory = 'C:\wc_data\2007-04-13_CWC_LFP_B\2007-4-13_13-45-34';
sess_data(cur_tot+3).directory = 'C:\wc_data\2007-04-16_CWC_LFP_C\2007-4-16_15-18-23';
sess_data(cur_tot+4).directory = 'C:\wc_data\2007-07-27_CWC_LFP_B'; %slender !no LF8!
sess_data(cur_tot+5).directory = 'C:\wc_data\2007-08-14_CWC_LFP'; %thick
sess_data(cur_tot+6).directory = 'C:\wc_data\2007-08-15_CWC_LFP';  %!no LF8!
sess_data(cur_tot+7).directory = 'C:\wc_data\2007-08-17_CWC_LFP_A\2007-8-17_16-33-25'; %!no LF8!
% sess_data(cur_tot+8).directory = 'C:\wc_data\2007-08-17-CWC_LFP_B'; %!no LF8!
sess_data(cur_tot+8).directory = 'C:\wc_data\2007-08-21_CWC_LFP\2007-8-21_14-4-42'; %!no LF8!

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

sess_data(cur_tot+1).directory = 'C:\wc_data\2006-02-28_TH_LFP';
sess_data(cur_tot+2).directory = 'C:\wc_data\2006-03-10_TH_LFP';
sess_data(cur_tot+3).directory = 'C:\wc_data\2007-01-15_TH_LFP\2007-1-15_19-2-10';
sess_data(cur_tot+4).directory = 'C:\wc_data\2007-01-16_TH_LFP\2007-1-16_16-07-47';
sess_data(cur_tot+5).directory = 'C:\wc_data\2007-01-19_TH_LFP\2007-1-19_19-42-24';

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

sess_data(cur_tot+1).directory = 'C:\wc_data\2006-12-17_CWC_LFP_B\2006-12-17_18-47-4';
sess_data(cur_tot+2).directory = 'C:\wc_data\2006-12-19_CWC_LFP_A\2006-12-19_18-36-3';
sess_data(cur_tot+3).directory = 'C:\wc_data\2006-12-19_CWC_LFP_B\2006-12-19_20-41-15';
sess_data(cur_tot+4).directory = 'C:\wc_data\2006-12-20_CWC_LFP\2006-12-20_18-14-19';

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

sess_data(cur_tot+1).directory = 'C:\wc_data\2007-05-30_CWC_LFP_B';

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

sess_data(cur_tot+1).directory = 'C:\wc_data\2006-05-05_BG_LFP_A';
sess_data(cur_tot+2).directory = 'C:\wc_data\2006-05-05_BG_LFP_B';

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

%% CA1 Pyramidal Neurons (N=13)
n_ca1pyr = 13;
cur_tot = cur_tot+n_bg;

sess_data(cur_tot+1).directory = 'C:\wc_data\2004-11-04_HC_LFP\2004-11-4_17-35-8'; 
sess_data(cur_tot+2).directory = 'C:\wc_data\2004_8_13_HC_LFP\2004-8-13_18-45-0';
sess_data(cur_tot+3).directory = 'C:\wc_data\2004-10-28_HC_LFP\2004-10-28_15-33-41';
sess_data(cur_tot+4).directory = 'C:\wc_data\2004-12-21_HC_LFP\2004-12-21_16-42-12';
sess_data(cur_tot+5).directory = 'C:\wc_data\2005-11-18_HC_LFP_B';
sess_data(cur_tot+6).directory = 'C:\wc_data\2005-11-21_HC_LFP';
sess_data(cur_tot+7).directory = 'C:\wc_data\2005-11-22_HC_LFP';
sess_data(cur_tot+8).directory = 'C:\wc_data\2005-11-23\parta';
sess_data(cur_tot+9).directory = 'C:\wc_data\2005-11-23\partb';
sess_data(cur_tot+10).directory = 'C:\wc_data\2005-11-23\partc';
sess_data(cur_tot+11).directory = 'C:\wc_data\2005-11-24\parta';
sess_data(cur_tot+12).directory = 'C:\wc_data\2006-01-13_HC_LFP';
sess_data(cur_tot+13).directory = 'C:\wc_data\2006-01-22_HC_LFP';
sess_data(cur_tot+14).directory = 'C:\wc_data\2006-01-30_HC_LFP';

sess_data(cur_tot+1).name = '2004-11-04';
sess_data(cur_tot+2).name = '2004_8_13';
sess_data(cur_tot+3).name = '2004-10-28';
sess_data(cur_tot+4).name = '2004-12-21';
sess_data(cur_tot+5).name = '2005-11-18';
sess_data(cur_tot+6).name = '2005-11-21';
sess_data(cur_tot+7).name = '2005-11-22';
sess_data(cur_tot+8).name = '2005-11-23_A';
sess_data(cur_tot+9).name = '2005-11-23_B';
sess_data(cur_tot+10).name = '2005-11-23_C';
sess_data(cur_tot+11).name = '2005-11-24_A';
sess_data(cur_tot+12).name = '2006-01-13';
sess_data(cur_tot+13).name = '2006-01-30';

for i = 1:n_ca1pyr
    sess_data(cur_tot+i).cell_type = 'pyramidal';
    sess_data(cur_tot+i).layer = 'U';
    sess_data(cur_tot+i).region = 'CA1';
end

%% CA1 InterNeurons (N=8)
n_ca1int = 8;
cur_tot = cur_tot+n_ca1pyr;

sess_data(cur_tot+1).directory = 'C:\wc_data\2005-12-02_HC_LFP'; %SLM int
sess_data(cur_tot+2).directory = 'C:\wc_data\2004-11-05_HC_LFP\2004-11-5_18-16-49';%Slm int
sess_data(cur_tot+3).directory = 'C:\wc_data\2005-11-18_HC_LFP_A'; %?
sess_data(cur_tot+4).directory = 'C:\wc_data\2005-11-23\partd';% SLM int
sess_data(cur_tot+5).directory = 'C:\wc_data\2005-11-24\partb';% SLM int
sess_data(cur_tot+6).directory = 'C:\wc_data\2005-11-30_HC_LFP';%SLM int
sess_data(cur_tot+7).directory = 'C:\wc_data\2006-05-24_HC_LFP'; %SLM int
sess_data(cur_tot+8).directory = 'C:\wc_data\2006-08-09_HC_LFP\2006-8-9_20-8-24'; %SLM int

sess_data(cur_tot+1).name = '2005-12-02';
sess_data(cur_tot+2).name = '2004-11-05';
sess_data(cur_tot+3).name = '2005-11-18';
sess_data(cur_tot+4).name = '2005-11-23';
sess_data(cur_tot+5).name = '2005-11-24';
sess_data(cur_tot+6).name = '2005-11-30';
sess_data(cur_tot+7).name = '2006-05-24';
sess_data(cur_tot+8).name = '2006-08-09';

for i = 1:n_ca1int
    sess_data(cur_tot+i).cell_type = 'interneuron';
    sess_data(cur_tot+i).layer = 'U';
    sess_data(cur_tot+i).region = 'CA1';
end

%% CA3 Pyramidal Neurons (N=10)
n_ca3pyr = 10;
cur_tot = cur_tot+n_ca1int;

sess_data(cur_tot+1).directory = 'C:\wc_data\2006-02-02_HC_LFP\2006-2-2_17-26-57';
sess_data(cur_tot+2).directory = 'C:\wc_data\2006-02-10_HC_LFP';
sess_data(cur_tot+3).directory = 'C:\wc_data\2006-02-15_HC_LFP\2006-2-15_19-55-30';
sess_data(cur_tot+4).directory = 'C:\wc_data\2006-03-23_HC_LFP\2006-3-23_18-29-6';
sess_data(cur_tot+5).directory = 'C:\wc_data\2006-04-13_HC_LFP_A\2006-4-13_20-28-58';
sess_data(cur_tot+6).directory = 'C:\wc_data\2006-04-13_HC_LFP_B\2006-4-14_0-12-37';
sess_data(cur_tot+7).directory = 'C:\wc_data\2006-04-14_HC_LFP_A\2006-4-14_18-35-24';
sess_data(cur_tot+8).directory = 'C:\wc_data\2006-04-14_HC_LFP_C\2006-4-14_22-22-51';
sess_data(cur_tot+9).directory = 'C:\wc_data\2006-11-19_HC_LFP_B\2006-11-19_19-34-24';
sess_data(cur_tot+10).directory = 'C:\wc_data\2006-11-22_HC_LFP\2006-11-22_16-03-13';

sess_data(cur_tot+1).name = '2006-02-02';
sess_data(cur_tot+2).name = '2006-02-10';
sess_data(cur_tot+3).name = '2006-02-15';
sess_data(cur_tot+4).name = '2006-03-23';
sess_data(cur_tot+5).name = '2006-04-13_A';
sess_data(cur_tot+6).name = '2006-04-13_B';
sess_data(cur_tot+7).name = '2006-04-14_A';
sess_data(cur_tot+8).name = '2006-04-14_C';
sess_data(cur_tot+9).name = '2006-11-19_B';
sess_data(cur_tot+10).name = '2006-11-22';

for i = 1:n_ca3pyr
    sess_data(cur_tot+i).cell_type = 'pyramidal';
    sess_data(cur_tot+i).layer = 'U';
    sess_data(cur_tot+i).region = 'CA3';
end

%% Dentate Gyrus Granule Cells (N=10)
n_dg_gran = 10;
cur_tot = cur_tot+n_ca3pyr;

sess_data(cur_tot+1).directory = 'C:\wc_data\2006-03-22_HC_LFP';
sess_data(cur_tot+2).directory = 'C:\wc_data\2006-02-03_HC_LFP_A';
sess_data(cur_tot+3).directory = 'C:\wc_data\2006-04-11_HC_LFP\2006-4-11_19-7-30';
sess_data(cur_tot+4).directory = 'C:\wc_data\2006-04-12_HC_LFP_A\2006-4-12_21-7-39';
sess_data(cur_tot+5).directory = 'C:\wc_data\2006-05-12_HC_LFP';
sess_data(cur_tot+6).directory = 'C:\wc_data\2006-08-29_HC_LFP\2006-8-29_17-46-35';
sess_data(cur_tot+7).directory = 'C:\wc_data\2006-09-25_HC_LFP';
sess_data(cur_tot+8).directory = 'C:\wc_data\2006-09-26_HC_LFP';
sess_data(cur_tot+9).directory = 'C:\wc_data\2006-11-18_HC_LFP';
sess_data(cur_tot+10).directory = 'C:\wc_data\2006-11-21_HC_LFP_A\2006-11-21_17-24-7';

sess_data(cur_tot+1).name = '2006-03-22';
sess_data(cur_tot+2).name = '2006-02-03';
sess_data(cur_tot+3).name = '2006-04-11';
sess_data(cur_tot+4).name = '2006-04-12_A';
sess_data(cur_tot+5).name = '2006-05-12';
sess_data(cur_tot+6).name = '2006-08-29';
sess_data(cur_tot+7).name = '2006-09-25';
sess_data(cur_tot+8).name = '2006-09-26';
sess_data(cur_tot+9).name = '2006-11-18';
sess_data(cur_tot+10).name = '2006-11-21_A';

for i = 1:n_dg_gran
    sess_data(cur_tot+i).cell_type = 'granule';
    sess_data(cur_tot+i).layer = 'U';
    sess_data(cur_tot+i).region = 'DG';
end


%% NEW L3MEC
cd C:\WC_Germany\persistent_downs\
load ./new_mec_dir.mat
n_new_mec = 25;
cur_tot = cur_tot + n_dg_gran;

for i = 1:n_new_mec
sess_data(cur_tot+i).directory = new_mec_dir{i};
sess_data(cur_tot+i).region = 'MEC';
sess_data(cur_tot+i).layer = '3';
end

%% NEW SVEN
cur_tot = cur_tot + n_new_mec;

cd C:\WC_Germany\persistent_downs\
load ./new_pdown_dir
cur_uset = find(new_pdown_use == 1);
n_new_new_mec = length(cur_uset);
for ii = 1:n_new_new_mec
    sess_data(cur_tot+ii).directory = new_pdown_dir{cur_uset(ii)};
    sess_data(cur_tot+ii).region = 'MEC';
    sess_data(cur_tot+ii).layer = '3';
    
    sess_ep(cur_tot+ii) = new_pdown_ep(cur_uset(ii));
    sess_dp(cur_tot+ii) = new_pdown_dp(cur_uset(ii));
end

%%
cd(strcat(drive_letter,':\WC_Germany\overall_EC\'))
mec = find_struct_field_vals(sess_data,'region','MEC');
lec = find_struct_field_vals(sess_data,'region','LEC');
layer3 = find_struct_field_vals(sess_data,'layer','3');
layer2 = find_struct_field_vals(sess_data,'layer','2');
layer23 = find_struct_field_vals(sess_data,'layer','2/3');
superficial = sort(unique([layer2 layer3 layer23]));
pyramidal = find_struct_field_vals(sess_data,'cell_type','pyramidal');
multipolar = find_struct_field_vals(sess_data,'cell_type','multipolar');
fan = find_struct_field_vals(sess_data,'cell_type','fan');
stellate = find_struct_field_vals(sess_data,'cell_type','stellate');
DG = find_struct_field_vals(sess_data,'region','DG');

sure = find_struct_field_vals(sess_data,'class_certainty',1);

l3mec = intersect(mec,layer3);
l3lec = intersect(lec,layer3);
l3mec_s = intersect(l3mec,sure);
l3lec_s = intersect(l3lec,sure);
l2mec = intersect(mec,layer2);
l2lec = intersect(lec,layer2);
l2mec_s = intersect(l2mec,sure);
l2lec_s = intersect(l2lec,sure);
pyramidal_s = intersect(pyramidal,sure);
l3mec_p = intersect(pyramidal_s,l3mec);
l3lec_p = intersect(pyramidal_s,l3lec);
sup_lec = intersect(superficial,lec);
sup_mec = intersect(superficial,mec);

cd C:\WC_Germany\persistent_downs\
save overall_EC_dir sess_data l3lec_p l3mec_p l3mec* l3lec* l2* mec lec pyramidal ...
    pyramidal_s stellate fan multipolar superficial layer* sure sup* sess_ep sess_dp