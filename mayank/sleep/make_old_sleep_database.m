clear all
addpath('~/James_scripts/NeuralynxMatlabImportExport_v501/')
addpath(genpath('~/James_scripts_old/Nlx2Mat_relDec11/'))
addpath('~/Analysis/Mayank/sleep/');
addpath('~/Analysis/Mayank/sven_thomas_combined/')

%MEC
% MP not really that stable. Maybe some decent cortical UDS but cant see clear persistence. Probably not analyzable
data(1).dir = '/Users/james/Data/Mayank/sleep/2012-7-5_Patch/2012-7-5_17-11-29';
data(1).ipsiLFPs = [2:8];
data(1).contraLFPs = [];
data(1).MPfile = 'CSC1.ncs';
data(1).ipsi_L = 0;
data(1).contra_L = 0;
data(1).heka_dir = '/Users/james/Data/Mayank/sleep/2012-7-5_Patch/2012_07_05_1.mat';
data(1).celltype = 'MECL3';
data(1).MPloc = 'MEC';
data(1).good_bounds = [0 192];

%MEC
% MP strange. LFP doesn't really show clear UDS either. Not usable
data(2).dir = '/Users/james/Data/Mayank/sleep/2012-7-8_2/2012-7-8_13-9-45';
data(2).ipsiLFPs = [2:8];
data(2).contraLFPs = [];
data(2).MPfile = 'CSC1.ncs';
data(2).ipsi_L = 0;
data(2).contra_L = 0;
data(2).heka_dir = '';
data(2).celltype = 'MECL3';
data(2).MPloc = 'MEC';
data(2).good_bounds = [0 Inf];


%MEC
%Wierd MP with 'incomplete' looking UP states. Decent cortical UDS. Maybe
%examples of pers downs? But they aren't really clear because of the wierd
%MP UDS
data(3).dir = '/Users/james/Data/Mayank/sleep/2012-7-10_Sleep_WC/2012-7-10_13-49-40';
data(3).ipsiLFPs = [2:8];
data(3).contraLFPs = [];
data(3).MPfile = 'CSC1.ncs';
data(3).ipsi_L = 0;
data(3).contra_L = 0;
data(3).heka_dir = '/Users/james/Data/Mayank/sleep/2012-7-10_Sleep_WC/2012_07_10_sleep_1.mat';
data(3).celltype = 'MECL3';
data(3).MPloc = 'MEC';
data(3).good_bounds = [0 490];


%MEC
%some pretty good cortical UDS. MP is pretty good too. Good examples of
%pers ups, and maybe a few pers downs.
data(4).dir = '/Users/james/Data/Mayank/sleep/2012-7-13#1/2012-7-13_13-31-49';
data(4).ipsiLFPs = [2:8];
data(4).contraLFPs = [];
data(4).MPfile = 'CSC1.ncs';
data(4).ipsi_L = 0;
data(4).contra_L = 0;
data(4).heka_dir = '/Users/james/Data/Mayank/sleep/2012-7-13#1/2012_07_13_sleep_3.mat';
data(4).celltype = 'MECL3';
data(4).MPloc = 'MEC';
data(4).good_bounds = [0 515];


%MEC
%not much clear UDS (a few good brief epochs). MP is pretty wierd though.
%Probably not usable. Not really clear examples of pers.
data(5).dir = '/Users/james/Data/Mayank/sleep/2012-7-13#2/2012-7-13_16-54-15';
data(5).ipsiLFPs = [2:8];
data(5).contraLFPs = [];
data(5).MPfile = 'CSC1.ncs';
data(5).ipsi_L = 0;
data(5).contra_L = 0;
data(5).heka_dir = '/Users/james/Data/Mayank/sleep/2012-7-13#2/2012_07_13_sleep_6.mat';
data(5).celltype = 'MECL3';
data(5).MPloc = 'MEC';
data(5).good_bounds = [0 Inf];


cd ~/Analysis/Mayank/sleep/
save sleep_dirs_old data
