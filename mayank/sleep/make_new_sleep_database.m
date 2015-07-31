clear all
addpath('~/James_scripts/NeuralynxMatlabImportExport_v501/')
addpath(genpath('~/James_scripts_old/Nlx2Mat_relDec11/'))
addpath('~/Analysis/Mayank/sleep/');

dd = 1;
%cortical
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-09-16_12-26-27';
data(dd).ipsiLFPs = []; %ispilateral LFP channels
data(dd).contraLFPs = [2 6 10 14 18 22 26 30]; %contralateral LFP channels
data(dd).MPfile = 'CSC33.ncs'; %name of MP ncs file
data(dd).ipsi_L = 0; %indicator if ispi files are named 'L'
data(dd).contra_L = 0; %indicator if conta files are named 'L'
data(dd).MPloc = 'Ctx'; %location of cell
data(dd).good_bounds = [0 Inf]; %bounds on which MP data is usable

dd = dd + 1;
%d2 cortical
%some decent epochs of cortical UDS. The MP is a bit wierd though and
%doesnt have too clear of UDS. Not sure if usable or not...
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-09-17_12-45-54';
data(dd).ipsiLFPs = [2 6 10 14 18 22 26 30];
data(dd).contraLFPs = [];
data(dd).MPfile = 'CSC33.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 0;
data(dd).MPloc = 'Ctx';
data(dd).good_bounds = [0 290];

dd = dd + 1;
%d3 cortical
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-09-30_14-18-08';
data(dd).ipsiLFPs = [];
data(dd).contraLFPs = [1:16];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 0;
data(dd).MPloc = 'Ctx';
data(dd).good_bounds = [0 124];

dd = dd + 1;
%d4 cortical
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-10-01_13-08-37';
data(dd).ipsiLFPs = [];
data(dd).contraLFPs = [1:16];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 0;
data(dd).MPloc = 'Ctx';
data(dd).good_bounds = [0 546];

dd = dd + 1;
%d5 cortical
%some pretty good UDS towards the beginning. Probably usable but not great
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-10-07_12-30-20';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [1:16];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 1;
data(dd).contra_L = 0;
data(dd).MPloc = 'Ctx';
data(dd).good_bounds = [0 800];

dd = dd + 1;
%d6 cortical
%great UDS, especially at the beginning. Both MP and LFP are great.
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-10-09_12-29-04';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [1:16];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'Ctx';
data(dd).good_bounds = [0 Inf];

dd = dd + 1;
%d7 MEC [probable L3MEC. can see large prolonged UPs at time with robust
%spk rates. also lower freq uds].
%not great UDS in LFP or MP. Some brief epochs with decent UDS. Might be
%able to see some examples of pers ups. Not pers downs.
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-10-14_14-04-11';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 0;
data(dd).MPloc = 'MEC';
data(dd).good_bounds = [0 609];

dd = dd + 1;
%d8 MEC [cell type somewhat unclear. I'd guess L2MEC based on small peak in
%theta band, and somewhat weaker US and low spk rates]
%not really clear UDS. I would say not usable, but there might be some
%instances that could be argued for pers states.
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-10-22/2014-10-22_11-57-18';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 0;
data(dd).MPloc = 'MEC';
data(dd).good_bounds = [0 417];

dd = dd + 1;
%d9 MEC [clear L3MEC. high rate, large US, low-freq]
%great MEC example. Long rec. Good cortical UDS towards the end. Very clear
%MEC MP UDS, with nice pups and pdowns.
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-10-31_13-43-29';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 0;
data(dd).MPloc = 'MEC';
data(dd).good_bounds = [0 Inf];

dd = dd + 1;
%d10 MEC. [probably L3MEC. high rate, large US, and low-freq]
%Short, not much cortical UDS. Probably not useful
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-11-04_12-27-44';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [13];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'MEC';
data(dd).good_bounds = [0 95];

dd = dd + 1;
%d11 MEC. [very likely L2 stellate. weaker 'humped' us, lower rate, and
%clear theta peak]
%Long enough, with some decent cortical UDS. Wierd MP, dominated by
%DS, with choppy ups. Maybe some pdowns.
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-11-04_15-31-37';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [13];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'MEC';
data(dd).good_bounds = [0 468];

dd = dd + 1;
%d12 MEC. [probably L2. Not that clear though. slight theta peak, and
%humped US]
%Not much clear UDS in either cortical LFP or MP
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-11-06_13-41-25';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [13];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'MEC';
data(dd).good_bounds = [0 400];

dd = dd + 1;
%d13 MEC. [not clear, but L3 if had to guess]
%Nice recording with clear UDS, and nice MP. But very short, so only a
%few nice examples.
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-11-06_15-20-38';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [13];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'MEC';
data(dd).good_bounds = [0 75];

dd = dd + 1;
%d14 MEC. [clear L3. large US, high rate, low-freq]
%Pretty good rec. Decent cortical UDS and MP. Some good pups and
%pdowns.
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-11-08_14-36-03';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [13];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'MEC';
data(dd).good_bounds = [0 510];

dd = dd + 1;
%d15 MEC. [pretty clear L2. Weak US, low rate, maybe a hint of theta]
%Decent rec. Some decent cortical UDS. MP is pretty DS heavy, with
%choppy UPs. Probably some pers DOWNS. Maybe some pers UPS. Kinda hard to
%tell though.
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-11-12_10-57-15';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [13];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'MEC';
data(dd).good_bounds = [0 242];

dd = dd + 1;
%d16 MEC. [unclear]. Sven says probably cortical neuron
%Cortical UDS is pretty good at times. The MP is not really stable,
%and doesn't have clear UDS. So, very likely not usable.
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-11-13_12-48-52';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [13];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'MEC';
data(dd).good_bounds = [0 310];

dd = dd + 1;
%d17 MEC. [either bad MP rec or maybe an L2 cell, not clear though]. Sven
%says probably cortical neuron
%Great cortical UDS. MP not really stable, or just not showing clear
%UDS. So, very likely not usable.
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-11-13_13-25-58';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [13];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'MEC';
data(dd).good_bounds = [0 540];

dd = dd + 1;
%d18 MEC. [probably L2, but not that clear. weak US, low rate, not much of a clear theta peak though]
%Nice cortical UDS. MP does not have very clear UDS, but occasionally
%its decently clear. Probably a few pdowns and pups here, but wont be the
%best examples. Possibly usable.
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-11-15_12-48-19';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [13];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'MEC';
data(dd).good_bounds = [0 780];

dd = dd + 1;
%d19 MEC. [likely an L2 because of weak US, somewhat sparser rate, and a
%hint of theta]
%Very nice long rec. Good cortical UDS. MEC MP pretty good, with some
%nice examples of Pup and Pdown at times. The UDS are a bit strange
%(changing very often), but should be usable. MP also appears to be
%degrading slowly towards the end.
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-11-15_13-08-31';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [13];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'MEC';
data(dd).good_bounds = [0 1094];

dd = dd + 1;
%d20 MEC. [almost certainly not L3. weak US, low rate, theta peak]
%Probably too short to use for much. Decent cortical UDS. MEC MP not
%too stable. Also has DS-heavy with choppy UPS. Probably not usable.
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-11-20_12-54-16';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [13];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'MEC';
data(dd).good_bounds = [0 226];

dd = dd + 1;
%d21 MEC. [obvious L3. huge US, high rate, low-freq]
%Nice recording. A bit on the short side. Cortical UDS is OK. MEC MP
%is fantastic! Very clear bimodal UDS. Nice strong PUps. Maybe some PDowns
%as well
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-11-20_13-09-58';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [13];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'MEC';
data(dd).good_bounds = [0 284];

dd = dd + 1;
%d22. [not L3. Weak US, low rate. not an obvious theta peak, so maybe
%cortical]
%Maybe MEC or cortical. Short choppy US. Some decent epochs of
%cortical UDS. Maybe some pers DS, but not the clearest examples.
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-11-29_13-12-01';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [13];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'unsure'; %prob cortical (POR or visual)
data(dd).good_bounds = [0 882];

dd = dd + 1;
%d23 MP has high rate and does not show clear UDS. Cortical LFP shows decent
%UDS. MIght be interpreted as lots of pers UPs, but that's not really what
%it looks like
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-12-09_15-04-21';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [1:16];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 1;
data(dd).contra_L = 0;
data(dd).MPloc = 'Ctx'; 
data(dd).good_bounds = [0 779];

dd = dd + 1;
%d24. Some epochs with decent UDS in MP and LFP. Seem pretty synched when
%both have clear UDS. Should have some usable epochs, though not the best
%examples.
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-12-09_15-29-56';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [1:16];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'Ctx'; 
data(dd).good_bounds = [0 579];

dd = dd + 1;
%d25. Too short of a rec to really be usable. 
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-12-10_11-43-50';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [1:16];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'Ctx'; 
data(dd).good_bounds = [0 144];

dd = dd + 1;
%d26. MP not good, with choppy unstable UDS. LFP also doesn't show
%clear/sustained UDS. Probably not usable.
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-12-11_11-53-25';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [1:16];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'Ctx'; 
data(dd).good_bounds = [0 532];

dd = dd + 1;
%d27. MP and LFP both have some OK UDS. Likely some usable epochs, but not
%very clear, so not good examples.
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-12-11_12-33-16';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [1:16];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'Ctx'; 
data(dd).good_bounds = [0 359];

dd = dd + 1;
%d28. UDS not very clear.
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-11-04_13-08-16';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [1:16];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'MEC'; 
data(dd).good_bounds = [0 155];

dd = dd + 1;
%d29. Pretty good UDS in LFP, and MP. Looks like some pers UPs. Not many
%(if any) clear pers downs
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-11-06_14-08-44';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [1:16];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'MEC'; 
data(dd).good_bounds = [0 235];

cd ~/Analysis/Mayank/sleep/
save sleep_dirs data
