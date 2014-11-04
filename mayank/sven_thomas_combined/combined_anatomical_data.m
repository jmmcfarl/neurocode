
clear anatomy
%% L3MEC Old
%depth  postrhinal   rostro-caudal   dorsoventral   lateralT  lateralS
%mediolateralHpc

anatomy = [...
    nan	nan	-5	-3.1	3.6	nan	nan
nan	nan	-4.9	-3.5	3.8	nan	nan
nan	nan	-5	-3	3.5	nan	nan
nan	nan	-5	-2.8	3.5	nan	nan
390	500	-5.1	-2.5	3.7	3.85	1
340	nan	-4.9	-2	3.6	3.85	1
250	200	-5	-2.3	3.8	4	0
210	600	-5.1	-2.5	3.5	3.55	3
240	600	-5	-2.5	3.2	3.4	4
310	500	-5	-2.8	3.3	3.55	3
280	550	-5.1	-2.5	3.1	3.25	5
360	550	-5	-2.4	3.9	3.85	1
340	450	-5	-2.2	3.8	4	0
370	350	-5.1	-2.5	3.3	3.55	3
420	50	-5.1	-2.4	3.4	3.7	2
310	250	-5	-2.4	3.7	3.85	1
350	330	-5	-2.2	3.3	3.4	4
430	450	-5	-2.7	3.5	3.85	1
430	600	-5.1	-2.7	3.7	4	0
200	50	-5.1	-2.3	3.8	4	0
250	75	-4.9	-2.7	4	4.15	-1
350	50	-5	-2.2	3.6	3.85	1
];

%% L3LEC
anatomy = [anatomy;
275	nan	-4.4	-3.9	4	nan	nan
280	nan	-3.8	-3.9	4.2	nan	nan
240	nan	-4.1	-3.7	4.2	nan	nan
260	nan	-4.5	-3.5	4	nan	nan
230	nan	-3.6	-4	4.1	nan	nan
260	nan	-3.8	-4.8	3.8	nan	nan
300	nan	-4.2	-4.4	3.9	nan	nan
440	nan	-4.3	-4.5	3.8	nan	nan
190	nan	-3.5	-4.1	4	nan	nan
220	nan	-3.7	-4.2	4.1	nan	nan
410	nan	-3.9	-4.5	3.7	nan	nan
300	nan	-4	-4.6	3.7	nan	nan
230	nan	-3.9	-4.6	3.7	nan	nan
];

%% L3MEC new
anatomy = [anatomy;
    300	300	-5	-2.3	3.6	3.55	3
300	300	-5	-2.5	3.8	3.85	1
275	175	-5	-2.8	3.8	3.7	2
210	525	-5	-3	3.6	3.7	2
250	350	-5	-2.6	3.6	3.7	2
225	350	-5	-3.2	3.7	3.85	1
350	350	-5	-2.8	3.7	4	0
175	425	-5	-2.9	3.8	3.85	1
450	125	-5.1	-3	3.7	4	0
425	300	-5	-2.5	3.6	3.85	1
350	650	-5	-3.3	3.3	3.55	3
400	850	-4.9	-3.5	3.4	3.55	3
500	550	-5	-3	3.5	3.7	2
400	100	-5	-2.5	3.8	3.85	1
350	500	-5.1	-3.4	3.7	3.85	1
600	600	-5	-3.1	3.8	3.85	1
500	350	-5.1	-3.4	3.9	4	0
325	125	-4.9	-2.4	3.9	4	0
450	400	-5.1	-3.6	3.9	4	0
];

%% L3MEC NP
anatomy = [anatomy;
    250	nan	-5.1	-3	3.4	nan nan
210	nan	-5	-2.5	3.2	nan nan
210	nan	-5.2	-2.3	3.2	nan nan
250	nan	-5	-2.2	3.4	nan nan
250	nan	-5	-2.7	4	4.15	-1
230	nan	-5.1	-3.5	3.8	4	0
350	nan	-5.1	-3	3.9	4	0
140	200	-4.9	-3	3.7	nan	nan
];

%% L3LEC NP
anatomy = [anatomy;
260	nan	-4.3	-4.1	3.9	nan	nan
340	nan	-4	-4.2	3.9	nan	nan
490	nan	-3.9	-4.8	3.2	nan	nan
300	nan	-4.2	-4.1	4	nan	nan
190	nan	-4	-4.4	4	nan	nan
];

%%
load ./combined_dir_nd.mat
uset = sort([l3mec l3lec l3mec_np l3lec_np]);
all_cells = 1:length(combined_dir);
l3mec = find(ismember(all_cells(uset),l3mec));
l3lec = find(ismember(all_cells(uset),l3lec));
l3mec_np = find(ismember(all_cells(uset),l3mec_np));
l3lec_np = find(ismember(all_cells(uset),l3lec_np));

load ./combined_core_analysis_fin_nd_np

%%
l3mec_np(l3mec_np==62) =[];
int = 62;
offset = 0.5;
scale = 1000;
 scatter3(anatomy(l3mec,3),anatomy(l3mec,4),anatomy(l3mec,5),scale*fract_rt2_ups(l3mec)+offset,'r','linewidth',2);
 hold on
 scatter3(anatomy(l3lec,3),anatomy(l3lec,4),anatomy(l3lec,5),scale*fract_rt2_ups(l3lec)+offset,'b','linewidth',2);
 scatter3(anatomy(l3mec_np,3),anatomy(l3mec_np,4),anatomy(l3mec_np,5),scale*fract_rt2_ups(l3mec_np)+offset,'k','linewidth',2);
  scatter3(anatomy(l3lec_np,3),anatomy(l3lec_np,4),anatomy(l3lec_np,5),scale*fract_rt2_ups(l3lec_np)+offset,'g','linewidth',2);
  scatter3(anatomy(int,3),anatomy(int,4),anatomy(int,5),scale*fract_rt2_ups(int)+offset,'c','linewidth',2);
  xlabel('Rostrocaudal')
  ylabel('Dorsoventral')
  zlabel('Mediolateral')
  