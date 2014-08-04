clear all; flen =6; tcell  = 26; 

%% natural stimli from recording 75
cd ~/Data/blanche/rec_75/matlabdata/   
load dstimpsRec75.mat; 
stf75=load('stimfiles75.mat'); 

allstims = dstimps75; 
nconds   = length(allstims)
cnames   = [stf75.stimfiles]; cellfun(@disp,cnames)

keys = {'af','wn','pn','pt','ps','ns'}
rids = cellfun(@(tkey)...
	find(strcmp(cellfun(@(x)x(1:2),cnames,'UniformOutput',0),tkey)),...
	keys,'UniformOutput',0); 

natstims76 = rids{6}(5:end); 
tns        = rids{6}(3:4); cellfun(@disp,cnames(tns))
tpn        = rids{3}(3:4); cellfun(@disp,cnames(tpn))
tps        = rids{5}(3:4); cellfun(@disp,cnames(tps))

% tconds = [tps tpn]; ctype ='pn' %% + pn 35/45
tconds = [tpn]; ctype ='pn' %% + pn 35/45
tconds

selstim   = []; 
% cur_means = [];
% cur_vars = [];
for icond=1:length(tconds); 
	tcond = tconds(icond)
    cur_stim = allstims{tcond};
    cur_stim = zscore(cur_stim); %z-score normalization
	selstim=[selstim;cur_stim];
%     cur_means = [cur_means; mean(cur_stim)];
%     cur_vars = [cur_vars; var(cur_stim)];
end; 

NT   = size(selstim,1); SDIM  = size(selstim,2); NeK   = flen*SDIM;
X    = zeros(NT-flen+1,NeK);
for i = flen:NT; X(i-flen+1,1:NeK) = reshape(selstim((i-flen+1):i,:),1,NeK); end
tfilename = sprintf('StimXmat_z-%s.mat',ctype)
save(tfilename,'-v7.3','X');

disp(' -- whitening'); drawnow;
[coefs,scores] = princomp(X); scorevars = var(scores);

tfilename = sprintf('PCScores_z-%s.mat',ctype)
save(tfilename,'-v7.3','coefs','scores','scorevars','flen','SDIM');



% % the spikes
% load spksegsRec75.mat; load spksegsRec76.mat; load stdparsRec75; 
% allspks  = [spksegs75,spksegs76]';  ncells = size(allspks,2)
% cellfun(@(x)size(x,1),allstims,'UniformOutput',0) 
% stimlen = dt*6000
% nspks = cellfun(@length,allspks); 
% 
% boxplot(nspks','Labels',cnames,'labelorientation','inline')
% 
% 
% tconds = [tns,natstims76]; ctype ='ns'  %% natstims76 + ns 35/45
% cellfun(@disp,cnames(tconds))
% 
% aselspks = cell(ncells,1); 
% for icell = 1:ncells; 
% 	selspks = []; 
% 	for icond=1:length(tconds); 
% 		tcond = tconds(icond)
% 		selspks=[selspks;((icond-1)*stimlen)+allspks{tcond,icell}];
% 	end;
% 	aselspks{icell} = selspks; 
% end
% ecdf(aselspks{1}); hold on; for icell =2:10; ecdf(aselspks{icell}); end; hold off; 
% 
% tfilename = sprintf('spks7576-%s.mat',ctype); save(tfilename,'-v7.3','aselspks');
% 
% 
% 
% tconds = [tpn,natstims76]; ctype ='pn' %% + pn 35/45
% aselspks = cell(ncells,1); 
% for icell = 1:ncells; 
% 	selspks = []; 
% 	for icond=1:length(tconds); 
% 		tcond = tconds(icond)
% 		selspks=[selspks;((icond-1)*stimlen)+allspks{tcond,icell}];end;
% 	aselspks{icell} = selspks; 
% end
% ecdf(aselspks{1}); hold on; for icell =2:10; ecdf(aselspks{icell}); end; hold off; 
% 
% tfilename = sprintf('spks7576-%s.mat',ctype); save(tfilename,'-v7.3','aselspks');
% 
% 
% 
% tconds = [tps,natstims76]; ctype ='ps' %% + ps 35/45
% aselspks = cell(ncells,1); 
% for icell = 1:ncells; 
% 	selspks = []; 
% 	for icond=1:length(tconds); 
% 		tcond = tconds(icond) 
% 		selspks=[selspks;((icond-1)*stimlen)+allspks{tcond,icell}];end;
% 	aselspks{icell} = selspks; 
% end
% ecdf(aselspks{1}); hold on; for icell =2:10; ecdf(aselspks{icell}); end; hold off; 
% 
% tfilename = sprintf('spks7576-%s.mat',ctype); save(tfilename,'-v7.3','aselspks');
% 

% %% alltogether now
% load spksegsRec75.mat; load spksegsRec76.mat; load stdparsRec75; 
% allspks  = [spksegs75,spksegs76]';  ncells = size(allspks,2)
% cellfun(@(x)size(x,1),allstims,'UniformOutput',0) 
% stimlen = dt*6000; nspks = cellfun(@length,allspks); 
% 
% tconds = [tpn,tps,tns,natstims76]; ctype='all'; 
% aselspks = cell(ncells,1); 
% for icell = 1:ncells; 
% 	selspks = []; 
% 	for icond=1:length(tconds); 
% 		tcond = tconds(icond)
% 		selspks=[selspks;((icond-1)*stimlen)+allspks{tcond,icell}];
% 	end;
% 	aselspks{icell} = selspks; 
% end
% ecdf(aselspks{1}); hold on; for icell =2:10; ecdf(aselspks{icell}); end; hold off; 
% 
% tfilename = sprintf('spks7576-%s.mat',ctype); save(tfilename,'-v7.3','aselspks');
