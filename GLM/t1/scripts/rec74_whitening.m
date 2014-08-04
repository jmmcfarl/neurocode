
clear all; flen =6; 

n_reps = 24;

%% natural stimli from recording 74
cd ~/Data/blanche/rec_74/matlabdata/   
load dstimpsRec74.mat; 
cd ~/Data/blanche/rec_74
stf74=load('stimfiles74.mat'); 

nconds   = length(dstimps74)
cnames   = [stf74.stimfiles]; cellfun(@disp,cnames)

keys = {'ar','wr','pr','nr'}
rids = cellfun(@(tkey)...
	find(strcmp(cellfun(@(x)x(1:2),cnames,'UniformOutput',0),tkey)),...
	keys,'UniformOutput',0); 

tns        = rids{4}(3:4); cellfun(@disp,cnames(tns))
tpn        = rids{3}(3:4); cellfun(@disp,cnames(tpn))

tconds = [tpn,tns]; ctype='all'; 
% tconds = [rids{3}]; ctype='pink'; 

selstim   = []; 
% cur_means = [];
% cur_vars = [];
for icond=1:length(tconds); 
	tcond = tconds(icond)
    cur_stim = dstimps74{tcond};
%     cur_stim = (cur_stim - mean(cur_stim(:)))/std(cur_stim(:)); %z-score normalization
	selstim=[selstim;repmat(cur_stim,n_reps,1)];
%     cur_means = [cur_means; mean(cur_stim)];
%     cur_vars = [cur_vars; var(cur_stim)];
end; 

NT   = size(selstim,1); SDIM  = size(selstim,2); NeK   = flen*SDIM;
X    = zeros(NT-flen+1,NeK);
for i = flen:NT; X(i-flen+1,1:NeK) = reshape(selstim((i-flen+1):i,:),1,NeK); end
cd ~/Data/blanche/rec_74/matlabdata/   
tfilename = sprintf('StimXmat_z-%s.mat',ctype)
save(tfilename,'-v7.3','X');

disp(' -- whitening'); drawnow;
[coefs,scores] = princomp(X); scorevars = var(scores);

tfilename = sprintf('PCScores_z-%s.mat',ctype)
save(tfilename,'-v7.3','coefs','scores','scorevars','flen','SDIM');


%%
% the spikes
cd ~/Data/blanche/rec_74/matlabdata/   
load spksegsRec74.mat; load ~/Data/blanche/rec_75/matlabdata/stdparsRec75; 
allspks  = spksegs74';  

% %merge cell 10 and 14 (really same cell) to create cell 27
% cell10_spks = allspks(:,10);
% cell14_spks = allspks(:,14);
% temp_comb = cell(size(allspks,1),1);
% for i = 1:size(allspks,1)
%    temp_comb{i} = unique([cell10_spks{i};cell14_spks{i}]); 
% end
% allspks = [allspks temp_comb];

ncells = size(allspks,2)
cellfun(@(x)size(x,1),dstimps74,'UniformOutput',0) 
stimlen = dt*6000
nspks = cellfun(@length,allspks); 

boxplot(nspks','Labels',cnames,'labelorientation','inline')

cellfun(@disp,cnames(tconds))

aselspks = cell(ncells,1); 
for icell = 1:ncells; 
	selspks = []; 
	for icond=1:length(tconds); 
		tcond = tconds(icond)
		selspks=[selspks;((icond-1)*stimlen)+allspks{tcond,icell}];
	end;
	aselspks{icell} = selspks; 
end
% ecdf(aselspks{1}); hold on; for icell =2:10; ecdf(aselspks{icell}); end; hold off; 
cd ~/Data/blanche/rec_74/matlabdata/   
tfilename = sprintf('spks74-%s.mat',ctype); save(tfilename,'-v7.3','aselspks');
