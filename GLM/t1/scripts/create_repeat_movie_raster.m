clear all
close all
datdir = '~/Data/blanche/rec_74';
cd(datdir)
cd matlabdata
load ./spksegsRec74
load ./spksRec74
%%
cd ~/James_scripts/GLM/t1/rec74_modfits/
dt = 0.05;
cell_range = 1:25
win_width = 0.25;
for ii = 1:16
cur_seq = ii;

spmat = spksegs74(cell_range,cur_seq);
range = [0 5];
% rep_bound_times = rep_bound_times(1:24,:);
mk_size = 1;
ncells = length(spmat);
mint = range(1);
maxt = range(2);

spacing = 1;
y = 0;

f1 = figure;set(f1,'Position',[500 1000 1000 1000]); 
set(f1,'PaperUnits','centimeters');
set(f1, 'PaperSize', [20 25]);
set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))])
sps = 1:9;
% subplot(3,3,sps(mod(sps,3)~=1))
hold on
if nargin < 5
    colors = zeros(ncells,3);
    colors([2:2:ncells],1) = 1;
    %     colors([2:2:ncells],1) = 1;
    %     colors([2:2:ncells],2) = 0.5;
    colors([1:2:ncells],3) = 1;
end

nreps = 24;
cur_rep_bound_times = rep_bound_times(((cur_seq-1)*24+1):cur_seq*24,:);
cur_rep_bound_times = cur_rep_bound_times - seq_bound_times(cur_seq,1);
for cc = 1:ncells
    for n = 1:nreps
        cur_spks = spmat{cc}(spmat{cc} > cur_rep_bound_times(n,1) & spmat{cc} < cur_rep_bound_times(n,2));
        cur_spks = cur_spks - cur_rep_bound_times(n,1);
        cur_spks = cur_spks(cur_spks > range(1) & cur_spks < range(2));
        plot(cur_spks,ones(size(cur_spks))*y,'*','color',colors(cc,:),'markersize',mk_size);
        y = y + spacing;
    end
    y = y + spacing;
    line([mint maxt],[y y],'color','k')
    y = y + spacing;
end

axis([mint maxt 0 (y+1)])
set(gca,'LineWidth',0.6,'FontSize',14,'YTickLabel',[])
xlabel('Time (s)','fontsize',16)
ylabel('Cell number','fontsize',16)
title(stimfiles{cur_seq})
print(stimfiles{cur_seq},'-dpng');
close all
end


%%
subplot(3,3,4); colormap(gray)
dsf = 8;
sub_dt = .02/dsf;

seq_beg = (cur_seq-1)*250*24;
cur_onsets = fonsets(seq_beg+1:(seq_beg+250),1);
cur_onsets = (cur_onsets-cur_onsets(1))/1000;
aviobj = avifile('example.avi','fps',15);
for i = 50:80
    subplot(3,3,4);
    imagesc(reshape(dstimps74{cur_seq}(i,:),32,32));
    subplot(3,3,sps(mod(sps,3)~=1))
    for j = 1:dsf
        xlim([-win_width win_width]/2+cur_onsets(i)+(j-1)*sub_dt);
        F = getframe(f1);
        aviobj = addframe(aviobj,F);
    end
end
close(f1)
aviobj = close(aviobj);