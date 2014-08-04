clear all
close all
cd ~/Data/bruce/2_27_12

%% load parallel data
load lemM232.Cells.rls.Pp.mat
par = AllCellRes;
%% load orthoganol data
load lemM232.Cells.rls.Op.mat
orth = AllCellRes;

%%
cd ~/James_scripts/bruce

for i = 1:10
    subplot(2,1,1)
    plot(par(i).xpos,par(i).means,orth(i).xpos,orth(i).means,'r')
    xlim([3.5 6])
    subplot(2,1,2)
    plot(par(i).ypos,par(i).means,orth(i).ypos,orth(i).means,'r')
    xlim([-3.5 -1])
i
pause(2)
end