function [fglm] = fitGLMweights_jmm(glm,ustim,spkbs)
% USAGE: [fglm] = fitGLMweights(glm,ustim,spkbs)
%   fits modul-weights + constant of a glm

nmods = length(glm.mods); 

%compile a matrix representing the output of each NL module
fX    = []; 
for i = 1:nmods
	fX = [fX procModul(glm.mods(i),ustim)];
end

%% without regularisation
fglm = fitWeights_jmm(glm,fX,spkbs);

end

