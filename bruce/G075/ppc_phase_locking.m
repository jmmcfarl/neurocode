function [ppc0] = ppc_phase_locking(phases)

phase_vecs = [cos(phases) sin(phases)];

phase_dist = pdist(phase_vecs,@phase_dist);
ppc0 = mean(phase_dist);

end

function D = phase_dist(XI,XJ)
D = XJ(:,1)*XI(1) + XJ(:,2)*XI(2);
end
