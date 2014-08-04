function [nodeMap,edgeMap] = UGM_makeMRFmaps(nNodes,edgeStruct,ising,tied)
% Assumes that all nodes have the same number of states

nStates = edgeStruct.nStates(1);
nEdges = edgeStruct.nEdges;

nodeMap = zeros(nNodes,nStates,'int32');
if tied
	for s = 1:nStates
		nodeMap(:,s) = s;
	end
else
	nodeMap(:) = 1:numel(nodeMap);
end
nNodeParams = max(nodeMap(:));

edgeMap = zeros(nStates,nStates,nEdges,'int32');
if tied
	switch ising
		case 1
			for s = 1:nStates
				edgeMap(s,s,:) = nNodeParams+1;
			end
		case 2
			for s = 1:nStates
				edgeMap(s,s,:) = nNodeParams+s;
			end
		case 0
			s = 1;
			for s1 = 1:nStates
				for s2 = 1:nStates
					edgeMap(s1,s2,:) = nNodeParams+s;
					s = s+1;
				end
			end
	end
else
	switch ising
		case 1
			for e = 1:nEdges
				for s = 1:nStates
					edgeMap(s,s,e) = nNodeParams+e;
				end
			end
		case 2
			se = 1;
			for e = 1:nEdges
				for s = 1:nStates
					edgeMap(s,s,e) = nNodeParams+se;
					se = se+1;
				end
			end
		case 0
			sse = 1;
			for e = 1:nEdges
				for s1 = 1:nStates
					for s2 = 1:nStates
						edgeMap(s1,s2,e) = nNodeParams+sse;
						sse = sse+1;
					end
				end
			end
	end
end