function [wtx, wrx, D] = getBeamformingWeights(hEst, nLayers, scOffset, noRBs)
[~,~,N,M] = size(hEst);
scNo = scOffset + 1;
% Get channel coefficients over one RB just to get beamforming weights
hEst = hEst(scNo:scNo+(12*noRBs-1),:,:,:);
H = permute(mean(reshape(hEst,[],N,M)),[2,3,1]);
[U,D,V] = svd(H);
wtx = V(:,1:nLayers).';
wrx = U(:,1:nLayers)';
end