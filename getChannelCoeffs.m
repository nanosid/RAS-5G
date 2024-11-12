function [H] = getChannelCoeffs(hEst, scOffset, noRBs)
[~,~,N,M] = size(hEst);
scNo = scOffset + 1;
% Get channel coefficients over one RB just to get beamforming weights
hEst = hEst(scNo:scNo+(12*noRBs-1),:,:,:);
H = permute(mean(reshape(hEst,[],N,M)),[2,3,1]);
end