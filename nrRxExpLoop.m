clear;
N = 5;
snr = zeros(N,1);
rssi = zeros(N,1);
var = 1;
% for var = 1:N
% 	nrRx;
% 	snr(var) = max(dmrsEst)
% 	rssi(var) = meas.RSSIPerAntenna
% end

%% tests
nrRxGetRSSI;
snr(var) = max(dmrsEst)
rssi(var) = meas.RSSIPerAntenna
var = var + 1;