function [maxErr] = estimateCDLChannel(paths)

% setup carrier and CSI-RS config
carrier = nrCarrierConfig;
carrier.NSizeGrid = 25;
carrier.SubcarrierSpacing = 15;
carrier.NSlot = 1;
carrier.NFrame = 0

csirs = nrCSIRSConfig;
csirs.CSIRSType = {'nzp','zp'};
csirs.CSIRSPeriod = {[5 1],[5 1]};
csirs.Density = {'one','one'};
csirs.RowNumber = [3 5];
csirs.SymbolLocations = {1,6};
csirs.SubcarrierLocations = {6,4};
csirs.NumRB = 25

powerCSIRS = 0;
disp(['CSI-RS power scaling: ' num2str(powerCSIRS) ' dB']);
sym = nrCSIRS(carrier,csirs);
csirsSym = sym*db2mag(powerCSIRS);
csirsInd = nrCSIRSIndices(carrier,csirs);
ports = max(csirs.NumCSIRSPorts);   % Number of antenna ports
txGrid = nrResourceGrid(carrier,ports);
txGrid(csirsInd) = csirsSym;
%plotGrid(size(txGrid),csirsInd,csirsSym);
[txWaveform,ofdmInfo] = nrOFDMModulate(carrier,txGrid);

% setup CDL-C channel
v = 0;										% UE velocity in kmph
fc = 5.925e9;								% Carrier frequency in Hz
c = physconst('lightspeed');				% Speed of light in m/s
fd = (v*1000/3600)/c*fc;					% UE max Doppler shift in Hz

% Configure basic channel model parameters
channel = nrCDLChannel;
%channel.DelayProfile = 'CDL-C';
%channel.DelaySpread = 10e-9;

% Configure custom delay profile
channel.DelayProfile = 'custom';
R = paths;										% # of paths
load('CustomDelayProfile.mat');
channel.PathDelays = channelInfo.PathDelays(1:R);
channel.AveragePathGains = channelInfo.AveragePathGains(1:R);
channel.AnglesAoA = channelInfo.AnglesAoA(1:R);
channel.AnglesAoD = channelInfo.AnglesAoD(1:R);
channel.AnglesZoA = channelInfo.AnglesZoA(1:R);
channel.AnglesZoD = channelInfo.AnglesZoD(1:R);

channel.CarrierFrequency = fc;
channel.MaximumDopplerShift = fd;

% Configure transmit and receive array layout in the form [M N P Mg Ng]
% M = # of rows of antenna in each array
% N = # of columns of antenna in each array
% Mg = # of rows of array panels
% Ng = # of columns of array panels
% P = polarization angles
channel.TransmitAntennaArray.Size = [1 4 1 1 1];	% Single Tx array with 2 antennas
channel.ReceiveAntennaArray.Size = [1 1 1 1 1];		% Single Rx array with 1 antenna

chInfo = info(channel);
maxChDelay = ceil(max(chInfo.PathDelays*channel.SampleRate)) + chInfo.ChannelFilterDelay;
txWaveform = [txWaveform; zeros(maxChDelay,size(txWaveform,2))];

[rxWaveform,pathGains] = channel(txWaveform);
pathFilters = getPathFilters(channel);
H_actual = nrPerfectChannelEstimate(carrier,pathGains,pathFilters);
SNRdB = 5;           % in dB
SNR = 10^(SNRdB/10);  % Linear value
N0 = 1/sqrt(2.0*R*double(ofdmInfo.Nfft)*SNR); % Noise variance
rng(0);
noise = N0*complex(randn(size(rxWaveform)),randn(size(rxWaveform)));
rxWaveform = rxWaveform + noise;

% Disable ZP-CSI-RS resource, not going to be used for timing and channel
% estimation
csirs.CSIRSPeriod = {[5 1],'off'};
% Generate reference symbols and apply power scaling
refSym = db2mag(powerCSIRS)*nrCSIRS(carrier,csirs);
% Generate reference indices
refInd = nrCSIRSIndices(carrier,csirs);
offset = nrTimingEstimate(carrier,rxWaveform,refInd,refSym)
rxWaveform = rxWaveform(1+offset:end,:);
rxGrid = nrOFDMDemodulate(carrier,rxWaveform); % Of size K-by-L-by-R

% compare estimated and actual channel
cdmLen = [2 1]; % Corresponds to CDMType = 'FD-CDM2'
[H_est,nVar] = nrChannelEstimate(carrier,rxGrid,refInd,refSym,'CDMLengths',cdmLen);
estSNR = -10*log10(nVar);
disp(['estimated SNR = ' num2str(estSNR) ' dB'])

% figure;
% % Plot the estimated channel
% subplot(1,2,1)
% imagesc(abs(H_est(:,:,1,1)));
% colorbar;
% title('Estimated Channel')
% axis xy;
% xlabel('OFDM Symbols');
% ylabel('Subcarriers');
% 
% % Plot the actual channel
% subplot(1,2,2)
% imagesc(abs(H_actual(:,:,1,1)));
% colorbar;
% title('Actual Channel')
% axis xy;
% xlabel('OFDM Symbols');
% ylabel('Subcarriers');

H_err = (H_est - H_actual(:,:,:,1:size(H_est,4)));
[minErr,maxErr] = bounds(abs(H_err),'all');
disp(['Absolute value of the channel estimation error is in the range of [' num2str(minErr) ', ' num2str(maxErr) ']'])
end