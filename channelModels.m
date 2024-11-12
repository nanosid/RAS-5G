%% Testing different CDL channel model delay profiles and number of paths

v = 0;										% UE velocity in kmph
fc = 5.925e9;								% Carrier frequency in Hz
c = physconst('lightspeed');				% Speed of light in m/s
fd = (v*1000/3600)/c*fc;					% UE max Doppler shift in Hz

% Configure basic channel model parameters
channel = nrCDLChannel;
channel.DelayProfile = 'CDL-C';
channel.DelaySpread = 10e-9;

% Configure custom delay profile
% channel.DelayProfile = 'custom';
% R = 8;										% # of paths
% load('CustomDelayProfile.mat');
% channel.PathDelays = channelInfo.PathDelays(1:R);
% channel.AveragePathGains = channelInfo.AveragePathGains(1:R);
% channel.AnglesAoA = channelInfo.AnglesAoA(1:R);
% channel.AnglesAoD = channelInfo.AnglesAoD(1:R);
% channel.AnglesZoA = channelInfo.AnglesZoA(1:R);
% channel.AnglesZoD = channelInfo.AnglesZoD(1:R);

channel.CarrierFrequency = fc;
channel.MaximumDopplerShift = fd;

% Configure transmit and receive array layout in the form [M N P Mg Ng]
% M = # of rows of antenna in each array
% N = # of columns of antenna in each array
% Mg = # of rows of array panels
% Ng = # of columns of array panels
% P = polarization angles
channel.TransmitAntennaArray.Size = [1 2 1 1 1];	% Single Tx array with 2 antennas
channel.ReceiveAntennaArray.Size = [1 1 1 1 1];		% Single Rx array with 1 antenna

displayChannel(channel,'LinkEnd','Tx');
displayChannel(channel,'LinkEnd','Rx');