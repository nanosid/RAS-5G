% Define basic parameters
fc = 5.925e9;
c = physconst("LightSpeed");

% Define gNB and RAS properties
gNBPos = [43.021952, -77.461666];				% Ionia, NY radio telescope
gNBAntSize = [4 4];
gNBAntDir = [0 0].';
RASPos = [42.929852, -77.500154];
RASAntSize = [1 1];
RASAntDir = [0 0].';
reflectionsOrder = 1;
SCS = 15;
NRB = 52;

% Define gNb and RAS objects
gNBSite = txsite("Name","Victor_gNB","Latitude",gNBPos(1),"Longitude",gNBPos(2),"AntennaAngle",gNBAntDir(1:2),...
	"AntennaHeight",32,"TransmitterFrequency",fc);
RASSite = rxsite("Name","Ionia_RAS","Latitude",RASPos(1),"Longitude",RASPos(2),"AntennaAngle",RASAntDir(1:2),...
	"AntennaHeight",5);

% Define Raytracing prop model and display plotted rays
pm = propagationModel("raytracing","Method","sbr","MaxNumReflections",reflectionsOrder);
RAS_rays = raytrace(gNBSite,RASSite,pm,"Type","pathloss");
%plot(rays{1})

% Obtain path gains and directions
% Normalize path delays to 0 seconds, first path occurs at t=0
pathToAs = [RAS_rays{1}.PropagationDelay] - min([RAS_rays{1}.PropagationDelay]);
avgPathGains = -[RAS_rays{1}.PathLoss];
pathAoDs = [RAS_rays{1}.AngleOfDeparture];
pathAoAs = [RAS_rays{1}.AngleOfArrival];
isLOS = any([RAS_rays{1}.LineOfSight]);

% Setup CDL channel model
RAS_channel = nrCDLChannel;
RAS_channel.DelayProfile = 'Custom';
RAS_channel.PathDelays = pathToAs;
RAS_channel.AveragePathGains = avgPathGains;
RAS_channel.AnglesAoD = pathAoDs(1,:);
RAS_channel.AnglesZoD = 90-pathAoDs(2,:);
RAS_channel.AnglesAoA = pathAoAs(1,:);
RAS_channel.AnglesZoA = 90-pathAoAs(2,:);
RAS_channel.HasLOSCluster = isLOS;
RAS_channel.CarrierFrequency = fc;
RAS_channel.NormalizePathGains = false;
RAS_channel.NormalizeChannelOutputs = false;

% Setup antenna array properties
lambda = c/fc;
RASArray = phased.NRRectangularPanelArray('Size',[RASAntSize(1:2) 1 1],'Spacing', [0.5*lambda*[1 1] 1 1]);
RASArray.ElementSet = {phased.IsotropicAntennaElement};   % isotropic antenna element
RAS_channel.ReceiveAntennaArray = RASArray;
RAS_channel.ReceiveArrayOrientation = [RASAntDir(1); (-1)*RASAntDir(2); 0];  % the (-1) converts elevation to downtilt

gNBArray = phased.NRRectangularPanelArray('Size',[gNBAntSize(1:2) 1 1],'Spacing', [0.5*lambda*[1 1] 1 1]);
gNBArray.ElementSet = {phased.IsotropicAntennaElement};
RAS_channel.TransmitAntennaArray = gNBArray;
RAS_channel.TransmitArrayOrientation = [gNBAntDir(1); (-1)*gNBAntDir(2); 0];   % the (-1) converts elevation to downtilt

% Obtain channel path gains
ofdmInfo = nrOFDMInfo(NRB,SCS);
RAS_channel.SampleRate = ofdmInfo.SampleRate;
%channel.ChannelFiltering = false;

% Design sample waveform
RAS_channelInfo = info(RAS_channel);
T = RAS_channel.SampleRate * 1e-3;
RAS_Nt = RAS_channelInfo.NumTransmitAntennas;
RAS_Nr = RAS_channelInfo.NumReceiveAntennas;
txWaveform = complex(randn(T,RAS_Nt),randn(T,RAS_Nt));
[RAS_rxWaveform,RAS_pathGains,RAS_sampleTimes] = RAS_channel(txWaveform);
noise = wgn(size(RAS_rxWaveform,1),size(RAS_rxWaveform,2),-137);
RAS_snr_default = snr(RAS_rxWaveform,noise);

% Setup sample UE to compare SNRs after beamforming and nullification
UEPos = [43.019708, -77.461838];
UEAntSize = [1 1];
UEAntDir = [0 0].';
UESite = rxsite("Name","Victor_UE","Latitude",UEPos(1),"Longitude",UEPos(2),"AntennaAngle",UEAntDir(1:2),...
	"AntennaHeight",2);
UE_rays = raytrace(gNBSite,UESite,pm,"Type","pathloss");
pathToAs = [UE_rays{1}.PropagationDelay] - min([UE_rays{1}.PropagationDelay]);
avgPathGains = -[UE_rays{1}.PathLoss];
pathAoDs = [UE_rays{1}.AngleOfDeparture];
pathAoAs = [UE_rays{1}.AngleOfArrival];
isLOS = any([UE_rays{1}.LineOfSight]);

UE_channel = nrCDLChannel;
UE_channel.DelayProfile = 'Custom';
UE_channel.PathDelays = pathToAs;
UE_channel.AveragePathGains = avgPathGains;
UE_channel.AnglesAoD = pathAoDs(1,:);
UE_channel.AnglesZoD = 90-pathAoDs(2,:);
UE_channel.AnglesAoA = pathAoAs(1,:);
UE_channel.AnglesZoA = 90-pathAoAs(2,:);
UE_channel.HasLOSCluster = isLOS;
UE_channel.CarrierFrequency = fc;
UE_channel.NormalizePathGains = false;
UE_channel.NormalizeChannelOutputs = false;
UE_channel.SampleRate = ofdmInfo.SampleRate;

UEArray = phased.NRRectangularPanelArray('Size',[UEAntSize(1:2) 1 1],'Spacing', [0.5*lambda*[1 1] 1 1]);
UEArray.ElementSet = {phased.IsotropicAntennaElement};   % isotropic antenna element
UE_channel.ReceiveAntennaArray = UEArray;
UE_channel.ReceiveArrayOrientation = [UEAntDir(1); (-1)*UEAntDir(2); 0];  % the (-1) converts elevation to downtilt
UE_channel.TransmitAntennaArray = gNBArray;
UE_channel.TransmitArrayOrientation = [gNBAntDir(1); (-1)*gNBAntDir(2); 0];   % the (-1) converts elevation to downtilt

UE_channelInfo = info(UE_channel);
T = UE_channel.SampleRate * 1e-3;
RAS_Nt = UE_channelInfo.NumTransmitAntennas;
[UE_rxWaveform,UE_pathGains,UE_sampleTimes] = UE_channel(txWaveform);
noise = wgn(size(UE_rxWaveform,1),size(UE_rxWaveform,2),-89);
UE_snr_default = snr(UE_rxWaveform,noise);

% Derive beamforming matrix to nullify RAS and focus on UE
% Get gNB-UE channel coefficients over all RBs and OFDM symbols
UE_pathFilters = getPathFilters(UE_channel);
nSlots = 0;
[UE_offset,~] = nrPerfectTimingEstimate(UE_pathGains, UE_pathFilters);
hest_gU = nrPerfectChannelEstimate(UE_pathGains, UE_pathFilters, NRB, SCS, nSlots, UE_offset, UE_sampleTimes);

% Get beamformer to maximize SNR at UE
nLayers = 1;
scOffset = 0;
noRBs = 1;
[w_gNB, w_UE, ~] = getBeamformingWeights(hest_gU, nLayers, scOffset, noRBs);

% Get gNB-RAS channel coefficients over all RBs and OFDM symbols
RAS_pathFilters = getPathFilters(RAS_channel);
[RAS_offset,~] = nrPerfectTimingEstimate(RAS_pathGains, RAS_pathFilters);
hest_gR = nrPerfectChannelEstimate(RAS_pathGains, RAS_pathFilters, NRB, SCS, nSlots, RAS_offset, RAS_sampleTimes);

% Now get beamformer which also nullifies RAS
hest_gR_temp = permute(mean(reshape(hest_gR,[],RAS_Nr,RAS_Nt)),[2,3,1]);
P = null(hest_gR_temp);
% h_gR = squeeze(mean(RAS_pathGains, 1));
% P = null(h_gR);
w_copt = P*P'*w_gNB.';

% SNR after beamforming
wf_bf = complex(randn(T,1),randn(T,1));
wf_bf = wf_bf*w_copt.';
[RAS_wf_bf,~,~] = RAS_channel(wf_bf);
RAS_snr_bf = snr(RAS_wf_bf,noise);
RAS_snr_gain = RAS_snr_bf - RAS_snr_default;
[UE_wf_bf,~,~] = UE_channel(wf_bf);
UE_snr_bf = snr(UE_wf_bf,noise);
UE_snr_gain = UE_snr_bf - UE_snr_default;

% Show beampattern of gNB
show(UESite);
show(gNBSite);
show(RASSite);
gNBSite.Antenna = clone(RAS_channel.TransmitAntennaArray); % need a clone, otherwise setting the Taper weights would affect the channel array
gNBSite.Antenna.Taper = w_copt;
pattern(gNBSite,fc,"Size",20000);