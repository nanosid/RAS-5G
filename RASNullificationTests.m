%Define basic parameters
fc = 5.925e9;
c = physconst("LightSpeed");
% i = 4
gNBLocs = [43.030833, -77.456944; 43.021952, -77.461666; 42.924344, -77.625622; 42.923448, -77.612715;	42.985532, -77.478548; ...
	42.973523, -77.482448; 42.963357, -77.486774; 42.955440, -77.488604; 42.941536, -77.493440; 42.930223, -77.501120];
UELocs = [43.029590, -77.447416; 43.019708, -77.461838; 42.926072, -77.624249; 42.921467, -77.610440; 42.986065, -77.480522; ...
	42.974998, -77.479186; 42.966559, -77.484972; 42.955329, -77.492702; 42.941755, -77.488698; 42.932849, -77.503660];
numants = [1 2 4 8 16 32];
N = 1000;
UE_snr_gain = zeros(length(gNBLocs),length(numants),N);
RAS_snr_gain = zeros(length(gNBLocs),length(numants),N);

% Loop over locations and num of antennas
for i=1:length(gNBLocs)
	for j=1:length(numants)
		parfor var=1:N
		% Define gNB and RAS properties
		gNBPos = gNBLocs(i,:);				
		gNBAntSize = [numants(j) 2];
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
		
		% Display gNB and RAS location on map
		%show(gNBSite);
		%show(RASSite);
		
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
		UEPos = UELocs(i,:);
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
		
		% Get gNB-RAS channel coefficients over all RBs and OFDM symbols
		RAS_pathFilters = getPathFilters(RAS_channel);
		[RAS_offset,~] = nrPerfectTimingEstimate(RAS_pathGains, RAS_pathFilters);
		hest_gR = nrPerfectChannelEstimate(RAS_pathGains, RAS_pathFilters, NRB, SCS, nSlots, RAS_offset, RAS_sampleTimes);
		
		% Get beamformer to maximize SNR at UE
		nLayers = 1;
		scOffset = 0;
		noRBs = 1;
		[w_gNB, w_UE, ~] = getBeamformingWeights(hest_gU, nLayers, scOffset, noRBs);
		
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
		RAS_snr_gain(i,j,var) = RAS_snr_bf - RAS_snr_default;
		[UE_wf_bf,~,~] = UE_channel(wf_bf);
		UE_snr_bf = snr(UE_wf_bf,noise);
		UE_snr_gain(i,j,var) = UE_snr_bf - UE_snr_default;
	end
end
end
% figure;
% surf(RAS_snr_gain);
% figure;
% surf(UE_snr_gain);