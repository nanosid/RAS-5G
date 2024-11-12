clear;
% Define constant parameters
c = physconst("LightSpeed");
gNBAntDir = [0 0].';
RASPos = [42.929852, -77.500154];
RASAntSize = [1 1];
RASAntDir = [0 0].';
UEAntSize = [1 1];
UEAntDir = [0 0].';
SCS = 15;
NRB = 52;
nSlots = 0;
nLayers = 1;
scOffset = 0;
noRBs = 1;
RASSite = rxsite("Name","Ionia_RAS","Latitude",RASPos(1),"Longitude",RASPos(2),"AntennaAngle",RASAntDir(1:2),"AntennaHeight",5);

% Define possible values for metrics
gNBLocs = [43.030833, -77.456944; 43.021952, -77.461666; 42.924344, -77.625622; 42.923448, -77.612715;	42.985532, -77.478548; ...
	42.973523, -77.482448; 42.963357, -77.486774; 42.955440, -77.488604; 42.941536, -77.493440; 42.930223, -77.501120];
numants = [2 4 8 16];
fcs = [739e6 1880e6 5925e6 7125e6];
numrefs = [1 2 3 4];
numUEs = [1 8 16 32];
gNBDist = zeros(1,10);
for i = 1:size(gNBLocs,1)
	[gNBDist(i),~] = distance(gNBLocs(i,1),gNBLocs(i,2),RASPos(1),RASPos(2));
end
gNBDist = deg2rad(gNBDist).*(earthRadius*1e-3);
N = 16;
UE_snr_default = zeros(length(gNBLocs),length(numants),length(fcs),length(numrefs),length(numUEs),N);
UE_snr_bf = zeros(size(UE_snr_default));
UE_snr_gain = zeros(size(UE_snr_bf));
RAS_snr_default = zeros(size(UE_snr_gain));
RAS_snr_bf = zeros(size(UE_snr_gain));
RAS_snr_gain = zeros(size(RAS_snr_default));
latency = zeros(size(RAS_snr_gain));

for i1 = 1:length(gNBLocs)
	for i2 = 1:length(numants)
		for i3 = 1:length(fcs)
			for i4 = 1:length(numrefs)
				for i5 = 1:length(numUEs)
					iterationID = sprintf("Iteration %d,%d,%d,%d,%d...\n",i1,i2,i3,i4,i5);
					disp(iterationID);
					parfor var = 1:N
						gNBPos = gNBLocs(i1,:);				
						gNBAntSize = [numants(i2) 2];
						fc = fcs(i3);
						reflectionsOrder = numrefs(i4);
						K = numUEs(i5);
						gNBSite = txsite("Name","Victor_gNB","Latitude",gNBPos(1),"Longitude",gNBPos(2),"AntennaAngle",...
							gNBAntDir(1:2),"AntennaHeight",32,"TransmitterFrequency",fc);
						
						% Define Raytracing prop model and display plotted rays
						pm = propagationModel("raytracing","Method","sbr","MaxNumReflections",reflectionsOrder);
						RAS_rays = raytrace(gNBSite,RASSite,pm,"Type","pathloss");
	
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
						RAS_channel.ReceiveArrayOrientation = [RASAntDir(1); (-1)*RASAntDir(2); 0]; 
						
						gNBArray = phased.NRRectangularPanelArray('Size',[gNBAntSize(1:2) 1 1],'Spacing', [0.5*lambda*[1 1] 1 1]);
						gNBArray.ElementSet = {phased.IsotropicAntennaElement};
						RAS_channel.TransmitAntennaArray = gNBArray;
						RAS_channel.TransmitArrayOrientation = [gNBAntDir(1); (-1)*gNBAntDir(2); 0];

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
						if mean(isnan(txWaveform),'all')
							nanIdx = isnan(txWaveform);
							txWaveform(nanIdx) = 0;
						end
						[RAS_rxWaveform,RAS_pathGains,RAS_sampleTimes] = RAS_channel(txWaveform);
						noise_RAS = wgn(size(RAS_rxWaveform,1),size(RAS_rxWaveform,2),-137);
						RAS_snr_default(i1,i2,i3,i4,i5,var) = snr(RAS_rxWaveform,noise_RAS);

						% Get gNB-RAS channel coefficients over all RBs and OFDM symbols
						RAS_pathFilters = getPathFilters(RAS_channel);
						[RAS_offset,~] = nrPerfectTimingEstimate(RAS_pathGains, RAS_pathFilters);
						hest_gR = nrPerfectChannelEstimate(RAS_pathGains, RAS_pathFilters, NRB, SCS, nSlots,...
							RAS_offset, RAS_sampleTimes);
						hest_gR_temp = permute(mean(reshape(hest_gR,[],RAS_Nr,RAS_Nt)),[2,3,1]);

						% Define UE sites
						% locDisp = abs(randn(K,2).*5e-3);
						% UEPos = locDisp + [ones(K,1).*gNBPos(1) ones(K,1).*gNBPos(2)];
						UENames = "UE"+string(1:K);
						UEPos = zeros(K,2);
						UESites = rxsite("Name",UENames,"Latitude",0,"Longitude",0,"AntennaAngle",UEAntDir(1:2),"AntennaHeight",2);

						% Now derive channel for K UEs and the beamforming weights for each UE
						hest_gU_all = zeros(K,RAS_Nt);
						w_gNB_all = zeros(K,RAS_Nt);
						noise_UE = zeros(size(RAS_rxWaveform,1),K);

						UE_channel = nrCDLChannel;
						for i=1:K
							pathToAs = [];
    						while(isempty(pathToAs))
        						locDisp = abs(randn(1,2).*5e-3);
        						UEPos(i,:) = locDisp + gNBPos;
        						UESites(i).Latitude = UEPos(i,1);
        						UESites(i).Longitude = UEPos(i,2);
        						UE_rays = raytrace(gNBSite,UESites(i),pm,"Type","pathloss");
        						pathToAs = [UE_rays{1}.PropagationDelay] - min([UE_rays{1}.PropagationDelay]);
    						end
% 							if isempty(pathToAs)
% 								continue;
% 							else
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
								
								UEArray = phased.NRRectangularPanelArray('Size',[UEAntSize(1:2) 1 1],'Spacing',...
									[0.5*lambda*[1 1] 1 1]);
								UEArray.ElementSet = {phased.IsotropicAntennaElement};   % isotropic antenna element
								UE_channel.ReceiveAntennaArray = UEArray;
								UE_channel.ReceiveArrayOrientation = [UEAntDir(1); (-1)*UEAntDir(2); 0];
								UE_channel.TransmitAntennaArray = gNBArray;
								UE_channel.TransmitArrayOrientation = [gNBAntDir(1); (-1)*gNBAntDir(2); 0];
								
								UE_channelInfo = info(UE_channel);
								T = UE_channel.SampleRate * 1e-3;
								UE_Nt = UE_channelInfo.NumTransmitAntennas;
								UE_Nr = UE_channelInfo.NumReceiveAntennas;
								[~,UE_pathGains,UE_sampleTimes] = UE_channel(txWaveform);
								UE_pathFilters = getPathFilters(UE_channel);
								[UE_offset,~] = nrPerfectTimingEstimate(UE_pathGains, UE_pathFilters);
								hest_gU = nrPerfectChannelEstimate(UE_pathGains, UE_pathFilters, NRB, SCS, nSlots,...
									UE_offset, UE_sampleTimes);
								hest_gU_all(i,:) = getChannelCoeffs(hest_gU, scOffset, noRBs);
								UE_rxWaveform = txWaveform*hest_gU_all(i,:)';
                                noise_UE(:,i) = wgn(size(UE_rxWaveform,1),size(UE_rxWaveform,2),-89);
								UE_snr_default(i1,i2,i3,i4,i5,var) = UE_snr_default(i1,i2,i3,i4,i5,var) + ...
									snr(UE_rxWaveform,noise_UE(:,i));
                                [w_gNB_all(i,:), ~, ~] = getBeamformingWeights(hest_gU, nLayers, scOffset, noRBs);
% 							end
						end
						UE_snr_default(i1,i2,i3,i4,i5,var) = UE_snr_default(i1,i2,i3,i4,i5,var)/K;
						tic
						hest_gU_all = hest_gU_all./norm(hest_gU_all);
						[U,D,V] = svd(hest_gU_all);
						w_gNB = V(:,1:nLayers).';
						P = null(hest_gR_temp);
						w_copt = P*P'*w_gNB.';
						latency(i1,i2,i3,i4,i5,var) = toc;

						wf = complex(randn(T,1),randn(T,1));
						wf_bf = wf*w_copt.';
						if mean(isnan(wf_bf),'all')
							nanIdx = isnan(wf_bf);
							wf_bf(nanIdx) = 0;
						end
						[RAS_wf_bf,~,~] = RAS_channel(wf_bf);
						RAS_snr_bf(i1,i2,i3,i4,i5,var) = snr(RAS_wf_bf,noise_RAS);
						RAS_snr_gain(i1,i2,i3,i4,i5,var) = RAS_snr_bf(i1,i2,i3,i4,i5,var) - RAS_snr_default(i1,i2,i3,i4,i5,var);

						% SNR at K UEs
						for i = 1:K
							if(UE_snr_default(i1,i2,i3,i4,i5,var) == 0)
								continue;
							else
								UE_wf_bf = wf_bf*hest_gU_all(i,:)';
								UE_snr_bf(i1,i2,i3,i4,i5,var) = UE_snr_bf(i1,i2,i3,i4,i5,var) + snr(UE_wf_bf,noise_UE(:,i));
							end
						end
						UE_snr_bf(i1,i2,i3,i4,i5,var) = UE_snr_bf(i1,i2,i3,i4,i5,var)/K;
						UE_snr_gain(i1,i2,i3,i4,i5,var) = (UE_snr_bf(i1,i2,i3,i4,i5,var) - UE_snr_default(i1,i2,i3,i4,i5,var))./5;
					end
				end
			end
		end
	end
end
str_N = sprintf("N_%d",N);
save(str_N);