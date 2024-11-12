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
numants = [1 2 4 8 16];
fc = 5925e6;
numrefs = [1 2 3 4];
numUEs = [1 8 16];
gNBDist = zeros(1,10);
for i = 1:size(gNBLocs,1)
	[gNBDist(i),~] = distance(gNBLocs(i,1),gNBLocs(i,2),RASPos(1),RASPos(2));
end
gNBDist = deg2rad(gNBDist).*(earthRadius*1e-3);
N = 10000;
latency_ants = zeros(length(numants),length(numUEs),N);
latency_refs = zeros(length(numrefs),length(numUEs),N);
latency_dist = zeros(length(gNBLocs),length(numUEs),N);

for rowidx = 1:3
	switch rowidx
		case 1
			L = length(numants);
		case 2
			L = length(numrefs);
		case 3
			L = length(gNBLocs);
	end
	for i1 = 1:L
		for i2 = 1:length(numUEs)
			iterationID = sprintf("Iteration %d,%d,%d...\n",rowidx,i1,i2);
			disp(iterationID);
			K = numUEs(i2);
			switch rowidx
				case 1
					gNBAntSize = [numants(i1) 2];
					reflectionsOrder = 3;
					gNBPos = gNBLocs(1,:);
				case 2
					gNBAntSize = [8 2];
					reflectionsOrder = numrefs(i1);
					gNBPos = gNBLocs(1,:);
				case 3
					gNBAntSize = [8 2];
					reflectionsOrder = 3;
					gNBPos = gNBLocs(i1,:);
			end
			M = prod(gNBAntSize);
			for var = 1:N
% 				gNBSite = txsite("Name","Victor_gNB","Latitude",gNBPos(1),"Longitude",gNBPos(2),"AntennaAngle",...
% 					gNBAntDir(1:2),"AntennaHeight",32,"TransmitterFrequency",fc);
% 				
% 				% Define Raytracing prop model and display plotted rays
% 				pm = propagationModel("raytracing","Method","sbr","MaxNumReflections",reflectionsOrder);
% 				RAS_rays = raytrace(gNBSite,RASSite,pm,"Type","pathloss");
% 	
% 				% Obtain channel path gains
% 				ofdmInfo = nrOFDMInfo(NRB,SCS);
% 				RAS_channel = getChannelObj(RAS_rays,fc,ofdmInfo.SampleRate);
% 				RAS_channel.SampleRate = ofdmInfo.SampleRate;
% 				%channel.ChannelFiltering = false;
% 	
% 				% Setup antenna array properties
% 				lambda = c/fc;
% 				RASArray = phased.NRRectangularPanelArray('Size',[RASAntSize(1:2) 1 1],'Spacing', [0.5*lambda*[1 1] 1 1]);
% 				RASArray.ElementSet = {phased.IsotropicAntennaElement};   % isotropic antenna element
% 				RAS_channel.ReceiveAntennaArray = RASArray;
% 				RAS_channel.ReceiveArrayOrientation = [RASAntDir(1); (-1)*RASAntDir(2); 0]; 
% 				
% 				gNBArray = phased.NRRectangularPanelArray('Size',[gNBAntSize(1:2) 1 1],'Spacing', [0.5*lambda*[1 1] 1 1]);
% 				gNBArray.ElementSet = {phased.IsotropicAntennaElement};
% 				RAS_channel.TransmitAntennaArray = gNBArray;
% 				RAS_channel.TransmitArrayOrientation = [gNBAntDir(1); (-1)*gNBAntDir(2); 0];
% 				
% 				% Design sample waveform
% 				RAS_channelInfo = info(RAS_channel);
% 				T = RAS_channel.SampleRate * 1e-3;
% 				RAS_Nt = RAS_channelInfo.NumTransmitAntennas;
% 				RAS_Nr = RAS_channelInfo.NumReceiveAntennas;
% 				txWaveform = complex(randn(T,RAS_Nt),randn(T,RAS_Nt));
% 				if mean(isnan(txWaveform),'all')
% 					nanIdx = isnan(txWaveform);
% 					txWaveform(nanIdx) = 0;
% 				end
% 				[RAS_rxWaveform,RAS_pathGains,RAS_sampleTimes] = RAS_channel(txWaveform);
% 				noise_RAS = wgn(size(RAS_rxWaveform,1),size(RAS_rxWaveform,2),-137);
% 	% 			RAS_snr_default(rowidx,i1,i2,var) = snr(RAS_rxWaveform,noise_RAS);
% 	
% 				% Get gNB-RAS channel coefficients over all RBs and OFDM symbols
% 				RAS_pathFilters = getPathFilters(RAS_channel);
% 				[RAS_offset,~] = nrPerfectTimingEstimate(RAS_pathGains, RAS_pathFilters);
% 				hest_gR = nrPerfectChannelEstimate(RAS_pathGains, RAS_pathFilters, NRB, SCS, nSlots,...
% 					RAS_offset, RAS_sampleTimes);
% 				hest_gR_temp = permute(mean(reshape(hest_gR,[],RAS_Nr,RAS_Nt)),[2,3,1]);
% 				hest_gR_temp = hest_gR_temp./norm(hest_gR_temp);
% 	
% 				% Define UE sites
% 				% locDisp = abs(randn(K,2).*5e-3);
% 				% UEPos = locDisp + [ones(K,1).*gNBPos(1) ones(K,1).*gNBPos(2)];
% 				K = numUEs(i2);
% 				UENames = "UE"+string(1:K);
% 				UEPos = zeros(K,2);
% 				UESites = rxsite("Name",UENames,"Latitude",0,"Longitude",0,"AntennaAngle",UEAntDir(1:2),"AntennaHeight",2);
% 	
% 				% Now derive channel for K UEs and the beamforming weights for each UE
% 				hest_gU_all = zeros(K,RAS_Nt);
% 				w_gNB_all = zeros(K,RAS_Nt);
% 				noise_UE = zeros(size(RAS_rxWaveform,1),K);
% 				UE_channel = nrCDLChannel;
% 				for i=1:K
% 					pathToAs = [];
% 					while(isempty(pathToAs))
% 						locDisp = abs(randn(1,2).*5e-3);
% 						UEPos(i,:) = locDisp + gNBPos;
% 						UESites(i).Latitude = UEPos(i,1);
% 						UESites(i).Longitude = UEPos(i,2);
% 						UE_rays = raytrace(gNBSite,UESites(i),pm,"Type","pathloss");
% 						pathToAs = [UE_rays{1}.PropagationDelay] - min([UE_rays{1}.PropagationDelay]);
% 					end
% 					UE_channel = getChannelObj(UE_rays,fc,ofdmInfo.SampleRate);
% 					UE_channel.SampleRate = ofdmInfo.SampleRate;
% 					
% 					UEArray = phased.NRRectangularPanelArray('Size',[UEAntSize(1:2) 1 1],'Spacing',...
% 						[0.5*lambda*[1 1] 1 1]);
% 					UEArray.ElementSet = {phased.IsotropicAntennaElement};   % isotropic antenna element
% 					UE_channel.ReceiveAntennaArray = UEArray;
% 					UE_channel.ReceiveArrayOrientation = [UEAntDir(1); (-1)*UEAntDir(2); 0];
% 					UE_channel.TransmitAntennaArray = gNBArray;
% 					UE_channel.TransmitArrayOrientation = [gNBAntDir(1); (-1)*gNBAntDir(2); 0];
% 					
% 					UE_channelInfo = info(UE_channel);
% 					T = UE_channel.SampleRate * 1e-3;
% 					UE_Nt = UE_channelInfo.NumTransmitAntennas;
% 					UE_Nr = UE_channelInfo.NumReceiveAntennas;
% 					[UE_rxWaveform,UE_pathGains,UE_sampleTimes] = UE_channel(txWaveform);
% 					UE_pathFilters = getPathFilters(UE_channel);
% 					[UE_offset,~] = nrPerfectTimingEstimate(UE_pathGains, UE_pathFilters);
% 					hest_gU = nrPerfectChannelEstimate(UE_pathGains, UE_pathFilters, NRB, SCS, nSlots,...
% 						UE_offset, UE_sampleTimes);
% 					hest_gU_all(i,:) = getChannelCoeffs(hest_gU, scOffset, noRBs);
% 					hest_gU_all(i,:) = hest_gU_all(i,:)./norm(hest_gU_all(i,:));
%                 	noise_UE(:,i) = wgn(size(UE_rxWaveform,1),size(UE_rxWaveform,2),-89);
% 				end
				hest_gU_all = randn(K,M,'like',1i);
				hest_gR_temp = randn(1,M,'like',1i);
				tic
				[U,~,V] = svd(hest_gU_all);
				w_gNB = V(:,1:nLayers).';
				P = null(hest_gR_temp);
				w_copt = P*P'*w_gNB.';
				latency = toc/1e-3;
				switch rowidx
					case 1
						latency_ants(i1, i2, var) = latency;
					case 2
						latency_refs(i1, i2, var) = latency;
					case 3
						latency_dist(i1, i2, var) = latency;
				end
			end
		end
	end
end

%% Plots
figure;
plot(1:length(numants),mean(latency_ants,3));
xlabel("No. of Antennas");
ylabel("Latency (ms)");
semiPlots;
leg = legend("K=1","K=8","K=16","location","northwest");
leg.ItemTokenSize = [10 18];
leg.Orientation = 'horizontal';
ylim([0 0.2]);
% yticks([0.2 0.4 0.6 0.8 1]);
xticks(1:length(numants));
xticklabels(numants.*2);

figure;
plot(numrefs,mean(latency_refs,3));
xlabel("Max No. of Reflections");
ylabel("Latency (ms)");
semiPlots;
xticks(numrefs)
leg = legend("K=1","K=8","K=16","location","northwest");
leg.ItemTokenSize = [10 18];
leg.Orientation = 'horizontal';
ylim([0 0.2]);
% yticks([0.2 0.4 0.6 0.8 1]);

figure;
plot(gNBDist,mean(latency_dist,3));
xlabel("gNB Distance (km)");
ylabel("Latency (ms)");
semiPlots;
leg = legend("K=1","K=8","K=16","location","northwest");
leg.ItemTokenSize = [10 18];
leg.Orientation = 'horizontal';
ylim([0 0.2]);
% yticks([0.2 0.4 0.6 0.8 1]);