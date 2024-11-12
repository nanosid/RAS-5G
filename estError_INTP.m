clear;
% Define constant parameters
c = physconst("LightSpeed");
e = 10^(-220/10);
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
gNBLocs = [42.930223, -77.501120];
numants = [4 8 16];
fcs = [739e6 1880e6 5925e6];
numUEs = [1 8 16];
gNBDist = distance(gNBLocs(1),gNBLocs(2),RASPos(1),RASPos(2));
gNBDist = deg2rad(gNBDist)*(earthRadius*1e-3);
N = 2;
mu = 0:1e-3:1e-2;
sigma = [1e-4, 5e-4, 1e-3];

RAS_null_heu = zeros(3,length(mu),length(sigma),length(numants),N);
RAS_null_heu_err = zeros(size(RAS_null_heu));
RAS_null_opt = zeros(size(RAS_null_heu));
RAS_null_opt_err = zeros(size(RAS_null_heu));

for rowidx = 1:3
for i1 = 1:length(mu)
	for i2 = 1:length(sigma)
		for i3 = 1:length(numants)
			iterationID = sprintf("Iteration %d,%d,%d,%d...\n",rowidx,i1,i2,i3);
			disp(iterationID);
			gNBPos = gNBLocs;				
			reflectionsOrder = 3;
			switch rowidx
				case 1
					gNBAntSize = [numants(i3) 2];
					fc = 5925e6;
					K = 16;
				case 2
					gNBAntSize = [8 2];
					fc = fcs(i3);
					K = 16;
				case 3
					gNBAntSize = [8 2];
					fc = 5925e6;
					K = numUEs(i3);
			end
			M = prod(gNBAntSize);
			for var = 1:N
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

				% Get gNB-RAS channel coefficients over all RBs and OFDM symbols
				RAS_pathFilters = getPathFilters(RAS_channel);
				[RAS_offset,~] = nrPerfectTimingEstimate(RAS_pathGains, RAS_pathFilters);
				hest_gR = nrPerfectChannelEstimate(RAS_pathGains, RAS_pathFilters, NRB, SCS, nSlots,...
					RAS_offset, RAS_sampleTimes);
				hest_gR_temp = permute(mean(reshape(hest_gR,[],RAS_Nr,RAS_Nt)),[2,3,1]);
				hest_gR_temp = hest_gR_temp./norm(hest_gR_temp);
				hest_gR_err = hest_gR_temp + mu(i1) + sigma(i2)*complex(randn(size(hest_gR_temp)),randn(size(hest_gR_temp)));

				% Define UE sites
				% locDisp = abs(randn(K,2).*5e-3);
				% UEPos = locDisp + [ones(K,1).*gNBPos(1) ones(K,1).*gNBPos(2)];
				UENames = "UE"+string(1:K);
				UEPos = zeros(K,2);
				UESites = rxsite("Name",UENames,"Latitude",0,"Longitude",0,"AntennaAngle",UEAntDir(1:2),"AntennaHeight",2);

				% Now derive channel for K UEs and the beamforming weights for each UE
				hest_gU_all = zeros(K,RAS_Nt);
				hest_gU_err = zeros(K,RAS_Nt);
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
					[UE_rxWaveform,UE_pathGains,UE_sampleTimes] = UE_channel(txWaveform);
					UE_pathFilters = getPathFilters(UE_channel);
					[UE_offset,~] = nrPerfectTimingEstimate(UE_pathGains, UE_pathFilters);
					hest_gU = nrPerfectChannelEstimate(UE_pathGains, UE_pathFilters, NRB, SCS, nSlots,...
						UE_offset, UE_sampleTimes);
					hest_gU_all(i,:) = getChannelCoeffs(hest_gU, scOffset, noRBs);
					hest_gU_all(i,:) = hest_gU_all(i,:)./norm(hest_gU_all(i,:));
					noise_UE(:,i) = wgn(size(UE_rxWaveform,1),size(UE_rxWaveform,2),-59);
				end
				% No error
				[~,~,V] = svd(hest_gU_all);
				w_gNB = V(:,1:nLayers).';
				P = null(hest_gR_temp);
				w_copt = P*P'*w_gNB.';

				% Error in h_gR
				P_err = null(hest_gR_err);
				w_copt_err_gR = P_err*P_err'*w_gNB.';

				% Opt no error
				cvx_begin
					variable w(1,M) complex
					minimize(norm(hest_gU_all*w' - ones(K,1)));
					subject to
						norm(hest_gR_temp*w') <= e;
				cvx_end
				w_opt = w;

				% Opt error in h_gR
				cvx_begin
					variable w(1,M) complex
					minimize(norm(hest_gU_all*w' - ones(K,1)));
					subject to
						norm(hest_gR_err*w') <= e;
				cvx_end
				w_opt_err_gR = w;

				RAS_null_heu_err(rowidx,i1,i2,i3,var) = 10*log10(norm(hest_gR_temp*w_copt_err_gR).^2);
				RAS_null_heu(rowidx,i1,i2,i3,var) = 10*log10(norm(hest_gR_temp*w_copt).^2);
				RAS_null_opt_err(rowidx,i1,i2,i3,var) = 10*log10(norm(hest_gR_temp*w_opt_err_gR').^2);
				RAS_null_opt(rowidx,i1,i2,i3,var) = 10*log10(norm(hest_gR_temp*w_opt').^2);
			end
		end
	end
	end
end

%% Plots
plot(mu,permute(mean(RAS_null_heu(1,:,:,1,:) - RAS_null_heu_err(1,:,:,1,:),5),[3 2 1 4]))
hold on
plot(mu,permute(mean(RAS_null_opt(1,:,:,1,:) - RAS_null_opt_err(1,:,:,1,:),5),[3 2 1 4]))
semiPlots
leg = legend("w_{heu}, \sigma^2=10^{-5}","w_{heu}, \sigma^2=5\times10^{-5}","w_{heu}, \sigma^2=10^{-4}",...
	"w_{opt}, \sigma^2=10^{-5}","w_{opt}, \sigma^2=5\times10^{-5}","w_{opt}, \sigma^2=10^{-4}");
xlabel("\mu");
ylabel("\Delta IP_{hgR}");
leg.ItemTokenSize = [10 18];
leg.FontSize = 7;
title("M=8, f_c=5.9, K=16");
