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
numants = [4 8 16];
fcs = [739e6 1880e6 5925e6];
numrefs = [1 2 3 4];
numUEs = [1 8 16];
gNBDist = zeros(1,10);
for i = 1:size(gNBLocs,1)
	[gNBDist(i),~] = distance(gNBLocs(i,1),gNBLocs(i,2),RASPos(1),RASPos(2));
end
gNBDist = deg2rad(gNBDist).*(earthRadius*1e-3);
N = 2;
UE_snr_default = zeros(3,length(gNBLocs),length(numants),N,numUEs(end));
UE_snr_bf = zeros(size(UE_snr_default));
% UE_snr_gain = zeros(size(UE_snr_bf));
RAS_null_default = zeros(3,length(gNBLocs),length(numants),N);
RAS_null_bf = zeros(size(RAS_null_default));
% RAS_snr_gain = zeros(size(UE_snr_default));
sumrate_default = zeros(size(RAS_null_default));
sumrate_bf = zeros(size(RAS_null_default));
latency = zeros(size(RAS_null_default));

for rowidx = 3
for i1 = 1:length(gNBLocs)
	for i2 = 1:length(numants)
		iterationID = sprintf("Iteration %d,%d,%d...\n",rowidx,i1,i2);
		disp(iterationID);
		gNBPos = gNBLocs(i1,:);				
		reflectionsOrder = 3;
		switch rowidx
			case 1
				gNBAntSize = [numants(i2) 2];
				fc = 5925e6;
				K = 16;
			case 2
				gNBAntSize = [8 2];
				fc = fcs(i2);
				K = 16;
			case 3
				gNBAntSize = [8 2];
				fc = 5925e6;
				K = numUEs(i2);
		end
		for var = 1:N
			gNBSite = txsite("Name","Victor_gNB","Latitude",gNBPos(1),"Longitude",gNBPos(2),"AntennaAngle",...
				gNBAntDir(1:2),"AntennaHeight",32,"TransmitterFrequency",fc);
			
			% Define Raytracing prop model and display plotted rays
			pm = propagationModel("raytracing","Method","sbr","MaxNumReflections",reflectionsOrder);
			RAS_rays = raytrace(gNBSite,RASSite,pm,"Type","pathloss");

			% Obtain channel path gains
			ofdmInfo = nrOFDMInfo(NRB,SCS);
			RAS_channel = getChannelObj(RAS_rays,fc,ofdmInfo.SampleRate);
			RAS_channel.SampleRate = ofdmInfo.SampleRate;
			%channel.ChannelFiltering = false;

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
% 			RAS_snr_default(rowidx,i1,i2,var) = snr(RAS_rxWaveform,noise_RAS);

			% Get gNB-RAS channel coefficients over all RBs and OFDM symbols
			RAS_pathFilters = getPathFilters(RAS_channel);
			[RAS_offset,~] = nrPerfectTimingEstimate(RAS_pathGains, RAS_pathFilters);
			hest_gR = nrPerfectChannelEstimate(RAS_pathGains, RAS_pathFilters, NRB, SCS, nSlots,...
				RAS_offset, RAS_sampleTimes);
			hest_gR_temp = permute(mean(reshape(hest_gR,[],RAS_Nr,RAS_Nt)),[2,3,1]);
			hest_gR_temp = hest_gR_temp./norm(hest_gR_temp);

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
				UE_channel = getChannelObj(UE_rays,fc,ofdmInfo.SampleRate);
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
                noise_UE(:,i) = wgn(size(UE_rxWaveform,1),size(UE_rxWaveform,2),-89);
			end
			tic
			[U,~,V] = svd(hest_gU_all);
			w_gNB = V(:,1:nLayers).';
			P = null(hest_gR_temp);
			w_copt = P*P'*w_gNB.';
			latency(rowidx,i1,i2,var) = toc/1e-3;

			wf = complex(randn(T,1),randn(T,1));
			for i = 1:K
				UE_wf_bf = (wf*w_gNB)*hest_gU_all(i,:)';
				UE_snr_default(rowidx,i1,i2,var,i) = snr(UE_wf_bf,noise_UE(:,i));
				if UE_snr_default(rowidx,i1,i2,var,i) > -1
					sumrate_default(rowidx,i1,i2,var) = sumrate_default(rowidx,i1,i2,var) + log2(1+UE_snr_default(rowidx,i1,i2,var,i));
				else
					sumrate_default(rowidx,i1,i2,var) = 0;
					break;
				end
				UE_wf_bf = (wf*w_copt.')*hest_gU_all(i,:)';
				UE_snr_bf(rowidx,i1,i2,var,i) = snr(UE_wf_bf,noise_UE(:,i));
				if UE_snr_bf(rowidx,i1,i2,var,i) > -1
					sumrate_bf(rowidx,i1,i2,var) = sumrate_bf(rowidx,i1,i2,var) + log2(1+UE_snr_bf(rowidx,i1,i2,var,i));
				else
					sumrate_bf(rowidx,i1,i2,var) = 0;
					break;
				end
			end
			RAS_null_default(rowidx,i1,i2,var) = 10*log10(norm(hest_gR_temp*w_gNB').^2);
			RAS_null_bf(rowidx,i1,i2,var) = 10*log10(norm(hest_gR_temp*w_copt).^2);

% 			wf_bf = wf*w_copt.';
% 			if mean(isnan(wf_bf),'all')
% 				nanIdx = isnan(wf_bf);
% 				wf_bf(nanIdx) = 0;
% 			end
% 			[RAS_wf_bf,~,~] = RAS_channel(wf_bf);
% 			RAS_snr_bf(rowidx,i1,i2,var) = snr(RAS_wf_bf,noise_RAS);
% 			RAS_snr_gain(rowidx,i1,i2,var) = RAS_snr_bf(rowidx,i1,i2,var) - RAS_snr_default(rowidx,i1,i2,var);
% 
% 			% SNR at K UEs
% 			for i = 1:K
% 				UE_wf_bf = wf_bf*hest_gU_all(i,:)';
% 				UE_snr_bf(rowidx,i1,i2,var) = UE_snr_bf(rowidx,i1,i2,var) + snr(UE_wf_bf,noise_UE(:,i));
% 			end
% 			UE_snr_bf(rowidx,i1,i2,var) = UE_snr_bf(rowidx,i1,i2,var)/K;
% 			UE_snr_gain(rowidx,i1,i2,var) = (UE_snr_bf(rowidx,i1,i2,var) - UE_snr_default(rowidx,i1,i2,var))./5;
		end
	end
end
end

%% Plots

%{
figure;
plot(gNBDist,permute(mean(RAS_null_default(1,:,1:3,:),4),[3 2 1]))
hold on;
plot(gNBDist,permute(mean(RAS_null_bf(1,:,1:3,:),4),[3 2 1]))
xlabel("gNB Distance (km)");
ylabel("Interference Power (dB)");
semiPlots;
leg = legend("M=4","M=8","M=16","location","best");
leg.ItemTokenSize = [10 18];
leg.Orientation = 'horizontal';
ylim([-400 0]);
yticks([-400 -300 -200 -100 0]);

figure;
plot(gNBDist,permute(mean(RAS_null_default(2,:,1:3,:),4),[3 2 1]))
hold on;
plot(gNBDist,permute(mean(RAS_null_bf(2,:,1:3,:),4),[3 2 1]))
xlabel("gNB Distance (km)");
ylabel("Interference Power (dB)");
semiPlots;
leg = legend("fc=0.7","fc=1.8","fc=5.9","location","best");
leg.ItemTokenSize = [10 18];
leg.Orientation = 'horizontal';
ylim([-400 0]);
yticks([-400 -300 -200 -100 0]);

figure;
plot(gNBDist,permute(mean(RAS_null_default(3,:,1:3,:),4),[3 2 1]))
hold on;
plot(gNBDist,permute(mean(RAS_null_bf(3,:,1:3,:),4),[3 2 1]))
xlabel("gNB Distance (km)");
ylabel("Interference Power (dB)");
semiPlots;
leg = legend("K=1","K=8","K=16","location","best");
leg.ItemTokenSize = [10 18];
leg.Orientation = 'horizontal';
ylim([-400 0]);
yticks([-400 -300 -200 -100 0]);

figure;
plot(gNBDist,permute(mean(sumrate_default(1,:,1:3,:),4),[3 2 1]))
hold on;
plot(gNBDist,permute(mean(sumrate_bf(1,:,1:3,:),4),[3 2 1]))
xlabel("gNB Distance (km)");
ylabel("Sum Rate (bits/s/Hz)");
semiPlots;
leg = legend("M=4","M=8","M=16","location","best");
leg.ItemTokenSize = [10 18];
leg.Orientation = 'horizontal';
ylim([0 100]);
yticks([0 25 50 75 100]);

figure;
plot(gNBDist,permute(mean(sumrate_default(2,:,1:3,:),4),[3 2 1]))
hold on;
plot(gNBDist,permute(mean(sumrate_bf(2,:,1:3,:),4),[3 2 1]))
xlabel("gNB Distance (km)");
ylabel("Sum Rate (bits/s/Hz)");
semiPlots;
leg = legend("fc=0.7","fc=1.8","fc=5.9","location","best");
leg.ItemTokenSize = [10 18];
leg.Orientation = 'horizontal';
ylim([0 100]);
yticks([0 25 50 75 100]);

figure;
plot(gNBDist,permute(mean(sumrate_default(3,:,1:3,:),4),[3 2 1]))
hold on;
plot(gNBDist,permute(mean(sumrate_bf(3,:,1:3,:),4),[3 2 1]))
xlabel("gNB Distance (km)");
ylabel("Sum Rate (bits/s/Hz)");
semiPlots;
leg = legend("K=1","K=8","K=16","location","best");
leg.ItemTokenSize = [10 18];
leg.Orientation = 'horizontal';
ylim([0 100]);
yticks([0 25 50 75 100]);

figure;
plot(gNBDist,permute(mean(latency(1,:,1:3,:),4),[3 2 1]))
xlabel("gNB Distance (km)");
ylabel("Latency (ms)");
semiPlots;
leg = legend("M=4","M=8","M=16","location","best");
leg.ItemTokenSize = [10 18];
leg.Orientation = 'horizontal';
ylim([0 1]);
yticks([0.2 0.4 0.6 0.8]);

figure;
plot(gNBDist,permute(mean(latency(2,:,1:3,:),4),[3 2 1]))
xlabel("gNB Distance (km)");
ylabel("Latency (ms)");
semiPlots;
leg = legend("fc=0.7","fc=1.8","fc=5.9","location","best");
leg.ItemTokenSize = [10 18];
leg.Orientation = 'horizontal';
ylim([0 1]);
yticks([0.2 0.4 0.6 0.8]);
%}

figure;
plot(gNBDist,permute(mean(latency(3,:,1:3,:),4),[3 2 1]))
xlabel("gNB Distance (km)");
ylabel("Latency (ms)");
semiPlots;
leg = legend("K=1","K=8","K=16","location","best");
leg.ItemTokenSize = [10 18];
leg.Orientation = 'horizontal';
ylim([0.1 0.8]);
yticks([0.2 0.4 0.6 0.8]);