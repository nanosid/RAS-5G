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
gNBLocs = [42.930223, -77.501120];
gNBDist = distance(gNBLocs(1,1),gNBLocs(1,2),RASPos(1),RASPos(2));
gNBDist = deg2rad(gNBDist).*(earthRadius*1e-3);
gNBPos = gNBLocs;
gNBAntSize = [8 2];
fc = 1880e6;
reflectionsOrder = 3;
NtestUEs = 10;
K = 8;
N = 128;
UE_snr_default = zeros(3,NtestUEs,N);
UE_snr_bf = zeros(size(UE_snr_default));
RAS_snr_default = zeros(size(UE_snr_default));
RAS_snr_bf = zeros(size(UE_snr_default));
latency = zeros(size(UE_snr_default));

gNBSite = txsite("Name","Victor_gNB","Latitude",gNBPos(1),"Longitude",gNBPos(2),"AntennaAngle",...
	gNBAntDir(1:2),"AntennaHeight",32,"TransmitterFrequency",fc);

% Define Raytracing prop model and display plotted rays
pm = propagationModel("raytracing","Method","sbr","MaxNumReflections",reflectionsOrder);
RAS_rays = raytrace(gNBSite,RASSite,pm,"Type","pathloss");

% Obtain path gains and directions
ofdmInfo = nrOFDMInfo(NRB,SCS);
RAS_channel = getChannelObj(RAS_rays,fc,ofdmInfo.SampleRate);

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
noise_RAS = wgn(size(RAS_rxWaveform,1),size(RAS_rxWaveform,2),0);
RAS_snr_default(:,:) = snr(RAS_rxWaveform,noise_RAS);

% Get gNB-RAS channel coefficients over all RBs and OFDM symbols
RAS_pathFilters = getPathFilters(RAS_channel);
[RAS_offset,~] = nrPerfectTimingEstimate(RAS_pathGains, RAS_pathFilters);
hest_gR = nrPerfectChannelEstimate(RAS_pathGains, RAS_pathFilters, NRB, SCS, nSlots, RAS_offset, RAS_sampleTimes);
hest_gR_temp = permute(mean(reshape(hest_gR,[],RAS_Nr,RAS_Nt)),[2,3,1]);
hest_gR_temp = hest_gR_temp./norm(hest_gR_temp);

testUE = rxsite("Name","testUE","Latitude",0,"Longitude",0,"AntennaAngle",UEAntDir(1:2),"AntennaHeight",2);

pathToAs = [];
while(isempty(pathToAs))	
	locDisp = randn(1,2).*1e-3;
	testUE.Latitude = locDisp(1) + gNBPos(1);
	testUE.Longitude = locDisp(2) + gNBPos(2);
	test_UE_rays = raytrace(gNBSite,testUE,pm,"Type","pathloss");
	pathToAs = [test_UE_rays{1}.PropagationDelay] - min([test_UE_rays{1}.PropagationDelay]);
end
testLat = testUE.Latitude;
testLong = testUE.Longitude;

for rowidx = 1:3
	for i1 = 1:NtestUEs
		iterationID = sprintf("Iteration %d,%d...\n",rowidx,i1);
		disp(iterationID);
		switch rowidx
			case 1
				testUE.Latitude = RASPos(1);
				testUE.Longitude = RASPos(2);
				testUE.AntennaHeight = 2.^i1;
			case 2
				testUE.Latitude = testLat;
				testUE.Longitude = testLong;
				testUE.AntennaHeight = 2.^i1;
			case 3
				testUE.Latitude = testLat + ((i1*10)/earthRadius)*(180/pi);
				testUE.Longitude = testLong;
				testUE.AntennaHeight = 30;
		end
		test_UE_rays = raytrace(gNBSite,testUE,pm,"Type","pathloss");
		UE_channel = getChannelObj(test_UE_rays,fc,ofdmInfo.SampleRate);

		UEArray = phased.NRRectangularPanelArray('Size',[UEAntSize(1:2) 1 1],'Spacing',...
			[0.5*lambda*[1 1] 1 1]);
		UEArray.ElementSet = {phased.IsotropicAntennaElement};   % isotropic antenna element
		UE_channel.ReceiveAntennaArray = UEArray;
		UE_channel.ReceiveArrayOrientation = [UEAntDir(1); (-1)*UEAntDir(2); 0];
		UE_channel.TransmitAntennaArray = gNBArray;
		UE_channel.TransmitArrayOrientation = [gNBAntDir(1); (-1)*gNBAntDir(2); 0];
		UE_channelInfo = info(UE_channel);
		UE_Nt = UE_channelInfo.NumTransmitAntennas;
		UE_Nr = UE_channelInfo.NumReceiveAntennas;
		[UE_rxWaveform,UE_pathGains,UE_sampleTimes] = UE_channel(txWaveform);
		UE_pathFilters = getPathFilters(UE_channel);
		[UE_offset,~] = nrPerfectTimingEstimate(UE_pathGains, UE_pathFilters);
		hest_gU = nrPerfectChannelEstimate(UE_pathGains, UE_pathFilters, NRB, SCS, nSlots,...
			UE_offset, UE_sampleTimes);
		testhest_gU = getChannelCoeffs(hest_gU, scOffset, noRBs);
		testhest_gU = testhest_gU./norm(testhest_gU);
		%UE_rxWaveform = txWaveform*testhest_gU';
        testnoise_UE = wgn(size(UE_rxWaveform,1),size(UE_rxWaveform,2),-39);
		%UE_snr_default(rowidx,i1) = snr(UE_rxWaveform,testnoise_UE);
        [testw_gNB, ~, ~] = getBeamformingWeights(hest_gU, nLayers, scOffset, noRBs);

		parfor var = 1:N
			hest_gU_all = zeros(K,RAS_Nt);
			w_gNB_all = zeros(K,RAS_Nt);
			noise_UE = zeros(size(RAS_rxWaveform,1),K);
			UENames = "UE"+string(1:K-1);
			UESites = rxsite("Name",UENames,"Latitude",0,"Longitude",0,"AntennaAngle",UEAntDir(1:2),"AntennaHeight",2);
			for i = 1:K-1
				pathToAs = [];
				while(isempty(pathToAs))
					locDisp = randn(1,2).*1e-3;
					UESites(i).Latitude = locDisp(1) + gNBPos(1);
					UESites(i).Longitude = locDisp(2) + gNBPos(2);
					UE_rays = raytrace(gNBSite,UESites(i),pm,"Type","pathloss");
					pathToAs = [UE_rays{1}.PropagationDelay] - min([UE_rays{1}.PropagationDelay]);
				end
				UE_channel = getChannelObj(UE_rays,fc,ofdmInfo.SampleRate);				
				UEArray = phased.NRRectangularPanelArray('Size',[UEAntSize(1:2) 1 1],'Spacing',[0.5*lambda*[1 1] 1 1]);
				UEArray.ElementSet = {phased.IsotropicAntennaElement};   % isotropic antenna element
				UE_channel.ReceiveAntennaArray = UEArray;
				UE_channel.ReceiveArrayOrientation = [UEAntDir(1); (-1)*UEAntDir(2); 0];
				UE_channel.TransmitAntennaArray = gNBArray;
				UE_channel.TransmitArrayOrientation = [gNBAntDir(1); (-1)*gNBAntDir(2); 0];
				UE_channelInfo = info(UE_channel);
				UE_Nt = UE_channelInfo.NumTransmitAntennas;
				UE_Nr = UE_channelInfo.NumReceiveAntennas;
				[UE_rxWaveform,UE_pathGains,UE_sampleTimes] = UE_channel(txWaveform);
				UE_pathFilters = getPathFilters(UE_channel);
				[UE_offset,~] = nrPerfectTimingEstimate(UE_pathGains, UE_pathFilters);
				hest_gU = nrPerfectChannelEstimate(UE_pathGains, UE_pathFilters, NRB, SCS, nSlots,...
				UE_offset, UE_sampleTimes);
				hest_gU_all(i,:) = getChannelCoeffs(hest_gU, scOffset, noRBs);
				hest_gU_all(i,:) = hest_gU_all(i,:)./norm(hest_gU_all(i,:));
			end
			tic
			hest_gU_all(K,:) = testhest_gU;
			[U,D,V] = svd(hest_gU_all);
			w_gNB = V(:,1:nLayers).';
			P = null(hest_gR_temp);
			w_copt = P*P'*w_gNB.';
			latency(rowidx,i1,var) = toc;
			
			wf = complex(randn(T,1),randn(T,1));
			wf_bf = wf*w_copt.';
			if mean(isnan(wf_bf),'all')
				nanIdx = isnan(wf_bf);
				wf_bf(nanIdx) = 0;
			end
% 			[RAS_wf_bf,~,~] = RAS_channel(wf_bf);
% 			RAS_snr_bf(rowidx,i1) = RAS_snr_bf(rowidx,i1) + snr(RAS_wf_bf,noise_RAS);
% 			UE_wf_bf = wf_bf*hest_gU_all(K,:)';
% 			UE_snr_bf(rowidx,i1) = UE_snr_bf(rowidx,i1) + snr(UE_wf_bf,testnoise_UE);
% 			UE_wf_bf = (wf*w_gNB)*hest_gU_all(K,:)';
% 			UE_snr_default(rowidx,i1) = UE_snr_default(rowidx,i1) + snr(UE_wf_bf,testnoise_UE);
			[RAS_wf_bf,~,~] = RAS_channel(wf_bf);
			RAS_snr_bf(rowidx,i1,var) = snr(RAS_wf_bf,noise_RAS);
			UE_wf_bf = wf_bf*hest_gU_all(K,:)';
			UE_snr_bf(rowidx,i1,var) = snr(UE_wf_bf,testnoise_UE);
			UE_wf_bf = (wf*w_gNB)*hest_gU_all(K,:)';
			UE_snr_default(rowidx,i1,var) = snr(UE_wf_bf,testnoise_UE);
		end
% 		RAS_snr_bf(rowidx,i1) = RAS_snr_bf(rowidx,i1)/N;
% 		UE_snr_default(rowidx,i1) = UE_snr_default(rowidx,i1)/N;
% 		UE_snr_bf(rowidx,i1) = UE_snr_bf(rowidx,i1)/N;
% 		latency(rowidx,i1) = latency(rowidx,i1)/N;
	end
end

%% Plots
plot(1:10,mean(UE_snr_bf,3));
hold on;
plot(1:10,mean(UE_snr_default,3));
semiPlots
xticks(1:10);
xlim([1 10]);
xlabel("UE Position Index");
% ylim([-10 30]);
% yticks(-10:10:30);
ylabel("Moving UE SNR (dB)");
leg = legend("Case 1","Case 2","Case 3","location","best");
leg.ItemTokenSize = [10 18];
leg.Orientation = 'horizontal';