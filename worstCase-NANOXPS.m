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
% NtestUEs = 10;
K = 8;
% N = 8;
UE_snr_default = zeros(1,K);
UE_snr_bf = zeros(1,K);
% latency = zeros(size(UE_snr_default));

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

testUENames = "testUE"+string(1:K);
testUEs = rxsite("Name",testUENames,"Latitude",RASPos(1),"Longitude",RASPos(2),"AntennaAngle",UEAntDir(1:2),"AntennaHeight",5);
xLocs = [1 2 2 1 -1 -2 -2 -1];
yLocs = [2 1 -1 -2 -2 -1 1 2];

hest_gU_all = zeros(K,RAS_Nt);
noise_UE = zeros(size(RAS_rxWaveform,1),K);
for i=1:K
	testUEs(i).Latitude = testUEs(i).Latitude + (yLocs(i)/earthRadius)*(180/pi);
	testUEs(i).Longitude = testUEs(i).Longitude + (xLocs(i)/earthRadius)*(180/pi);
	UE_rays = raytrace(gNBSite,testUEs(i),pm,"Type","pathloss");
	pathToAs = [UE_rays{1}.PropagationDelay] - min([UE_rays{1}.PropagationDelay]);
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
	noise_UE(:,i) = wgn(size(UE_rxWaveform,1),size(UE_rxWaveform,2),-59);
end

[U,D,V] = svd(hest_gU_all);
w_gNB = V(:,1:nLayers).';
P = null(hest_gR_temp);
w_copt = P*P'*w_gNB.';

wf = complex(randn(T,1),randn(T,1));
sumrate_default = 0;
sumrate_bf = 0;
for i = 1:K
	UE_wf_bf = (wf*w_gNB)*hest_gU_all(i,:)';
	UE_snr_default(i) = snr(UE_wf_bf,noise_UE(:,i));
	if UE_snr_default(i) > -1
		sumrate_default = sumrate_default + log2(1+UE_snr_default(i));
	else
		sumrate_default = 0;
		break;
	end
	UE_wf_bf = (wf*w_copt.')*hest_gU_all(i,:)';
	UE_snr_bf(i) = snr(UE_wf_bf,noise_UE(:,i));
	if UE_snr_bf(i) > -1
		sumrate_bf = sumrate_bf + log2(1+UE_snr_bf(i));
	else
		sumrate_bf = 0;
		break;
	end
end
RAS_null_default = 10*log10(norm(hest_gR_temp*w_gNB').^2)
RAS_null_bf = 10*log10(norm(hest_gR_temp*w_copt).^2)
sumrate_default
sumrate_bf

% Display beampattern of gNB
show(gNBSite); 
show(RASSite);
show(testUEs);
gNBSite.Antenna = clone(RAS_channel.TransmitAntennaArray);
gNBSite.Antenna.Taper = w_copt;
pattern(gNBSite,fc,"Size",250);

%{
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
				pathToAs = [];
				while(isempty(pathToAs))
					locDisp = randn(1,2).*1e-3;
					testUE.Latitude = locDisp(1) + gNBPos(1);
					testUE.Longitude = locDisp(2) + gNBPos(2);
					UE_rays = raytrace(gNBSite,testUE,pm,"Type","pathloss");
					pathToAs = [UE_rays{1}.PropagationDelay] - min([UE_rays{1}.PropagationDelay]);
				end
				testUE.AntennaHeight = 2.^i1;
			case 3
				testUE.AntennaHeight = 30;
				locDisp = randn(1,2).*1e-3;
				testUE.Latitude = gNBPos(1) + locDisp(1) + ((i1*10)/earthRadius)*(180/pi);
				testUE.Longitude = gNBPos(2) + locDisp(2);
		end
		UE_rays = raytrace(gNBSite,testUE,pm,"Type","pathloss");
		UE_channel = getChannelObj(UE_rays,fc,ofdmInfo.SampleRate);

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
		%UE_rxWaveform = txWaveform*testhest_gU';
        testnoise_UE = wgn(size(UE_rxWaveform,1),size(UE_rxWaveform,2),-89);
		UE_snr_default(rowidx,i1) = snr(UE_rxWaveform,testnoise_UE);
        [testw_gNB, ~, ~] = getBeamformingWeights(hest_gU, nLayers, scOffset, noRBs);

		for var = 1:N
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
			end
			tic
			hest_gU_all(K,:) = testhest_gU;
			hest_gU_all = hest_gU_all./norm(hest_gU_all);
			[U,D,V] = svd(hest_gU_all);
			w_gNB = V(:,1:nLayers).';
			P = null(hest_gR_temp);
			w_copt = P*P'*w_gNB.';
			latency(rowidx,i1) = latency(rowidx,i1) + toc;
			
			wf = complex(randn(T,1),randn(T,1));
			wf_bf = wf*w_copt.';
			if mean(isnan(wf_bf),'all')
				nanIdx = isnan(wf_bf);
				wf_bf(nanIdx) = 0;
			end
			[RAS_wf_bf,~,~] = RAS_channel(wf_bf);
			RAS_snr_bf(rowidx,i1) = RAS_snr_bf(rowidx,i1) + snr(RAS_wf_bf,noise_RAS);
			UE_wf_bf = wf_bf*hest_gU_all(K,:)';
			UE_snr_bf(rowidx,i1) = UE_snr_bf(rowidx,i1) + snr(UE_wf_bf,testnoise_UE);
		end
		RAS_snr_bf(rowidx,i1) = RAS_snr_bf(rowidx,i1)/N;
		UE_snr_bf(rowidx,i1) = UE_snr_bf(rowidx,i1)/N;
		latency(rowidx,i1) = latency(rowidx,i1)/N;
	end
end
%}