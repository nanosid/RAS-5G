clear;
% Define constant parameters
c = physconst("LightSpeed");
e = 10^(-220/10);
gNBAntDir = [0 0].';
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

gNBAntSize = [8 2];
M = prod(gNBAntSize);
fc = 1880e6;
reflectionsOrder = 1;
K = 8;
lambda = c/fc;
L1 = 1400;
ofdmInfo = nrOFDMInfo(NRB,SCS);

% viewer = siteviewer("Basemap","openstreetmap","Buildings","RIT_kimballLoop.osm");
gNBPos = [43.083695, -77.680426];
UEPos = [43.083376, -77.680509];
RASPos = [43.083378, -77.680697];

gNBSite = txsite("Name","GCI_gNB","Latitude",gNBPos(1),"Longitude",gNBPos(2),"AntennaAngle",gNBAntDir(1:2),"AntennaHeight",10,"TransmitterFrequency",fc);
RASSite = rxsite("Name","GCI_RAS","Latitude",RASPos(1),"Longitude",RASPos(2),"AntennaAngle",RASAntDir(1:2),"AntennaHeight",1);
UESite = rxsite("Name","GCI_UE","Latitude",UEPos(1),"Longitude",UEPos(2),"AntennaAngle",UEAntDir(1:2),"AntennaHeight",1);

pm = propagationModel("raytracing","Method","sbr","MaxNumReflections",1);
T = ofdmInfo.SampleRate * 1e-3;
txWaveform = complex(randn(T,M),randn(T,M));
if mean(isnan(txWaveform),'all')
	nanIdx = isnan(txWaveform);
	txWaveform(nanIdx) = 0;
end

RAS_rays = raytrace(gNBSite,RASSite,pm,"Type","pathloss");
RAS_channel = getChannelObj(RAS_rays,fc,ofdmInfo.SampleRate);
RASArray = phased.NRRectangularPanelArray('Size',[RASAntSize(1:2) 1 1],'Spacing', [0.5*lambda*[1 1] 1 1]);
RASArray.ElementSet = {phased.IsotropicAntennaElement};   % isotropic antenna element
RAS_channel.ReceiveAntennaArray = RASArray;
RAS_channel.ReceiveArrayOrientation = [RASAntDir(1); (-1)*RASAntDir(2); 0]; 
gNBArray = phased.NRRectangularPanelArray('Size',[gNBAntSize(1:2) 1 1],'Spacing', [0.5*lambda*[1 1] 1 1]);
gNBArray.ElementSet = {phased.IsotropicAntennaElement};
RAS_channel.TransmitAntennaArray = gNBArray;
RAS_channel.TransmitArrayOrientation = [gNBAntDir(1); (-1)*gNBAntDir(2); 0];

[RAS_rxWaveform,RAS_pathGains,RAS_sampleTimes] = RAS_channel(txWaveform);
RAS_pathFilters = getPathFilters(RAS_channel);
[RAS_offset,~] = nrPerfectTimingEstimate(RAS_pathGains, RAS_pathFilters);
hest_gR = nrPerfectChannelEstimate(RAS_pathGains, RAS_pathFilters, NRB, SCS, nSlots, RAS_offset, RAS_sampleTimes);
hest_gR = getChannelCoeffs(hest_gR, scOffset, noRBs);

%% Erroneous channels
hest_err = hest_gR;
parfor i1=1:L1
	hest_err = hest_err + wgn(size(hest_gR,1),size(hest_gR,2),-350,'complex');
end
hest_gR_temp = hest_err;

%% Channels of a volume
%{
width = -0.5:0.1:0.5;
L2 = length(width);
hest_temp_all = zeros(L2^3,M);
for i1 = 1:L2
	for i2 = 1:L2
		for i3 = 1:L2
			iterationID = sprintf("Iteration %d %d %d...\n",i1,i2,i3);
			disp(iterationID);
			tempSite = rxsite("Name","GCI_RAS_temp","Latitude",RASPos(1),"Longitude",RASPos(2),"AntennaAngle",RASAntDir(1:2),"AntennaHeight",1);
			tempSite.Latitude = tempSite.Latitude + (width(i1)/earthRadius)*(180/pi);
			tempSite.Longitude = tempSite.Longitude + (width(i2)/earthRadius)*(180/pi);
			tempSite.AntennaHeight = tempSite.AntennaHeight + width(i3);
			temp_rays = raytrace(gNBSite,tempSite,pm,"Type","pathloss");
			temp_channel = getChannelObj(temp_rays,fc,ofdmInfo.SampleRate);

			tempArray = phased.NRRectangularPanelArray('Size',[RASAntSize(1:2) 1 1],'Spacing', [0.5*lambda*[1 1] 1 1]);
			tempArray.ElementSet = {phased.IsotropicAntennaElement};   % isotropic antenna element
			temp_channel.ReceiveAntennaArray = tempArray;
			temp_channel.ReceiveArrayOrientation = [RASAntDir(1); (-1)*RASAntDir(2); 0]; 
			
			gNBArray = phased.NRRectangularPanelArray('Size',[gNBAntSize(1:2) 1 1],'Spacing', [0.5*lambda*[1 1] 1 1]);
			gNBArray.ElementSet = {phased.IsotropicAntennaElement};
			temp_channel.TransmitAntennaArray = gNBArray;
			temp_channel.TransmitArrayOrientation = [gNBAntDir(1); (-1)*gNBAntDir(2); 0];

			[~,temp_pathGains,temp_sampleTimes] = temp_channel(txWaveform);
			temp_pathFilters = getPathFilters(temp_channel);
			[temp_offset,~] = nrPerfectTimingEstimate(temp_pathGains, temp_pathFilters);
			hest_temp = nrPerfectChannelEstimate(temp_pathGains, temp_pathFilters, NRB, SCS, nSlots, temp_offset, temp_sampleTimes);
			hest_temp_all((i1-1)*L2^2 + (i2-1)*L2 + i3,:) = getChannelCoeffs(hest_temp, scOffset, noRBs);
		end
	end
end
%}
%% Mean and Variance of raytraced channels
%{
hest_test_all = zeros(L1,M);
parfor i1 = 1:L1
	test_rays = raytrace(gNBSite,RASSite,pm,"Type","pathloss");
	test_channel = getChannelObj(test_rays,fc,ofdmInfo.SampleRate);				
	testArray = phased.NRRectangularPanelArray('Size',[UEAntSize(1:2) 1 1],'Spacing',[0.5*lambda*[1 1] 1 1]);
	testArray.ElementSet = {phased.IsotropicAntennaElement};   % isotropic antenna element
	gNBArray = phased.NRRectangularPanelArray('Size',[gNBAntSize(1:2) 1 1],'Spacing', [0.5*lambda*[1 1] 1 1]);
	gNBArray.ElementSet = {phased.IsotropicAntennaElement};
	test_channel.ReceiveAntennaArray = testArray;
	test_channel.ReceiveArrayOrientation = [UEAntDir(1); (-1)*UEAntDir(2); 0];
	test_channel.TransmitAntennaArray = gNBArray;
	test_channel.TransmitArrayOrientation = [gNBAntDir(1); (-1)*gNBAntDir(2); 0];
	txWaveform = complex(randn(T,M),randn(T,M));
	[~,test_pathGains,test_sampleTimes] = test_channel(txWaveform);
	test_pathFilters = getPathFilters(test_channel);
	[test_offset,~] = nrPerfectTimingEstimate(test_pathGains, test_pathFilters);
	hest_test = nrPerfectChannelEstimate(test_pathGains, test_pathFilters, NRB, SCS, nSlots,...
	test_offset, test_sampleTimes);
	hest_test_all(i1,:) = getChannelCoeffs(hest_test, scOffset, noRBs);
end

mu = mean(real(hest_test_all),1) + mean(imag(hest_test_all),1)*1i;
sigma2 = var(real(hest_test_all),0,1) + var(imag(hest_test_all),0,1)*1i;
%}

%% Set up UEs

UENames = "UE"+string(1:K);
UESites = rxsite("Name",UENames,"Latitude",UEPos(1),"Longitude",UEPos(2),"AntennaAngle",UEAntDir(1:2),"AntennaHeight",2);

% Now derive channel for K UEs and the beamforming weights for each UE
hest_gU_all = zeros(K,M);
noise_UE = zeros(size(RAS_rxWaveform,1),K);

UE_channel = nrCDLChannel;
for i=1:K
	pathToAs = [];
	while(isempty(pathToAs))
		locDisp = abs(randn(1,2).*5e-3);
		UESites(i).Latitude = UEPos(1) + locDisp(1);
		UESites(i).Longitude = UEPos(2) + locDisp(2);
		UE_rays = raytrace(gNBSite,UESites(i),pm,"Type","pathloss");
		pathToAs = [UE_rays{1}.PropagationDelay] - min([UE_rays{1}.PropagationDelay]);
	end
	UE_channel = getChannelObj(UE_rays,fc,ofdmInfo.SampleRate);
	UEArray = phased.NRRectangularPanelArray('Size',[UEAntSize(1:2) 1 1],'Spacing',...
		[0.5*lambda*[1 1] 1 1]);
	UEArray.ElementSet = {phased.IsotropicAntennaElement};   % isotropic antenna element
	UE_channel.ReceiveAntennaArray = UEArray;
	UE_channel.ReceiveArrayOrientation = [UEAntDir(1); (-1)*UEAntDir(2); 0];
	UE_channel.TransmitAntennaArray = gNBArray;
	UE_channel.TransmitArrayOrientation = [gNBAntDir(1); (-1)*gNBAntDir(2); 0];
	
	UE_channelInfo = info(UE_channel);
	[UE_rxWaveform,UE_pathGains,UE_sampleTimes] = UE_channel(txWaveform);
	UE_pathFilters = getPathFilters(UE_channel);
	[UE_offset,~] = nrPerfectTimingEstimate(UE_pathGains, UE_pathFilters);
	hest_gU = nrPerfectChannelEstimate(UE_pathGains, UE_pathFilters, NRB, SCS, nSlots,...
		UE_offset, UE_sampleTimes);
	hest_gU_all(i,:) = getChannelCoeffs(hest_gU, scOffset, noRBs);
	% hest_gU_all(i,:) = hest_gU_all(i,:)./norm(hest_gU_all(i,:));
    noise_UE(:,i) = wgn(size(UE_rxWaveform,1),size(UE_rxWaveform,2),-89);
end

%% Beamforming
[~,~,V] = svd(hest_gU_all);
w_gNB = V(:,1:nLayers).';
P = null(hest_gR_temp);
w_heu = P*P'*w_gNB.';

intPower = 10*log10(norm(hest_gR*w_heu).^2)