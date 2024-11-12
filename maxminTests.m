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
gNBDist = distance(gNBLocs(1,1),gNBLocs(1,2),RASPos(1),RASPos(2));
gNBDist = deg2rad(gNBDist).*(earthRadius*1e-3);
gNBPos = gNBLocs;
gNBAntSize = [8 2];
M = prod(gNBAntSize);
fc = 1880e6;
reflectionsOrder = 3;
% NtestUEs = 10;
K = 8;
% N = 8;
UE_snr_default = zeros(1,K);
UE_snr_heu = zeros(1,K);
UE_snr_opt = zeros(1,K);

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

h_k = hest_gU_all.';
H = zeros(M,M,K);
for i = 1:K
	H(:,:,i) = h_k(:,i)*h_k(:,i)';
end
hgR = hest_gR_temp.';
HgR = hgR*hgR';
cvx_begin
	variable tau nonnegative
	variable X(M,M) hermitian
	variable s(K,1)
	minimize(-tau)
	subject to
		norm(diag(X*HgR)) <= e
		for i = 1:K
			trace(X*H(:,:,i)) - s(i) == tau
			s(i) >= 0
		end
		trace(X) <= 1
		X == hermitian_semidefinite(M)
cvx_end

L = 10000;
x = zeros(M,L);
for i = 1:L
	x(:,i) = randn(M,1,'like',1i);
	[U,D,~] = eig(X);
	x(:,i) = U*sqrt(D)*x(:,i);
	x(:,i) = x(:,i)./norm(x(:,i));
	scale = e./(norm(x(:,i)'*hgR).^2);
	x(:,i) = x(:,i).*(scale.^(0.5));
end
[~,idx] = max(min(abs(x'*h_k),[],2));
w_opt = x(:,idx).';

wf = complex(randn(T,1),randn(T,1));
sumrate_default = 0;
sumrate_heu = 0;
sumrate_opt = 0;

for i = 1:K
	UE_wf_bf = (wf*w_gNB)*hest_gU_all(i,:)';
	UE_snr_default(i) = snr(UE_wf_bf,noise_UE(:,i));
	UE_wf_bf = (wf*w_copt.')*hest_gU_all(i,:)';
	UE_snr_heu(i) = snr(UE_wf_bf,noise_UE(:,i));
	UE_wf_bf = (wf*(w_opt./(scale.^0.25)))*hest_gU_all(i,:)';
	UE_snr_opt(i) = snr(UE_wf_bf,noise_UE(:,i));
end
if mean(UE_snr_default(i)) > -1
	sumrate_default = sum(log2(1+UE_snr_default(i)));
else
	sumrate_default = sum(log2(1+10.^(UE_snr_default(i)./10)));
end
if mean(UE_snr_heu(i)) > -1
	sumrate_heu = sum(log2(1+UE_snr_heu(i)));
else
	sumrate_heu = sum(log2(1+10.^(UE_snr_heu(i)./10)));
end
if mean(UE_snr_opt(i)) > -1
	sumrate_opt = sum(log2(1+UE_snr_opt(i)));
else
	sumrate_opt = sum(log2(1+10.^(UE_snr_opt(i)./10)));
end
RAS_null_default = 10*log10(norm(hest_gR_temp*w_gNB').^2)
RAS_null_heu = 10*log10(norm(hest_gR_temp*w_copt).^2)
RAS_null_opt = 10*log10(norm(hest_gR_temp*w_opt').^2)
sumrate_default
sumrate_heu
sumrate_opt

%% Evaluate REM around RAS
X = -200:1:200;	
Y = -200:1:200;
remPow_gNB = zeros(length(X),length(Y));
remPow_heu = zeros(size(remPow_gNB));
remPow_opt = zeros(size(remPow_gNB));
for i=1:length(X)
	parfor j=1:length(Y)
		iterationID = sprintf("Iteration %d %d...\n",i,j);
		disp(iterationID);
		remLoc = rxsite("Name","remLoc","Latitude",RASPos(1),"Longitude",RASPos(2),"AntennaAngle",UEAntDir(1:2),"AntennaHeight",5);
		remLoc.Latitude = remLoc.Latitude + (Y(j)/earthRadius)*(180/pi);
		remLoc.Longitude = remLoc.Longitude + (X(i)/earthRadius)*(180/pi);
		UE_rays = raytrace(gNBSite,remLoc,pm,"Type","pathloss");
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
		rem_hest = getChannelCoeffs(hest_gU, scOffset, noRBs);
		rem_hest = rem_hest./norm(rem_hest);
		if(X(i) == 0 && Y(j) == 0)
			rem_hest = hest_gR_temp;
		end
		remPow_gNB(i,j) = 10*log10(norm(rem_hest*w_gNB').^2);
		remPow_heu(i,j) = 10*log10(norm(rem_hest*w_copt).^2);
		remPow_opt(i,j) = 10*log10(norm(rem_hest*w_opt').^2);
	end
end

%% Plot REM

% Default beamformer
figure;
s = surface(X,Y,remPow_gNB);
s.EdgeColor = 'none';
colorbar;
ylabel('Latitude (m)');
xlabel('Longitude (m)');
hold on
plot(0,0,'or');
text(0,10,'RAS');
text(41,-116,'gNB');
plot(41,-106,'ob');

% Heuristic Solution
figure;
s = surface(X,Y,remPow_heu);
s.EdgeColor = 'none';
colorbar;
ylabel('Latitude (m)');
xlabel('Longitude (m)');
hold on
plot(0,0,'or');
text(0,10,'RAS');
text(41,-116,'gNB');
plot(41,-106,'ob');

% Max-min solution
figure;
s = surface(X,Y,remPow_opt);
s.EdgeColor = 'none';
colorbar;
ylabel('Latitude (m)');
xlabel('Longitude (m)');
hold on
plot(0,0,'or');
text(0,10,'RAS');
text(41,-116,'gNB');
plot(41,-106,'ob');