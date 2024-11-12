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

gNBAntSize = [1 2];
M = prod(gNBAntSize);
fc = 5350e6;
reflectionsOrder = 1;
K = 1;

viewer = siteviewer("Basemap","openstreetmap","Buildings","RIT_kimballLoop.osm");
gNBPos = [43.083695, -77.680426];
UEPos = [43.083376, -77.680509];
RASPos = [43.083378, -77.680697];

gNBSite = txsite("Name","GCI_gNB","Latitude",gNBPos(1),"Longitude",gNBPos(2),"AntennaAngle",gNBAntDir(1:2),"AntennaHeight",10,"TransmitterFrequency",fc);
RASSite = rxsite("Name","GCI_RAS","Latitude",RASPos(1),"Longitude",RASPos(2),"AntennaAngle",RASAntDir(1:2),"AntennaHeight",1);
UESite = rxsite("Name","GCI_UE","Latitude",UEPos(1),"Longitude",UEPos(2),"AntennaAngle",UEAntDir(1:2),"AntennaHeight",1);

gNBSite.show();
UESite.show();
RASSite.show();

%% Plot rays
pm = propagationModel("raytracing","Method","sbr","MaxNumReflections",reflectionsOrder);
UE_rays = raytrace(gNBSite,UESite,pm,"Type","pathloss");
RAS_rays = raytrace(gNBSite,RASSite,pm,"Type","pathloss");

plot(UE_rays{1});
plot(RAS_rays{1});

%% Nullify RAS
ofdmInfo = nrOFDMInfo(NRB,SCS);
RAS_channel = getChannelObj(RAS_rays,fc,ofdmInfo.SampleRate);

lambda = c/fc;
RASArray = phased.NRRectangularPanelArray('Size',[RASAntSize(1:2) 1 1],'Spacing', [0.5*lambda*[1 1] 1 1]);
RASArray.ElementSet = {phased.IsotropicAntennaElement};   % isotropic antenna element
RAS_channel.ReceiveAntennaArray = RASArray;
RAS_channel.ReceiveArrayOrientation = [RASAntDir(1); (-1)*RASAntDir(2); 0]; 

gNBArray = phased.NRRectangularPanelArray('Size',[gNBAntSize(1:2) 1 1],'Spacing', [0.1*[1 1] 1 1]);
gNBArray.ElementSet = {phased.IsotropicAntennaElement};
RAS_channel.TransmitAntennaArray = gNBArray;
RAS_channel.TransmitArrayOrientation = [gNBAntDir(1); (-1)*gNBAntDir(2); 0];

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
hest_gR = nrPerfectChannelEstimate(RAS_pathGains, RAS_pathFilters, NRB, SCS, nSlots, RAS_offset, RAS_sampleTimes);
hest_gR_temp = permute(mean(reshape(hest_gR,[],RAS_Nr,RAS_Nt)),[2,3,1]);
% hest_gR_temp = hest_gR_temp./norm(hest_gR_temp);

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
hest_gU = nrPerfectChannelEstimate(UE_pathGains, UE_pathFilters, NRB, SCS, nSlots, UE_offset, UE_sampleTimes);
hest_gU = getChannelCoeffs(hest_gU, scOffset, noRBs);
noise_UE = wgn(size(UE_rxWaveform,1),size(UE_rxWaveform,2),-59);

[~,~,V] = svd(hest_gU);
w_gNB = V(:,1:nLayers).';
P = null(hest_gR_temp);
w_heu = P*P'*w_gNB.';

h_k = hest_gU_all.';
H = zeros(M,M,K(rowidx));
for i = 1:K(rowidx)
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
		for i = 1:K(rowidx)
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

IP_def = 10*log10(norm(hest_gR_temp*w_gNB').^2)
IP_heu = 10*log10(norm(hest_gR_temp*w_heu).^2)
IP_opt = 10*log10(norm(hest_gR_temp*w_opt').^2)

wf = complex(randn(T,1),randn(T,1));
SR_def = log2(1+snr((wf*w_gNB)*hest_gU.',noise_UE))
SR_heu = log2(1+snr((wf*w_heu.')*hest_gU.',noise_UE))
SR_opt = log2(1+snr((wf*(w_opt./(scale.^0.25)))*hest_gU.',noise_UE))