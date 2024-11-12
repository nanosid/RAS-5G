clear;
% Define basic parameters
fc = 5.925e9;
c = physconst("LightSpeed");
lambda = c/fc;

% Define gNB and RAS sites
gNBPos = [42.930223, -77.501120];				% Ionia, NY radio telescope
gNBAntSize = [4 4];
gNBAntDir = [0 0].';
RASPos = [42.929852, -77.500154];
RASAntSize = [1 1];
RASAntDir = [0 0].';
reflectionsOrder = 4;
SCS = 15;
NRB = 52;
nSlots = 0;
nLayers = 1;
scOffset = 0;
noRBs = 1;
gNBSite = txsite("Name","Victor_gNB","Latitude",gNBPos(1),"Longitude",gNBPos(2),"AntennaAngle",gNBAntDir(1:2),...
	"AntennaHeight",32,"TransmitterFrequency",fc);
RASSite = rxsite("Name","Ionia_RAS","Latitude",RASPos(1),"Longitude",RASPos(2),"AntennaAngle",RASAntDir(1:2),...
	"AntennaHeight",5);

% Define channel models for RAS
pm = propagationModel("raytracing","Method","sbr","MaxNumReflections",reflectionsOrder);
RAS_rays = raytrace(gNBSite,RASSite,pm,"Type","pathloss");
pathToAs = [RAS_rays{1}.PropagationDelay] - min([RAS_rays{1}.PropagationDelay]);
avgPathGains = -[RAS_rays{1}.PathLoss];
pathAoDs = [RAS_rays{1}.AngleOfDeparture];
pathAoAs = [RAS_rays{1}.AngleOfArrival];
isLOS = any([RAS_rays{1}.LineOfSight]);

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

RASArray = phased.NRRectangularPanelArray('Size',[RASAntSize(1:2) 1 1],'Spacing', [0.5*lambda*[1 1] 1 1]);
RASArray.ElementSet = {phased.IsotropicAntennaElement};   % isotropic antenna element
RAS_channel.ReceiveAntennaArray = RASArray;
RAS_channel.ReceiveArrayOrientation = [RASAntDir(1); (-1)*RASAntDir(2); 0];  % the (-1) converts elevation to downtilt

gNBArray = phased.NRRectangularPanelArray('Size',[gNBAntSize(1:2) 1 1],'Spacing', [0.5*lambda*[1 1] 1 1]);
gNBArray.ElementSet = {phased.IsotropicAntennaElement};
RAS_channel.TransmitAntennaArray = gNBArray;
RAS_channel.TransmitArrayOrientation = [gNBAntDir(1); (-1)*gNBAntDir(2); 0];   % the (-1) converts elevation to downtilt

ofdmInfo = nrOFDMInfo(NRB,SCS);
RAS_channel.SampleRate = ofdmInfo.SampleRate;
%channel.ChannelFiltering = false;

% Design sample waveform for RAS
RAS_channelInfo = info(RAS_channel);
T = RAS_channel.SampleRate * 1e-3;
RAS_Nt = RAS_channelInfo.NumTransmitAntennas;
RAS_Nr = RAS_channelInfo.NumReceiveAntennas;
txWaveform = complex(randn(T,RAS_Nt),randn(T,RAS_Nt));
[RAS_rxWaveform,RAS_pathGains,RAS_sampleTimes] = RAS_channel(txWaveform);
noise_RAS = wgn(size(RAS_rxWaveform,1),size(RAS_rxWaveform,2),-137);
RAS_snr_default = snr(RAS_rxWaveform,noise_RAS);

% Get gNB-RAS channel coefficients over all RBs and OFDM symbols
RAS_pathFilters = getPathFilters(RAS_channel);
[RAS_offset,~] = nrPerfectTimingEstimate(RAS_pathGains, RAS_pathFilters);
hest_gR = nrPerfectChannelEstimate(RAS_pathGains, RAS_pathFilters, NRB, SCS, nSlots, RAS_offset, RAS_sampleTimes);
hest_gR_temp = permute(mean(reshape(hest_gR,[],RAS_Nr,RAS_Nt)),[2,3,1]);

% Define UE sites
% UEPos = [43.019708, -77.461838];
K = 8;
UEAntSize = [1 1];
UEAntDir = [0 0].';
% locDisp = abs(randn(K,2).*5e-3);
% UEPos = locDisp + [ones(K,1).*gNBPos(1) ones(K,1).*gNBPos(2)];
UENames = "UE"+string(1:K);
UESites = rxsite("Name",UENames,"Latitude",0,"Longitude",0,"AntennaAngle",UEAntDir(1:2),"AntennaHeight",2);
UEPos = zeros(K,2);

% Now derive channel for K UEs and the beamforming weights for each UE
UE_snr_default = zeros(1,K);
UE_snr_bf = zeros(1,K);
hest_gU_all = zeros(K,RAS_Nt);
w_gNB_all = zeros(K,RAS_Nt);
noise_UE = zeros(size(RAS_rxWaveform,1),K);
% UE_wf_bf = zeros(size(noise_UE));
% wf = complex(randn(T,1),randn(T,1));

for i=1:K
	pathToAs = [];
    while(isempty(pathToAs))
        locDisp = randn(1,2).*5e-3;
        UEPos(i,:) = locDisp + gNBPos;
        % UEName = "UE"+string(i);
        % UESite = rxsite("Name",UEName,"Latitude",UEPos(i,1),"Longitude",UEPos(i,2),"AntennaAngle",UEAntDir(1:2),"AntennaHeight",2);
        UESites(i).Latitude = UEPos(i,1);
        UESites(i).Longitude = UEPos(i,2);
        UE_rays = raytrace(gNBSite,UESites(i),pm,"Type","pathloss");
        pathToAs = [UE_rays{1}.PropagationDelay] - min([UE_rays{1}.PropagationDelay]);
    end
    % UE_rays = raytrace(gNBSite,UESites(i),pm,"Type","pathloss");
	% pathToAs = [UE_rays{1}.PropagationDelay] - min([UE_rays{1}.PropagationDelay]);
	%if isempty(pathToAs)
		%continue;
	%else
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
		UE_Nt = UE_channelInfo.NumTransmitAntennas;
		UE_Nr = UE_channelInfo.NumReceiveAntennas;
		[~,UE_pathGains,UE_sampleTimes] = UE_channel(txWaveform);
		UE_pathFilters = getPathFilters(UE_channel);
		[UE_offset,~] = nrPerfectTimingEstimate(UE_pathGains, UE_pathFilters);
		hest_gU = nrPerfectChannelEstimate(UE_pathGains, UE_pathFilters, NRB, SCS, nSlots, UE_offset, UE_sampleTimes);
		hest_gU_all(i,:) = getChannelCoeffs(hest_gU, scOffset, noRBs);
        UE_rxWaveform = txWaveform*hest_gU_all(i,:)';
		noise_UE(:,i) = wgn(size(UE_rxWaveform,1),size(UE_rxWaveform,2),-89);
		UE_snr_default(i) = snr(UE_rxWaveform,noise_UE(:,i));
        [w_gNB_all(i,:), ~, ~] = getBeamformingWeights(hest_gU, nLayers, scOffset, noRBs);
% 		wf_bf = wf*w_gNB_all(i,:);
% 		[UE_wf_bf(:,i),~,~] = UE_channel(wf_bf);
% 		UE_snr_bf(i) = snr(UE_wf_bf(:,i),noise_UE(:,i));
	%end
end
% UE_snr_gain = UE_snr_bf - UE_snr_default;
% [~,idx] = max(UE_snr_gain);
% w_gNB = w_gNB_all(idx,:);

% Derive beamformer for K UEs
tic
hest_gU_all = hest_gU_all./norm(hest_gU_all);
[U,D,V] = svd(hest_gU_all);

% A = zeros(RAS_Nt,1);
% A(1:K) = diag(D);
% w_gNB = (V\A).';
w_gNB = V(:,1:nLayers).';

% A = diag(D);
% w_gNB = ones(1,RAS_Nt);
% del = ones(1,RAS_Nt);
% eta = 0.1;
% count = 0;
% test = -1;
% while(norm(del)*eta > 1e-1)
% 	for i = 1:K
% 		del = hest_gU_all(i,:).*(w_gNB*hest_gU_all(i,:)') - hest_gU_all(i,:).*conj(A(i));
% 		norm(del)
% 		w_gNB = w_gNB - eta.*del;
% 	end
% 	count = count + 1;
% end

% W = eye(RAS_Nt);
% for i = 1:K
% 	[~,~,Q] = svd(hest_gU_all(i,:));
% 	W = W*Q;
% end
% w_gNB = W(:,1:nLayers).';

% UE_snr_gain_mean = zeros(1,K);
% wf = complex(randn(T,1),randn(T,1));
% for j = 1:K
% 	wf_bf = wf*w_gNB_all(j,:);
% 	for i = 1:K
% 		if(UE_snr_default(i) == 0)
% 			continue;
% 		else
% 			[UE_wf_bf,~,~] = UE_channel(wf_bf);
% 			UE_snr_bf(i) = snr(UE_wf_bf,noise_UE(:,i));
% 		end
% 	end
% 	UE_snr_gain_mean(j) = mean(UE_snr_bf - UE_snr_default);
% end
% [~,idx] = max(UE_snr_gain_mean);
% w_gNB = w_gNB_all(idx,:);

% Now get beamformer which also nullifies RAS
P = null(hest_gR_temp);
w_copt = P*P'*w_gNB.';
toc
% w_copt = w_copt./norm(w_copt);

% SNR after beamforming
wf = complex(randn(T,1),randn(T,1));
wf_bf = wf*w_copt.';
[RAS_wf_bf,~,~] = RAS_channel(wf_bf);
RAS_snr_bf = snr(RAS_wf_bf,noise_RAS);
RAS_snr_gain = RAS_snr_bf - RAS_snr_default;

% SNR at K UEs
for i = 1:K
	if(UE_snr_default(i) == 0)
		continue;
	else
		UE_wf_bf = wf_bf*hest_gU_all(i,:)';
		UE_snr_bf(i) = snr(UE_wf_bf,noise_UE(:,i));
	end
end
UE_snr_gain = (UE_snr_bf - UE_snr_default)./5;

% Display beampattern of gNB
show(gNBSite); 
show(RASSite);
show(UESites);
gNBSite.Antenna = clone(RAS_channel.TransmitAntennaArray);
gNBSite.Antenna.Taper = w_copt;
pattern(gNBSite,fc,"Size",1000);