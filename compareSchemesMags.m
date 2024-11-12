clear;
% Define constant parameters
c = physconst("LightSpeed");
e = 10^(-220/10);
pb = [1 5 10];
hmag = [1e-4 1e-3 1e-2 1e-1 1];
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
% gNBDist = distance(gNBLocs(1,1),gNBLocs(1,2),RASPos(1),RASPos(2));
% gNBDist = deg2rad(gNBDist).*(earthRadius*1e-3);
gNBPos = gNBLocs;
gNBAntSize = [8 2];
M = prod(gNBAntSize);
fc = 1880e6;
reflectionsOrder = 3;
K = 8;
N = 20;
UE_snr_default = zeros(length(hmag),length(pb),K,N);
UE_snr_apr = zeros(size(UE_snr_default));
UE_snr_cf = zeros(size(UE_snr_default));
UE_snr_opt = zeros(size(UE_snr_default));
latency = zeros(4,length(hmag),length(pb),N);
nullification = zeros(size(latency));
sumrate = zeros(size(latency));

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

% Get gNB-RAS channel coefficients over all RBs and OFDM symbols
RAS_pathFilters = getPathFilters(RAS_channel);
[RAS_offset,~] = nrPerfectTimingEstimate(RAS_pathGains, RAS_pathFilters);
hest_gR = nrPerfectChannelEstimate(RAS_pathGains, RAS_pathFilters, NRB, SCS, nSlots, RAS_offset, RAS_sampleTimes);
hest_gR_temp = permute(mean(reshape(hest_gR,[],RAS_Nr,RAS_Nt)),[2,3,1]);
hest_gR_temp = hest_gR_temp./norm(hest_gR_temp);

for var = 1:N
for i1=1:length(hmag)
for i2=1:length(pb)
	% hest_gR_temp = hest_gR_temp./norm(hest_gR_temp).*hmag(i1);
	UENames = "UE"+string(1:K);
	UESites = rxsite("Name",UENames,"Latitude",0,"Longitude",0,"AntennaAngle",UEAntDir(1:2),"AntennaHeight",2);
	hest_gU_all = zeros(K,RAS_Nt);
	noise_UE = zeros(size(RAS_rxWaveform,1),K);
	for i=1:K
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
		hest_gU = nrPerfectChannelEstimate(UE_pathGains, UE_pathFilters, NRB, SCS, nSlots, UE_offset, UE_sampleTimes);
		hest_gU_all(i,:) = getChannelCoeffs(hest_gU, scOffset, noRBs);
		hest_gU_all(i,:) = hest_gU_all(i,:)./norm(hest_gU_all(i,:)).*hmag(i1);
		noise_UE(:,i) = wgn(size(UE_rxWaveform,1),size(UE_rxWaveform,2),-69);
	end
	
	% Approximate solution
	tic
	[~,~,V] = svd(hest_gU_all);
	w_gNB = V(:,1:nLayers).';
	latency(1,i1,i2,var) = toc;
	
	tic
	P = null(hest_gR_temp);
	w_apr = (P*P')*w_gNB.';
	latency(2,i1,i2,var) = toc;

    % CVX for no null
    tic
    cvx_begin
        variable w(1,M) complex
        minimize(norm(hest_gU_all*w' - pb(i2).*ones(K,1)));
		subject to
			norm(w) <= pb(i2);
    cvx_end
	w_cf = w;
	latency(3,i1,i2,var) = toc;

	% CVX with RAS null
	tic
	cvx_begin
		variable w(1,M) complex
		minimize(norm(hest_gU_all*w' - pb(i2).*ones(K,1)));
		subject to
			norm(hest_gR_temp*w') <= e*10^(-pb(i2)/10);
			norm(w) <= pb(i2);
	cvx_end
	w_opt = w;
	latency(4,i1,i2,var) = toc;

	wf = complex(randn(T,1),randn(T,1));	
	for i = 1:K
		UE_wf_bf = (wf*w_gNB)*hest_gU_all(i,:)';
		UE_snr_default(i1,i2,i,var) = snr(UE_wf_bf,noise_UE(:,i))+pb(i2);
		UE_wf_bf = (wf*w_apr.')*hest_gU_all(i,:)';
		UE_snr_apr(i1,i2,i,var) = snr(UE_wf_bf,noise_UE(:,i))+pb(i2);
		UE_wf_bf = (wf*w_cf)*hest_gU_all(i,:)';
		UE_snr_cf(i1,i2,i,var) = snr(UE_wf_bf,noise_UE(:,i));
		UE_wf_bf = (wf*w_opt)*hest_gU_all(i,:)';
		UE_snr_opt(i1,i2,i,var) = snr(UE_wf_bf,noise_UE(:,i));
	end
	nullification(1,i1,i2,var) = 10*log10(norm(hest_gR_temp*w_gNB').^2)+pb(i2);
	nullification(2,i1,i2,var) = 10*log10(norm(hest_gR_temp*w_apr).^2)+pb(i2);
	nullification(3,i1,i2,var) = 10*log10(norm(hest_gR_temp*w_cf').^2);
	nullification(4,i1,i2,var) = 10*log10(norm(hest_gR_temp*w_opt').^2);
	if mean(UE_snr_default(i1,i2,:,var)) > -1
		sumrate(1,i1,i2,var) = sum(log2(1+UE_snr_default(i1,i2,:,var)));
	else
		sumrate(1,i1,i2,var) = sum(log2(1+10.^(UE_snr_default(i1,i2,:,var)./10)));
	end
	if mean(UE_snr_apr(i1,i2,:,var)) > -1
		sumrate(2,i1,i2,var) = sum(log2(1+UE_snr_apr(i1,i2,:,var)));
	else
		sumrate(2,i1,i2,var) = sum(log2(1+10.^(UE_snr_apr(i1,i2,:,var)./10)));
	end
	if mean(UE_snr_cf(i1,i2,:,var)) > -1
		sumrate(3,i1,i2,var) = sum(log2(1+UE_snr_cf(i1,i2,:,var)));
	else
		sumrate(3,i1,i2,var) = sum(log2(1+10.^(UE_snr_cf(i1,i2,:,var)./10)));
	end
	if mean(UE_snr_opt(i1,i2,:,var)) > -1
		sumrate(4,i1,i2,var) = sum(log2(1+UE_snr_opt(i1,i2,:,var)));
	else
		sumrate(4,i1,i2,var) = sum(log2(1+10.^(UE_snr_opt(i1,i2,:,var)./10)));
	end
	iterationID = sprintf("Iteration %d %d %d...\n",var, i1, i2);
	disp(iterationID);
end
end
end
null_plot = mean(nullification,4);
sumrate_plot = abs(mean(sumrate,4));
latency_plot = mean(latency,4)./(1e-3);

%% Plots

% pb = 1
figure;
plot(1:length(hmag),null_plot(:,:,1));
comparePlots;
xticks(1:length(hmag));
xticklabels(["10^{-4}" "10^{-3}" "10^{-2}" "10^{-1}" "1"]);
xlabel("Channel Magnitude");
ylabel("Interference Power (dBW)");
yticks([-400 -300 -200 -100 0]);
ylim([-400 50]);
title("Power Budget, \rho_b = 1dBW");

figure;
plot(1:length(hmag),sumrate_plot(:,:,1));
comparePlots;
xticks(1:length(hmag));
xticklabels(["10^{-4}" "10^{-3}" "10^{-2}" "10^{-1}" "1"]);
xlabel("Channel Magnitude");
ylabel("Sum Rate (bits/s/Hz)");
yticks([0 25 50 75 100]);
ylim([0 75]);
title("Power Budget, \rho_b = 1dBW");

% pb = 5
figure;
plot(1:length(hmag),null_plot(:,:,2));
comparePlots;
xticks(1:length(hmag));
xticklabels(["10^{-4}" "10^{-3}" "10^{-2}" "10^{-1}" "1"]);
xlabel("Channel Magnitude");
ylabel("Interference Power (dBW)");
yticks([-400 -300 -200 -100 0]);
ylim([-400 50]);
title("Power Budget, \rho_b = 5dBW");

figure;
plot(1:length(hmag),sumrate_plot(:,:,2));
comparePlots;
xticks(1:length(hmag));
xticklabels(["10^{-4}" "10^{-3}" "10^{-2}" "10^{-1}" "1"]);
xlabel("Channel Magnitude");
ylabel("Sum Rate (bits/s/Hz)");
yticks([0 25 50 75 100]);
ylim([0 75]);
title("Power Budget, \rho_b = 5dBW");

% pb = 10
figure;
plot(1:length(hmag),null_plot(:,:,3));
comparePlots;
xticks(1:length(hmag));
xticklabels(["10^{-4}" "10^{-3}" "10^{-2}" "10^{-1}" "1"]);
xlabel("Channel Magnitude");
ylabel("Interference Power (dBW)");
yticks([-400 -300 -200 -100 0]);
ylim([-400 50]);
title("Power Budget, \rho_b = 10dBW");

figure;
plot(1:length(hmag),sumrate_plot(:,:,3));
comparePlots;
xticks(1:length(hmag));
xticklabels(["10^{-4}" "10^{-3}" "10^{-2}" "10^{-1}" "1"]);
xlabel("Channel Magnitude");
ylabel("Sum Rate (bits/s/Hz)");
yticks([0 25 50 75 100]);
ylim([0 75]);
title("Power Budget, \rho_b = 10dBW");