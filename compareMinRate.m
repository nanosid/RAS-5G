clear;
% Define constant parameters
c = physconst("LightSpeed");
e = 10^(-220/10);
pb = 1;
hmag = 1;
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
K = [1 2 4 8 16];
% K = [1 4 8];
N = 100;
UE_snr_default = zeros(length(K),K(end));
UE_snr_apr = zeros(size(UE_snr_default));
UE_snr_cf = zeros(size(UE_snr_default));
UE_snr_opt = zeros(size(UE_snr_default));
latency = zeros(4,length(K),N);
nullification = zeros(4,length(K),N);
sumrate = zeros(4,length(K),N);
minrate = zeros(4,length(K),N);

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
hest_gR_temp = hest_gR_temp./norm(hest_gR_temp).*hmag;

for var = 1:N
UE_snr_default = zeros(length(K),K(end));
UE_snr_apr = zeros(size(UE_snr_default));
UE_snr_cf = zeros(size(UE_snr_default));
UE_snr_opt = zeros(size(UE_snr_default));
for rowidx=1:length(K)
	iterationID = sprintf("Iteration %d %d...\n",var, rowidx);
	disp(iterationID);
	UENames = "UE"+string(1:K(rowidx));
	UESites = rxsite("Name",UENames,"Latitude",0,"Longitude",0,"AntennaAngle",UEAntDir(1:2),"AntennaHeight",2);
	hest_gU_all = zeros(K(rowidx),RAS_Nt);
	noise_UE = zeros(size(RAS_rxWaveform,1),K(rowidx));
	sumk = zeros(RAS_Nt,RAS_Nt);
	prodk = eye(RAS_Nt,RAS_Nt);
	for i=1:K(rowidx)
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
		hest_gU_all(i,:) = hest_gU_all(i,:)./norm(hest_gU_all(i,:)).*hmag;
		sumk = sumk + hest_gU_all(i,:)'*hest_gU_all(i,:);
		prodk = prodk*(hest_gU_all(i,:)'*hest_gU_all(i,:));
		noise_UE(:,i) = wgn(size(UE_rxWaveform,1),size(UE_rxWaveform,2),-59);
	end
	
	% Approximate solution
	tic
	%hest_gU_all = hest_gU_all./norm(hest_gU_all);
	[U,~,V] = svd(hest_gU_all);
	w_gNB = V(:,1:nLayers).';
	latency(1,rowidx,var) = toc;
	
	P = null(hest_gR_temp);
	w_apr = P*P'*w_gNB.';
	latency(2,rowidx,var) = toc;

	% iSGD
% 	A = diag(D);
% 	w_opt = ones(1,RAS_Nt);
% 	del = ones(1,RAS_Nt);
% 	eta = 0.01;
% 	count = 0;
% 	test = -1;
% 	tic
% 	while(norm(del)*eta > 1e-3)
% 		for i = 1:K(rowidx)
% 			del = hest_gU_all(i,:).*(w_opt*hest_gU_all(i,:)') - hest_gU_all(i,:).*conj(A(i));
% 			norm(del)
% 			w_opt = w_opt - eta.*del;
% 		end
% 		count = count + 1;
% 	end
% 	latency(2,rowidx,var) = toc;

	% Closed form solution
	% tic;
	% [V,D] = eig(prodk);
	% %prodk = V*diag(nthroot(diag(D),K))/V;
	% G = hest_gR_temp/(sumk+prodk);
	% w_cf = (e^(1/2)/(hest_gR_temp*G')).*G;
	% latency(3,rowidx,var) = toc;

    % CVX for no null
    % tic
    % cvx_begin
    %     variable w(1,M) complex
    %     minimize(norm(hest_gU_all*w' - pb.*ones(K(rowidx),1)));
	% 	subject to
	% 		norm(w) <= pb;
    % cvx_end
	% w_cf = w;
	% latency(3,rowidx,var) = toc;

	% CVX with RAS null
	% tic
	% cvx_begin
	% 	variable w(1,M) complex
	% 	minimize(norm(hest_gU_all*w' - pb.*ones(K(rowidx),1)));
	% 	subject to
	% 		norm(hest_gR_temp*w') <= e*10^(-pb/10);
	% 		norm(w) <= pb;
	% cvx_end
	% w_opt = w;
	% latency(4,rowidx,var) = toc;

	h_k = hest_gU_all.';
	H = zeros(M,M,K(rowidx));
	for i = 1:K(rowidx)
		H(:,:,i) = h_k(:,i)*h_k(:,i)';
	end
	hgR = hest_gR_temp.';
	HgR = hgR*hgR';

	% Max-min fair with no RAS null
	tic	
	cvx_begin
		variable tau nonnegative
		variable X(M,M) hermitian
		variable s(K,1)
		minimize(-tau)
		subject to
			% norm(diag(X*HgR)) <= e
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
		% scale = e./(norm(x(:,i)'*hgR).^2);
		% x(:,i) = x(:,i).*(scale.^(0.25));
	end
	[~,idx] = max(min(abs(x'*h_k),[],2));
	w_cf = x(:,idx).';
	latency(3,rowidx,var) = toc;

	% Max-min fair with RAS null
	tic
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
	latency(4,rowidx,var) = toc;

	wf = complex(randn(T,1),randn(T,1));
	
	for i = 1:K(rowidx)
		UE_wf_bf = (wf*w_gNB)*hest_gU_all(i,:)';
		UE_snr_default(rowidx,i) = snr(UE_wf_bf,noise_UE(:,i));
		UE_wf_bf = (wf*w_apr.')*hest_gU_all(i,:)';
		UE_snr_apr(rowidx,i) = snr(UE_wf_bf,noise_UE(:,i));
		UE_wf_bf = (wf*w_cf)*hest_gU_all(i,:)';
		UE_snr_cf(rowidx,i) = snr(UE_wf_bf,noise_UE(:,i));
		UE_wf_bf = (wf*(w_opt./(scale.^0.25)))*hest_gU_all(i,:)';
		UE_snr_opt(rowidx,i) = snr(UE_wf_bf,noise_UE(:,i));
	end
	nullification(1,rowidx,var) = 10*log10(norm(hest_gR_temp*w_gNB').^2);
	nullification(2,rowidx,var) = 10*log10(norm(hest_gR_temp*w_apr).^2);
	nullification(3,rowidx,var) = 10*log10(norm(hest_gR_temp*w_cf').^2);
	nullification(4,rowidx,var) = 10*log10(norm(hest_gR_temp*w_opt').^2);
	if mean(UE_snr_default(rowidx,:)) > -1
		sumrate(1,rowidx,var) = sum(log2(1+UE_snr_default(rowidx,:)));
		minrate(1,rowidx,var) = min(log2(1+UE_snr_default(rowidx,1:K(rowidx))));
	else
		sumrate(1,rowidx,var) = sum(log2(1+(10.^(UE_snr_default(rowidx,:)))./10));
		minrate(1,rowidx,var) = min(log2(1+(10.^(UE_snr_default(rowidx,1:K(rowidx))))./10));
	end
	if mean(UE_snr_apr(rowidx,:)) > -1
		sumrate(2,rowidx,var) = sum(log2(1+UE_snr_apr(rowidx,:)));
		minrate(2,rowidx,var) = min(log2(1+UE_snr_apr(rowidx,1:K(rowidx))));
	else
		sumrate(2,rowidx,var) = sum(log2(1+(10.^(UE_snr_apr(rowidx,:)))./10));
		minrate(2,rowidx,var) = min(log2(1+(10.^(UE_snr_apr(rowidx,1:K(rowidx))))./10));
	end
	if mean(UE_snr_cf(rowidx,:)) > -1
		sumrate(3,rowidx,var) = sum(log2(1+UE_snr_cf(rowidx,:)));
		minrate(3,rowidx,var) = min(log2(1+UE_snr_cf(rowidx,1:K(rowidx))));
	else
		sumrate(3,rowidx,var) = sum(log2(1+(10.^(UE_snr_cf(rowidx,:)))./10));
		minrate(3,rowidx,var) = min(log2(1+(10.^(UE_snr_cf(rowidx,1:K(rowidx))))./10));
	end
	if mean(UE_snr_opt(rowidx,:)) > -1
		sumrate(4,rowidx,var) = sum(log2(1+UE_snr_opt(rowidx,:)));
		minrate(4,rowidx,var) = min(log2(1+UE_snr_opt(rowidx,1:K(rowidx))));
	else
		sumrate(4,rowidx,var) = sum(log2(1+(10.^(UE_snr_opt(rowidx,:)))./10));
		minrate(4,rowidx,var) = min(log2(1+(10.^(UE_snr_opt(rowidx,1:K(rowidx))))./10));
	end
end
end
null_plot = mean((nullification(:,:,:)),3);
sumrate_plot = abs(mean(sumrate(:,:,:),3));
latency_plot = mean(latency(:,:,:),3)./(1e-3);
minrate_plot = abs(mean(minrate(:,:,:),3));
save("schemeComparison");
%% Plots
figure;
plot(1:length(K),null_plot);
ylabel("Interference Power (dBW)");
comparePlots;
f.Position = [5 5 2 2];
a.FontSize = 9;
leg.ItemTokenSize = [10 18];
leg.Orientation = 'horizontal';
leg.FontSize = 5;
leg.NumColumns = 2;
yticks([-400 -300 -200 -100 0]);
xticks(1:5);
xticklabels(["1" "2" "4" "8" "16"]);
% xticks(1:3);
% xticklabels(["1" "4" "8"]);

figure;
plot(1:length(K),sumrate_plot);
ylabel("Sum Rate (bits/s/Hz)");
comparePlots;
f.Position = [5 5 2 2];
a.FontSize = 9;
leg.ItemTokenSize = [10 18];
leg.Orientation = "horizontal";
leg.FontSize = 5;
leg.NumColumns = 1;
yticks([0 25 50 75 100]);
xticks(1:5);
xticklabels(["1" "2" "4" "8" "16"]);
% xticks(1:3);
% xticklabels(["1" "4" "8"]);

figure;
plot(1:length(K),minrate_plot);
ylabel("Min Rate (bits/s/Hz)");
comparePlots;
f.Position = [5 5 2 2];
a.FontSize = 9;
leg.ItemTokenSize = [10 18];
leg.Orientation = "horizontal";
leg.FontSize = 5;
leg.NumColumns = 2;
ylim([0 10]);
% yticks([0 25 50 75 100]);
xticks(1:5);
xticklabels(["1" "2" "4" "8" "16"]);
% xticks(1:3);
% xticklabels(["1" "4" "8"]);

figure;
semilogy(1:length(K),latency_plot);
ylabel("Latency (ms)");
comparePlots;
f.Position = [5 5 2 2];
a.FontSize = 9;
leg.ItemTokenSize = [10 18];
leg.Orientation = 'horizontal';
leg.FontSize = 5;
leg.NumColumns = 2;
yticks([1e-2 1e-1 1 1e1 1e2 1e3]);
ylim([1e-2 1e3]);
xticks(1:5);
xticklabels(["1" "2" "4" "8" "16"]);
% xticks(1:3);
% xticklabels(["1" "4" "8"]);

% Display beampattern of gNB
% show(gNBSite); 
% show(RASSite);
% show(UESites);
% gNBSite.Antenna = clone(RAS_channel.TransmitAntennaArray);
% gNBSite.Antenna.Taper = w_apr;
% pattern(gNBSite,fc,"Size",250);