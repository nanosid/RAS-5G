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
ofdmInfo = nrOFDMInfo(NRB,SCS);
RASSite = rxsite("Name","Ionia_RAS","Latitude",RASPos(1),"Longitude",RASPos(2),"AntennaAngle",RASAntDir(1:2),"AntennaHeight",5);

% Define possible values for metrics
gNBLocs = [42.930223, -77.501120];
numants = [8 16];
fcs = [739e6 1880e6 5925e6];
numUEs = [1 4 8];
gNBDist = distance(gNBLocs(1),gNBLocs(2),RASPos(1),RASPos(2));
gNBDist = deg2rad(gNBDist)*(earthRadius*1e-3);
N = 100;
sigma = 0:-5:-50;

RAS_null_heu = zeros(3,length(sigma),length(numants),N);
RAS_null_heu_err = zeros(size(RAS_null_heu));
RAS_null_heu_err_k = zeros(size(RAS_null_heu));
RAS_null_opt = zeros(size(RAS_null_heu));
RAS_null_opt_err = zeros(size(RAS_null_heu));
RAS_null_opt_err_k = zeros(size(RAS_null_heu));
RAS_null_def = zeros(size(RAS_null_heu));
RAS_null_def_err = zeros(size(RAS_null_heu));
RAS_null_def_err_k = zeros(size(RAS_null_heu));

UE_snr_heu = zeros(3,length(sigma),length(numants),N,numUEs(end));
UE_snr_heu_err = zeros(size(UE_snr_heu));
UE_snr_heu_err_k = zeros(size(UE_snr_heu));
UE_snr_opt = zeros(size(UE_snr_heu));
UE_snr_opt_err = zeros(size(UE_snr_heu));
UE_snr_opt_err_k = zeros(size(UE_snr_heu));
UE_snr_def = zeros(size(UE_snr_heu));
UE_snr_def_err = zeros(size(UE_snr_heu));
UE_snr_def_err_k = zeros(size(UE_snr_heu));

sumrate_heu = zeros(size(RAS_null_heu));
sumrate_heu_err = zeros(size(RAS_null_heu));
sumrate_heu_err_k = zeros(size(RAS_null_heu));
sumrate_opt = zeros(size(RAS_null_heu));
sumrate_opt_err = zeros(size(RAS_null_heu));
sumrate_opt_err_k = zeros(size(RAS_null_heu));
sumrate_def = zeros(size(RAS_null_heu));
sumrate_def_err = zeros(size(RAS_null_heu));
sumrate_def_err_k = zeros(size(RAS_null_heu));

for rowidx = 1
	for i1 = 1:length(sigma)
		for i2 = 1:length(numants)
			gNBPos = gNBLocs;
			reflectionsOrder = 3;
			switch rowidx
				case 1
					gNBAntSize = [numants(i2) 2];
					fc = 5925e6;
					K = 8;
				case 2
					gNBAntSize = [8 2];
					fc = fcs(i2);
					K = 8;
				case 3
					gNBAntSize = [8 2];
					fc = 5925e6;
					K = numUEs(i2);
			end
			M = prod(gNBAntSize);
			for var = 1:N
				gNBSite = txsite("Name","Victor_gNB","Latitude",gNBPos(1),"Longitude",gNBPos(2),"AntennaAngle",...
					gNBAntDir(1:2),"AntennaHeight",32,"TransmitterFrequency",fc);
			
				% Define Raytracing prop model and display plotted rays
				pm = propagationModel("raytracing","Method","sbr","MaxNumReflections",reflectionsOrder);
				RAS_rays = raytrace(gNBSite,RASSite,pm,"Type","pathloss");
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
				hest_gR = nrPerfectChannelEstimate(RAS_pathGains, RAS_pathFilters, NRB, SCS, nSlots, RAS_offset, RAS_sampleTimes);
				hest_gR_temp = getChannelCoeffs(hest_gR,scOffset,noRBs);
				hest_gR_temp = hest_gR_temp./norm(hest_gR_temp);
				% hest_gR_err = hest_gR_temp + wgn(size(hest_gR_temp,1),size(hest_gR_temp,2),sigma(i1),'complex');
				amp = abs(hest_gR_temp) + ((10^(2*sigma(i1)/10).*randn(size(hest_gR_temp))));
				theta = angle(hest_gR_temp) + ((10^(2*sigma(i1)/10).*randn(size(hest_gR_temp))));
				% hest_gR_err = amp.*exp(1i.*(theta+((10^(sigma(i1)/10).*randn(size(theta))))));
				% hest_gR_err = (amp+((10^(sigma(i1)/10).*randn(size(amp))))).*exp(1i.*theta);
				hest_gR_err = amp.*exp(1i.*theta);

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
					UE_channel = getChannelObj(UE_rays,fc,ofdmInfo.SampleRate);
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
					% hest_gU_err(i,:) = hest_gU_all(i,:) + wgn(size(hest_gU_all(i,:),1),size(hest_gU_all(i,:),2),sigma(i1),'complex');
                	amp_U = abs(hest_gU_all(i,:)) + ((10^(2*sigma(i1)/10).*randn(size(hest_gU_all(i,:)))));
					theta_U = angle(hest_gU_all(i,:)) + ((10^(2*sigma(i1)/10).*randn(size(hest_gU_all(i,:)))));
					hest_gU_err(i,:) = amp_U.*exp(1i.*theta_U);
					noise_UE(:,i) = wgn(size(UE_rxWaveform,1),size(UE_rxWaveform,2),-59);
				end
				% No error
				[~,~,V] = svd(hest_gU_all);
				w_gNB = V(:,1:nLayers).';
				P = null(hest_gR_temp);
				w_copt = P*P'*w_gNB.';
				
				% Error in h_k
				[~,~,V] = svd(hest_gU_err);
				w_gNB_err_k = V(:,1:nLayers).';
				w_copt_err_k = P*P'*w_gNB_err_k.';

				% Error in h_gR
				P_err = null(hest_gR_err);
				w_copt_err_gR = P_err*P_err'*w_gNB.';

				% Opt no error
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

				% Opt error in h_gR
				h_k = hest_gU_all.';
				H = zeros(M,M,K);
				for i = 1:K
					H(:,:,i) = h_k(:,i)*h_k(:,i)';
				end
				hgR = hest_gR_err.';
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
				w_opt_err_gR = x(:,idx).';

				% Opt error in h_k
				h_k = hest_gU_err.';
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
				w_opt_err_k = x(:,idx).';

				wf = complex(randn(T,1),randn(T,1));
				for i = 1:K
					% Error in h_k
					UE_wf_bf = (wf*w_copt_err_k.')*hest_gU_all(i,:)';
					UE_snr_heu_err_k(rowidx,i1,i2,var,i) = snr(UE_wf_bf,noise_UE(:,i));
					if UE_snr_heu_err_k(rowidx,i1,i2,var,i) > -1
						sumrate_heu_err_k(rowidx,i1,i2,var) = sumrate_heu_err_k(rowidx,i1,i2,var) + log2(1+UE_snr_heu_err_k(rowidx,i1,i2,var,i));
					end

					% Error in h_gR
					UE_wf_bf = (wf*w_copt_err_gR.')*hest_gU_all(i,:)';
					UE_snr_heu_err(rowidx,i1,i2,var,i) = snr(UE_wf_bf,noise_UE(:,i));
					if UE_snr_heu_err(rowidx,i1,i2,var,i) > -1
						sumrate_heu_err(rowidx,i1,i2,var) = sumrate_heu_err(rowidx,i1,i2,var) + log2(1+UE_snr_heu_err(rowidx,i1,i2,var,i));
					end
					
					% No error
					UE_wf_bf = (wf*w_copt.')*hest_gU_all(i,:)';
					UE_snr_heu(rowidx,i1,i2,var,i) = snr(UE_wf_bf,noise_UE(:,i));
					if UE_snr_heu(rowidx,i1,i2,var,i) > -1
						sumrate_heu(rowidx,i1,i2,var) = sumrate_heu(rowidx,i1,i2,var) + log2(1+UE_snr_heu(rowidx,i1,i2,var,i));
					end

					% Opt error in h_k
					UE_wf_bf = (wf*(w_opt_err_k./(scale.^0.25)))*hest_gU_all(i,:)';
					UE_snr_opt_err_k(rowidx,i1,i2,var,i) = snr(UE_wf_bf,noise_UE(:,i));
					if UE_snr_opt_err_k(rowidx,i1,i2,var,i) > -1
						sumrate_opt_err_k(rowidx,i1,i2,var) = sumrate_opt_err_k(rowidx,i1,i2,var) + log2(1+UE_snr_opt_err_k(rowidx,i1,i2,var,i));
					end

					% Opt error in h_gR
					UE_wf_bf = (wf*(w_opt_err_gR./(scale.^0.25)))*hest_gU_all(i,:)';
					UE_snr_opt_err(rowidx,i1,i2,var,i) = snr(UE_wf_bf,noise_UE(:,i));
					if UE_snr_opt_err(rowidx,i1,i2,var,i) > -1
						sumrate_opt_err(rowidx,i1,i2,var) = sumrate_opt_err(rowidx,i1,i2,var) + log2(1+UE_snr_opt_err(rowidx,i1,i2,var,i));
					end
					
					% Opt no error
					UE_wf_bf = (wf*(w_opt./(scale.^0.25)))*hest_gU_all(i,:)';
					UE_snr_opt(rowidx,i1,i2,var,i) = snr(UE_wf_bf,noise_UE(:,i));
					if UE_snr_opt(rowidx,i1,i2,var,i) > -1
						sumrate_opt(rowidx,i1,i2,var) = sumrate_opt(rowidx,i1,i2,var) + log2(1+UE_snr_opt(rowidx,i1,i2,var,i));
					end

					% Def error in h_k
					UE_wf_bf = (wf*w_gNB_err_k)*hest_gU_all(i,:)';
					UE_snr_def_err_k(rowidx,i1,i2,var,i) = snr(UE_wf_bf,noise_UE(:,i));
					if UE_snr_def_err_k(rowidx,i1,i2,var,i) > -1
						sumrate_def_err_k(rowidx,i1,i2,var) = sumrate_def_err_k(rowidx,i1,i2,var) + log2(1+UE_snr_def_err_k(rowidx,i1,i2,var,i));
					end

					% Def error in h_gR
					UE_wf_bf = (wf*w_gNB)*hest_gU_all(i,:)';
					UE_snr_def_err(rowidx,i1,i2,var,i) = snr(UE_wf_bf,noise_UE(:,i));
					if UE_snr_def_err(rowidx,i1,i2,var,i) > -1
						sumrate_def_err(rowidx,i1,i2,var) = sumrate_def_err(rowidx,i1,i2,var) + log2(1+UE_snr_def_err(rowidx,i1,i2,var,i));
					end
					
					% Def no error
					UE_wf_bf = (wf*w_gNB)*hest_gU_all(i,:)';
					UE_snr_def(rowidx,i1,i2,var,i) = snr(UE_wf_bf,noise_UE(:,i));
					if UE_snr_def(rowidx,i1,i2,var,i) > -1
						sumrate_def(rowidx,i1,i2,var) = sumrate_def(rowidx,i1,i2,var) + log2(1+UE_snr_def(rowidx,i1,i2,var,i));
					end
				end
				RAS_null_heu_err_k(rowidx,i1,i2,var) = 10*log10(norm(hest_gR_temp*w_copt_err_k).^2);
				RAS_null_heu_err(rowidx,i1,i2,var) = 10*log10(norm(hest_gR_temp*w_copt_err_gR).^2);
				RAS_null_heu(rowidx,i1,i2,var) = 10*log10(norm(hest_gR_temp*w_copt).^2);
				RAS_null_opt_err_k(rowidx,i1,i2,var) = 10*log10(norm(hest_gR_temp*w_opt_err_k').^2);
				RAS_null_opt_err(rowidx,i1,i2,var) = 10*log10(norm(hest_gR_temp*w_opt_err_gR').^2);
				RAS_null_opt(rowidx,i1,i2,var) = 10*log10(norm(hest_gR_temp*w_opt').^2);
				RAS_null_def_err_k(rowidx,i1,i2,var) = 10*log10(norm(hest_gR_temp*w_gNB_err_k').^2);
				RAS_null_def_err(rowidx,i1,i2,var) = 10*log10(norm(hest_gR_temp*w_gNB').^2);
				RAS_null_def(rowidx,i1,i2,var) = 10*log10(norm(hest_gR_temp*w_gNB').^2);
				iterationID = sprintf("Iteration %d,%d,%d,%d...\n",rowidx,i1,i2,var);
				disp(iterationID);
			end
		end
	end
end

idx = RAS_null_heu == -Inf;
RAS_null_heu(idx) = -350;
idx = RAS_null_opt == -Inf;
RAS_null_opt(idx) = -350;
save("estErrorFull_100_3");

%% Plots

% figure;
% plot(sigma,permute(mean(RAS_null_def_err(1,:,1,:) - RAS_null_def(1,:,1,:),4),[3 2 1]));
% hold on
% plot(sigma,permute(mean(RAS_null_opt_err(1,:,1,:) - RAS_null_opt(1,:,1,:),4),[3 2 1]));
% plot(sigma,permute(mean(RAS_null_heu_err(1,:,1,:) - RAS_null_heu(1,:,1,:),4),[3 2 1]));
% plot(sigma,permute(mean(RAS_null_def_err(1,:,2,:) - RAS_null_def(1,:,2,:),4),[3 2 1]));
% plot(sigma,permute(mean(RAS_null_opt_err(1,:,2,:) - RAS_null_opt(1,:,2,:),4),[3 2 1]));
% plot(sigma,permute(mean(RAS_null_heu_err(1,:,2,:) - RAS_null_heu(1,:,2,:),4),[3 2 1]));
% semiPlots
% f.Position = [4 4 3 2];
% xlabel('Error Variance');
% xticklabels(flip(["1","10^{-1}","10^{-2}","10^{-3}","10^{-4}","10^{-5}"]));
% ylabel('Change in Int. Pow. (dBW)');
% ylim([-50 350]);
% leg = legend("No Null, M=8","w_{mmn}, M=8","w_{heu}, M=8","No Null, M=16","w_{mmn}, M=16","w_{heu}, M=16",'location','northwest');
% leg.ItemTokenSize = [10 18];
% leg.FontSize = 5;
% leg.NumColumns = 2;
% title('Error in gNB-RAS Channel');
% 
% figure;
% plot(sigma,permute(mean(RAS_null_def_err_k(1,:,1,:) - RAS_null_def(1,:,1,:),4),[3 2 1]));
% hold on
% plot(sigma,permute(mean(RAS_null_opt_err_k(1,:,1,:) - RAS_null_opt(1,:,1,:),4),[3 2 1]));
% plot(sigma,permute(mean(RAS_null_heu_err_k(1,:,1,:) - RAS_null_heu(1,:,1,:),4),[3 2 1]));
% plot(sigma,permute(mean(RAS_null_def_err_k(1,:,2,:) - RAS_null_def(1,:,2,:),4),[3 2 1]));
% plot(sigma,permute(mean(RAS_null_opt_err_k(1,:,2,:) - RAS_null_opt(1,:,2,:),4),[3 2 1]));
% plot(sigma,permute(mean(RAS_null_heu_err_k(1,:,2,:) - RAS_null_heu(1,:,2,:),4),[3 2 1]));
% semiPlots
% f.Position = [4 4 3 2];
% xlabel('Error Variance');
% xticklabels(flip(["1","10^{-1}","10^{-2}","10^{-3}","10^{-4}","10^{-5}"]));
% ylabel('Change in Int. Pow. (dBW)');
% ylim([-50 350]);
% leg = legend("No Null, M=8","w_{mmn}, M=8","w_{heu}, M=8","No Null, M=16","w_{mmn}, M=16","w_{heu}, M=16",'location','northwest');
% leg.ItemTokenSize = [10 18];
% leg.FontSize = 5;
% leg.NumColumns = 2;
% title('Error in gNB-UE Channels');
% 
% figure;
% plot(sigma,permute(mean(sumrate_def_err(1,:,1,:) - sumrate_def(1,:,1,:),4),[3 2 1]));
% hold on
% plot(sigma,permute(mean(sumrate_opt_err(1,:,1,:) - sumrate_opt(1,:,1,:),4),[3 2 1]));
% plot(sigma,permute(mean(sumrate_heu_err(1,:,1,:) - sumrate_heu(1,:,1,:),4),[3 2 1]));
% plot(sigma,permute(mean(sumrate_def_err(1,:,2,:) - sumrate_def(1,:,2,:),4),[3 2 1]));
% plot(sigma,permute(mean(sumrate_opt_err(1,:,2,:) - sumrate_opt(1,:,2,:),4),[3 2 1]));
% plot(sigma,permute(mean(sumrate_heu_err(1,:,2,:) - sumrate_heu(1,:,2,:),4),[3 2 1]));
% semiPlots
% f.Position = [4 4 3 2];
% xlabel('Error Variance');
% xticklabels(flip(["1","10^{-1}","10^{-2}","10^{-3}","10^{-4}","10^{-5}"]));
% ylabel('Change in Sum Rate (bits/s/Hz)');
% ylim([-15 15]);
% yticks(-15:5:15);
% leg = legend("No Null, M=8","w_{mmn}, M=8","w_{heu}, M=8","No Null, M=16","w_{mmn}, M=16","w_{heu}, M=16",'location','northwest');
% leg.ItemTokenSize = [10 18];
% leg.FontSize = 5;
% leg.NumColumns = 2;
% title('Error in gNB-RAS Channel');
% 
% figure;
% plot(sigma,permute(mean(sumrate_def_err_k(1,:,1,:) - sumrate_def(1,:,1,:),4),[3 2 1]));
% hold on
% plot(sigma,permute(mean(sumrate_opt_err_k(1,:,1,:) - sumrate_opt(1,:,1,:),4),[3 2 1]));
% plot(sigma,permute(mean(sumrate_heu_err_k(1,:,1,:) - sumrate_heu(1,:,1,:),4),[3 2 1]));
% plot(sigma,permute(mean(sumrate_def_err_k(1,:,2,:) - sumrate_def(1,:,2,:),4),[3 2 1]));
% plot(sigma,permute(mean(sumrate_opt_err_k(1,:,2,:) - sumrate_opt(1,:,2,:),4),[3 2 1]));
% plot(sigma,permute(mean(sumrate_heu_err_k(1,:,2,:) - sumrate_heu(1,:,2,:),4),[3 2 1]));
% semiPlots
% f.Position = [4 4 3 2];
% xlabel('Error Variance');
% xticklabels(flip(["1","10^{-1}","10^{-2}","10^{-3}","10^{-4}","10^{-5}"]));
% ylabel('Change in Sum Rate (bits/s/Hz)');
% ylim([-15 15]);
% yticks(-15:5:15);
% leg = legend("No Null, M=8","w_{mmn}, M=8","w_{heu}, M=8","No Null, M=16","w_{mmn}, M=16","w_{heu}, M=16",'location','northwest');
% leg.ItemTokenSize = [10 18];
% leg.FontSize = 5;
% leg.NumColumns = 2;
% title('Error in gNB-UE Channels');

%% Better plots

figure;
plot(sigma,permute(mean(RAS_null_def_err(1,:,1,:),4),[3 2 1]));
hold on
plot(sigma,permute(mean(RAS_null_opt_err(1,:,1,:),4),[3 2 1]));
plot(sigma,permute(mean(RAS_null_heu_err(1,:,1,:),4),[3 2 1]));
plot(sigma,permute(mean(RAS_null_def(1,:,1,:),4),[3 2 1]));
plot(sigma,permute(mean(RAS_null_opt(1,:,1,:),4),[3 2 1]));
plot(sigma,permute(mean(RAS_null_heu(1,:,1,:),4),[3 2 1]));
semiPlots
f.Position = [4 4 3 1.5];
xlabel('Error Variance');
xticklabels(flip(["1","10^{-1}","10^{-2}","10^{-3}","10^{-4}","10^{-5}"]));
ylabel('Interference Power (dBW)');
ylim([-350 0]);
leg = legend("No Null","w_{mmn}","w_{heu}",'location','northwest');
leg.ItemTokenSize = [10 18];
leg.FontSize = 7;
leg.NumColumns = 3;
% title('Error in gNB-RAS Channel');

figure;
plot(sigma,permute(mean(RAS_null_def_err_k(1,:,1,:),4),[3 2 1]));
hold on
plot(sigma,permute(mean(RAS_null_opt_err_k(1,:,1,:),4),[3 2 1]));
plot(sigma,permute(mean(RAS_null_heu_err_k(1,:,1,:),4),[3 2 1]));
plot(sigma,permute(mean(RAS_null_def(1,:,1,:),4),[3 2 1]));
plot(sigma,permute(mean(RAS_null_opt(1,:,1,:),4),[3 2 1]));
plot(sigma,permute(mean(RAS_null_heu(1,:,1,:),4),[3 2 1]));
semiPlots
f.Position = [4 4 3 1.5];
xlabel('Error Variance');
xticklabels(flip(["1","10^{-1}","10^{-2}","10^{-3}","10^{-4}","10^{-5}"]));
ylabel('Interference Power (dBW)');
ylim([-350 0]);
leg = legend("No Null","w_{mmn}","w_{heu}",'location','northwest');
leg.ItemTokenSize = [10 18];
leg.FontSize = 7;
leg.NumColumns = 1;
% title('Error in gNB-RAS Channel');

figure;
plot(sigma,permute(mean(sumrate_def_err(1,:,1,:),4),[3 2 1]));
hold on
plot(sigma,permute(mean(sumrate_opt_err(1,:,1,:),4),[3 2 1]));
plot(sigma,permute(mean(sumrate_heu_err(1,:,1,:),4),[3 2 1]));
plot(sigma,permute(mean(sumrate_def(1,:,1,:),4),[3 2 1]));
plot(sigma,permute(mean(sumrate_opt(1,:,1,:),4),[3 2 1]));
plot(sigma,permute(mean(sumrate_heu(1,:,1,:),4),[3 2 1]));
semiPlots
f.Position = [4 4 3 1.5];
xlabel('Error Variance');
xticklabels(flip(["1","10^{-1}","10^{-2}","10^{-3}","10^{-4}","10^{-5}"]));
ylabel('Sum Rate (bits/s/Hz)');
ylim([0 50]);
% yticks(-15:5:15);
leg = legend("No Null","w_{mmn}","w_{heu}",'location','northwest');
leg.ItemTokenSize = [10 18];
leg.FontSize = 7;
leg.NumColumns = 1;
% title('Error in gNB-RAS Channel');

figure;
plot(sigma,permute(mean(sumrate_def_err_k(1,:,1,:),4),[3 2 1]));
hold on
plot(sigma,permute(mean(sumrate_opt_err_k(1,:,1,:),4),[3 2 1]));
plot(sigma,permute(mean(sumrate_heu_err_k(1,:,1,:),4),[3 2 1]));
plot(sigma,permute(mean(sumrate_def(1,:,1,:),4),[3 2 1]));
plot(sigma,permute(mean(sumrate_opt(1,:,1,:),4),[3 2 1]));
plot(sigma,permute(mean(sumrate_heu(1,:,1,:),4),[3 2 1]));
semiPlots
f.Position = [4 4 3 1.5];
xlabel('Error Variance');
xticklabels(flip(["1","10^{-1}","10^{-2}","10^{-3}","10^{-4}","10^{-5}"]));
ylabel('Sum Rate (bits/s/Hz)');
ylim([0 50]);
% yticks(-15:5:15);
leg = legend("No Null","w_{mmn}","w_{heu}",'location','northwest');
leg.ItemTokenSize = [10 18];
leg.FontSize = 7;
leg.NumColumns = 1;
% title('Error in gNB-RAS Channel');