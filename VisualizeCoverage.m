%Create antenna object
f0 = 5.925e9;
ant = dipole('Length',1);
% f = figure;
% show(ant)
% view(-220,30)
% 
% %Calculate electric field
% z = -15;
% x = (-250:4:250)*1e3;
% y = (-100:4:400)*1e3;
% 
% [X,Y] = meshgrid(x,y);
% numpoints = length(x)*length(y);
% points = [X(:)'; Y(:)'; z*ones(1,numel(X))];
% 
% E = EHfields(ant,f0,points); % Units: V/m
% 
% Emag = zeros(1,numpoints);
% for m=1:numpoints
%     Emag(m) = norm(E(:,m)/sqrt(2));
% end
% Emag = 20*log10(reshape(Emag,length(y),length(x))); % Units: dBV/m
% Emag = Emag + 120; % Units: dBuV/m
% 
% d_min = min(Emag(:));
% d_max = max(Emag(:));
% del = (d_max-d_min)/12;
% d_vec = round(d_min:del:d_max);
% 
% if isvalid(f)
%     close(f)
% end
% figure
% contourf(X*1e-3,Y*1e-3,Emag,d_vec,'showtext','on')
% title('Field Strength (dB\muV/m) on flat Earth (1V tx)')
% xlabel('lateral (km)')
% ylabel('boresight (km)')
% c = colorbar; 
% set(get(c,'title'),'string','dB\muV/m')

%Setup gNB site
lat =  38.416865;
lon = -79.982815;
h = 7;
az = 10;
xyrot = wrapTo180(az - 90);

tx = txsite('Name','Snow Shoe', ...
    'Latitude',lat, ...
    'Longitude',lon, ...
    'Antenna',ant, ...
    'AntennaHeight',h, ...
    'AntennaAngle',xyrot, ...
    'TransmitterFrequency',f0, ...
	'TransmitterPower', 120);

% %Calculate Tx power
% Z = impedance(tx.Antenna,tx.TransmitterFrequency);
% If = feedCurrent(tx.Antenna,tx.TransmitterFrequency);
% Irms = norm(If)/sqrt(2);
% Ptx = real(Z)*(Irms)^2;
% tx.TransmitterPower = Ptx;

% Launch Site Viewer with no terrain
viewer = siteviewer("Terrain", "none");
viewer.Basemap = 'topographic'; 
sigStrengths = [-100 -80 -60 -40 -20 0];
coverage(tx,'freespace', ...
    'SignalStrengths',sigStrengths, ...
    'Colormap','parula', ...
    'ColorLimits',[-100 0])