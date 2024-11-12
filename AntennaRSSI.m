%Create antenna object
f0 = 5.925e9;
ant = dipole('Length',1);

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
	'TransmitterPower', 20);

%Calculate Tx power
% Z = impedance(tx.Antenna,tx.TransmitterFrequency);
% If = feedCurrent(tx.Antenna,tx.TransmitterFrequency);
% Irms = norm(If)/sqrt(2);
% Ptx = real(Z)*(Irms)^2;
% tx.TransmitterPower = Ptx;

rx = rxsite('Name','Green Bank 100m', ...
        'Latitude', 38.433006, ...
        'Longitude', -79.839753, ...
        'ReceiverSensitivity', -90);

pm = propagationModel("longley-rice");
ssTwoRef = sigstrength(rx,tx,pm)

