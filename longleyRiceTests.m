%Load antenna object
load("dipole4x2.mat");

%Setup telescope site
lat2 = 42.929852;
lon2 = -77.500154;
% lat2 = 38.433006;
% lon2 = -79.839753;
rx = rxsite('Name','Green Bank 100m', ...
        'Latitude', lat2, ...
        'Longitude', lon2, ...
        'ReceiverSensitivity', -90);

%Setup gNB site
f0 = 5.925e9;
lat =  43.0308333;
lon = -77.4569444;
% lat = 38.6469444;
% lon = -79.5420833;
h = 32;
az = azimuth(lat, lon, lat2, lon2)+10;
xyrot = wrapTo180(az);

tx = txsite('Name','Snow Shoe', ...
    'Latitude',lat, ...
    'Longitude',lon, ...
    'Antenna',ant, ...
    'AntennaHeight',h, ...
    'AntennaAngle',xyrot, ...
    'TransmitterFrequency',f0, ...
	'TransmitterPower', 0.25);

pm = propagationModel("longley-rice");
coverage(tx,rx,pm,"SignalStrengths",-100:0);
ss = sigstrength(rx,tx,pm)