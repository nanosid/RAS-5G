clear;
% Define constant parameters
c = physconst("LightSpeed");
gNBAntDir = [0 0].';
gNBAntSize = [8 2];
gNBPos = [42.930223, -77.501120];
fc = 60e9;
gNBSite = txsite("Name","Victor_gNB","Latitude",gNBPos(1),"Longitude",gNBPos(2),"AntennaAngle",...
	gNBAntDir(1:2),"AntennaHeight",32,"TransmitterFrequency",fc);
lambda = c/fc;
gNBArray = phased.NRRectangularPanelArray('Size',[gNBAntSize(1:2) 1 1],'Spacing', [0.5*lambda*[1 1] 1 1]);
gNBArray.ElementSet = {phased.NRAntennaElement};
gNBSite.Antenna = gNBArray;
show(gNBSite);
pattern(gNBSite,fc,"Size",250);