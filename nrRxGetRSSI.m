carrier = nrCarrierConfig;
carrier.NSizeGrid = 26;
carrier.SubcarrierSpacing = 15;
carrier.NSlot = 1;
carrier.NFrame = 0;

csirs = nrCSIRSConfig;
csirs.CSIRSType = {'nzp'};
csirs.CSIRSPeriod = {[5 1]};
csirs.Density = {'one'};
csirs.RowNumber = 2;
csirs.SymbolLocations = {1};
csirs.SubcarrierLocations = {6};
csirs.NumRB = 26;

powerCSIRS = 0;
refSym = db2mag(powerCSIRS)*nrCSIRS(carrier,csirs);
refInd = nrCSIRSIndices(carrier,csirs);
offset = nrTimingEstimate(carrier,rxWaveform,refInd,refSym);
csirswf = rxWaveform(offset-7680:1+offset,:);

rxGrid = nrOFDMDemodulate(carrier,csirswf);
meas = nrCSIRSMeasurements(carrier,csirs,rxGrid);
meas.RSSIPerAntenna = meas.RSSIPerAntenna - 160;