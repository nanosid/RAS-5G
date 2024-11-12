function [channelObj] = getChannelObj(raysObj,fc,sr)
pathToAs = [raysObj{1}.PropagationDelay] - min([raysObj{1}.PropagationDelay]);
avgPathGains = -[raysObj{1}.PathLoss];
pathAoDs = [raysObj{1}.AngleOfDeparture];
pathAoAs = [raysObj{1}.AngleOfArrival];
isLOS = any([raysObj{1}.LineOfSight]);
channelObj = nrCDLChannel;
channelObj.DelayProfile = 'Custom';
channelObj.PathDelays = pathToAs;
channelObj.AveragePathGains = avgPathGains;
channelObj.AnglesAoD = pathAoDs(1,:);
channelObj.AnglesZoD = 90-pathAoDs(2,:);
channelObj.AnglesAoA = pathAoAs(1,:);
channelObj.AnglesZoA = 90-pathAoAs(2,:);
channelObj.HasLOSCluster = isLOS;
channelObj.CarrierFrequency = fc;
channelObj.NormalizePathGains = false;
channelObj.NormalizeChannelOutputs = false;
channelObj.SampleRate = sr;
end

