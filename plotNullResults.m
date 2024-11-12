load('UE_snr_gain_1000.mat');
load('RAS_snr_gain_1000.mat');
RAS_null = -flip(mean(RAS_snr_gain,3));
UE_null = -flip(mean(UE_snr_gain,3));
dist = 1:10;
lineSpec = ['-*r'; '-og'; '-xb'; '-+c'; '-^m'; '-vy'];
figure;
for i=1:size(RAS_null,2)
	plot(dist,RAS_null(:,i),lineSpec(i,:));
	hold on;
end
grid on;
legend("M=2","M=4","M=8","M=16","M=32","M=64");
ylim([25 40]);
xlim([1 10]);
xlabel('gNB-RAS distance (km)');
ylabel('Nullification (dB)');
figure;
for i=1:size(UE_null,2)
	plot(dist,UE_null(:,i),lineSpec(i,:));
	hold on;
end
grid on;
legend("M=2","M=4","M=8","M=16","M=32","M=64");
ylim([0 35]);
xlim([1 10]);
xlabel('gNB-RAS distance (km)');
ylabel('Nullification (dB)');