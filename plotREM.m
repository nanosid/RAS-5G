clear;

interIdx = 1:10:401;
load('REMTest_2.mat');
%% w_heu REM
figure;
s = surf(X(interIdx),Y(interIdx),remPow_heu(interIdx,interIdx),'FaceColor','interp');
figProps(s);

%% w_gNB REM
figure;
s = surf(X(interIdx),Y(interIdx),remPow_gNB(interIdx,interIdx),'FaceColor','interp');
figProps(s);

%% w_opt REM
figure;
s = surf(X(interIdx),Y(interIdx),remPow_opt(interIdx,interIdx),'FaceColor','interp');
figProps(s);

%% Local Functions
function [] = figProps(s)
view([0 90]);
s.EdgeColor = 'none';
c = colorbar;
c.Label.String = 'Interference Power (dBW)';
f = gcf;
a = gca;
f.Units = 'inches';
f.Position = [6 6 2.5 1.5];
a.FontName = 'Times New Roman';
a.FontSizeMode = 'manual';
a.FontSize = 9;
a.LabelFontSizeMultiplier = 1;
grid on;
a.GridAlpha = 0.4;
xlabel('Longitude Distance (m)');
ylabel('Latitude Distance (m)');
t = text(-35,35,"RAS");
t.FontName = 'Times New Roman';
t = text(40,-155,"gNB");
t.FontName = 'Times New Roman';
hold on;
plot(40,-110,'ob');
plot(0,0,'or');
end