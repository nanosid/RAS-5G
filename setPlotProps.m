%% Size and position
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

%% Line Modification
lineCount = size(a.Children,1);
gNBLegend = ["gNB=11km","gNB=9km","gNB<1km"];
%antennaLegend = ["M=4","M=8","M=16","M=32"];
antennaLegend = ["M=4","M=8","M=16","M=4 null","M=8 null","M=16 null"];
%fcLegend = ["fc=0.7GHz","fc=1.8GHz","fc=5.9GHz","fc=7GHz"];
fcLegend = ["fc=0.7GHz","fc=1.8GHz","fc=5.9GHz","fc=0.7GHz null","fc=1.8GHz null","fc=5.9GHz null"];
refsLegend = ["R=1","R=2","R=3","R=4"];
%UELegend = ["K=1","K=8","K=16","K=32"];
UELegend = ["K=1","K=8","K=16","K=1 null","K=8 null","K=16 null"];
lineSpec = ['-*r'; '-og'; '-xb'; '-^m'];
%lineSpec = ['-*r'; '-og'; '-xb'];
lineColor = [1 0 0; 0 1 0; 0 0 1; 1 0 1];
%lineColor = [1 0 0; 0 1 0; 0 0 1];
for i = 1:lineCount
	a.Children(i).Marker = lineSpec(mod(i,4)+1,2);
	%a.Children(i).Marker = lineSpec(mod(i,3)+1,2);
	a.Children(i).Color = lineColor(mod(i,4)+1,:);
	%a.Children(i).Color = lineColor(mod(i,3)+1,:);
	if i < 5
		a.Children(i).LineStyle = "-";
	end
end