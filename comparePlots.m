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

leg = legend("No Null","w_{heu}","w_{mm}","w_{mmn}","location","best");
leg.ItemTokenSize = [10 18];
leg.Orientation = 'horizontal';
xlabel("No. of UEs");

lineCount = size(a.Children,1);
lineSpec = ['-^k'; '-*r'; '-og'; '-xb'];
lineColor = [0 0 0; 1 0 0; 0 0.7 0; 0 0 1];
for i = 1:lineCount
	a.Children(5-i).Marker = lineSpec(i,2);
	a.Children(5-i).Color = lineColor(i,:);
end