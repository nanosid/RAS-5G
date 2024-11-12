f = gcf;
a = gca;
f.Units = 'inches';
f.Position = [6 6 3 1.5];
a.FontName = 'Times New Roman';
a.FontSizeMode = 'manual';
a.FontSize = 9;
a.LabelFontSizeMultiplier = 1;
grid on;
a.GridAlpha = 0.4;

lineCount = size(a.Children,1);
lineSpec = ['-*r'; '-og'; '-xb'];
lineColor = [1 0 0; 0 0.7 0; 0 0 1];
for i = 1:lineCount
	a.Children(i).LineStyle = "--";
	a.Children(i).Marker = lineSpec(mod(i,3)+1,2);
	a.Children(i).Color = lineColor(mod(i,3)+1,:);
	if i < 4
		a.Children(i).LineStyle = "-";
	end
end