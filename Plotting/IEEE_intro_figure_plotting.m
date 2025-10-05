

fig_num = 1211313;
DSC_Mod_Plotting_from_Data(fig_num, '017', 2, 33, [], [] ,[], 210)

fig = figure(fig_num);
hold on
scatter3(NaN, NaN, NaN, MarkerEdgeColor="black", MarkerFaceColor=[0.7,0.7,0.7], Marker="square", SizeData=200);
plot(NaN, NaN, Color="red");
hold off

xlabel("x [AU]", FontWeight="bold")
ylabel("y [AU]", FontWeight="bold")

legend_size = 16;
[lgd, icons] = legend("Sun", "Earth", "Mars", "Relay Satellites", "Downlink Path", FontSize = legend_size+2, FontName="Times New Roman");
change_Legend_Icon_Sizes(icons, [legend_size, legend_size, legend_size, legend_size])
set(gca, FontSize=16, FontName="Times New Roman");
set(lgd, FontSize=legend_size+2, LineWidth=0.5);
axis off
set(icons, FontSize=legend_size+1, FontName="Times New Roman");




function [] = change_Legend_Icon_Sizes(icons, sizes)
    N = length(sizes);
    marks = findobj(icons, 'Type', 'patch');

    for i = 1:N
        marks(i).MarkerSize = sizes(i);
        marks(i).LineWidth = 1;
    end

end