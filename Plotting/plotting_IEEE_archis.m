% Plotting IEEE Aerospace Paper Architecture Examples

width = 600;
height = 600;

Nf = 3;
Nl = 100;
run_idx = '018';
legend_size = 16;

etR_years = (1:3:3*length(bandwidths))/(365.25);

filename = "../Data/run"+run_idx+"/Nf"+int2str(Nf)+"Nl"+int2str(Nl)+".mat";
load(filename, "XSats_i", "bandwidths", "NSats", "X_opt")

figure_num = 2323;
P2_name = 'Mars';
[F1, p2_handle] = plot_DSC_Modular_Structure(XSats_i, run_idx, figure_num, P2_name, 125);
set(F1, 'Position', [F1.Position(1) F1.Position(2) width height])
[lgd, icons] = legend("Sun", "Earth", "Mars", "HRN Satellites");
change_Legend_Icon_Sizes(icons, [legend_size, legend_size, legend_size, legend_size])
set(gca, FontSize=16, FontName="Times New Roman");
set(lgd, FontSize=legend_size+2, LineWidth=0.5);
uistack(p2_handle, "up")
set(icons, FontSize=legend_size+1, FontName="Times New Roman");


%%

F2 = figure(figure_num+10000);
set(F2, 'Position', [F2.Position(1) F2.Position(2) width height])
plot(etR_years,bandwidths, LineWidth=1.5)
hold on
yline(mean(bandwidths), LineWidth=2.5, LineStyle="--")
yline(min(bandwidths), Color="red", LineWidth=2.5, LineStyle="-.")
hold off
xlabel("Simulation time [years]", FontWeight="bold")
ylabel("Mars-Earth Data Rate [Mbps]", FontWeight="bold")
legend("Network Data Rate", "Mean Data Rate", "Min. Data Rate")
set(gca, FontSize=16, FontName="Times New Roman")
ylim([0 50*ceil(max(bandwidths)/50)])
xlim([0, max(etR_years)])


% Display key metrics
% 
% disp(NSats)
% disp(X_opt)

function [] = change_Legend_Icon_Sizes(icons, sizes)
    N = length(sizes);
    marks = findobj(icons, 'Type', 'patch');

    for i = 1:N
        marks(i).MarkerSize = sizes(i);
        marks(i).LineWidth = 1;
    end

end
