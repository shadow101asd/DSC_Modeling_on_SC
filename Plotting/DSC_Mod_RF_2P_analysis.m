%% Collect Data from DSC_Modular_RF_2P Runs

run_idx = '018';
Nf_range = 1:3;
Nl_range = 1:125;

Bs = NaN([length(Nl_range), length(Nf_range)]);
minBs = Bs;
for Nf = Nf_range
    for Nl = Nl_range
        filename = "../Data/run"+run_idx+"/Nf"+int2str(Nf)+"Nl"+int2str(Nl)+".mat";
        
        if isfile(filename)
            % Extract and save data of interest
            Bs(Nl, Nf) = load(filename, "B").B;
            data = load(filename, "bandwidths");
            if isfield(data, "bandwidths")
                bandwidths = data.bandwidths;
                minBs(Nl,Nf) = min(bandwidths);
            else
                minBs(Nl,Nf) = NaN;
            end
        else
            warning("Missing file: " + filename)
        end
    end
end

%% Plotting

plot_range = Nl_range(1:end);
w = 750; % Width
h = 350; % Height

f1 = figure(str2double(run_idx));
set(f1, 'Position', [f1.Position(1) f1.Position(2) w h])
clf(str2double(run_idx))

h1 = plot(plot_range, Bs(plot_range,1), LineWidth=2);
hold on
h2 = plot(plot_range, Bs(plot_range,2), LineWidth=2);
h3 = plot(plot_range, Bs(plot_range,3), LineWidth=2);
hold off

lgd = [];
for Nf = Nf_range
    lgd = [lgd, "Nf = " + int2str(Nf)];
end
legend(lgd, Location='northwest', AutoUpdate='off')

uistack(h2,'top')
uistack(h1,'top')

xlim([min(Nl_range), max(Nl_range)])
ylim([0 ceil(max(max(Bs))/100)*100])
xlabel("Number of Starship Launches N_l", FontWeight="bold")
ylabel("Average Mars-Earth Data Rate [Mbps]", FontWeight="bold")
% title("Run " + run_idx + ": Ongoing")
set(gca, FontName = "Times New Roman", FontSize=14)
set(gca, 'Box','off');

f2 = figure(100+str2double(run_idx));
set(f2, 'Position', [f2.Position(1) f2.Position(2) w (h/2+30)])
clf(100+str2double(run_idx))

plot(plot_range, minBs(plot_range,:), LineWidth=2)

lgd = [];
for Nf = Nf_range
    lgd = [lgd, "Nf = " + int2str(Nf)];
end
legend(lgd, Location='northwest')
xlim([min(Nl_range), max(Nl_range)])
ylim([0 ceil(max(max(minBs))/100)*100])
xlabel("Number of Starship Launches N_l", FontWeight="bold")
ylabel("Min. Data Rate [Mbps]", FontWeight="bold")
% title("Run " + run_idx + ": Ongoing")
set(gca, FontName = "Times New Roman", FontSize=14)
set(gca, 'Box','off');


%% This part is manual for now

% Arch. I
Nl_R1 = 1:40;
B1 = Bs(Nl_R1,1);

% Arch II
Nl_R2 = 41:max(Nl_range);
B2 = Bs(Nl_R2,1);

% Arch III
Nl_R3 = 12:58;
B3 = Bs(Nl_R3,2);

% Arch IV
Nl_R4 = 59:max(Nl_range);
B4 = Bs(Nl_R4,2);

% Arch V
Nl_R5 = 43:max(Nl_range);
B5 = Bs(Nl_R5,3);


f3 = figure(200+str2double(run_idx));
set(f3, 'Position', [f3.Position(1) f3.Position(2) w h])
clf(200+str2double(run_idx))


xlabel("Number of Starship Launches N_l", FontWeight="bold")
ylabel("Average Mars-Earth Data Rate [Mbps]", FontWeight="bold")
set(gca, FontName = "Times New Roman", FontSize=14)

hold on
p1 = plot(Nl_R1, B1,LineWidth=2);
labelCurves(p1, "I", 0.55, -1)

p2 = plot(Nl_R2, B2,LineWidth=2);
labelCurves(p2, "II", 0.5, 1)

p3 = plot(Nl_R3, B3,LineWidth=2);
labelCurves(p3, "III", 0.5, 1)

p4 = plot(Nl_R4, B4,LineWidth=2);
labelCurves(p4, "IV", 0.5, 1)

p5 = plot(Nl_R5, B5,LineWidth=2);
labelCurves(p5, "V", 0.4, 1)

hold off

xlim([min(Nl_range), max(Nl_range)])
ylim([0 ceil(max(max(Bs))/100)*100])
set(gca, 'Box','off');


%% Functions


function labelCurves(h, labels, frac, over_under)
% h: line handles from plot(...)
% labels: string/cellstr (or [] to use DisplayName)
% frac: 0..1 position along each curve (default 0.95)

if nargin < 3 || isempty(frac), frac = 0.95; end
if nargin < 2 || isempty(labels)
    labels = arrayfun(@(hk) string(hk.DisplayName), h);
end
ax = ancestor(h(1),'axes');
yl = ylim(ax); 
pad = 0.1*(yl(2)-yl(1));
pad = max(pad, 50);

for k = 1:numel(h)
    x = h(k).XData; y = h(k).YData;
    m = isfinite(x) & isfinite(y); x = x(m); y = y(m);
    i = max(2, round(frac*numel(x)));
    t = text(x(i), y(i)+ over_under*pad, " "+labels(k), ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'Color', h(k).Color, ...     % <-- label color = line color
        'FontWeight','bold', ...
        'FontUnits','points', ...
        'FontName', "Times New Roman", ...
        'FontSize', 20, ... % 'EdgeColor', "black", ... % 'LineWidth', 1, ... % 'BackgroundColor', lighten(h(k).Color, 0.9), ...
        'Margin',2, ...
        'Interpreter','none', ...
        'Clipping','on');
end

end
