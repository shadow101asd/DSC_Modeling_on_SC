% Plotting an example FC for the Parametrization section

% Add paths

addpath("../../Run/") % Run folder
addpath("../") % Main Plotting scripts

% Constants

AU = 149600000; % km
muSu = 1.327124400419393e+11;

% Parameters

Nws = [1 1 1 2 3 4 5 1 5  1  1  3];
NMs = [1 2 3 3 2 5 4 5 1 -1 -2 -2];
K = length(Nws);

aF = AU;
eF = 0.35;
w0F = pi/2;
f0F = 0;

etR = 0;
res_outline = 1000;

X_outlines = zeros(6,K,res_outline-1);

for i = 1:K
    X_outlines(:,i,:) = generate2DP_Flower(aF,eF,w0F,f0F,Nws(i),-NMs(i),etR,muSu,res_outline-1);
end
X_outlines(:,:,end+1) = X_outlines(:,:,1); % close loop visually

%% Plotting

lim = 1.05*aF*(1+eF)/AU;
n = 44552448;
h = 300; % pixels
w = 300; % pixels
screen_h = 700;
screen_w = 1500; % pixels


fs = [];
for i = 1:K
    fs(i) = figure(n+i);
    clf(fs(i))
    scatter(0,0,SizeData=150, MarkerEdgeColor="black", MarkerFaceColor="yellow", Marker="pentagram", LineWidth=0.75);
    hold on
    plot(permute(X_outlines(1,i,:)/AU, [3,1,2]),permute(X_outlines(2,i,:)/AU, [3,1,2]), ...
        LineStyle="-", LineWidth=2, Color="black");
    viscircles([0,0], lim, Color="white")
    hold off
    axis equal
    axis off
    xlim([-lim lim])
    ylim([-lim lim])
    set(fs(i), Position=[(mod(w*(i-1),screen_w)) (screen_h-h*floor(w*(i-1)/screen_w)) w h])
    % hold on
    % annotation('rectangle',[0 0 1 1],'Color','white');
    % hold off
end


%% Saving

for i = 1:K
    filename = "./Figures/appendix_figs/FC_apdx_" + int2str(i) + ".eps";
    exportgraphics(figure(fs(i)), filename, BackgroundColor="white")
end