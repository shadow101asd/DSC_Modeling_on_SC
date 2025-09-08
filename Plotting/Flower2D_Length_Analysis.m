
% Parameters

mu_dummy = 1;
etR = 0;

W_lim = 100;
Ws = -W_lim:W_lim;
F_lim = 100;
Fs = -F_lim:F_lim;

resolution = 1e4;
a_norm = 1;
eF = 0.5;
Ki = [a_norm, eF, 0, 0, pi/2, 0];

Ls = zeros([2*W_lim+1 2*F_lim+1]);

%% Run
for nW = 1:length(Ws)
    W = Ws(nW);
    parfor nF = 1:length(Fs)
        F = Fs(nF);
        XSats = NSATSpropagateFromKeplerians_2DFlower(Ki,resolution,W,F,etR,mu_dummy);
        XSats(:,:,resolution+1) = XSats(:,:,1); % Close loop for plotting and analysis
        Ls(nW,nF) = sum(sqrt(sum(diff(XSats(1:2,1,:),2).^2, 2)));

        % Display progress
        % disp("Analyses run: " + num2str((nW-1)*length(Fs) + nF) + "/" + num2str(length(Ws)*length(Fs)))
    end
    disp("Analyses run: " + num2str(nW) + "/" + num2str(length(Ws)))
end

%% Plotting

[X,Y] = meshgrid(Ws, Fs);
x = [min(X) max(X)];
y = [min(Y) max(Y)];

figure(12121)
% contourf(X,Y,Ls)
h1 = imagesc(x,y,Ls);
set(gca, 'YDir','normal');
c = colorbar;
ylabel(c,'Constellation Path Length [NDU]', FontSize=12);
xlabel("W")
ylabel("F")
