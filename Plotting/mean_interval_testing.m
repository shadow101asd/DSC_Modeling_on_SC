
[~,~,~,Dpaths_fast,Dequivs_fast] = networkAnalysis_FAST(XEa,XMa,"Earth","Mars",XSats_opt, 1.2);

dts = 1:50;
means = zeros([1 length(dts)]);

for i = 1:length(dts)
    dt = dts(i);
    means(i) = mean(Dequivs_fast(1:dt:end))/AU;
end

figure(2)
plot(dts, (means-means(1))/means(1) * 100)
xlabel("interval skip")
ylabel("% change in mean")