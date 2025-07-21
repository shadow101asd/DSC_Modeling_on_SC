%% Collect Data from DSC_Modular_RF_2P Runs

run_idx = '014';
Nf_range = 1:2;
Nl_range = 1:100;

Bs = NaN([length(Nl_range), length(Nf_range)]);
for Nf = Nf_range
    for Nl = Nl_range
        filename = "../Data/run"+run_idx+"/Nf"+int2str(Nf)+"Nl"+int2str(Nl)+".mat";
        
        if isfile(filename)
            % Extract and save data of interest
            Bs(Nl, Nf) = load(filename, "B").B;
        else
            warning("Missing file: " + filename)
        end
    end
end

%% Plotting

figure(str2double(run_idx))
plot(Nl_range, Bs)

lgd = [];
for Nf = Nf_range
    lgd = [lgd, "Nf = " + int2str(Nf)];
end
legend(lgd)
ylim([0 ceil(max(max(Bs))/100)*100])
xlabel("# of SS Launches")
ylabel("Datarate [Mbps]")
title("Run " + run_idx + ": Ongoing")