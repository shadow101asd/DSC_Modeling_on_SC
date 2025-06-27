%% Collect Data from DSC_Modular_RF_2P Runs

run_idx = '010';
Nf_range = 1:4;
Nl_range = 1:60;

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

figure(12)
plot(Nl_range, Bs)
legend("Nf = 1", "Nf = 2", "Nf = 3", "Nf = 4")