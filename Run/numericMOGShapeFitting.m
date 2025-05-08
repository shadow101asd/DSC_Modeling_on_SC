function [fit_rho, fit_theta, gof_rho, gof_theta, fake_t] = numericMOGShapeFitting(eMOG, fourierDegree, resolution)
    
    if nargin == 2
        resolution = 1e4;
    end
    % Form MOG
    [theta, rho] = getPolarRepresentationOfMOG(1e8, eMOG, 1, resolution); % rho is in arbitrary units

    % Let's normalize rho such that the min distance is 1
    rho = rho/min(rho);

    % get rid of discontinuity in theta
    [~, idx] = max(theta);
    theta(idx+1:end) = theta(idx+1:end) + 2*pi;
    
    % Fitting
    fit_type = "fourier" + int2str(fourierDegree);
    fake_t = linspace(0,1,length(rho))';

    [fit_rho, gof_rho] = fit(fake_t, rho, fit_type);

    [fit_theta_temp, gof_theta] = fit(fake_t, theta - linspace(0,2*pi,length(theta))', fit_type);

    % Incorporate linear term into the fit

    new_form = formula(fit_theta_temp) + " + z*x";
    % Formatting
    new_form = regexprep(new_form,'\n+','');
    new_form = regexprep(new_form,'  ',' ');
    new_form = regexprep(new_form,'  ',' ');
    new_form = regexprep(new_form,'  ',' ');

    % Get alphabetical order
    cnames = coeffnames(fit_theta_temp);
    [~, I] = sort(cnames);
    
    coeffs = coeffvalues(fit_theta_temp);
    coeffs = coeffs(I);
    new_coeffs = [coeffs, 2*pi];
    new_coeffs = num2cell(new_coeffs);

    fit_theta = cfit(fittype(new_form), new_coeffs{:});
end