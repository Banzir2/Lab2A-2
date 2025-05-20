close all; clear;
C=33*1e-9;
R = [11, 9.65, 8.42, 7.43, 6.82, 5.55, 4.31, 3.219, 1.825, 1.005] * 1e3;
file_names = {'110ohm', '96ohm', '84ohm', '74ohm', '66ohm', '55ohm', '43ohm', '32ohm', '18ohm', '10ohm'};

% Preallocate storage

avg_time_constants = zeros(1, length(file_names));
avg_time_error=zeros(1, length(file_names));
for k = 1:length(file_names)
     folder = file_names{k};

    % Read files
    file1 = fullfile("PartA", folder, "meas1.csv");
    file2 = fullfile("PartA", folder, "meas2.csv");
    data1 = readmatrix(file1);
    data2 = readmatrix(file2);
    
    % Extract time and voltage
    t1 = data1(:, 1); v1 = data1(:, 2);
    t2 = data2(:, 1); v2 = data2(:, 2); scal2=data2(:, 3);


    % Select time window for charging phase
    idx_i = find(t1 > 0, 1);
    idx_f = find(t1 > 0.0035, 1);
    t1_new = t1(idx_i:idx_f);
    t2_new = t2(idx_i:idx_f);
    v1_new = v1(idx_i:idx_f);
    v2_new = v2(idx_i:idx_f);
    scal2_new=scal2(idx_i:idx_f);

   % Shift time so it starts at 0
    t2_shifted = t2_new - t2_new(1);
    
    % Correct model parameters
    v_max = max(v2_new);
    v_min = min(v2_new);
    a = v_max - v_min;  % Voltage span
    c = v_min;          % Starting voltage
    
    % Define model: a*(1 - exp(-t/b)) + c
    fitModel = fittype(@(b, x) a * (1 - exp(-x ./ b)) + c, ...
                       'independent', 'x', 'coefficients', 'b');
    
    % Fit it
    [fitResult, gof] = fit(t2_shifted, v2_new, fitModel);
    
    % Plot
    figure;
    plot(t1_new, v1_new, 'b', 'DisplayName', 'Signal generator Voltage');
    hold on;

    %this is the error form the scope
    v_2_error=scope_err(scal2_new);
    plot(t2_new, v2_new, 'r', 'DisplayName', 'Capacitor Voltage');
    
    % Fit curve
    t_fit = linspace(0, max(t2_shifted), 500);
    v_fit = a * (1 - exp(-t_fit ./ fitResult.b)) + c;
    plot(t_fit + t2_new(1), v_fit, 'k--', 'LineWidth', 2, ...
         'DisplayName', 'Fit (V = a(1 - e^{-t/b}) + c)');
    
    title(sprintf('Voltage vs Time (Charging) for R = %.2f KΩ', R(k)/1000));
    xlabel('Time [s]');
    ylabel('Voltage [V]');
    legend('Location', 'best');
    grid on;
    
    % Show result
    % Confidence interval to get error in b (95%)
    ci = confint(fitResult, 0.95);  
    b_error = (ci(2) - ci(1)) / 2;
    b_charge=num2str(fitResult.b);
    disp(['Fitted time constant b = ', num2str(fitResult.b), ' ± ', num2str(b_error), ' s']);

    % Calculate residuals

    v_fit_full = a * (1 - exp(-t2_shifted ./ fitResult.b)) + c;
    residuals = v2_new - v_fit_full;
    
    % Estimate uncertainty: assume constant error, use std of residuals
    sigma = v_2_error;
    
    % Reduced chi-squared
    N = length(v2_new);
    p = 1; % only 'b' is fitted
    chi_squared = sum((residuals ./ sigma).^2);
    chi_squared_reduced = chi_squared / (N - p);
    
    % Display
    disp(['Reduced Chi-Squared: ', num2str(chi_squared_reduced)]);




        % Select time window for discharging phase
    idx_discharge_start = find(t1 > 0.0043, 1);
    idx_discharge_end   = find(t1 > 0.008, 1);

    t1_discharge = t1(idx_discharge_start:idx_discharge_end);
    t2_discharge = t2(idx_discharge_start:idx_discharge_end);
    v1_discharge = v1(idx_discharge_start:idx_discharge_end);
    v2_discharge = v2(idx_discharge_start:idx_discharge_end);
    scale2_discharge = scal2(idx_discharge_start:idx_discharge_end);

    % Shift time so it starts at 0
    t2_discharge_shifted = t2_discharge - t2_discharge(1);

    % Extract voltage bounds
    v_max_discharge = max(v2_discharge);
    v_min_discharge = min(v2_discharge);
    a_discharge = v_max_discharge - v_min_discharge;
    c_discharge = v_min_discharge;

    % Define exponential decay model: a*exp(-t/b) + c
    fitModel_discharge = fittype(@(b, x) a_discharge * exp(-x ./ b) + c_discharge, ...
                                 'independent', 'x', 'coefficients', 'b');

    % Fit the model
    [fitResult_discharge, gof_discharge] = fit(t2_discharge_shifted, v2_discharge, fitModel_discharge);

    % Plot
    figure;
    plot(t1_discharge, v1_discharge, 'b', 'DisplayName', 'Signal generator Voltage');
    hold on;

    v2_error_discharge = scope_err(scale2_discharge);
    plot(t2_discharge, v2_discharge, 'r', 'DisplayName', 'Capacitor Voltage');

    % Plot fit
    t_fit_discharge = linspace(0, max(t2_discharge_shifted), 500);
    v_fit_discharge = a_discharge * exp(-t_fit_discharge ./ fitResult_discharge.b) + c_discharge;
    plot(t_fit_discharge + t2_discharge(1), v_fit_discharge, 'k--', 'LineWidth', 2, ...
         'DisplayName', 'Fit (V = a·e^{-t/b} + c)');

    title(sprintf('Voltage vs Time (Discharging) for R = %.2f KΩ', R(k)/1000));
    xlabel('Time [s]');
    ylabel('Voltage [V]');
    legend('Location', 'best');
    grid on;

    % Show fit result with error
    ci_discharge = confint(fitResult_discharge, 0.95);  
    b_error_discharge = (ci_discharge(2) - ci_discharge(1)) / 2;
    b_discharge= num2str(fitResult_discharge.b);
    disp(['Fitted time constant b (discharging) = ', num2str(fitResult_discharge.b), ...
          ' ± ', num2str(b_error_discharge), ' s']);

    % Residuals and chi-squared
    v_fit_full_discharge = a_discharge * exp(-t2_discharge_shifted ./ fitResult_discharge.b) + c_discharge;
    residuals_discharge = v2_discharge - v_fit_full_discharge;

    sigma_discharge = v2_error_discharge;
    N_discharge = length(v2_discharge);
    p_discharge = 1;

    chi_squared_discharge = sum((residuals_discharge ./ sigma_discharge).^2);
    chi_squared_reduced_discharge = chi_squared_discharge / (N_discharge - p_discharge);

    disp(['Reduced Chi-Squared (discharging): ', num2str(chi_squared_reduced_discharge)]);
     b_math= fitResult.b;
     b_discharge_math=fitResult_discharge.b;
     avg_b = (b_math + b_discharge_math) / 2;
     avg_time_constants(k)=avg_b;
     avg_time_error(k)=(b_error+b_error_discharge)/2;

end
% Final plot: decay time vs resistance
figure;
hold on;

% Plot error bars
errorbar(R, avg_time_constants, avg_time_error, avg_time_error, '.', ...
    'MarkerSize', 20, 'LineWidth', 1.5, 'DisplayName', 'Data with error bars');

% Linear fit: y = a*x (no intercept)
fit_linear = fit(R(:), avg_time_constants(:), 'a*x', ...
                 'Weights', 1./(avg_time_error(:)).^2);  % weighted least squares

% Evaluate fitted line
R_fit = linspace(min(R), max(R), 500);
tau_fit = fit_linear.a * R_fit;

% Plot the fit
plot(R_fit, tau_fit, 'r--', 'LineWidth', 2, 'DisplayName', sprintf('Fit: y = %.3g·x', fit_linear.a));

% Add labels and title
xlabel('Resistance [\Omega]', 'FontSize', 16);
ylabel('Decay Time [s]', 'FontSize', 16);
title('Decay Time vs Resistance with Linear Fit', 'FontSize', 16);
legend('Location', 'northwest');
grid on;

% Improve axis appearance
ax = gca;
ax.FontSize = 14;

% --- Chi-squared calculation for the fit y = a*x
y_model = fit_linear.a * R;
residuals = avg_time_constants - y_model;
chi_squared = sum((residuals ./ avg_time_error).^2);
reduced_chi_squared = chi_squared / (length(R) - 1);  % 1 parameter: 'a'

% --- Get fit parameter and its error
ci = confint(fit_linear, 0.95); % 95% confidence interval
a_val = fit_linear.a;
a_err = (ci(2) - ci(1)) / 2;

% Display in Command Window
disp(['Fitted slope a = ', num2str(a_val), ' ± ', num2str(a_err), ' s/Ω']);
disp(['Reduced Chi-Squared = ', num2str(reduced_chi_squared)]);

