close all; clear;
C = 33e-9; % Capacitance in Farads
R = [11, 9.65, 8.42, 7.43, 6.82, 5.55, 4.31, 3.219, 1.825, 1.005] * 1e3;
file_names = {'10ohm', '110ohm', '18ohm', '32ohm', '43ohm', '55ohm', '66ohm', '74ohm', '84ohm', '96ohm'};

% Preallocate storage
avg_time_constants = zeros(1, length(file_names));
avg_chi_squareds = zeros(1, length(file_names));

for k = 1:length(file_names)
    folder = file_names{k};

    % Read files
    file1 = fullfile("PartA", folder, "meas1.csv");
    file2 = fullfile("PartA", folder, "meas2.csv");
    data1 = readmatrix(file1);
    data2 = readmatrix(file2);

    t1 = data1(:, 1); v1 = data1(:, 2);
    t2 = data2(:, 1); v2 = data2(:, 2); scal2 = data2(:, 3);

    % Charging phase
    idx_i = find(t1 > 0, 1);
    idx_f = find(t1 > 0.0035, 1);
    t2_new = t2(idx_i:idx_f);
    v2_new = v2(idx_i:idx_f);
    scal2_new = scal2(idx_i:idx_f);
    t2_shifted = t2_new - t2_new(1);

    v_max = max(v2_new); v_min = min(v2_new);
    a = v_max - v_min; c = v_min;

    fitModel = fittype(@(b, x) a * (1 - exp(-x ./ b)) + c, ...
        'independent', 'x', 'coefficients', 'b');

    [fitResult, ~] = fit(t2_shifted, v2_new, fitModel);
    b_charge = fitResult.b;

    % Charging residuals
    v_fit_full = a * (1 - exp(-t2_shifted ./ b_charge)) + c;
    residuals = v2_new - v_fit_full;
    sigma = scope_err(scal2_new);
    chi_sq_charge = sum((residuals ./ sigma).^2) / (length(v2_new) - 1);

    % Discharging phase
    idx_discharge_start = find(t1 > 0.0043, 1);
    idx_discharge_end   = find(t1 > 0.008, 1);
    t2_discharge = t2(idx_discharge_start:idx_discharge_end);
    v2_discharge = v2(idx_discharge_start:idx_discharge_end);
    scale2_discharge = scal2(idx_discharge_start:idx_discharge_end);
    t2_discharge_shifted = t2_discharge - t2_discharge(1);

    v_max_d = max(v2_discharge); v_min_d = min(v2_discharge);
    a_d = v_max_d - v_min_d; c_d = v_min_d;

    fitModel_d = fittype(@(b, x) a_d * exp(-x ./ b) + c_d, ...
        'independent', 'x', 'coefficients', 'b');

    [fitResult_d, ~] = fit(t2_discharge_shifted, v2_discharge, fitModel_d);
    b_discharge = fitResult_d.b;

    % Discharging residuals
    v_fit_full_d = a_d * exp(-t2_discharge_shifted ./ b_discharge) + c_d;
    residuals_d = v2_discharge - v_fit_full_d;
    sigma_d = scope_err(scale2_discharge);
    chi_sq_discharge = sum((residuals_d ./ sigma_d).^2) / (length(v2_discharge) - 1);

    % Average time constant and residual
    avg_b = (b_charge + b_discharge) / 2;
    avg_chi_sq = (chi_sq_charge + chi_sq_discharge) / 2;

    avg_time_constants(k) = avg_b;
    avg_chi_squareds(k) = avg_chi_sq;
end

% Plot average time constant vs residual
figure;
plot(avg_time_constants, R, 'o-', 'LineWidth', 2);
xlabel('Reduced Chi-Squared (Average)');
ylabel('Average Time Constant [s]');
title('Average Time Constant vs. Residual');
grid on;
