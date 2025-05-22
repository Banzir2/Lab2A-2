close all; clear;

freq_10=[32.98, 36.5, 30.05];
Amplitude_10ohm=[18.6,6.24,6.24];
Amplitude_otot10ohm=[1.1, 1.26, 1.28];

freq_100=[32.97,37.23,29.08];
Amplitude_otot100ohm=[1.14,1.26, 1.29 ];
Amplitude_100ohm=[14.8, 5.04, 4.88];

freq_1k=[33.2,44.49, 13.93];
Amplitude_otot1kohm=[1.26,1.3, 1.3];
Amplitude_1kohm=[4.56, 1.52,1.52];

freq_220=[33.07,38.45, 27.74];
Amplitude_otot10kohm=[1.18,1.28, 1.29];
Amplitude_10kohm=[11.6, 3.8,3.84];

freq_420=[35.11, 37.37, 39.25, 43.18, 46.38, 33.16, 31.29, 28.03, 24.34, 20.18, 83.19, 168.4, 10.86, 5.215]*1e3;

file_names = ["35k", "37k", "39k", "43k", "46k", "33k", "31k", "28k", "24k", "20k", "83k", "168k", "10k", "5k"];

ampl_diff = zeros(size(file_names));
phase_diff = zeros(size(file_names));
err_amp = zeros(size(file_names));
err_phase = zeros(size(file_names));
for i = 1:size(file_names, 2)
    data1 = readmatrix(strcat('470ohm/', file_names(i), '/meas1.csv'));
    data2 = readmatrix(strcat('470ohm/', file_names(i), '/meas2.csv'));

    % Model function to fit: y = a * exp(b * x)
    freq = freq_420(i);
    model = @(params, x) params(1) * cos(2*pi*freq * x + params(2));

    lb = [0, -pi];
    ub = [25, pi];

    [paramsFitted1, ~, residual, ~, ~, ~, jacobian] = lsqcurvefit(model, [0.4, 0], data1(:, 1), data1(:, 2), lb, ub);
    var_res = sum(residual.^2) / (length(residual) - length(paramsFitted1));
    cov_matrix = var_res * inv(jacobian' * jacobian);
    std_errors_1 = full(sqrt(diag(cov_matrix)));
    [paramsFitted2, ~, residual, ~, ~, ~, jacobian] = lsqcurvefit(model, [0.02, 0], data2(:, 1), data2(:, 2), [0, -pi], ub);
    var_res = sum(residual.^2) / (length(residual) - length(paramsFitted2));
    cov_matrix = var_res * inv(jacobian' * jacobian);
    std_errors_2 = full(sqrt(diag(cov_matrix)));

    ampl_diff(i) = abs(paramsFitted2(1) / paramsFitted1(1));
    phase_diff(i) = mod(paramsFitted1(2) - paramsFitted2(2), pi);
    err_amp(i) = scope_err(data1(1, 3)) + scope_err(data2(1, 3));
    err_phase(i) = std_errors_1(2) + std_errors_2(2);
end

clear alpha beta cov_matrix delta dx Fs gamma i jacobian k L lb model N P paramsFitted1 paramsFitted2 refinedFreq reifnedIndex residual std_errors_1 std_errors_2 ub var_res w xdata Y y_detrended y_windowed ydata B0 data1 data2

function r = weighted_residuals(model, x, xdata, ydata, weights)
    r = weights .* (model(x, xdata) - ydata);
end

xhisquare = zeros(1, 2);

figure; hold on; % Amp diff
errorbar(freq_420, ampl_diff, err_amp, '.', 'MarkerSize', 25);

modelFun = @(params, x) 1 ./ (sqrt(params(1)*(x.^2) + (params(2)*(x.^2) - 1).^2));
startPoint = [(740*1e-9)^2, (23e-3 * 1e-9)];
xData = 2*pi*freq_420(:);     % Ensure column vector
yData = ampl_diff(:);
lb = [0, 0];
ub = [Inf, Inf];
options = optimoptions('lsqcurvefit', 'TolFun', 1e-8);
params_fit = lsqcurvefit(modelFun, startPoint, xData, yData, lb, ub, options);

scatter(freq_420, modelFun(params_fit, 2*pi*freq_420));
diff2 = (ampl_diff - modelFun(params_fit, 2*pi*freq_420)).^2;
xhisquare(1) = sum(diff2 ./ (err_amp.^2 + (pi*freq_420*1e-6).^2)) / 12;
xlabel('Frequency [Hz]', 'FontSize', 16);
ylabel('Amplitude ratio', 'FontSize', 16);
title('Amplitude ratio by input frequency', 'FontSize', 16);
ax = gca;
ax.FontSize = 14;

figure; hold on; % Phase diff
errorbar(freq_420, phase_diff - pi/2, err_phase * 30, '.', 'MarkerSize', 25);
f = fit(freq_420.', phase_diff.' - pi/2, 'atan(a*x+b/x)', 'StartPoint', [25e-3 / 740, -1/(740*1e-9)]);
scatter(freq_420, f(freq_420))
diff2 = (phase_diff - f(freq_420).' - pi/2).^2;
xhisquare(2) = sum(diff2 ./ (err_phase.^2 + (pi*freq_420*1e-6).^2)) / 10;
xlabel('Frequency [Hz]', 'FontSize', 16);
ylabel('Phase diff [rad]', 'FontSize', 16);
title('Phase diff by input frequency', 'FontSize', 16);
ax = gca;
ax.FontSize = 14;