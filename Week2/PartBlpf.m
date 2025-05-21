close all; clear;
C_33n=37.45*1e-9;
L=2.152/1000;
R_100k=100.1*1e3;
R_1k=999;
freq_lpf=[2.290,10.47,43.92,106.8,267.4,733.7,1.430*1e3,5.123*1e3,15.85*1e3, 54.35*1e3,106.4*1e3, 191.3*1e3,230.0*1e3];
file_names = ["2-H", "10-H", "43-H", "106-H", "267-H", "733-H", "1.4k-H", "5k-H", "15k-H", "54k-H", "106k-H", "191k-H", "230k-H"];

ampl_diff = zeros(size(file_names));
phase_diff = zeros(size(file_names));
err_amp = zeros(size(file_names));
err_phase = zeros(size(file_names));
for i = 1:size(file_names, 2)
    data1 = readmatrix(strcat('PartB-LowPass/', file_names(i), '/meas1.csv'));
    data2 = readmatrix(strcat('PartB-LowPass/', file_names(i), '/meas2.csv'));
    
    % Model function to fit: y = a * exp(b * x)
    model = @(params, x) params(1) * cos(2*pi*freq_lpf(i) * x + params(2));

    lb = [0, -pi];
    ub = [Inf, pi];

    [paramsFitted1, ~, residual, ~, ~, ~, jacobian] = lsqcurvefit(model, [0.4, 0], data1(:, 1), data1(:, 2), lb, ub);
    var_res = sum(residual.^2) / (length(residual) - length(paramsFitted1));
    cov_matrix = var_res * inv(jacobian' * jacobian);
    std_errors_1 = full(sqrt(diag(cov_matrix)));
    [paramsFitted2, ~, residual, ~, ~, ~, jacobian] = lsqcurvefit(model, [0.02, 0], data2(:, 1), data2(:, 2), [0, -pi], ub);
    var_res = sum(residual.^2) / (length(residual) - length(paramsFitted2));
    cov_matrix = var_res * inv(jacobian' * jacobian);
    std_errors_2 = full(sqrt(diag(cov_matrix)));

    ampl_diff(i) = abs(paramsFitted2(1) / paramsFitted1(1));
    phase_diff(i) = mod(paramsFitted2(2) - paramsFitted1(2), pi);
    err_amp(i) = scope_err(data1(1, 3)) + scope_err(data2(1, 3));
    err_phase(i) = std_errors_1(2) + std_errors_2(2);
end

for i = 1:size(phase_diff, 2)
    if phase_diff(i) > 1
        phase_diff(i) = phase_diff(i) - pi;
    end
end

clear alpha beta C_33n cov_matrix delta dx Fs gamma i jacobian k lb model N P paramsFitted1 paramsFitted2 refinedFreq reifnedIndex residual std_errors_1 std_errors_2 ub var_res w xdata Y y_detrended y_windowed ydata B0 data1 data2

function r = weighted_residuals(model, x, xdata, ydata, weights)
    r = weights .* (model(x, xdata) - ydata);
end

xhisquare = zeros(1, 2);
startpoint = [R_1k / (2*pi*L), 0];
model_amp = @(params, x) (1 ./ sqrt((x / params(1)).^2 + 1));
model_phase = @(params, x) -atan(x / params(1));
lb = [0, -Inf];
ub = [Inf, Inf];

fun = @(params, x) weighted_residuals(model_amp, params, freq_lpf, ampl_diff, 1 ./ sqrt(err_amp.^2 + (pi*freq_lpf*1e-6)));
[paramsFittedAmp, ~, residual, ~, ~, ~, jacobian] = lsqcurvefit(fun, startpoint, freq_lpf, zeros(size(ampl_diff)), lb, ub);
var_res = sum(residual.^2) / (length(residual) - length(paramsFittedAmp));
cov_matrix = var_res * inv(jacobian' * jacobian);
std_errors_amp = full(sqrt(diag(cov_matrix)));

fun = @(params, x) weighted_residuals(model_phase, params, freq_lpf, phase_diff, 1 ./ sqrt(err_phase.^2 + (pi*freq_lpf*1e-4)));
[paramsFittedPhase, ~, residual, ~, ~, ~, jacobian] = lsqcurvefit(fun, startpoint, freq_lpf, zeros(size(phase_diff)), lb, ub);
var_res = sum(residual.^2) / (length(residual) - length(paramsFittedPhase));
cov_matrix = var_res * inv(jacobian' * jacobian);
std_errors_phase = full(sqrt(diag(cov_matrix)));

figure; hold on; % Amp ratio
errorbar(log(freq_lpf), ampl_diff, err_amp * 4, '.', 'MarkerSize', 25);
diff2 = (ampl_diff - model_amp(paramsFittedAmp, freq_lpf)).^2;
xhisquare(1) = sum(diff2 ./ (err_amp.^2 + (pi*freq_lpf*1e-6).^2)) / 10;
%plot(log(freq_lpf), model_amp(paramsFittedAmp, freq_lpf));
xlabel('Log(Frequency)', 'FontSize', 16);
ylabel('Amplitude ratio', 'FontSize', 16);
title('Amplitude ratio by input frequency 1K\Omega', 'FontSize', 16);
ax = gca;
ax.FontSize = 14;

figure; hold on; % Phase diff
errorbar(log(freq_lpf), phase_diff, err_phase * 15, '.', 'MarkerSize', 25);
diff2 = (phase_diff - model_phase(paramsFittedPhase, freq_lpf)).^2;
xhisquare(2) = sum(diff2 ./ (err_phase.^2 + (pi*freq_lpf*1e-4).^2)) / 10;
%plot(log(freq_lpf), model_phase(paramsFittedPhase, freq_lpf));
xlabel('Log(Frequency)', 'FontSize', 16);
ylabel('Phase diff [rad]', 'FontSize', 16);
title('Phase diff by input frequency 1K\Omega', 'FontSize', 16);
ax = gca;
ax.FontSize = 14;

abs(paramsFittedAmp(1) - startpoint(1)) / std_errors_amp(1)
abs(paramsFittedPhase(1) - startpoint(1)) / std_errors_phase(1)