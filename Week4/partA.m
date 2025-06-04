close all; clear;

%real vlaues
R1=219.65;
R2=218.15;
Ctag=[944.1e-12,2.026e-9 ,94.76e-9];
C1=206.5e-12;
C2=202.2e-12;
L1=121.78e-3;
L2=98.52e-3;

Z1 = @(params, f) params(1) + (1i*2*pi*f*params(2)) + 1 ./ (1i*2*pi*f*params(3)) + 1 ./ (1i*2*pi*f*Ctag(1));
Z2 = @(params, f) params(1) + (1i*2*pi*f*params(2)) + 1 ./ (1i*2*pi*f*params(3)) + 1 ./ (1i*2*pi*f*Ctag(1));
Zm = @(f) 1 ./ (1i*2*pi*f*Ctag(1));
model_R1 = @(params, f) params(1) * abs((params(4)*Z2(params(5:7), f)) ./ (Z1(params(1:3), f).*Z2(params(5:7), f) - Zm(f).^2));
model_R2 = @(params, f) params(5) * abs((params(4)*Zm(f)) ./ (Z1(params(1:3), f).*Z2(params(5:7), f) - Zm(f).^2));

%1nF
f=[39.817,41.142, 42.11,44.07, 47.09, 38.98, 37.134, 34.928, 32.07, 33.62, 30.425, 20.00, 56.72, 105.63, 4.931]*1e3;
V1=[1.22,0.84, 0.612, 0.388, 0.253,1.08, 0.26, 0.508, 0.584, 1.01, 0.354, 0.092, 0.126, 0.0268, 0.0156];
V2=[0.812,0.448, 0.282, 0.136, 0.06, 0.92, 0.592, 0.82, 0.38, 1.04, 0.176, 0.024, 0.016, 0.068, 0.0038];

f_space = linspace(min(f), max(f), 100000);

[f, newIdx] = sort(f);
V1 = V1(newIdx);
V2 = V2(newIdx);

figure; hold on;
% plot(f, V1);
% plot(f, V2);

startPoint = [R1, L1, C1, 1.3, R2, L2, C2];
xData = f(:);     % Ensure column vector
yData = V1(:);
lb = zeros(size(startPoint));
ub = ones(size(startPoint)) * Inf;
[params_fit, ~, residual, ~, ~, ~, jacobian] = lsqcurvefit(model_R1, startPoint, xData, yData, lb, ub);
var_res = sum(residual.^2) / (length(residual) - length(params_fit));
cov_matrix = var_res * inv(jacobian' * jacobian);
std_errors = full(sqrt(diag(cov_matrix)));

plot(f_space, model_R1(params_fit, f_space));
plot(xData, yData);

%2.2nF
f=[33.65,34.63, 36.71, 37.645, 35.62, 34.11, 33.08, 36.00,37.17 , 43.36, 26.88,16.666, 52.998 ]*1e3;
V1=[1.1, 0.760, 1.58, 1.24, 1, 0.92, 1.06, 1.3, 1.46, 0.308, 0.22, 0.0748, 0.0138];
V2=[1,0.96, 0.94, 0.62, 1, 1, 0.72, 1.02, 0.8, 0.056, 0.048, 0.0108, 0.014];

[f, newIdx] = sort(f);
V1 = V1(newIdx);
V2 = V2(newIdx);

figure; hold on;
plot(f, V1);
plot(f, V2);

%100nF
f=[ 33.33, 33.11, 33.2, 33.4 ,33.54, 33.6, 33.71, 33.8, 34.11, 34.3, 32.97,33.065, 34.57, 34.9, 36.22, 40.16, 32.01, 31.26];
V1=[ 2.4, 2.06,2.18, 2.58, 2.94, 3.00, 3.2, 3.28, 3.14, 2.8, 1.84, 2, 2.4, 1.92, 0.96, 0.372, 1, 0.728];
V2=[ 0.08, 0.08,0.08, 0.08, 0.2, 0.2, 0.24, 0.24, 0.24, 0.24, 0.08, 0.08, 0.24, 0.12, 0.04 , 0.016,0.04, 0.028];

[f, newIdx] = sort(f);
V1 = V1(newIdx);
V2 = V2(newIdx);

figure; hold on;
plot(f, V1);
plot(f, V2);
