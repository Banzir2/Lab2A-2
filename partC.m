close all; clear;

RminPot = 0.2;
RmaxPot = 98.6;
R100 = 99.1;
R220 = 219.1;
R470 = 464;

Vb = 1.181;
Vb_error=Hioki(1.181,0,0);
V=[4, 269.7 , 438, 554, 638, 703, 703, 753, 795, 844, 872, 895, 918,957 , 986, 1027, 1042, 1064, 1083, 1095]/1000;
I=[17.53, 13.59, 11.07, 9.36, 8.1, 7.14, 7.08,6.33, 5.73, 5.34, 4.9, 4.53, 4.19,3.56, 3.11, 2.46, 2.24, 1.88, 1.57, 1.39 ]/1000;
R=[0.2, 19.8, 39.6, 59.3, 78.9, 98.6, 99.4, 119.2, 139, 158.7, 178.5, 197.7, 219.1, 268.8, 317.6, 415, 464, 563, 694, 780];
dR=Beckman(R,0);
dV=Hioki(V,0,0);
dI=Beckman(I,0);


% Calculate theoretical voltage
V_theory = I .* R;

% Plot
figure;
errorbar(R, V,dV,dV,dR,dR ,'bo', 'LineWidth', 1.5, 'DisplayName', 'Measured V','MarkerSize',5);
hold on;
plot(R, V_theory, 'ro', 'LineWidth', 1.5, 'DisplayName', 'Theoretical V = I·R', 'MarkerSize', 5);

xlabel('R [\Omega]', 'FontSize', 12);
ylabel('V [V]', 'FontSize', 12);
title('Voltage vs Resistance', 'FontSize', 14);
legend('Location', 'best');
grid on;
hold on;
% Fitting V = a*R / (R + b)
fitType = fittype('a*x / (x + b)', 'independent', 'x', 'coefficients', {'a', 'b'});
[fitResult, gof] = fit(R.', V.', fitType, 'StartPoint', [1, 1]);
coeffs = coeffvalues(fitResult);         % Extract a and b
confInt = confint(fitResult);            % 95% confidence intervals

% Calculate standard errors (half the width of confidence interval)
a = coeffs(1);
b = coeffs(2);
a_err = (confInt(2,1) - confInt(1,1)) / 2;
b_err = (confInt(2,2) - confInt(1,2)) / 2;

[f,Y,Z]=fit(R.', V.', fitType);
diff = V - f(R).';
xhisquare = sum((diff.^2) ./ ((dR.*I).^2 + dV.^2)) / 18;
% Add fit curve
R_fit = linspace(min(R), max(R), 500);
V_fit = fitResult.a * R_fit ./ (R_fit + fitResult.b);
plot(R_fit, V_fit, 'k--', 'LineWidth', 1,'HandleVisibility','off');

xlabel('R [Omega]', 'FontSize', 12);
ylabel('V [V]', 'FontSize', 12);
title('Voltage vs Resistance', 'FontSize', 14);
legend('Location', 'best');
grid on;
% Plot2
diffV=Vb-V;
figure;
[f, X, O] = fit(I.', diffV.', 'p1*x');
diff = diffV - f(I).';
negXhisquare = sum((diff.^2) ./ ((dI.*R).^2 + (dV+Vb_error).^2)) / 19;
p1 = f.p1;  % Slope
%p2 = f.p2;  % Intercept
%ci = confint(f);           % 2x2 matrix: [lower; upper]
%p1_error = (ci(2,1) - ci(1,1)) / 2;
%p2_error = (ci(2,2) - ci(1,2)) / 2;
ci = confint(f);                % 2x1 matrix: [lower; upper]
p1_error = (ci(2,1) - ci(1,1)) / 2;
plot(f);
hold on;
errorbar(I, diffV,dV+Vb_error,dV+Vb_error,dI,dI ,'bo', 'LineWidth', 1.5, 'DisplayName', 'Measured V','MarkerSize',5);
xlabel('I [A]', 'FontSize', 12);
ylabel('V [V]', 'FontSize', 12);
title('difference in voltage vs current', 'FontSize', 14);
legend('Location', 'best');
grid on;
hold on;
