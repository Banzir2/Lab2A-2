close all; clear;

R22Mblack = 25.09 * 1e6;
R22Mred = 26.37 * 1e6;
R100K = 100 * 1e3;
Vb = 8.98;
R = [R22Mblack, R22Mred, R22Mblack + R22Mred];

V1 = [37.2, 36.1, 17.7] / 1000;
V1b = 106.5 / 1000;
Vxb = 8.88;

r = (R100K*44e6*Vxb) / (V1b*44e6-R100K*Vxb);

[dV, dR] = Hioki(V1, R, 0);
figure; hold on;
errorbar(R, V1, 3*dV, 3*dV, 12*dR, 12*dR, '.', 'MarkerSize', 20);
legend('Data', 'FontSize', 16);
xlabel('Resistance [ohm]', 'FontSize', 16);
ylabel('Voltage [volt]', 'FontSize', 16);
title('R1 voltage by Rx resistance', 'FontSize', 16);
ax = gca;
ax.FontSize = 14;