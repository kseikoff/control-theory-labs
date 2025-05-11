% plant parameters
A=[0 1; -1 2];
B=[1 2; -3 3];
C=[2 1; 3 -2];
D=[0 0; 0 0];

% A eig
A_eig = eig(A)

% W(s)
sys = ss(A, B, C, D);
W_s = tf(sys)

% zeros
zeros = zero(sys)

% poles
poles = pole(sys)

% Jordan matrix
[P, J] = jordan(A)
B_J = inv(P) * B
C_J = C * P

% out
U = [B A*B];
U_out = [C*U D]
rank(U_out)

% time
t = 0:0.01:5;

% w(t)
w1 = -exp(t) - 12*t.*exp(t);
w2 = 7*exp(t) + 3*t.*exp(t);
w3 = 9*exp(t) - 4*t.*exp(t);
w4 = t.*exp(t);

% h(t)
h1 = -12*t.*exp(t) + 11*exp(t) - 11;
h2 = 3*t.*exp(t) + 4*exp(t) - 4;
h3 = -4*t.*exp(t) + 13*exp(t) - 13;
h4 = t.*exp(t) - exp(t) + 1;

% w(t) renders
figure;
subplot(2,2,1)
plot(t, w1, 'b', 'LineWidth', 1.5); grid on;
xlabel('Time (s)', 'Interpreter','latex', 'FontSize', 12);
ylabel('$w_1(t)$', 'Interpreter','latex', 'FontSize', 12);
title('Weight Function $w_1(t)$', 'Interpreter','latex', 'FontSize', 14);

subplot(2,2,2)
plot(t, w2, 'g', 'LineWidth', 1.5); grid on;
xlabel('Time (s)', 'Interpreter','latex', 'FontSize', 12);
ylabel('$w_2(t)$', 'Interpreter','latex', 'FontSize', 12);
title('Weight Function $w_2(t)$', 'Interpreter','latex', 'FontSize', 14);

subplot(2,2,3)
plot(t, w3, 'r', 'LineWidth', 1.5); grid on;
xlabel('Time (s)', 'Interpreter','latex', 'FontSize', 12);
ylabel('$w_3(t)$', 'Interpreter','latex', 'FontSize', 12);
title('Weight Function $w_3(t)$', 'Interpreter','latex', 'FontSize', 14);

subplot(2,2,4)
plot(t, w4, 'm', 'LineWidth', 1.5); grid on;
xlabel('Time (s)', 'Interpreter','latex', 'FontSize', 12);
ylabel('$w_4(t)$', 'Interpreter','latex', 'FontSize', 12);
title('Weight Function $w_4(t)$', 'Interpreter','latex', 'FontSize', 14);

% h(t) renders
figure;
subplot(2,2,1)
plot(t, h1, 'b', 'LineWidth', 1.5); grid on;
xlabel('Time (s)', 'Interpreter','latex', 'FontSize', 12);
ylabel('$h_1(t)$', 'Interpreter','latex', 'FontSize', 12);
title('Step Response $h_1(t)$', 'Interpreter','latex', 'FontSize', 14);

subplot(2,2,2)
plot(t, h2, 'g', 'LineWidth', 1.5); grid on;
xlabel('Time (s)', 'Interpreter','latex', 'FontSize', 12);
ylabel('$h_2(t)$', 'Interpreter','latex', 'FontSize', 12);
title('Step Response $h_2(t)$', 'Interpreter','latex', 'FontSize', 14);

subplot(2,2,3)
plot(t, h3, 'r', 'LineWidth', 1.5); grid on;
xlabel('Time (s)', 'Interpreter','latex', 'FontSize', 12);
ylabel('$h_3(t)$', 'Interpreter','latex', 'FontSize', 12);
title('Step Response $h_3(t)$', 'Interpreter','latex', 'FontSize', 14);

subplot(2,2,4)
plot(t, h4, 'm', 'LineWidth', 1.5); grid on;
xlabel('Time (s)', 'Interpreter','latex', 'FontSize', 12);
ylabel('$h_4(t)$', 'Interpreter','latex', 'FontSize', 12);
title('Step Response $h_4(t)$', 'Interpreter','latex', 'FontSize', 14);

% freq
omega = logspace(-2, 2, 1000);

% A(w)
A1 = sqrt(omega.^6 + 123*omega.^4 + 243*omega.^2 + 121) ./ (omega.^2 + 1).^2;
A2 = sqrt(49*omega.^6 + 114*omega.^4 + 81*omega.^2 + 16) ./ (omega.^2 + 1).^2;
A3 = sqrt(81*omega.^6 + 331*omega.^4 + 419*omega.^2 + 169) ./ (omega.^2 + 1).^2;
A4 = 1 ./ (omega.^2 + 1);

% A(w) renders
figure;

subplot(2,2,1)
semilogx(omega, A1, 'b', 'LineWidth', 1.5); grid on;
xlabel('Frequency $\omega$ (rad/s)', 'Interpreter','latex', 'FontSize', 12);
ylabel('$A_1(\omega)$', 'Interpreter','latex', 'FontSize', 12);
title('Amplitude Response $A_1(\omega)$', 'Interpreter','latex', 'FontSize', 14);

subplot(2,2,2)
semilogx(omega, A2, 'g', 'LineWidth', 1.5); grid on;
xlabel('Frequency $\omega$ (rad/s)', 'Interpreter','latex', 'FontSize', 12);
ylabel('$A_2(\omega)$', 'Interpreter','latex', 'FontSize', 12);
title('Amplitude Response $A_2(\omega)$', 'Interpreter','latex', 'FontSize', 14);

subplot(2,2,3)
semilogx(omega, A3, 'r', 'LineWidth', 1.5); grid on;
xlabel('Frequency $\omega$ (rad/s)', 'Interpreter','latex', 'FontSize', 12);
ylabel('$A_3(\omega)$', 'Interpreter','latex', 'FontSize', 12);
title('Amplitude Response $A_3(\omega)$', 'Interpreter','latex', 'FontSize', 14);

subplot(2,2,4)
semilogx(omega, A4, 'm', 'LineWidth', 1.5); grid on;
xlabel('Frequency $\omega$ (rad/s)', 'Interpreter','latex', 'FontSize', 12);
ylabel('$A_4(\omega)$', 'Interpreter','latex', 'FontSize', 12);
title('Amplitude Response $A_4(\omega)$', 'Interpreter','latex', 'FontSize', 14);

sgtitle('Amplitude-Frequency Characteristics', 'Interpreter','latex', 'FontSize', 16);

% L(w)
L1 = 10 * log10(omega.^6 + 123*omega.^4 + 243*omega.^2 + 121) - 40 * log10(omega.^2 + 1);
L2 = 10 * log10(49*omega.^6 + 114*omega.^4 + 81*omega.^2 + 16) - 40 * log10(omega.^2 + 1);
L3 = 10 * log10(81*omega.^6 + 331*omega.^4 + 419*omega.^2 + 169) - 40 * log10(omega.^2 + 1);
L4 = -20 * log10(omega.^2 + 1);

% L(w) renders
figure;

subplot(2,2,1)
semilogx(omega, L1, 'b', 'LineWidth', 1.5); grid on;
xlabel('Frequency $\omega$ (rad/s)', 'Interpreter','latex', 'FontSize', 12);
ylabel('$L_1(\omega)$ (dB)', 'Interpreter','latex', 'FontSize', 12);
title('Logarithmic Amplitude Response $L_1(\omega)$', 'Interpreter','latex', 'FontSize', 14);

subplot(2,2,2)
semilogx(omega, L2, 'g', 'LineWidth', 1.5); grid on;
xlabel('Frequency $\omega$ (rad/s)', 'Interpreter','latex', 'FontSize', 12);
ylabel('$L_2(\omega)$ (dB)', 'Interpreter','latex', 'FontSize', 12);
title('Logarithmic Amplitude Response $L_2(\omega)$', 'Interpreter','latex', 'FontSize', 14);

subplot(2,2,3)
semilogx(omega, L3, 'r', 'LineWidth', 1.5); grid on;
xlabel('Frequency $\omega$ (rad/s)', 'Interpreter','latex', 'FontSize', 12);
ylabel('$L_3(\omega)$ (dB)', 'Interpreter','latex', 'FontSize', 12);
title('Logarithmic Amplitude Response $L_3(\omega)$', 'Interpreter','latex', 'FontSize', 14);

subplot(2,2,4)
semilogx(omega, L4, 'm', 'LineWidth', 1.5); grid on;
xlabel('Frequency $\omega$ (rad/s)', 'Interpreter','latex', 'FontSize', 12);
ylabel('$L_4(\omega)$ (dB)', 'Interpreter','latex', 'FontSize', 12);
title('Logarithmic Amplitude Response $L_4(\omega)$', 'Interpreter','latex', 'FontSize', 14);

sgtitle('Logarithmic Amplitude-Frequency Characteristics', 'Interpreter','latex', 'FontSize', 16);

% phi(w)
phi1 = atan2(omega.^3 - 23*omega, 13*omega.^2 - 11);
phi2 = atan2(-7*omega.^3 + 15*omega, 18*omega.^2 - 4);
phi3 = atan2(-9*omega.^3 + 35*omega, 31*omega.^2 - 13);
phi4 = atan2(2*omega, -omega.^2 + 1);

% rad to deg
phi1 = rad2deg(phi1);
phi2 = rad2deg(phi2);
phi3 = rad2deg(phi3);
phi4 = rad2deg(phi4);

% phi(w)
figure;

subplot(2,2,1)
plot(omega, phi1, 'b', 'LineWidth', 1.5); grid on;
xlabel('Frequency $\omega$ (rad/s)', 'Interpreter','latex', 'FontSize', 12);
ylabel('$\varphi_1(\omega) (^\circ)$', 'Interpreter','latex', 'FontSize', 12);
title('Phase-Frequency Characteristic $\varphi_1(\omega)$', 'Interpreter','latex', 'FontSize', 14);

subplot(2,2,2)
plot(omega, phi2, 'g', 'LineWidth', 1.5); grid on;
xlabel('Frequency $\omega$ (rad/s)', 'Interpreter','latex', 'FontSize', 12);
ylabel('$\varphi_2(\omega) (^\circ)$', 'Interpreter','latex', 'FontSize', 12);
title('Phase-Frequency Characteristic $\varphi_2(\omega)$', 'Interpreter','latex', 'FontSize', 14);

subplot(2,2,3)
plot(omega, phi3, 'r', 'LineWidth', 1.5); grid on;
xlabel('Frequency $\omega$ (rad/s)', 'Interpreter','latex', 'FontSize', 12);
ylabel('$\varphi_3(\omega) (^\circ)$', 'Interpreter','latex', 'FontSize', 12);
title('Phase-Frequency Characteristic $\varphi_3(\omega)$', 'Interpreter','latex', 'FontSize', 14);

subplot(2,2,4)
plot(omega, phi4, 'm', 'LineWidth', 1.5); grid on;
xlabel('Frequency $\omega$ (rad/s)', 'Interpreter','latex', 'FontSize', 12);
ylabel('$\varphi_4(\omega) (^\circ)$', 'Interpreter','latex', 'FontSize', 12);
title('Phase-Frequency Characteristic $\varphi_4(\omega)$', 'Interpreter','latex', 'FontSize', 14);

sgtitle('Phase-Frequency Characteristics', 'Interpreter','latex', 'FontSize', 16);

% log phi(w)
figure;

subplot(2,2,1)
semilogx(omega, phi1, 'b', 'LineWidth', 1.5); grid on;
xlabel('Frequency $\omega$ (rad/s)', 'Interpreter','latex', 'FontSize', 12);
ylabel('$\varphi_1(\omega) (^\circ)$', 'Interpreter','latex', 'FontSize', 12);
title('Logarithmic Phase Characteristic $\varphi_1(\omega)$', 'Interpreter','latex', 'FontSize', 14);

subplot(2,2,2)
semilogx(omega, phi2, 'g', 'LineWidth', 1.5); grid on;
xlabel('Frequency $\omega$ (rad/s)', 'Interpreter','latex', 'FontSize', 12);
ylabel('$\varphi_2(\omega) (^\circ)$', 'Interpreter','latex', 'FontSize', 12);
title('Logarithmic Phase Characteristic $\varphi_2(\omega)$', 'Interpreter','latex', 'FontSize', 14);

subplot(2,2,3)
semilogx(omega, phi3, 'r', 'LineWidth', 1.5); grid on;
xlabel('Frequency $\omega$ (rad/s)', 'Interpreter','latex', 'FontSize', 12);
ylabel('$\varphi_3(\omega) (^\circ)$', 'Interpreter','latex', 'FontSize', 12);
title('Logarithmic Phase Characteristic $\varphi_3(\omega)$', 'Interpreter','latex', 'FontSize', 14);

subplot(2,2,4)
semilogx(omega, phi4, 'm', 'LineWidth', 1.5); grid on;
xlabel('Frequency $\omega$ (rad/s)', 'Interpreter','latex', 'FontSize', 12);
ylabel('$\varphi_4(\omega) (^\circ)$', 'Interpreter','latex', 'FontSize', 12);
title('Logarithmic Phase Characteristic $\varphi_4(\omega)$', 'Interpreter','latex', 'FontSize', 14);

sgtitle('Logarithmic Phase-Frequency Characteristics', 'Interpreter','latex', 'FontSize', 16);