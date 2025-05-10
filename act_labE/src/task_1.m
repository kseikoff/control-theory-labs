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

% time for renders
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
xlabel('Time (s)', 'FontSize', 12); ylabel('w_1(t)', 'FontSize', 12);
title('Weight Function w_1(t)', 'FontSize', 14);

subplot(2,2,2)
plot(t, w2, 'g', 'LineWidth', 1.5); grid on;
xlabel('Time (s)', 'FontSize', 12); ylabel('w_2(t)', 'FontSize', 12);
title('Weight Function w_2(t)', 'FontSize', 14);

subplot(2,2,3)
plot(t, w3, 'r', 'LineWidth', 1.5); grid on;
xlabel('Time (s)', 'FontSize', 12); ylabel('w_3(t)', 'FontSize', 12);
title('Weight Function w_3(t)', 'FontSize', 14);

subplot(2,2,4)
plot(t, w4, 'm', 'LineWidth', 1.5); grid on;
xlabel('Time (s)', 'FontSize', 12); ylabel('w_4(t)', 'FontSize', 12);
title('Weight Function w_4(t)', 'FontSize', 14);

% h(t) renders
figure;
subplot(2,2,1)
plot(t, h1, 'b', 'LineWidth', 1.5); grid on;
xlabel('Time (s)', 'FontSize', 12); ylabel('h_1(t)', 'FontSize', 12);
title('Step Response h_1(t)', 'FontSize', 14);

subplot(2,2,2)
plot(t, h2, 'g', 'LineWidth', 1.5); grid on;
xlabel('Time (s)', 'FontSize', 12); ylabel('h_2(t)', 'FontSize', 12);
title('Step Response h_2(t)', 'FontSize', 14);

subplot(2,2,3)
plot(t, h3, 'r', 'LineWidth', 1.5); grid on;
xlabel('Time (s)', 'FontSize', 12); ylabel('h_3(t)', 'FontSize', 12);
title('Step Response h_3(t)', 'FontSize', 14);

subplot(2,2,4)
plot(t, h4, 'm', 'LineWidth', 1.5); grid on;
xlabel('Time (s)', 'FontSize', 12); ylabel('h_4(t)', 'FontSize', 12);
title('Step Response h_4(t)', 'FontSize', 14);