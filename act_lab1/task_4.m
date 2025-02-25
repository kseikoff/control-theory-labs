% input data
A = [-16 -27 7; 6 9 -4; -5 -11 0];
C = [0 -5 -5];
f_t = @(t) -9 * exp(-4 * t) * cos(t) + 9 * exp(-4 * t) * sin(t);
t1 = 3;

% observability matrix
V = [C; C*A; C*A^2];
r = rank(V);
det_V = det(V);
disp(A^2);
disp(V);
disp(r);
disp(det_V);

% A matrix eigenvalues
A_e = eig(A);
disp(A_e);

% Houtus matrices
H_1 = [A-A_e(1)*eye(3); C];
r_1 = rank(H_1);
disp(H_1);
disp(r_1);

H_2 = [A-A_e(2)*eye(3); C];
r_2 = rank(H_2);
disp(H_2);
disp(r_2);

H_3 = [A-A_e(3)*eye(3); C];
r_3 = rank(H_3);
disp(H_3);
disp(r_3);

% Jordan matrix
[P, J] = jordan(A);
P1(:,1) = P(:,1);
P1(:,2) = imag(P(:,2));
P1(:,3) = real(P(:,3));
P1_inv = P1^-1; 
J_re = P1_inv * A * P1;
C_jre = C * P1;
disp(P1);
disp(P1_inv);
disp(J_re);
disp(C_jre);

% gramian
integrand = @(t) expm(A' * t) * (C' * C) * expm(A * t);
Q_t1 = integral(@(t) integrand(t), 0, t1, 'ArrayValued', true);
disp(Q_t1);

% gramian eigenvalues
e = eig(Q_t1);
disp(e);

% initial conditions x(0)
integrand_x0 = @(t) expm(A' * t) * C' * f_t(t);
X_int = integral(@(t) integrand_x0(t), 0, t1, 'ArrayValued', true);
Q_t1_pinv = pinv(Q_t1);
x0 = Q_t1_pinv * X_int;
disp(Q_t1_pinv);
disp(x0);

% system modeling
x_t = @(t) expm(A * t) * x0; 
y_t = @(t) C * x_t(t);

time = linspace(0, t1, 1000);
y_arr = arrayfun(y_t, time);
f_arr = arrayfun(f_t, time);

figure;
plot(time, y_arr, 'b', 'LineWidth', 1.5); hold on;
plot(time, f_arr, 'r--', 'LineWidth', 1.5);
xlabel('t');
ylabel('y(t)');
legend('y(t)', 'f(t)');
title('y(t) and f(t) comparison');
grid on;

err = y_arr - f_arr;

figure;
plot(time, err, 'k', 'LineWidth', 1.5);
xlabel('t');
ylabel('e(t)');
title('Error');
grid on;

% other initial conditions
y_wanted = V * x0;
disp(y_wanted);

x0_1 = [-5.4; 1.8; 0]; % x_3 = 0
x0_2 = [-3.4; 0.8; 1]; % x_3 = 1
x0_3 = [-7.4; 2.8; -1]; % x_3 = -1

[t, x00] = ode45(@(t, x) A * x, time, x0);
[t, x1] = ode45(@(t, x) A * x, time, x0_1);
[t, x2] = ode45(@(t, x) A * x, time, x0_2);
[t, x3] = ode45(@(t, x) A * x, time, x0_3);

y00 = C * x0;
y1 = C * x1';
y2 = C * x2';
y3 = C * x3';

% x_t
figure;

subplot(3, 1, 1);
plot(t, x1(:,1), 'r', 'LineWidth', 1.5); hold on;
plot(t, x2(:,1), 'g--', 'LineWidth', 1.5);
plot(t, x3(:,1), 'b:', 'LineWidth', 1.5);
plot(t, x00(:,1), 'k:', 'LineWidth', 1.5);
scatter(0, x0_1(1), 'ro', 'filled');
scatter(0, x0_2(1), 'go', 'filled');
scatter(0, x0_3(1), 'bo', 'filled');
scatter(0, x0(1), 'ko', 'filled');
xlabel('t');
ylabel('x_1(t)');
title('x(t): x_1(t)');
legend('x_1(0)', 'x_2(0)', 'x_3(0)', 'x_0(0)', 'Location', 'northwest');
grid on;

subplot(3, 1, 2);
plot(t, x1(:,2), 'r', 'LineWidth', 1.5); hold on;
plot(t, x2(:,2), 'g--', 'LineWidth', 1.5);
plot(t, x3(:,2), 'b:', 'LineWidth', 1.5);
plot(t, x00(:,2), 'k:', 'LineWidth', 1.5);
scatter(0, x0_1(2), 'ro', 'filled');
scatter(0, x0_2(2), 'go', 'filled');
scatter(0, x0_3(2), 'bo', 'filled');
scatter(0, x0(2), 'ko', 'filled');
xlabel('t');
ylabel('x_2(t)');
title('x(t): x_2(t)');
legend('x_1(0)', 'x_2(0)', 'x_3(0)', 'x_0(0)', 'Location', 'southwest');
grid on;

subplot(3, 1, 3);
plot(t, x1(:,3), 'r', 'LineWidth', 1.5); hold on;
plot(t, x2(:,3), 'g--', 'LineWidth', 1.5);
plot(t, x3(:,3), 'b:', 'LineWidth', 1.5);
plot(t, x00(:,3), 'k:', 'LineWidth', 1.5);
scatter(0, x0_1(3), 'ro', 'filled');
scatter(0, x0_2(3), 'go', 'filled');
scatter(0, x0_3(3), 'bo', 'filled');
scatter(0, x0(3), 'ko', 'filled');
xlabel('t');
ylabel('x_3(t)');
title('x(t): x_3(t)');
legend('x_1(0)', 'x_2(0)', 'x_3(0)', 'x_0(0)', 'Location', 'northwest');
grid on;

% y_t
figure;
plot(t, y1, 'r', 'LineWidth', 1.5); hold on;
plot(t, y2, 'g--', 'LineWidth', 1.5);
plot(t, y3, 'b:', 'LineWidth', 1.5);
plot(t, y00, 'k:', 'LineWidth', 1.5);
xlabel('t');
ylabel('y(t)');
title('y(t) for x_i(0)');
legend('x_1(0)', 'x_2(0)', 'x_3(0)', 'x_0(0)', 'Location', 'southeast');
grid on;