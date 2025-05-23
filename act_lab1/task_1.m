% input data
A = [1 -2 3; 2 -3 2; -2 1 -4];
B = [-3; -1; 3];
x1 = [4; 3; -3];
t1 = 3;

% controllability matrix
U = [B A*B A*A*B];
r = rank(U);
disp(U);
disp(r);

% A matrix eigenvalues
A_e = eig(A);
disp(A_e);

% Houtus matrices
H_1 = [A-A_e(1)*eye(3) B];
r_1 = rank(H_1);
disp(H_1);
disp(r_1);

H_2 = [A-A_e(2)*eye(3) B];
r_2 = rank(H_2);
disp(H_2);
disp(r_2);

H_3 = [A-A_e(3)*eye(3) B];
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
B_jre = P1_inv * B;
disp(P1);
disp(P1_inv);
disp(J_re);
disp(B_jre);

% gramian
integrand = @(t) expm(A * t) * (B * B') * expm(A' * t);
P_t1 = integral(@(t) integrand(t), 0, t1, 'ArrayValued', true);
disp(P_t1);

% gramian eigenvalues
e = eig(P_t1);
disp(e);

% u_t
u_t = @(t) B' * expm(A' * (t1 - t)) * inv(P_t1) * x1;

% u_t modeling
time = linspace(0, t1, 1000);

control = arrayfun(@(t) u_t(t), time, 'UniformOutput', false);
control = cell2mat(control);

figure;
plot(time, control, 'LineWidth', 1.5);
xlabel('t');
ylabel('u(t)');
title('Control');
grid on;

% x_t
dxdt = @(t, x) A * x + B * u_t(t);

% x_t modeling
[t, x] = ode45(dxdt, [0 t1], [0; 0; 0]);

figure;
plot(t, x(:,1), 'r', 'LineWidth', 1.5); hold on;
plot(t, x(:,2), 'g', 'LineWidth', 1.5);
plot(t, x(:,3), 'b', 'LineWidth', 1.5);
scatter(t1, x1(1), 'ro', 'filled');
scatter(t1, x1(2), 'go', 'filled');
scatter(t1, x1(3), 'bo', 'filled');
yline(x1(1), '--', 'Color', [0.5 0.5 0.5]); % x_1
yline(x1(2), '--', 'Color', [0.5 0.5 0.5]); % x_2
yline(x1(3), '--', 'Color', [0.5 0.5 0.5]); % x_3
xlabel('t');
ylabel('x(t)');
legend('x_1', 'x_2', 'x_3');
title('System state');
grid on;
