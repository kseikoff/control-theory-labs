% input data
A = [1 -2 3; 2 -3 2; -2 1 -4];
B = [3; 1; -1];
x1p = [4; 3; -3];
x1pp = [3; 3; -2];
t1 = 3;

% controllability matrix
U = [B A*B A*A*B];
r = rank(U);
detU = det(U);
disp(U);
disp(r);
disp(detU)

% check x1p in subspace
U_x1p = [U x1p];
r_x1p = rank(U_x1p);
disp(U_x1p);
disp(r_x1p);

% check x1pp in subspace
U_x1pp = [U x1pp];
r_x1pp = rank(U_x1pp);
disp(U_x1pp);
disp(r_x1pp);

if r_x1p == r
    x1 = x1p;
else
    x1 = x1pp;
end
disp(x1)

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
P1(:,1) = real(P(:,1));
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
u_t = @(t) B' * expm(A' * (t1 - t)) * pinv(P_t1) * x1;
disp(pinv(P_t1));

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
plot(t, x, 'LineWidth', 1.5);
yline(x1(1), '--', 'Color', [0.5 0.5 0.5]); % x_1
yline(x1(2), '--', 'Color', [0.5 0.5 0.5]); % x_2
yline(x1(3), '--', 'Color', [0.5 0.5 0.5]); % x_3
xlabel('t');
ylabel('x(t)');
legend('x_1', 'x_2', 'x_3');
title('System state');
grid on;
