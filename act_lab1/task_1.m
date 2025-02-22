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

% gramian
integrand = @(t) expm(A * t) * (B * B') * expm(A' * t);
P_t1 = integral(@(t) integrand(t), 0, t1, 'ArrayValued', true);
disp(P_t1);

% gramian eigenvalues
e = eig(P_t1);
disp(e);