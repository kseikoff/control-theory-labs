A = [1 -2 3; 2 -3 2; -2 1 -4];
B = [-3; -1; 3];
t1 = 3;

% Gramian
integrand = @(t) expm(A * t) * (B * B') * expm(A' * t);
P_t1 = integral(@(t) integrand(t), 0, t1, 'ArrayValued', true);
disp(P_t1);