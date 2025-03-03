% input data
A = [5 2 7; 2 1 2; -2 -3 -4];
B = [3; 1; -1];

% A matrix eigenvalues
A_e = eig(A);
disp(A_e);

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

% Г matrices
g_1 = [0 1 0; 0 0 1; -8 -12 -6];
y_1 = [1 0 0];
check_1 = [y_1; y_1*g_1; y_1*g_1^2];
rank_check_1 = rank(check_1);
disp(check_1);
disp(rank_check_1);

g_2 = [0 1 0; 0 0 1; -8000 -4440 -222];
y_2 = [1 0 0];
check_2 = [y_2; y_2*g_2; y_2*g_2^2];
rank_check_2 = rank(check_2);
disp(check_2);
disp(rank_check_2);

g_3 = [0 1 0; 0 0 1; -80 -48 -6];
y_3 = [1 0 0];
check_3 = [y_3; y_3*g_3; y_3*g_3^2];
rank_check_3 = rank(check_3);
disp(check_3);
disp(rank_check_3);