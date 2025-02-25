% input data
A = [1 -2 3; 2 -3 2; -2 1 -4];
B = [3; 1; -1];
C = [0 1 1; 0 -4 2];
D = [0; 0];

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

% output controllability matrix
U = [B A*B A^2*B];
U_out = [C*U D];
r_U_out = rank(U_out);
disp(U);
disp(U_out);
disp(r_U_out);

% new communication matrix
D_new = [1; 0];
U_out_new = [C*U D_new];
r_U_out_new = rank(U_out_new);
disp(U_out_new);
disp(r_U_out_new);