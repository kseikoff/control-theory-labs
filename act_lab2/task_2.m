% input data
A = [0 1 0 1;
     -26 -7 20 -11;
     0 1 -1 2;
     16 4 -14 8];
C = [-1 0 1 -1];

% A matrix eigenvalues
A_e = eig(A);
disp(A_e);

% Jordan matrix
[P, J] = jordan(A);
P_re(:,1) = real(P(:,1));
P_re(:,2) = imag(P(:,2));
P_re(:,3) = real(P(:,3));
P_re(:,4) = imag(P(:,4));
P_re_inv = P_re^-1; 
J_re = P_re_inv * A * P_re;
C_jre = C * P_re;
disp(P_re);
disp(P_re_inv);
disp(J_re);
disp(C_jre);