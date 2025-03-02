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
