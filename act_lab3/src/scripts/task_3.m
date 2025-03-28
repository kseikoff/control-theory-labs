% plant parameters
A = [5 2 7; 2 1 2; -2 -3 -4];
B = [3; 1; -1];

% A matrix eigenvalues
A_e = eig(A)

% Jordan matrix
[P, J] = jordan(A);
Pre(:,1) = P(:,1);
Pre(:,2) = imag(P(:,2));
Pre(:,3) = real(P(:,3))
Pre_inv = Pre^-1
J_re = Pre_inv * A * Pre
B_jre = Pre_inv * B