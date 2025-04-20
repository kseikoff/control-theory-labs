% plant parameters
A = [2 0 -4 2;
     0 2 -2 4;
    -4 -2 2 0;
     2 4 0 2];
B= [2 0;
    4 0;
    6 0;
    8 0];
C= [-2 2 2 2;
     2 0 0 2];
D=[0 3;
   0 1];

% A matrix eigenvalues
eig_A = eig(A)

% Jordan matrix
[P_J, A_J] = jordan(A)
B_J = P_J^-1 * B
C_J = C*P_J