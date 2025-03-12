% input data
A = [2 0 -4 2;
     0 2 -2 4;
    -4 -2 2 0;
     2 4 0 2];
B = [2;4;6;8];
C = [0 0 1 0;
     1 0 0 0];
D = [3;1];

% Jordan matrix
[P, J] = jordan(A);
P_inv = inv(P);
disp(P);
disp(J);
disp(P_inv);

B_new = P_inv*B;
C_new = C*P;
disp(B_new);
disp(C_new);