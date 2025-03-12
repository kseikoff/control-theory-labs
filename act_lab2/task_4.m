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

% G 2x2, Q 2x4, Y 2x2
G = [-2 1; 0 -2];
Y = [1 0; 0 1];
U = [Y G*Y];
disp(U);
disp(rank(U));

% observer synthesis
cvx_begin sdp
variable Q(2,4)
G*Q-Q*A == Y*C;
cvx_end
disp(Q);