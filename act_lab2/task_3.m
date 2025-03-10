% input data
A = [2 0 -4 2;
     0 2 -2 4;
    -4 -2 2 0;
     2 4 0 2];
B = [2;4;6;8];
C = [-2 2 2 2;
     2 0 0 2];
D = [3;1];

% A eigenvalues
disp(eig(A));

% controlability
U = [B A*B A^2*B A^3*B];
disp(rank(U));

% observability
V = [C;C*A;C*A^2;C*A^3];
disp(rank(V));

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

% G control
G_K = [0 1 0 0;
       0 0 1 0;
       0 0 0 1;
       -1 -4 -6 -4];
Y_K = [1 0 0 0];
O = [Y_K; Y_K*G_K; Y_K*G_K^2; Y_K*G_K^3];
disp(O);
disp(rank(O));

% regulator synthesis
cvx_begin sdp
variable P(4,4)
A*P-P*G_K == B*Y_K;
cvx_end

K = -Y_K*inv(P);
disp(P);
disp(K);

% A+BK eigenvalues
ABK = A+B*K;
disp(ABK);
disp(eig(ABK));

% G observe
G_L = [0 1 0 0;
       0 0 1 0;
       0 0 0 1;
       -108 -135 -63 -13];
Y_L = [1 0;
       0 1;
       0 0;
       0 0];
U = [Y_L G_L*Y_L G_L^2*Y_L G_L^3*Y_L];
disp(U);
disp(rank(U));

% observer synthesis
cvx_begin sdp
variable Q(4,4)
G_L*Q-Q*A == Y_L*C;
cvx_end

L = inv(Q)*Y_L;
disp(Q);
disp(L);

% A+LC eigenvalues
ALC = A+L*C;
disp(ALC);
disp(eig(ALC));