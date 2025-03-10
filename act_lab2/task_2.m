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

% G matrices
Y = [1; 0; 0; 0];

G1 = [0 1 0 0; 0 0 1 0; 0 0 0 1; -16 -32 -24 -8];
U1 = [Y G1*Y G1^2*Y G1^3*Y];
disp(U1);
disp(rank(U1));
disp(eig(G1));

G2 = [0 1 0 0; 0 0 1 0; 0 0 0 1; -16000000 -8888000 -448440 -2222];
U2 = [Y G2*Y G2^2*Y G2^3*Y];
disp(U2);
disp(rank(U2));
disp(eig(G2));

G3 = [0 1 0 0; 0 0 1 0; 0 0 0 1; -260 -132 -49 -8];
U3 = [Y G3*Y G3^2*Y G3^3*Y];
disp(U3);
disp(rank(U3));
disp(eig(G3));

% observer synthesis
cvx_begin sdp
variable Q1(4,4)
variable Q2(4,4)
variable Q3(4,4)
G1*Q1-Q1*A == Y*C;
G2*Q2-Q2*A == Y*C;
G3*Q3-Q3*A == Y*C;
cvx_end

L1 = inv(Q1)*Y;
disp(Q1);
disp(L1);

L2 = inv(Q2)*Y;
disp(Q2);
disp(L2);

L3 = inv(Q3)*Y;
disp(Q3);
disp(L3);

% A+LC eigenvalues
ALC1 = A+L1*C;
ALC2 = A+L2*C;
ALC3 = A+L3*C;

ALC1_eig = eig(ALC1);
ALC2_eig = eig(ALC2);
ALC3_eig = eig(ALC3);

disp(ALC1);
disp(ALC1_eig);

disp(ALC2);
disp(ALC2_eig);

disp(ALC3);
disp(ALC3_eig);