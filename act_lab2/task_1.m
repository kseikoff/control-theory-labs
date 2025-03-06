% input data
A = [5 2 7; 2 1 2; -2 -3 -4];
B = [3; 1; -1];

% A matrix eigenvalues
A_e = eig(A);
disp(A_e);

% Jordan matrix
[P, J] = jordan(A);
Pre(:,1) = P(:,1);
Pre(:,2) = imag(P(:,2));
Pre(:,3) = real(P(:,3));
Pre_inv = Pre^-1; 
J_re = Pre_inv * A * Pre;
B_jre = Pre_inv * B;
disp(Pre);
disp(Pre_inv);
disp(J_re);
disp(B_jre);

% G matrices
G1 = [0 1 0; 0 0 1; -8 -12 -6];
Y1 = [1 0 0];
O1 = [Y1; Y1*G1; Y1*G1^2];
rank_O1 = rank(O1);
disp(O1);
disp(rank_O1);

G2 = [0 1 0; 0 0 1; -8000 -4440 -222];
Y2 = [1 0 0];
O2 = [Y2; Y2*G2; Y2*G2^2];
rank_O2 = rank(O2);
disp(O2);
disp(rank_O2);

G3 = [0 1 0; 0 0 1; -80 -48 -6];
Y3 = [1 0 0];
O3 = [Y3; Y3*G3; Y3*G3^2];
rank_O3 = rank(O3);
disp(O3);
disp(rank_O3);

% regulator synthesis
cvx_begin sdp
variable P1(3,3)
variable P2(3,3)
variable P3(3,3)
A*P1-P1*G1 == B*Y1;
A*P2-P2*G2 == B*Y2;
A*P3-P3*G3 == B*Y3;
cvx_end

K1 = -Y1*inv(P1);
disp(P1);
disp(K1);

K2 = -Y2*inv(P2);
disp(P2);
disp(K2);

K3 = -Y3*inv(P3);
disp(P3);
disp(K3);

% A+BK eigenvalues
ABK1 = A+B*K1;
ABK2 = A+B*K2;
ABK3 = A+B*K3;

ABK1_eig = eig(ABK1);
ABK2_eig = eig(ABK2);
ABK3_eig = eig(ABK3);

disp(ABK1);
disp(ABK1_eig);

disp(ABK2);
disp(ABK2_eig);

disp(ABK3);
disp(ABK3_eig);

