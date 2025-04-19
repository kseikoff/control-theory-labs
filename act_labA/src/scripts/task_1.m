% plant parameters
A = [5 2 7; 2 1 2; -2 -3 -4];
B = [3; 1; -1];
x0 = [1;1;1];

% A matrix eigenvalues
eig_A = eig(A)

% Jordan matrix
[P_J, A_J] = jordan(A);
P_Jre(:,1) = P_J(:,1);
P_Jre(:,2) = imag(P_J(:,2));
P_Jre(:,3) = real(P_J(:,3));
A_Jre = P_Jre^-1 * A * P_Jre
B_jre = P_Jre^-1 * B

% solving Riccati
% aftermath
% Q=I,R=1: K=[-1.7473 -5.4627 -1.4004]; e={-1.4413,-3.8630,-2} J_min=9.7962
% Q=4I,R=1: K=[-2.0000 -6.9630 -0.9630]; e={-1,-7,-2} J_min=14.7764
% Q=I,R=4: K=[-1.6433 -4.9797 -1.5454]; e={-2,-2.1821+-0.6216i} J_min=33.8518
% Q=4I,R=4: K=[-1.7473 -5.4627 -1.4004]; e={-1.4413,-3.8630,-2} J_min=39.1846
a = 4;
Q= a*eye(3);
R = a*1;

% a1
[P,K,e]=icare(A,B,Q,R);
P
K=-inv(R)*B'*P
e=eig(A+B*K)

% quality functional
J_min = x0'*P*x0