% plant parameters
A = [0 1 0 1;
    -26 -7 20 -11;
     0 1 -1 2;
     16 4 -14 8];
C = [-1 0 1 -1];

% A matrix eigenvalues
eig_A = eig(A)

% Jordan matrix
[P_J, A_J] = jordan(A);
P_Jre(:,1) = real(P_J(:,1));
P_Jre(:,2) = imag(P_J(:,2));
P_Jre(:,3) = real(P_J(:,3));
P_Jre(:,4) = imag(P_J(:,4));
A_Jre = P_Jre^-1 * A * P_Jre
C_jre = C*P_Jre

% solving Riccati
% aftermath
% Q=I,R=1: L=[0.5233;-7.6440;0.0221;6.4594]; e={-0.3454+-1.4689i,-3.1349+-3.3599i}
% Q=100I,R=1: L=[7.6954;-15.6987;-4.5805;15.1421]; e={-0.3419+-1.4702i,-13.3670+-5.6243i}
% Q=I,R=100: L=[0.0079;-1.7287;0.5099;2.0960]; e={-0.4006+-1.1676i,-0.3965+-1.9848i}
% Q=100I,R=100: L=[0.5233;-7.6440;0.0221;6.4594]; e={-0.3454+-1.4689i,-3.1349+-3.3599i}
a = 100;
Q= a*eye(4);
R = a*1;

[P,L,e]=icare(A',C',Q,R);
P
L=-P*C'*R^-1
e=eig(A+L*C)