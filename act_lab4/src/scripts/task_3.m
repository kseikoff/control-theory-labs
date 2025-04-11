% plant parameters
A = [5 2 7;
     2 1 2;
    -2 -3 -4];
B = [3; 1; -1];
Bf = [-4 0 0 -1;
      0 0 0 0;
      4 0 0 0];

G = [25 6 -20 11;
     14 3 -10 4;
     40 11 -31 17;
     6 4 -4 3];

Cz = [-2 1 -1];
Dz = [-20 -6 16 -9];

C = [2 0 3];
D = [8 2 -6 4];

% matrices
bigC = [C D];
bigA = [A Bf;
      zeros(4,3) G];

% eig bigA
eig_bigA = eig(bigA)

% re jordan
[P, J] = jordan(bigA);
P_re(:,1) = P(:,1);
P_re(:,2) = imag(P(:,2));
P_re(:,3) = real(P(:,3));
P_re(:,4) = imag(P(:,4));
P_re(:,5) = real(P(:,5));
P_re(:,6) = imag(P(:,6));
P_re(:,7) = real(P(:,7))
Pre_inv = P_re^-1
bigA_J_re = Pre_inv * bigA * P_re
C_J_re = bigC * P_re

% solving Riccati
Q = 0;
v = 2;
R = 1;
a = 2;

Aa = A + eye(3) * a;
[Pk,K,e]=icare(Aa,sqrt(2)*B,Q,R);
K1=-inv(R)*B'*Pk
eK1=eig(A+B*K1)