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

% find L: solving LMI with control constraint
al = 1;
xw0=[0;0;0;1;1;1;1];
xw0_est=[0;0;0;0;0;0;0];
xwe0=xw0-xw0_est;

cvx_begin sdp
variable Q(7,7)
variable Y(7,1)
variable mumu
minimize mumu
Q>0.0001*eye(7);
bigA'*Q + Q*bigA+ 2*al*Q + bigC'*Y' + Y*bigC <= 0;
[Q xwe0;
 xwe0' 1]>0;
[Q Y;
 Y' mumu]>0;
cvx_end
mu = sqrt(mumu)
bigL=inv(Q)*Y
e_bigL=eig(bigA+bigL*bigC)

% K2_CD regulator synthesis
cvx_begin sdp
variable P1(3,4)
variable Y1(1,4)
P1*G-A*P1 == B*Y1+Bf;
C*P1 + D == 0;
cvx_end

K2_CD = Y1-K1*P1

% K2_CzDz regulator synthesis
cvx_begin sdp
variable P2(3,4)
variable Y2(1,4)
P2*G-A*P2 == B*Y2+Bf;
Cz*P2 + Dz == 0;
cvx_end

K2_CzDz = Y2-K1*P2

L1 = [bigL(1); bigL(2); bigL(3)]
L2 = [bigL(4); bigL(5); bigL(6); bigL(7)]

% the principle of the internal model
eig_G = eig(G)
bar_A_CD = [A+B*K1+L1*C B*K2_CD+Bf+L1*D;
            L2*C G+L2*D]
bar_A_CzDz = [A+B*K1+L1*C B*K2_CzDz+Bf+L1*D;
              L2*C G+L2*D]
eig_bar_A_CD = eig(bar_A_CD)
eig_bar_A_CzDz = eig(bar_A_CzDz)