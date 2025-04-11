% plant parameters
A = [5 2 7;
     2 1 2;
    -2 -3 -4];
B = [3; 1; -1];
Bg = 0;

G = [25 6 -20 11;
     14 3 -10 4;
     40 11 -31 17;
     6 4 -4 3];

Cz = [-2 1 -1];
Dz = [-20 -6 16 -9];

% solving Riccati
Q = 0;
v = 2;
R = 1;
a = 2;

Aa = A + eye(3) * a;
[Pk,K,e]=icare(Aa,sqrt(2)*B,Q,R);
K1=-inv(R)*B'*Pk
eK1=eig(A+B*K1)

% K2 regulator synthesis
cvx_begin sdp
variable P(3,4)
variable Y(1,4)
P*G-A*P == B*Y+Bg;
Cz*P + Dz == 0;
cvx_end

K2 = Y-K1*P