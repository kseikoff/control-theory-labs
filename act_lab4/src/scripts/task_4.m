% plant parameters
A= [0 1; 0 0];
B=[0;1];
C=[1 0];
G = [0 1 0 0 0 0;
    -1 0 0 0 0 0;
     0 0 0 3 0 0;
     0 0 -3 0 0 0;
     0 0 0 0 0 5;
     0 0 0 0 -5 0];
Cz = [-1 0];
Dz = [1 0 1 0 1 0];
Bg = 0;

[P1, J] = jordan(A)

% K1 regulator synthesis
Q = 0;
v = 2;
R = 1;
a = 2;

Aa = A + eye(2) * a;
[Pk,K,e]=icare(Aa,sqrt(2)*B,Q,R);
K1=-inv(R)*B'*Pk
eK1=eig(A+B*K1)

% K2
cvx_begin sdp
variable P(2,6)
variable Y(1,6)
P*G-A*P == B*Y+Bg;
Cz*P + Dz == 0;
cvx_end

K2 = Y-K1*P
