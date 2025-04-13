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
D = 0;

% G eigenvalues
G_eig = eig(G)

% A eigenvalues
A_eig = eig(A)

% Jordan matrix
[P1, J] = jordan(A);
Pjre(:,1) = P1(:,1);
Pjre(:,2) = imag(P1(:,2));
Pjre(:,3) = real(P1(:,3))
Pjre_inv = Pjre^-1
Aj_re = Pjre_inv * A * Pjre
B_jre = Pjre_inv * B

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
P*G-A*P == B*Y+Bf;
Cz*P + D == 0;
cvx_end

K2 = Y-K1*P