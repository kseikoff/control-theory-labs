% plant parameters
A = [5 2 7;
     2 1 2;
    -2 -3 -4];
B = [3; 1; -1];
Bf = [-4 0 0 -1;
      0 0 0 0;
      4 0 0 0];
x0 = [0; 0; 0];

G = [25 6 -20 11;
     14 3 -10 4;
     40 11 -31 17;
     6 4 -4 3];
wf0 = [1;1;1;1];

Cz = [-2 1 -1];
D = 0;

% G eigenvalues
G_eig = eig(G)

% A eigenvalues
A_eig = eig(A)

% Jordan matrix
[Aj, J] = jordan(A);
Ajre(:,1) = Aj(:,1);
Ajre(:,2) = imag(Aj(:,2));
Ajre(:,3) = real(Aj(:,3))
Ajre_inv = Ajre^-1
J_re = Ajre_inv * A * Ajre
B_jre = Ajre_inv * B

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