% plant parameters
A=[5 2 7;
   2 1 2;
  -2 -3 -4];
B=[3;1;-1];
Bf=[-4 -1;
     0 0;
     4 0];
C=[2 0 3];
D=2;
Df=[8 3];

Gf=[25 6 -20 11;
    14 3 -10 4;
    40 11 -31 17;
    6 4 -4 3];
Yf=[8 2 -6 4;
   -20 -6 16 -9];

Gg = [0 1 0;
     -1 0 0;
      0 0 0];
Yg=[4 0 -1];
wg0 = [0;1;1];

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

% G eigenvalues
Gf_eig = eig(Gf)
Gg_eig = eig(Gg)

% solving Riccati: feedback comp
Q = eye(3);
v = 2;
R = 1;
a = 2;

Aa = A + eye(3) * (a-0.00000000000001);
[Pk,K,e]=icare(Aa,sqrt(v)*B,Q,R);
K=-inv(R)*B'*Pk
eK=eig(A+B*K)

% check Frankis-Davison: Kg
Gg_eig(1)
check_Kg1 = [A-eye(3)*Gg_eig(1) B; C D]
rank(check_Kg1)

Gg_eig(2)
check_Kg2 = [A-eye(3)*Gg_eig(2) B; C D]
rank(check_Kg2)

Gg_eig(3)
check_Kg3 = [A-eye(3)*Gg_eig(3) B; C D]
rank(check_Kg3)

% solving Frankis-Davison: Kg
cvx_begin sdp
variable Pg(3,3)
variable Kg(1,3)
Pg*Gg-(A+B*K)*Pg == B*Kg;
(C+D*K)*Pg+D*Kg == Yg;
cvx_end

Pg=Pg
Kg=Kg

% check Frankis-Davison: Kf
Gf_eig(1)
check_Kf1 = [A-eye(3)*Gf_eig(1) B; C D]
rank(check_Kf1)

Gf_eig(2)
check_Kf2 = [A-eye(3)*Gf_eig(2) B; C D]
rank(check_Kf2)

Gf_eig(3)
check_Kf3 = [A-eye(3)*Gf_eig(3) B; C D]
rank(check_Kf3)

Gf_eig(4)
check_Kf4 = [A-eye(3)*Gf_eig(4) B; C D]
rank(check_Kf4)

% solving Frankis-Davison: Kf
cvx_begin sdp
variable Pf(3,4)
variable Kf(1,4)
Pf*Gf-(A+B*K)*Pf-Bf*Yf == B*Kf;
(C+D*K)*Pf+D*Kf == -Df*Yf;
cvx_end

Pf=Pf
Kf=Kf