% plant parameters
A=[0 1; -1 2];
B=[1 2; -3 3];
C=[2 1; 3 -2];
D=[-2 0; 0 1];

Bf=[1 2; -1 3];
Df=[1 0; 0 -1];
Cz=[4 2; -1 1];
Dz=[2 0; 0 1];

G=[-2.5 1;
   -1 -2.5];
Y=[1 0; 0 1];

Gw=[0 2 0 0 0 0;
   -2 0 0 0 0 0;
    0 0 0 3 0 0;
    0 0 -3 0 0 0;
    0 0 0 0 0 4;
    0 0 0 0 -4 0];
Y1=[0 0 0 9 0 0;
    3 0 0 0 0 0];
Y2=[6 0 0 0 0 0;
    0 0 0 8 0 0];
Yg=[0 0 0 0 3 0;
    0 0 0 0 0 6];

% out y
U = [B A*B];
U_out_y = [C*U D]
rank(U_out_y)

% out z
U_out_z = [Cz*U Dz]
rank(U_out_z)

% Jordan matrix
[P, J] = jordan(A)
B_J = inv(P) * B
C_J_y = C * P
C_J_z = Cz * P

% W_y(s)
sys_y = ss(A, B, C, D);
W_y_s = tf(sys_y)

% W_z(s)
sys_z = ss(A, B, Cz, Dz);
W_z_s = tf(sys_z)

% zeros
zeros_y = zero(sys_y)
zeros_z = zero(sys_z)

% check O(Y,G)
O_Y_G = [Y; Y*G]
rank_O_Y_G=rank(O_Y_G)

% K regulator synthesis
cvx_begin sdp
variable P(2,2)
A*P-P*G == B*Y;
cvx_end

K = -Y*inv(P)
ApBK=A+B*K;
eig_ApBK=eig(ApBK)

% check rank for K_w
CzpDzK=Cz+Dz*K;
eig_G=eig(G);
mat1 = [ApBK-eye(2)*eig_G(1) B; CzpDzK Dz]
mat2 = [ApBK-eye(2)*eig_G(2) B; CzpDzK Dz]
check_mat1=rank(mat1)
check_mat2=rank(mat2)

% solving Frankis-Davison: Kw
cvx_begin sdp
variable Pw(2,6)
variable Kw(2,6)
Pw*Gw-(A+B*K)*Pw-Bf*Y1 == B*Kw;
(Cz+Dz*K)*Pw+Dz*Kw == Yg;
cvx_end

Kw=Kw

% observer
null = [0 0; 0 0; 0 0; 0 0; 0 0; 0 0];
barA=[Gw null; Bf*Y1 A]
barB=[null; B]
barC=[Df*Y2-Yg C]

Ql = eye(8);
Rl = eye(2);

[Pl,barL,e]=icare(barA',barC',Ql,Rl);
barL=-Pl*barC'*Rl^-1
el=eig(barA+barL*barC)

barK=[Kw K]

barL'