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

K=[2.1111 -13.4448 1.6787];
Kg=[-0.0932 18.6951 -8.1152];
Kf=[-725.9021 -225.1491 586.1685 -359.3897];

G=[-1 0 0;
    0 -5 0;
    0 0 -10];
Y=[1; 1; 1];

% find Q
cvx_begin sdp
variable Q(3,3)
Q*Gg-G*Q == Y*Yg;
cvx_end 

Q=Q
invQ=inv(Q)

null1=[0 0 0;
       0 0 0;
       0 0 0;
       0 0 0];
null2=[0;0;0;0];
barA = [Gf null1;
        Bf*Yf A]
barB = [null2;B]
barC=[Df*Yf C]

% solving Riccati
QL = eye(7);
v = 1;
R = 1;

[P,L,e]=icare(barA',barC',QL,R);
L=-P*barC'*R^-1