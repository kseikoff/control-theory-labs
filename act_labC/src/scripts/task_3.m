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
Q =[2.0000 -2.0000 -1.0000;
    0.7692 -0.1538 -0.2000;
    0.3960 -0.0396 -0.1000];
invQ=inv(Q)

G=[-0.5 0 0 0;
    0 -1 0 0;
    0 0 -1.5 0;
    0 0 0 -2];
Y=[1 0 1 0;
   0 1 0 1];
newY=[1;1;1;1];

% find QL
cvx_begin sdp
variable QL(4,4)
QL*G-Gf'*QL == Yf'*Y;
cvx_end

QL=QL
L=-Y*QL^-1;
L=L'

% find barC
cvx_begin sdp
variable barC(2,3)
barC*Bf == eye(2);
cvx_end

barC=barC

F=Gf-L*Yf

% find Qout
cvx_begin sdp
variable Qout(4,4)
Qout*Gf-G*Qout == newY*Df*Yf;
cvx_end

Qout=Qout
invQout=inv(Qout)