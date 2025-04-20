% plant parameters
A = [2 0 -4 2;
     0 2 -2 4;
    -4 -2 2 0;
     2 4 0 2];
B= [2 0;
    4 0;
    6 0;
    8 0];
C= [-2 2 2 2;
     2 0 0 2];
D=[0 3;
   0 1];

% A matrix eigenvalues
eig_A = eig(A)

% Jordan matrix
[P_J, A_J] = jordan(A)
B_J = P_J^-1 * B
C_J = C*P_J

% regulator
Qr= eye(4);
Rr = 10*eye(2);

[Pr,K,e]=icare(A,B,Qr,Rr);
K=-inv(Rr)*B'*Pr
ek=eig(A+B*K)

% observer
Ql= 2*eye(4);
Rl = eye(2);

[Pl,L,e]=icare(A',C',Ql,Rl);
L=-Pl*C'*Rl^-1
el=eig(A+L*C)