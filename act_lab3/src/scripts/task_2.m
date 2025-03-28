% plant parameters
A=[2 0 -4 2;
   0 2 -2 4;
  -4 -2 2 0;
   2 4 0 2];
B=[2; 4; 6; 8];
C=[-2 2 2 2;
   2 0 0 2];

% A matrix eigenvalues
A_e = eig(A)

% Jordan decomposition
[PA, JA] = jordan(A)
PA_inv = inv(PA)
BA = PA_inv * B
CA = C * PA

% Desired decay rate
a1 = 4;
a2 = 1;

% case 1: ak == al
ak = a1;
al = a1;

% solving Riccati
Q = 0;
v = 2;
R = 1;

% find K
Aak = A + eye(4) * (ak-0.0000000001);
[P,K,e]=icare(Aak,sqrt(2)*B,Q,R);
K_case1=-inv(R)*B'*P
eK_case1=eig(A+B*K_case1)

% find L
x0=[1;1;1;1];
x0_est=[0;0;0;0];
e0=x0-x0_est;

% solving LMI with control constraint
% mumu 2x2, not scalar anymore
% minimizing matrix by its norm
cvx_begin sdp
variable Q(4,4)
variable Y(4,2)
variable mumu(2,2)
minimize norm(mumu)
Q>0.0001*eye(4);
A'*Q + Q*A+ 2*al*Q + C'*Y' + Y*C <= 0;
[Q e0;
 e0' 1]>0;
[Q Y;
 Y' mumu]>0;
cvx_end
L_case1=inv(Q)*Y
eL_case1=eig(A+L_case1*C)

% case 2: ak > al
ak = a1;
al = a2;

% K found case 1
K_case2 = K_case1
eK_case2 = eK_case1

% find L
cvx_begin sdp
variable Q(4,4)
variable Y(4,2)
variable mumu(2,2)
minimize norm(mumu)
Q>0.0001*eye(4);
A'*Q + Q*A+ 2*al*Q + C'*Y' + Y*C <= 0;
[Q e0;
 e0' 1]>0;
[Q Y;
 Y' mumu]>0;
cvx_end
L_case2=inv(Q)*Y
eL_case2=eig(A+L_case2*C)

% case 3: ak < al
ak = a2;
al = a1;

% find K
Q = eye(4)*0.000001;
Aak = A + eye(4) * ak;
[P,K,e]=icare(Aak,sqrt(2)*B,Q,R);
K_case3=-inv(R)*B'*P
eK_case3=eig(A+B*K_case3)

% L found case 1
L_case3 = L_case1
eL_case3 = eL_case1
