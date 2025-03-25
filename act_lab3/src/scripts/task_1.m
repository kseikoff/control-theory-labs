% plant parameters
A = [5 2 7; 2 1 2; -2 -3 -4];
B = [3; 1; -1];

% A matrix eigenvalues
A_e = eig(A)

% Jordan matrix
[P, J] = jordan(A);
Pre(:,1) = P(:,1);
Pre(:,2) = imag(P(:,2));
Pre(:,3) = real(P(:,3))
Pre_inv = Pre^-1
J_re = Pre_inv * A * Pre
B_jre = Pre_inv * B

% Desired decay rate
a1 = 2;
a2 = 0.1;

% solving LMI no restrictions on control
cvx_begin sdp
% a1
variable P1(3,3) symmetric
variable Y1(1,3)
P1 > 0.0001*eye(3);
P1*A' + A*P1 + 2*a1*P1 + Y1'*B'+ B*Y1 <= 0;
cvx_end

cvx_begin sdp
% a2
variable P2(3,3) symmetric
variable Y2(1,3)
P2 > 0.0001*eye(3);
P2*A' + A*P2 + 2*a2*P2 + Y2'*B'+ B*Y2 <= 0;
cvx_end

K1_a1 = Y1*inv(P1)
K1_a2 = Y2*inv(P2)

% A+BK1_ai eigenvalues
ABK1_a1 = A+B*K1_a1;
ABK1_a2 = A+B*K1_a2;
eig(ABK1_a1)
eig(ABK1_a2)

% solving LMI with control constraint
x0 = [1; 1; 1];

% a1
cvx_begin sdp
variable P12(3,3) symmetric
variable Y12(1,3)
variable mumu_a1
minimize mumu_a1
P12 > 0.0001*eye(3);
P12*A' + A*P12 + 2*a1*P12 + Y12'*B'+ B*Y12 <= 0;
[P12 x0;
 x0' 1] > 0;
[P12 Y12';
 Y12 mumu_a1] > 0;
cvx_end

cvx_begin sdp
% a2
variable P22(3,3) symmetric
variable Y22(1,3)
variable mumu_a2
minimize mumu_a2
P22 > 0.0001*eye(3);
P22*A' + A*P22 + 2*a2*P22 + Y22'*B'+ B*Y22 <= 0;
[P22 x0;
 x0' 1] > 0;
[P22 Y22';
 Y22 mumu_a2] > 0;
cvx_end

mu_a1 = sqrt(mumu_a1)
mu_a2 = sqrt(mumu_a2)

K2_a1 = Y12*inv(P12)
K2_a2 = Y22*inv(P22)

% A+BK2_ai eigenvalues
ABK2_a1 = A+B*K2_a1;
ABK2_a2 = A+B*K2_a2;
eig(ABK2_a1)
eig(ABK2_a2)
