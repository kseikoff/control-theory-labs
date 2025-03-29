% plant parameters
A = [5 2 7; 2 1 2; -2 -3 -4];
B = [3; 1; -1];

% A matrix eigenvalues
A_e = eig(A)

% Jordan matrix
[P, J] = jordan(A);
Pre(:,2) = imag(P(:,2));
Pre(:,3) = real(P(:,3))
Pre_inv = Pre^-1
Jre = Pre_inv * A * Pre
Bjre = Pre_inv * B

% truncation
Jre_less = Jre;
Jre_less(1,:)=[];
Jre_less(:,1)=[]
Bjre_less = Bjre;
Bjre_less(1,:)=[]

% plant parameters b,k,r
b = -2;
k = 4;
r = abs(b)/k;

% aftermath
% Q=I, R=1
% K_1 = [-1.8110966785064287110220433888041, -4.3160849593283698027780387510567, -1.8110966785064287110220433888041]
% Q=I, R=0
% K_2 = [-1.8117543951362541390851582228351, -4.3177192340462210264038924401546, -1.8117543951362541390851582228351]
% Q=0, R=1
% K_3 = [-1.85, -4.3, -1.85]
% Q=0, R=0
% K_4 = [-1.8534107402031930333817126269956, -4.0261248185776487663280116110305, -1.8534107402031930333817126269956]

% Riccati
Q = 0;
R = 0;

syms P_ [2 2]
K_ = -(inv(R+Bjre_less'*P_*Bjre_less)*Bjre_less'*P_*(Jre_less-b*eye(2)));
eqs = (Jre_less+Bjre_less*K_-b*eye(2))'*P_*(Jre_less+Bjre_less*K_-b*eye(2))-r^2*P_==-Q;

s = vpasolve(eqs, [P_], Random=true);
P_=[s.P_1_1 s.P_1_2;
    s.P_2_1 s.P_2_2];
K=[0 -(inv(R+Bjre_less'*P_*Bjre_less)*Bjre_less'*P_*(Jre_less-b*eye(2)))]*Pre^-1
e=eig(A+B*K)

% eigen check
th = 0:pi/50:2*pi;
xunit = r * cos(th)+ b;
yunit = r * sin(th);
plot(xunit, yunit);
hold on
plot(real(e),imag(e),"x")
axis equal
grid on
xlabel("Re")
ylabel("Im")
hold off
