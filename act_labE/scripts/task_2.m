% plant parameters
A=[0 1; -1 2];
B=[1 2; -3 3];
C=[2 1; 3 -2];
D=[-2 0; 0 1];

Bf=[1 2; -1 3];
Df=[1 0; 0 -1];
Cz=[4 2; -1 1];
Dz=[2 0; 0 1];

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