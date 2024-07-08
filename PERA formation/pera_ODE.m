function dx = pera_ODE(t,x)
%% PERA end-effector formation control
%% 2021.0109

%% SYSTEM PARAMETERS
g_C = 9.81;
d3 = 0.32; d5 = 0.28; a7 = 0.2;
d = [0 0 d3 0 d5 0 0];
m = [0; 0; 0; 2.9; 0; 0.8; 0.2]; 
d_1=d(1); d_2=d(2); d_3=d(3); d_4=d(4); d_5=d(5); d_6=d(6); d_7=d(7);
m_1=m(1); m_2=m(2); m_3=m(3); m_4=m(4); m_5=m(5); m_6=m(6); m_7=m(7);

theta = [0.32; 0.28; 0.2];
a = [g_C*(m_4+m_5+m_6+m_7)*d_3; g_C*(m_6+m_7)*d_5];


%% STATES
dx = zeros(96,1);
q1 = x(1:7); dq1 = x(8:14); % agent 1
q2 = x(15:21); dq2 = x(22:28); % agent 2
q3 = x(29:35); dq3 = x(36:42); % agent 3
q4 = x(43:49); dq4 = x(50:56); % agent 4
theta_hat_all = x(57:68);
theta1_hat = theta_hat_all(1:3);
theta2_hat = theta_hat_all(4:6);
theta3_hat = theta_hat_all(7:9);
theta4_hat = theta_hat_all(10:12);


%% Desired shape
d12 = 0.4;
d13 = 0.4;
d32 = 0.4;
d41 = 0.4;
d42 = 0.4;
d43 = 0.4;


%% Base
base1 = [0, 0, 0, 0];  
base2 = [0.5,  0.0, 0,  0];  
base3 = [0.5,  0.5, 0,  0];
base4 = [0.0,  0.5, 0,  0];


%% End-effector position
x1 = pera_HT(q1,base1);
x2 = pera_HT(q2,base2);
x3 = pera_HT(q3,base3);
x4 = pera_HT(q4,base4);

%% Controller
z12 = x1 - x2; e12 = norm(z12) - norm(d12);
z13 = x1 - x3; e13 = norm(z13) - norm(d13);
z32 = x3 - x2; e32 = norm(z32) - norm(d32);
z41 = x4 - x1; e41 = norm(z41) - norm(d41);
z42 = x4 - x2; e42 = norm(z42) - norm(d42);
z43 = x4 - x3; e43 = norm(z43) - norm(d43);
e_all = [e12; e13; e32; e41; e42; e43];


e_agent1 = z12/norm(z12)*e12 + z13/norm(z13)*e13 - z41/norm(z41)*e41;
e_agent2 = - z12/norm(z12)*e12 - z32/norm(z32)*e32 - z42/norm(z42)*e42;
e_agent3 = z32/norm(z32)*e32 - z13/norm(z13)*e13 - z43/norm(z43)*e43;
e_agent4 = z41/norm(z41)*e41 + z42/norm(z42)*e42 + z43/norm(z43)*e43;

B = [1  1  0 -1  0  0;
     -1  0 -1  0 -1  0;
      0 -1  1  0  0 -1;
      0  0  0  1  1  1];
B_bar = kron(B,eye(3));

Dz = blkdiag(z12/norm(z12), z13/norm(z13), z32/norm(z32), ...
    z41/norm(z41), z42/norm(z42), z43/norm(z43));

F = [z12 z13 zeros(3,1) -z41 zeros(3,1) zeros(3,1);
      -z12 zeros(3,1) -z32 zeros(3,1) -z42 zeros(3,1);
      zeros(3,1) -z13 z32 zeros(3,1) -z42 zeros(3,1);
      zeros(3,1) zeros(3,1) zeros(3,1) z41 z42 z43]; % 12*6

J1 = pera_J(q1,theta1_hat);
J2 = pera_J(q2,theta2_hat);
J3 = pera_J(q3,theta3_hat);
J4 = pera_J(q4,theta4_hat);
J = [J1, zeros(3,7), zeros(3,7), zeros(3,7); 
    zeros(3,7), J2, zeros(3,7), zeros(3,7);
    zeros(3,7), zeros(3,7), J3, zeros(3,7);
    zeros(3,7), zeros(3,7), zeros(3,7), J4]; %% 12 * 28
rank(J'*B_bar*Dz); size(J'*B_bar*Dz);


Kp = 120; Kd = 10; KI = 1; alpha=0.01; lambda = 1;

u1 = -Kp*J1'*e_agent1 - Kd*dq1 + KI*x(69:75);
u2 = -Kp*J2'*e_agent2 - Kd*dq2 + KI*x(76:82);
u3 = -Kp*J3'*e_agent3 - Kd*dq3 + KI*x(83:89);
u4 = -Kp*J4'*e_agent4 - Kd*dq4 + KI*x(90:96);


%% Agents
agent1 = pera_EL(q1,dq1,u1);
agent2 = pera_EL(q2,dq2,u2);
agent3 = pera_EL(q3,dq3,u3);
agent4 = pera_EL(q4,dq4,u4);

%% ODE functions
dx(1:14) = agent1;
dx(15:28) = agent2;
dx(29:42) = agent3;
dx(43:56) = agent4;
dx(57:59) = -lambda*pera_Znew(q1,e_agent1)'*( alpha*pera_Znew(q1,e_agent1)*theta1_hat - dq1 );
dx(60:62) = -lambda*pera_Znew(q2,e_agent2)'*( alpha*pera_Znew(q2,e_agent2)*theta2_hat - dq2 );
dx(63:65) = -lambda*pera_Znew(q3,e_agent3)'*( alpha*pera_Znew(q3,e_agent3)*theta3_hat - dq3 );
dx(66:68) = -lambda*pera_Znew(q4,e_agent4)'*( alpha*pera_Znew(q4,e_agent4)*theta4_hat - dq4 );
dx(69:75) = -KI*x(69:75) + u1;
dx(76:82) = -KI*x(76:82) + u2;
dx(83:89) = -KI*x(83:89) + u3;
dx(90:96) = -KI*x(90:96) + u4;
 

