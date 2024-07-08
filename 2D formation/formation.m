function dx = formation(t,x)
% % 20210601
% % 4 two-link manipulators
% % adaptive parameter estimation

dx = zeros(24,1);
%% ========================= Reference formation =========================
x1ss = [2.25; 1.3];
x2ss = [2.45; 1.3];
x3ss = [2.45; 1.5];
x4ss = [2.25; 1.5];

z12ss = 2;
z23ss = 2;
z34ss = 2;
z41ss = 2;
z13ss = 2*sqrt(2);

d1 = norm(z12ss);
d2 = norm(z23ss);
d3 = norm(z34ss);
d4 = norm(z41ss);
d5 = norm(z13ss);

%% ========================= system paramters =========================
a1 = [3.9000; 0.7500; 1.1250; 23.5200; 7.3500];
a2 = [3.9000; 0.7500; 1.1250; 23.5200; 7.3500];
a3 = [3.9000; 0.7500; 1.1250; 23.5200; 7.3500];
a4 = [3.9000; 0.7500; 1.1250; 23.5200; 7.3500];

L1 = [1.5; 1.5];
L2 = [1.5; 1.5];
L3 = [1.5; 1.5];
L4 = [1.5; 1.5];

base1 = [0; 0];
base2 = [6; 0];
base3 = [6; 6];
base4 = [0; 6];


%% ========================= States =========================
q1=x(1:2); dq1=x(3:4); 
q2=x(5:6); dq2=x(7:8); 
q3=x(9:10); dq3=x(11:12); 
q4=x(13:14); dq4=x(15:16); 
ahat1=x(17:18); 
ahat2=x(19:20); 
ahat3=x(21:22); 
ahat4=x(23:24);

[x1, v1] = XJ(L1,q1,dq1,base1); 
[x2, v2] = XJ(L2,q2,dq2,base2);
[x3, v3] = XJ(L3,q3,dq3,base3);
[x4, v4] = XJ(L4,q4,dq4,base4);

%% ========================= measurments =========================
z12 = x1 - x2;
z23 = x2 - x3;
z34 = x3 - x4;
z41 = x4 - x1;
z13 = x1 - x3;

e12 = norm(z12) - norm(z12ss);
e23 = norm(z23) - norm(z23ss);
e34 = norm(z34) - norm(z34ss);
e41 = norm(z41) - norm(z41ss);
e13 = norm(z13) - norm(z13ss);

%% ========================= gradient =========================
grad1 = z12/norm(z12)*e12 - z41/norm(z41)*e41 + z13/norm(z13)*e13;
grad2 = z23/norm(z23)*e23 - z12/norm(z12)*e12;
grad3 = z34/norm(z34)*e34 - z23/norm(z23)*e23 - z13/norm(z13)*e13;
grad4 = z41/norm(z41)*e41 - z34/norm(z34)*e34;


%% ========================= joint =========================
Z1 = Z_new(q1,grad1);
Z2 = Z_new(q2,grad1);
Z3 = Z_new(q3,grad1);
Z4 = Z_new(q4,grad1);


%% ========================= Jacobian =========================
Jhat1 = J(ahat1,q1);
Jhat2 = J(ahat1,q2);
Jhat3 = J(ahat1,q3);
Jhat4 = J(ahat1,q4);


%% ========================= joint-space =========================
kp = 800.0; kd = 180.0; alpha = 0.02;
kp1 = 100.0; kd1 = 60.0; alpha = 0.02;
u1 = -kp*Jhat1'*grad1 - kd*dq1;
u2 = -kp*Jhat2'*grad2 - kd*dq2;
u3 = -kp*Jhat3'*grad3 - kd*dq3;
u4 = -kp*Jhat4'*grad4 - kd*dq4;

%% ============================= Plant =========================
agent1 = EL2(a1,q1,dq1,u1);
agent2 = EL2(a2,q2,dq2,u2);
agent3 = EL2(a3,q3,dq3,u3);
agent4 = EL2(a4,q4,dq4,u4);


%% ========================== ODE Functions =========================
dx(1:4) = agent1;
dx(5:8) = agent2;
dx(9:12)= agent3;
dx(13:16)= agent4;

dx(17:18) = -Z1'*(-alpha*Z1*ahat1 - dq1);
dx(19:20) = -Z2'*(-alpha*Z2*ahat2 - dq2);
dx(21:22) = -Z3'*(-alpha*Z3*ahat3 - dq3);
dx(23:24) = -Z4'*(-alpha*Z4*ahat4 - dq4);

