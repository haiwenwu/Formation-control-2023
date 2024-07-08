function output = pera_EL(position, velocity, input)
%% pera dynamic model, i.e., EL equation
%% 7dof

q = position;
dq = velocity;
q_1 = q(1); q_2 = q(2); q_3 = q(3); q_4 = q(4); q_5 = q(5); q_6 = q(6); q_7 = q(7);
dq_1 = dq(1); dq_2 = dq(2); dq_3 = dq(3); dq_4 = dq(4); dq_5 = dq(5); dq_6 = dq(6); dq_7 = dq(7);

%% Parameters
d3 = 0.32; d5 = 0.28; a7 = 0.2;
d = [0 0 d3 0 d5 0 0];
I = 10.*[0.0286; 0.0033; 0.0053; 9.1e-4; 7.3e-4; 5e-4; 5e-4];
m = [0; 0; 0; 2.9; 0; 0.8; 0.2]; %% guess from gravity


d_1=d(1); d_2=d(2); d_3=d(3); d_4=d(4); d_5=d(5); d_6=d(6); d_7=d(7);
I_1=I(1); I_2=I(2); I_3=I(3); I_4=I(4); I_5=I(5); I_6=I(6); I_7=I(7);
m_1=m(1); m_2=m(2); m_3=m(3); m_4=m(4); m_5=m(5); m_6=m(6); m_7=m(7);

%% INERTIA (MASS) MATRIX
m11 = I_1+I_2+I_3+I_4+I_5+I_6+I_7+d_3^2*(m_4+m_5+m_6+m_7)+d_5^2*(m_6+m_7)...
    - d_3^2*(m_4+m_5+m_6+m_7)*cos(q_2)^2 - d_5^2*(m_6+m_7)*cos(q_3)^2 + 2*d_3*d_5*(m_6+m_7)*cos(q_4)...
    + d_5^2*m_6*cos(q_2)^2*cos(q_3)^2 - d_5^2*m_6*cos(q_2)^2*cos(q_4)^2 + d_5^2*m_7*cos(q_2)^2*cos(q_3)^2 ...
    + d_5^2*(m_6+m_7)*cos(q_3)^2*cos(q_4)^2 - d_5^2*m_7*cos(q_2)^2*cos(q_4)^2 ...
    - d_5^2*(m_6+m_7)*cos(q_2)^2*cos(q_3)^2*cos(q_4)^2 ...
    - 2*d_3*d_5*(m_6+m_7)*cos(q_2)^2*cos(q_4)^2 ...
    + 2*d_3*d_5*(m_6+m_7)*cos(q_2)*cos(q_3)*sin(q_2)*sin(q_4) ...
    + 2*d_5^2*(m_6+m_7)*cos(q_2)*cos(q_3)*cos(q_4)*sin(q_2)*sin(q_4);

m12 = - d_5*sin(q_3)*(m_6+m_7)*(d_3*cos(q_2)*sin(q_4) - d_5*cos(q_3)*sin(q_2)...
    + d_5*cos(q_2)*cos(q_4)*sin(q_4) + d_5*cos(q_3)*cos(q_4)^2*sin(q_2));
m13 = (I_3+I_4+I_5+I_6+I_7 + d_5^2*m_6+d_5^2*m_7)*cos(q_2)...
    - d_5^2*(m_6+m_7)*cos(q_2)*cos(q_4)^2 ...
    + d_5^2*(m_6+m_7)*cos(q_3)*cos(q_4)*sin(q_2)*sin(q_4) ...
    + d_3*d_5*(m_6+m_7)*cos(q_3)*sin(q_2)*sin(q_4);
m14 = sin(q_2)*sin(q_3)*(I_4+I_5+I_6+I_7+d_5^2*m_6+d_5^2*m_7 + d_3*d_5*m_6*cos(q_4) + d_3*d_5*m_7*cos(q_4));
m15 = (cos(q_2)*cos(q_4) - cos(q_3)*sin(q_2)*sin(q_4))*(I_5*I_6*I_7);
m16 = -(I_6+I_7)*(cos(q_5)*sin(q_2)*sin(q_3) + cos(q_2)*sin(q_4)*sin(q_5) + cos(q_3)*cos(q_4)*sin(q_2)*sin(q_5));
m17 = I_7*cos(q_2)*cos(q_4)*cos(q_6) - I_7*cos(q_3)*cos(q_6)*sin(q_2)*sin(q_4) ...
    + I_7*cos(q_2)*cos(q_5)*sin(q_4)*sin(q_6) - I_7*sin(q_2)*sin(q_3)*sin(q_5)*sin(q_6) ...
    + I_7*cos(q_3)*cos(q_4)*cos(q_5)*sin(q_2)*sin(q_6);

m21 = - d_5*sin(q_3)*(m_6+m_7)*(d_3*cos(q_2)*sin(q_4) - d_5*cos(q_3)*sin(q_2) ...
    + d_5*cos(q_2)*cos(q_4)*sin(q_4) + d_5*cos(q_3)*cos(q_4)^2*sin(q_2));
m22 = I_2+I_3+I_4+I_5+I_6+I_7+d_3^2*(m_4+m_5+m_6+m_7) ...
    + d_5^2*m_6*cos(q_3) + d_5^2*m_6*cos(q_4) ...
    + d_5^2*m_7*cos(q_3)^2 + d_5^2*m_7*cos(q_4)^2 + 2*d_3*d_5*m_6*cos(q_4) ...
    + 2*d_3*d_5*m_7*cos(q_4) - d_5^2*m_6*cos(q_3)^2*cos(q_4)^2 - d_5^2*m_7*cos(q_3)^2*cos(q_4)^2;
m23 = - d_5*sin(q_3)*sin(q_4)*(m_6+m_7)*(d_3+d_5*cos(q_4));
m24 = cos(q_3)*(I_4+I_5+I_6+I_7+d_5^2*m_6+d_5^2*m_7+d_3*d_5*m_6*cos(q_4) + d_3*d_5*m_7*cos(q_4));
m25 = sin(q_3)*sin(q_4)*(I_5+I_6+I_7);
m26 = -(cos(q_3)*cos(q_5) - cos(q_4)*sin(q_3)*sin(q_5))*(I_6+I_7);
m27 = I_7*cos(q_6)*sin(q_3)*sin(q_4) - I_7*cos(q_3)*sin(q_5)*sin(q_6) - I_7*cos(q_4)*cos(q_5)*sin(q_3)*sin(q_6);

m31 = (I_3+I_4+I_5+I_6+I_7+d_5^2*m_6+d_5^2*m_7)*cos(q_2) ...
    - d_5^2*(m_6*m_7)*cos(q_2)*cos(q_4)^2 ...
    + d_5^2*(m_6+m_7)*cos(q_3)*cos(q_4)*sin(q_2)*sin(q_4) ...
    + d_3*d_5*(m_6+m_7)*cos(q_3)*sin(q_2)*sin(q_4);
m32 = - d_5*sin(q_3)*sin(q_4)*(m_6+m_7)*(d_3+d_5*cos(q_4));
m33 = I_3+I_4+I_5+I_6+I_7 + d_5^2*m_6*sin(q_4)^2 + d_5^2*m_7*sin(q_4)^2;
m34 = 0;
m35 = (I_5+I_6+I_7)*cos(q_4);
m36 = - sin(q_4)*sin(q_5)*(I_6+I_7);
m37 = I_7*cos(q_4)*cos(q_6) + I_7*cos(q_5)*sin(q_4)*sin(q_6);

m41 = sin(q_2)*sin(q_3)*(I_4+I_5+I_6+I_7 + d_5^2*m_6+d_5^2*m_7 + d_3*d_5*m_6*cos(q_4) + d_3*d_5*m_7*cos(q_4));
m42 = cos(q_3)*(I_4+I_5+I_6+I_7 + d_5^2*m_6+d_5^2*m_7 + d_3*d_5*m_6*cos(q_4) + d_3*d_5*m_7*cos(q_4));
m43 = 0;
m44 = I_4+I_5+I_6+I_7 + d_5^2*m_6+d_5^2*m_7;
m45 = 0;
m46 = -cos(q_5)*(I_6+I_7);
m47 = - I_7*sin(q_5)*sin(q_6);

m51 = (cos(q_2)*cos(q_4) - cos(q_3)*sin(q_2)*sin(q_4))*(I_5+I_6+I_7);
m52 = sin(q_3)*sin(q_4)*(I_5+I_6+I_7);
m53 = cos(q_4)*(I_5+I_6+I_7);
m54 = 0;
m55 = I_5+I_6+I_7;
m56 = 0;
m57 = I_7*cos(q_6);

m61 = - (I_6+I_7)*(cos(q_5)*sin(q_2)*sin(q_3) + cos(q_2)*sin(q_4)*sin(q_5) + cos(q_3)*cos(q_4)*sin(q_2)*sin(q_5));
m62 = - (cos(q_3)*cos(q_5) - cos(q_4)*sin(q_3)*sin(q_5))*(I_6+I_7);
m63 = - sin(q_4)*sin(q_5)*(I_6+I_7);
m64 = - cos(q_5)*(I_6+I_7);
m65 = 0;
m66 = I_6+I_7;
m67 = 0;

m71 = I_7*cos(q_2)*cos(q_4)*cos(q_6) - I_7*cos(q_3)*cos(q_6)*sin(q_2)*sin(q_4) ...
    + I_7*cos(q_2)*cos(q_5)*sin(q_4)*sin(q_6) - I_7*sin(q_2)*sin(q_3)*sin(q_5)*sin(q_6) ...
    + I_7*cos(q_3)*cos(q_4)*cos(q_5)*sin(q_2)*sin(q_6);
m72 = I_7*cos(q_6)*sin(q_3)*sin(q_4) - I_7*cos(q_3)*sin(q_5)*sin(q_6) ...
    - I_7*cos(q_4)*cos(q_5)*sin(q_3)*sin(q_6);
m73 = I_7*cos(q_4)*cos(q_6) + I_7*cos(q_5)*sin(q_4)*sin(q_6);
m74 = - I_7*sin(q_5)*sin(q_6);
m75 = I_7*cos(q_6);
m76 = 0;
m77 = I_7;

M = [m11 m12 m13 m14 m15 m16 m17; ...
    m21 m22 m23 m24 m25 m26 m27; ...
    m31 m32 m33 m34 m35 m36 m37; ...
    m41 m42 m43 m44 m45 m46 m47; ...
    m51 m52 m53 m54 m55 m56 m57; ...
    m61 m62 m63 m64 m65 m66 m67; ...
    m71 m72 m73 m74 m75 m76 m77];
 M = eye(7);

%% CORIOLIS and CENTRIFUGAL MATRIX
C = pera_C(q,dq,I);


%% GRAVITY VECTOR
g_C = 9.81;
g1 = g_C*(m_6+m_7)*(d_5*(sin(q_4)*(cos(q_1)*sin(q_3) + cos(q_2)*cos(q_3)*sin(q_1)) + cos(q_4)*sin(q_1)*sin(q_2)) ...
    + d_3*sin(q_1)*sin(q_2)) + g_C*d_3*(m_4+m_5)*sin(q_1)*sin(q_2);
g2 = - g_C*(m_6+m_7)*(d_5*(cos(q_1)*cos(q_2)*cos(q_4) - cos(q_1)*cos(q_3)*sin(q_2)*sin(q_4)) + d_3*cos(q_1)*cos(q_2)) ...
    - g_C*d_3*(m_4+m_5)*cos(q_1)*cos(q_2);
g3 = g_C*d_5*(m_6+m_7)*sin(q_4)*(cos(q_3)*sin(q_1) + cos(q_1)*cos(q_2)*sin(q_3));
g4 = g_C*d_5*(m_6+m_7)*(cos(q_4)*(sin(q_1)*sin(q_3) - cos(q_1)*cos(q_2)*cos(q_3)) + cos(q_1)*sin(q_2)*sin(q_4));
g5 = 0;
g6 = 0;
g7 = 0;

g = [g1; g2; g3; g4; g5; g6; g7];

% a = [g_C*(m_4+m_5+m_6+m_7)*d_3; g_C*(m_6+m_7)*d_5];
% g = pera_Y(q)*a;


output = [dq;    inv(M)*(input - C*dq - g)];

end