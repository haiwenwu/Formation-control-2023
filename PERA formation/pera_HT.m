function position = pera_HT(q,base)
%% Compute Cartesian coordinate
%% 7dof

%% States
x1=base(1); x2=base(2); x3=base(3); z=base(4); 

%% System parameters
theta = [0.32; 0.28; 0.2];
theta_1=theta(1); theta_2=theta(2); theta_3=theta(3);
q_1 = q(1); q_2 = q(2); q_3 = q(3); q_4 = q(4); q_5 = q(5); q_6 = q(6); q_7 = q(7);

pos_11 = - cos(q_1)*sin(q_2); %% 0.32
pos_12 = sin(q_1)*sin(q_3)*sin(q_4) - cos(q_1)*cos(q_4)*sin(q_2) - cos(q_1)*cos(q_2)*cos(q_3)*sin(q_4); %% 0.28
pos_13 = - cos(q_3)*cos(q_5)*sin(q_1)*sin(q_7) - cos(q_1)*cos(q_2)*cos(q_5)*sin(q_3)*sin(q_7) + cos(q_1)*cos(q_4)*cos(q_7)*sin(q_2)*sin(q_6) - cos(q_3)*cos(q_6)*cos(q_7)*sin(q_1)*sin(q_5) + cos(q_1)*sin(q_2)*sin(q_4)*sin(q_5)*sin(q_7) + cos(q_4)*sin(q_1)*sin(q_3)*sin(q_5)*sin(q_7) - cos(q_7)*sin(q_1)*sin(q_3)*sin(q_4)*sin(q_6) - cos(q_1)*cos(q_2)*cos(q_3)*cos(q_4)*sin(q_5)*sin(q_7) + cos(q_1)*cos(q_2)*cos(q_3)*cos(q_7)*sin(q_4)*sin(q_6) - cos(q_1)*cos(q_2)*cos(q_6)*cos(q_7)*sin(q_3)*sin(q_5) - cos(q_1)*cos(q_5)*cos(q_6)*cos(q_7)*sin(q_2)*sin(q_4) - cos(q_4)*cos(q_5)*cos(q_6)*cos(q_7)*sin(q_1)*sin(q_3) + cos(q_1)*cos(q_2)*cos(q_3)*cos(q_4)*cos(q_5)*cos(q_6)*cos(q_7);

pos_21 = - sin(q_1)*sin(q_2); %% 0.32
pos_22 = - cos(q_4)*sin(q_1)*sin(q_2) - cos(q_1)*sin(q_3)*sin(q_4) - cos(q_2)*cos(q_3)*sin(q_1)*sin(q_4); %% 0.28
pos_23 = cos(q_1)*cos(q_3)*cos(q_5)*sin(q_7) + cos(q_1)*cos(q_3)*cos(q_6)*cos(q_7)*sin(q_5) - cos(q_2)*cos(q_5)*sin(q_1)*sin(q_3)*sin(q_7) - cos(q_1)*cos(q_4)*sin(q_3)*sin(q_5)*sin(q_7) + cos(q_4)*cos(q_7)*sin(q_1)*sin(q_2)*sin(q_6) + cos(q_1)*cos(q_7)*sin(q_3)*sin(q_4)*sin(q_6) + sin(q_1)*sin(q_2)*sin(q_4)*sin(q_5)*sin(q_7) + cos(q_1)*cos(q_4)*cos(q_5)*cos(q_6)*cos(q_7)*sin(q_3) - cos(q_2)*cos(q_3)*cos(q_4)*sin(q_1)*sin(q_5)*sin(q_7) + cos(q_2)*cos(q_3)*cos(q_7)*sin(q_1)*sin(q_4)*sin(q_6) - cos(q_2)*cos(q_6)*cos(q_7)*sin(q_1)*sin(q_3)*sin(q_5) - cos(q_5)*cos(q_6)*cos(q_7)*sin(q_1)*sin(q_2)*sin(q_4) + cos(q_2)*cos(q_3)*cos(q_4)*cos(q_5)*cos(q_6)*cos(q_7)*sin(q_1); 

pos_31 = cos(q_2); %% 0.32
pos_32 = cos(q_2)*cos(q_4) - cos(q_3)*sin(q_2)*sin(q_4); %% 0.28
pos_33 = - cos(q_2)*cos(q_4)*cos(q_7)*sin(q_6) - cos(q_5)*sin(q_2)*sin(q_3)*sin(q_7) - cos(q_2)*sin(q_4)*sin(q_5)*sin(q_7) + cos(q_2)*cos(q_5)*cos(q_6)*cos(q_7)*sin(q_4) - cos(q_3)*cos(q_4)*sin(q_2)*sin(q_5)*sin(q_7) + cos(q_3)*cos(q_7)*sin(q_2)*sin(q_4)*sin(q_6) - cos(q_6)*cos(q_7)*sin(q_2)*sin(q_3)*sin(q_5) + cos(q_3)*cos(q_4)*cos(q_5)*cos(q_6)*cos(q_7)*sin(q_2);

pos = [pos_11*theta_1 + pos_12*theta_2 + pos_13*theta_3;
    pos_21*theta_1 + pos_22*theta_2 + pos_23*theta_3;
    pos_31*theta_1 + pos_32*theta_2 + pos_33*theta_3];

%% Base Homogeneous Transformation
base_R = [cos(z) -sin(z) 0;
          sin(z)  cos(z) 0;
               0       0 1];
base_T = [x1; x2; x3];
base_HT = [base_R, base_T; zeros(1,3), 1];

%% End-effector position
p = base_HT*[pos; 1];
position = p(1:3);

end