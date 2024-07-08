function output = EL2(theta,position,velocity,torque)
%% subfunction: 2-DOF EL in state space

theta_1 = theta(1); 
theta_2 = theta(2);
theta_3 = theta(3); 
theta_4 = theta(4); 
theta_5 = theta(5);

q = position;
q_1 = q(1);
q_2 = q(2);
dq = velocity;
dq_1 = dq(1);
dq_2 = dq(2);

M = [theta_1 + 2*theta_3*cos(q_2), theta_2 + theta_3*cos(q_2); 
     theta_2 + theta_3*cos(q_2), theta_2];

C = [-theta_3*sin(q_2)*dq_2, -theta_3*sin(q_2)*(dq_1 + dq_2); 
      theta_3*sin(q_2)*dq_1, 0];
  
g = [theta_4*cos(q_1) + theta_5*cos(q_1 + q_2);
     theta_5*cos(q_1+q_2)];

output = [velocity;
          M^(-1)*(torque - C*dq - 0*g)];
end