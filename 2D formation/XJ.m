function [X, vel] = XJ(L,q,dq,base)
%% subfunction: kinematics, task-space position
%% L: length (system parameter)
%% q: joint angle
%% base: position of base

%% states
L1 = L(1); L2 = L(2);
q1 = q(1); q2 = q(2);

P = [L1*cos(q1)+L2*cos(q1+q2);
     L1*sin(q1)+L2*sin(q1+q2)];
 
X = P + base;

J = [-L1*sin(q1)-L2*sin(q1+q2), -L2*sin(q1+q2);
      L1*cos(q1)+L2*cos(q1+q2),  L2*cos(q1+q2)];
vel = J*dq;
       
end