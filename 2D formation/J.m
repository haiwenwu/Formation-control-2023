function output = J(L,q)
%% subfunction: Jacobian
%% L: length (system parameter)
%% q: joint angle

%% states
L1 = L(1); L2 = L(2);
q1 = q(1); q2 = q(2);

output = [-L1*sin(q1)-L2*sin(q1+q2), -L2*sin(q1+q2);
           L1*cos(q1)+L2*cos(q1+q2),  L2*cos(q1+q2)];
end