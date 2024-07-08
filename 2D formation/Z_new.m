function output = Z_new(q,zeta)
%% subfunction: generate regression matrix, 2×2
%% J'(q)zeta = Z(q,zeta)\theta

%% states
q1 = q(1); q2 = q(2);
zeta1 = zeta(1); zeta2 = zeta(2);

output = [-sin(q1)*zeta1+cos(q1)*zeta2, -sin(q1+q2)*zeta1+cos(q1+q2)*zeta2;
           0,  -sin(q1+q2)*zeta1+cos(q1+q2)*zeta2];
 
end