function output = dJhat(Lhat,dLhat,q,dq)
%% subfunction: \dot{hat{J}}(q)
%% Lhat: length estimate
%% dLhat: derivative of length estimate
%% q: joint angle
%% dq: joint velocity

%% states
Lhat1 = Lhat(1); Lhat2 = Lhat(2);
dLhat1 = dLhat(1); dLhat2 = Lhat(2);
q1 = q(1); q2 = q(2);
dq1 = dq(1); dq2 = dq(2);

output = [-dLhat1*sin(q1)-dLhat2*sin(q1+q2), -dLhat2*sin(q1+q2); ...
           dLhat1*cos(q1)+dLhat2*cos(q1+q2),  dLhat2*cos(q1+q2)] ...
       + ...
         [-Lhat1*cos(q1)*dq1-Lhat2*cos(q1+q2)*(dq1+dq2), -Lhat2*cos(q1+q2)*(dq1+dq2); ...
          -Lhat1*sin(q1)*dq1-Lhat2*sin(q1+q2)*(dq1+dq2), -Lhat2*sin(q1+q2)*(dq1+dq2)];
end