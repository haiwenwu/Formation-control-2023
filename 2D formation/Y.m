function output = Y(position,velocity,x,y)
%% subfunction: generate regresion matrix, 2¡Á5
%% Y(u1,u2,u3,y4) = Y(position,velocity,acc,velocity) = Y(q,dq,x,y)
%% Mx + Cy + g = tau

q_em = position; % u1
q_em_1 = q_em(1);
q_em_2 = q_em(2);
dq = velocity; % u2
dq_1 = dq(1);
dq_2 = dq(2);
x_1 = x(1); % u3
x_2 = x(2);
y_1 = y(1); % u4
y_2 = y(2);

Y13 = (2*x_1 + x_2)*cos(q_em_2) - (dq_2*y_1 + dq_1*y_2 + dq_2*y_2)*sin(q_em_2);
Y23 = x_1*cos(q_em_2) + dq_1*y_1*sin(q_em_2);

output = [x_1,       x_2, Y13, cos(q_em_1), cos(q_em_1 + q_em_2);
            0, x_1 + x_2, Y23,           0, cos(q_em_1 + q_em_2)];
end