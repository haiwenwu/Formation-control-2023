function Y = pera_Y(q)

q_1=q(1); q_2=q(2); q_3=q(3); q_4=q(4); q_5=q(5); q_6=q(6); q_7=q(7);


Y11 = sin(q_1)*sin(q_2);
Y12 = (sin(q_4)*(cos(q_1)*sin(q_3) + cos(q_2)*cos(q_3)*sin(q_1)) + cos(q_4)*sin(q_1)*sin(q_2));

Y21 = -cos(q_1)*cos(q_2);
Y22 = -(cos(q_1)*cos(q_2)*cos(q_4) - cos(q_1)*cos(q_3)*sin(q_2)*sin(q_4));

Y32 = sin(q_4)*(cos(q_3)*sin(q_1) + cos(q_1)*cos(q_2)*sin(q_3));

Y42 = (cos(q_4)*(sin(q_1)*sin(q_3) - cos(q_1)*cos(q_2)*cos(q_3)) + cos(q_1)*sin(q_2)*sin(q_4));


Y = [Y11 Y12;
    Y21 Y22;
      0 Y32;
      0 Y42;
      0   0;
      0   0;
      0   0]; %% 7*2

end