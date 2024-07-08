%% PERA end-effector formation control
%% 2021.0109
clc; clear all; close all;

q1_ini = [0; 0; 0; 0; 0; -pi/2; 0;     0; 0; 0; 0; 0; 0; 0];
q2_ini = q1_ini; 
q3_ini = q1_ini; 
q4_ini = q1_ini; 
theta_hat_ini = kron(ones(4,1),[0.3; 0.2; 0.1]);
x_ini = [q1_ini; q2_ini; q3_ini; q4_ini; theta_hat_ini; zeros(28,1)];
[t,x] = ode15s('pera_ODE',[0 50], x_ini);

%% Desired shape
theta = [0.32; 0.28; 0.2];

d12 = 0.4;
d13 = 0.4;
d32 = 0.4;
d41 = 0.4;
d42 = 0.4;
d43 = 0.4;

%% Base
base1 = [0, 0, 0, 0];  
base2 = [0.5,  0.0, 0,  0];   
base3 = [0.5,  0.5, 0,  0];
base4 = [0.0,  0.5, 0,  0];


%% Data
pos1 = zeros(3,length(t));
pos2 = zeros(3,length(t));
pos3 = zeros(3,length(t));
pos4 = zeros(3,length(t));
% pos5 = zeros(3,length(t));
rank_N = zeros(1,length(t));
pos12 = zeros(length(t),1);
pos13 = zeros(length(t),1);
pos32 = zeros(length(t),1);
pos41 = zeros(length(t),1);
pos42 = zeros(length(t),1);
pos43 = zeros(length(t),1);

B = [1  1  0 -1  0  0;
    -1  0 -1  0 -1  0;
     0 -1  1  0  0 -1;
     0  0  0  1  1  1];
B_bar = kron(B,eye(3)); % 12*18


J = zeros(12,28,length(t));
Dz = zeros(18,6,length(t));
JTDz = zeros(28,6,length(t));
rank_all = zeros(3,length(t));
tau = zeros(28,length(t));
u = zeros(28,length(t));
e_all = zeros(6,length(t));

for i=1:length(t)
    pos1(1:3,i) = pera_HT(x(i,1:7),base1);
    pos2(1:3,i) = pera_HT(x(i,15:21),base2);
    pos3(1:3,i) = pera_HT(x(i,29:35),base3);
    pos4(1:3,i) = pera_HT(x(i,43:49),base4);
   
    
    pos12(i) = norm(pos1(1:3,i) - pos2(1:3,i));
    pos13(i) = norm(pos1(1:3,i) - pos3(1:3,i));
    pos32(i) = norm(pos3(1:3,i) - pos2(1:3,i));
    pos41(i) = norm(pos4(1:3,i) - pos1(1:3,i));
    pos42(i) = norm(pos4(1:3,i) - pos2(1:3,i));
    pos43(i) = norm(pos4(1:3,i) - pos3(1:3,i));
    
    theta_hat_all = x(i,71:85);
    theta1_hat = theta_hat_all(1:3);
    theta2_hat = theta_hat_all(4:6);
    theta3_hat = theta_hat_all(7:9);
    theta4_hat = theta_hat_all(10:12);
    %     rank_N(i) = rank(pera_J(x(i,1:7),theta1_hat)) + rank(pera_J(x(i,15:21),theta2_hat)) + rank(pera_J(x(i,29:35),theta3_hat)) + rank(pera_J(x(i,43:49),theta4_hat));
    
    J(:,:,i) = blkdiag(pera_J(x(i,1:7),theta1_hat),pera_J(x(i,15:21),theta2_hat),pera_J(x(i,29:35),theta3_hat),pera_J(x(i,43:49),theta4_hat));
    Dz(:,:,i) = blkdiag((pos1(1:3,i) - pos2(1:3,i))/pos12(i), ...
                        (pos1(1:3,i) - pos3(1:3,i))/pos13(i), ...
                        (pos3(1:3,i) - pos2(1:3,i))/pos32(i), ...
                        (pos4(1:3,i) - pos1(1:3,i))/pos41(i), ... 
                        (pos4(1:3,i) - pos2(1:3,i))/pos42(i), ...
                        (pos4(1:3,i) - pos3(1:3,i))/pos43(i));
    JTDz(:,:,i) = J(:,:,i)'*B_bar*Dz(:,:,i);
    rank_all(:,i) = [rank(J(:,:,i)); rank(Dz(:,:,i)); rank(JTDz(:,:,i))];
    tau(:,i) = JTDz(:,:,i)*[pos12(i)-norm(d12), pos13(i)-norm(d13), pos32(i)-norm(d32), ...
        pos41(i)-norm(d41), pos42(i)-norm(d42), pos43(i)-norm(d43)]';
    Kp = 30; Kd = 30;
    u(:,i) = -Kp*tau(:,i) - Kd*[x(i,8:14), x(i,22:28), x(i,36:42), x(i,50:56)]';
    
    JJT(:,:,i) = J(:,:,i)*J(:,:,i)';
end

%clf;

figure
subplot(2,1,1)
plot(t,x(:,1:7),'linewidth',1.2);
hold on
grid on
box off
plot(t,x(:,15:21),'--','linewidth',1.2);
plot(t,x(:,29:35),'-.','linewidth',1.2);
plot(t,x(:,43:49),':','linewidth',1.2);
ylabel('Joint position (rad)','FontSize',12,'Interpreter','latex');

subplot(2,1,2)
plot(t,x(:,8:14),'linewidth',1.2);
hold on
grid on
box off
plot(t,x(:,22:28),'--','linewidth',1.2);
plot(t,x(:,36:42),'-.','linewidth',1.2);
plot(t,x(:,50:56),':','linewidth',1.2);
xlabel('Time (s)','FontSize',14,'Interpreter','latex');
ylabel('Joint velocity (rad/s)','FontSize',12,'Interpreter','latex');


figure
plot(t,pos12(:) - norm(d12), 'linewidth',1.2);
hold on
grid on
plot(t,pos13(:) - norm(d13), 'linewidth',1.2);
plot(t,pos32(:) - norm(d32), 'linewidth',1.2);
plot(t,pos41(:) - norm(d41), 'linewidth',1.2);
plot(t,pos42(:) - norm(d42), 'linewidth',1.2);
plot(t,pos43(:) - norm(d43), 'linewidth',1.2);
legend3 = legend('$\|z_{1}\| - \|z_{1}^{\star}\|$','$\|z_{2}\| - \|z_{2}^{\star}\|$','$\|z_{3}\| - \|z_{3}^{\star}\|$','$\|z_{4}\| - \|z_{4}^{\star}\|$','$\|z_{5}\| - \|z_{5}^{\star}\|$','$\|z_{6}\| - \|z_{6}^{\star}\|$');
set(legend3,'Interpreter','latex','FontSize',10,'EdgeColor',[0.8 0.8 0.8]);
xlabel('Time (s)','FontSize',14,'Interpreter','latex');
ylabel('Error (m)','FontSize',14,'Interpreter','latex');


figure
plot3(pos1(1,:),pos1(2,:),pos1(3,:),'linewidth',1.2);
hold on
grid on
plot3(pos2(1,:),pos2(2,:),pos2(3,:),'linewidth',1.2);
plot3(pos3(1,:),pos3(2,:),pos3(3,:),'linewidth',1.2);
plot3(pos4(1,:),pos4(2,:),pos4(3,:),'linewidth',1.2);
line([pos1(1,end) pos2(1,end)], [pos1(2,end) pos2(2,end)], [pos1(3,end) pos2(3,end)], 'Color','[0.5 0.5 0.5]','LineStyle','--', 'linewidth',1.2);
line([pos1(1,end) pos3(1,end)], [pos1(2,end) pos3(2,end)], [pos1(3,end) pos3(3,end)], 'Color','[0.5 0.5 0.5]','LineStyle','--', 'linewidth',1.2);
line([pos3(1,end) pos2(1,end)], [pos3(2,end) pos2(2,end)], [pos3(3,end) pos2(3,end)], 'Color','[0.5 0.5 0.5]','LineStyle','--', 'linewidth',1.2);
line([pos4(1,end) pos1(1,end)], [pos4(2,end) pos1(2,end)], [pos4(3,end) pos1(3,end)], 'Color','[0.5 0.5 0.5]','LineStyle','--', 'linewidth',1.2);
line([pos4(1,end) pos2(1,end)], [pos4(2,end) pos2(2,end)], [pos4(3,end) pos2(3,end)], 'Color','[0.5 0.5 0.5]','LineStyle','--', 'linewidth',1.2);
line([pos4(1,end) pos3(1,end)], [pos4(2,end) pos3(2,end)], [pos4(3,end) pos3(3,end)], 'Color','[0.5 0.5 0.5]','LineStyle','--', 'linewidth',1.2);
plot3(pos1(1,1),pos1(2,1),pos1(3,1),'x', 'Color','[0 0.4470 0.7410]', 'linewidth',1.2); % start point
plot3(pos2(1,1),pos2(2,1),pos2(3,1),'x', 'Color','[0.8500 0.3250 0.0980]', 'linewidth',1.2); % start point
plot3(pos3(1,1),pos3(2,1),pos3(3,1),'x', 'Color','[0.9290 0.6940 0.1250]', 'linewidth',1.2); % start point
plot3(pos4(1,1),pos4(2,1),pos4(3,1),'x', 'Color','[0.4940 0.1840 0.5560]', 'linewidth',1.2); % start point
plot3(pos1(1,end),pos1(2,end),pos1(3,end),'o', 'Color','[0 0.4470 0.7410]', 'linewidth',1.2); % start point
plot3(pos2(1,end),pos2(2,end),pos2(3,end),'o', 'Color','[0.8500 0.3250 0.0980]', 'linewidth',1.2); % start point
plot3(pos3(1,end),pos3(2,end),pos3(3,end),'o', 'Color','[0.9290 0.6940 0.1250]', 'linewidth',1.2); % start point
plot3(pos4(1,end),pos4(2,end),pos4(3,end),'o', 'Color','[0.4940 0.1840 0.5560]', 'linewidth',1.2); % start point
legend4 = legend('$x_1$','$x_2$','$x_3$','$x_4$');
set(legend4,...
    'Position',[0.723256546632676 0.72904761432466 0.136785715375628 0.211190480913435],...
    'Interpreter','latex',...
    'EdgeColor',[0.8 0.8 0.8]);
xlabel('X (m)','FontSize',14,'Interpreter','latex');
ylabel('Y (m)','FontSize',14,'Interpreter','latex');
zlabel('Z (m)','FontSize',14,'Interpreter','latex');

%% plot parameter estimate
figure
plot(t,x(:,57:68), 'linewidth',1.2);
hold on
grid on
xlabel('Time (s)','FontSize',14,'Interpreter','latex');
ylabel('$\hat{a}$','FontSize',14,'Interpreter','latex');
