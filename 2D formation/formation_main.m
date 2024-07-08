% Simulation of 4 two-link manipulators
% tested @MATLAB 2020a

clc; clear all; close all;

ini1 = [0 pi/3 0 0];
ini2 = [pi/2 pi/3 0 0];
ini3 = [pi pi/3 0 0];
ini4 = [3*pi/2 pi/3 0 0];
[t,x] = ode15s('formation', [0 30], [ini1, ini2, ini3, ini4, 2 2 2 2 2 2 2 2]);

%% Reference formation
x1ss = [2.25; 1.3];
x2ss = [2.45; 1.3];
x3ss = [2.45; 1.5];
x4ss = [2.25; 1.5];

z12ss = 2;
z23ss = 2;
z34ss = 2;
z41ss = 2;
z13ss = 2*sqrt(2);


d1 = norm(z12ss);
d2 = norm(z23ss);
d3 = norm(z34ss);
d4 = norm(z41ss);
d5 = norm(z13ss);

L1 = [1.5; 1.5];
L2 = [1.5; 1.5];
L3 = [1.5; 1.5];
L4 = [1.5; 1.5];

base1 = [0; 0];
base2 = [6; 0];
base3 = [6; 6];
base4 = [0; 6];

pos = [];
vel = [];
ed = [];

for i=1:length(t)
    q1 = x(i,1:2)';
    dq1 = x(i,3:4)';
    [x1, v1] = XJ(L1,q1,dq1,base1); %x1 = X(L1,q1,base1); v1 = J(L1,q1)*dq1;
    q2 = x(i,5:6)';
    dq2 = x(i,7:8)';
    [x2, v2] = XJ(L2,q2,dq2,base2); 
    q3 = x(i,9:10)';
    dq3 = x(i,11:12)';
    [x3, v3] = XJ(L3,q3,dq3,base3); 
    q4 = x(i,13:14)';
    dq4 = x(i,15:16)';
    [x4, v4] = XJ(L4,q4,dq4,base4); 
    pos = [pos; x1', x2', x3', x4'];
    vel = [vel; v1', v2', v3', v4'];
    ed = [ed; norm(x1-x2), norm(x2-x3), norm(x3-x4), norm(x4-x1), norm(x1-x3)];
end


%% plot end-effector position and velocity
figure
subplot(2,1,1)
plot(t,pos(:,1),'Color',[0 0.4470 0.7410],'linewidth',1.2);
hold on
grid on
box off
plot(t,pos(:,2), '--','Color',[0 0.4470 0.7410],'linewidth',1.2);
plot(t,pos(:,3), 'Color',[0.8500 0.3250 0.0980],'linewidth',1.2);
plot(t,pos(:,4), '--','Color',[0.8500 0.3250 0.0980],'linewidth',1.2);
plot(t,pos(:,5), 'Color',[0.9290 0.6940 0.1250],'linewidth',1.2);
plot(t,pos(:,6), '--','Color',[0.9290 0.6940 0.1250],'linewidth',1.2);
plot(t,pos(:,7), 'Color',[0.4940 0.1840 0.5560],'linewidth',1.2);
plot(t,pos(:,8), '--','Color',[0.4940 0.1840 0.5560],'linewidth',1.2);
legend1 = legend('$x_{11}$','$x_{12}$','$x_{21}$','$x_{22}$','$x_{31}$','$x_{32}$','$x_{41}$','$x_{42}$');
set(legend1,...
    'Position',[0.439849738160025 0.877672295185193 0.466232095151713 0.072428573989868],...
    'Orientation','horizontal',...
    'NumColumns',4,...
    'Interpreter','latex',...
    'EdgeColor',[0.8 0.8 0.8]);
% xlabel('Time (s)');
ylabel('End-effector position (m)','FontSize',10,'Interpreter','latex');
% xlim([0 50]);

subplot(2,1,2)
plot(t,vel(:,1),'Color',[0 0.4470 0.7410],'linewidth',1.2);
hold on
grid on
box off
plot(t,vel(:,2), '--','Color',[0 0.4470 0.7410],'linewidth',1.2);
plot(t,vel(:,3), 'Color',[0.8500 0.3250 0.0980],'linewidth',1.2);
plot(t,vel(:,4), '--','Color',[0.8500 0.3250 0.0980],'linewidth',1.2);
plot(t,vel(:,5), 'Color',[0.9290 0.6940 0.1250],'linewidth',1.2);
plot(t,vel(:,6), '--','Color',[0.9290 0.6940 0.1250],'linewidth',1.2);
plot(t,vel(:,7), 'Color',[0.4940 0.1840 0.5560],'linewidth',1.2);
plot(t,vel(:,8), '--','Color',[0.4940 0.1840 0.5560],'linewidth',1.2);
xlabel('Time (s)','FontSize',14,'Interpreter','latex');
ylabel('End-effector velocity (m/s)','FontSize',10,'Interpreter','latex');
% xlim([0 50]);

%% plot joint position and velocity
figure
subplot(2,1,1)
plot(t,x(:,1),'Color',[0 0.4470 0.7410],'linewidth',1.2);
hold on
grid on
box off
plot(t,x(:,2), '--','Color',[0 0.4470 0.7410],'linewidth',1.2);
plot(t,x(:,5), 'Color',[0.8500 0.3250 0.0980],'linewidth',1.2);
plot(t,x(:,6), '--','Color',[0.8500 0.3250 0.0980],'linewidth',1.2);
plot(t,x(:,9), 'Color',[0.9290 0.6940 0.1250],'linewidth',1.2);
plot(t,x(:,10), '--','Color',[0.9290 0.6940 0.1250],'linewidth',1.2);
plot(t,x(:,13), 'Color',[0.4940 0.1840 0.5560],'linewidth',1.2);
plot(t,x(:,14), '--','Color',[0.4940 0.1840 0.5560],'linewidth',1.2);
legend2 = legend('$q_{11}$','$q_{12}$','$q_{21}$','$q_{22}$','$q_{31}$','$q_{32}$','$q_{41}$','$q_{42}$');
set(legend2,...
    'Position',[0.444675580996204 0.892910390423288 0.456580409479355 0.072428573989868],...
    'Orientation','horizontal',...
    'NumColumns',4,...
    'Interpreter','latex',...
    'EdgeColor',[0.8 0.8 0.8]);
% xlabel('Time (s)');
ylabel('Joint position (rad)','FontSize',12,'Interpreter','latex');
% xlim([0 50]);

subplot(2,1,2)
plot(t,x(:,3),'Color',[0 0.4470 0.7410],'linewidth',1.2);
hold on
grid on
box off
plot(t,x(:,4), '--','Color',[0 0.4470 0.7410],'linewidth',1.2);
plot(t,x(:,7), 'Color',[0.8500 0.3250 0.0980],'linewidth',1.2);
plot(t,x(:,8), '--','Color',[0.8500 0.3250 0.0980],'linewidth',1.2);
plot(t,x(:,11), 'Color',[0.9290 0.6940 0.1250],'linewidth',1.2);
plot(t,x(:,12), '--','Color',[0.9290 0.6940 0.1250],'linewidth',1.2);
plot(t,x(:,15), 'Color',[0.4940 0.1840 0.5560],'linewidth',1.2);
plot(t,x(:,16), '--','Color',[0.4940 0.1840 0.5560],'linewidth',1.2);
xlabel('Time (s)','FontSize',14,'Interpreter','latex');
ylabel('Joint velocity (rad/s)','FontSize',12,'Interpreter','latex');
% xlim([0 50]);

%% plot inner distance error
figure
plot(t, ed(:,1)-d1*ones(length(t),1), 'linewidth',1.2);
hold on
grid on
plot(t, ed(:,2)-d2*ones(length(t),1), 'linewidth',1.2);
plot(t, ed(:,3)-d3*ones(length(t),1), 'linewidth',1.2);
plot(t, ed(:,4)-d4*ones(length(t),1), 'linewidth',1.2);
plot(t, ed(:,5)-d5*ones(length(t),1), 'linewidth',1.2);

legend3 = legend('$\|z_{1}\| - \|z_{1}^{\star}\|$','$\|z_{2}\| - \|z_{2}^{\star}\|$','$\|z_{3}\| - \|z_{3}^{\star}\|$','$\|z_{4}\| - \|z_{4}^{\star}\|$','$\|z_{5}\| - \|z_{5}^{\star}\|$');
set(legend3,'Interpreter','latex','FontSize',10,'EdgeColor',[0.8 0.8 0.8]);
xlabel('Time (s)','FontSize',14,'Interpreter','latex');
ylabel('Error (m)','FontSize',14,'Interpreter','latex');

%% plot X-Y
figure
plot(pos(:,1),pos(:,2), 'linewidth',1.2);
hold on; grid on;
plot(pos(:,3),pos(:,4), 'linewidth',1.2);
plot(pos(:,5),pos(:,6), 'linewidth',1.2);
plot(pos(:,7),pos(:,8), 'linewidth',1.2);
% plot(pos(1,1),pos(1,2), 'x','linewidth',1.2);
line([pos(end,1) pos(end,3)], [pos(end,2) pos(end,4)], 'Color','[0.5 0.5 0.5]','LineStyle','--', 'linewidth',1.2);
line([pos(end,3) pos(end,5)], [pos(end,4) pos(end,6)], 'Color','[0.5 0.5 0.5]','LineStyle','--', 'linewidth',1.2);
line([pos(end,5) pos(end,7)], [pos(end,6) pos(end,8)], 'Color','[0.5 0.5 0.5]','LineStyle','--', 'linewidth',1.2);
line([pos(end,7) pos(end,1)], [pos(end,8) pos(end,2)], 'Color','[0.5 0.5 0.5]','LineStyle','--', 'linewidth',1.2);
line([pos(end,1) pos(end,5)], [pos(end,2) pos(end,6)], 'Color','[0.5 0.5 0.5]','LineStyle','--', 'linewidth',1.2);
plot(pos(1,1),pos(1,2),'x', 'Color','[0 0.4470 0.7410]', 'linewidth',1.2); % start point
plot(pos(1,3),pos(1,4),'x', 'Color','[0.8500 0.3250 0.0980]', 'linewidth',1.2); % start point
plot(pos(1,5),pos(1,6),'x', 'Color','[0.9290 0.6940 0.1250]', 'linewidth',1.2); % start point
plot(pos(1,7),pos(1,8),'x', 'Color','[0.4940 0.1840 0.5560]', 'linewidth',1.2); % start point
plot(pos(end,1),pos(end,2),'o', 'Color','[0 0.4470 0.7410]', 'linewidth',1.2); % end point
plot(pos(end,3),pos(end,4),'o', 'Color','[0.8500 0.3250 0.0980]', 'linewidth',1.2); % end point
plot(pos(end,5),pos(end,6),'o', 'Color','[0.9290 0.6940 0.1250]', 'linewidth',1.2); % end point
plot(pos(end,7),pos(end,8),'o', 'Color','[0.4940 0.1840 0.5560]', 'linewidth',1.2); % end point
legend4 = legend('$x_1$','$x_2$','$x_3$','$x_4$');
set(legend4,...
    'Position',[0.774776799750116 0.660555552285815 0.111889866916551 0.150238098507836],...
    'Interpreter','latex',...
    'EdgeColor',[0.8 0.8 0.8]);
xlabel('X (m)','FontSize',14,'Interpreter','latex');
ylabel('Y (m)','FontSize',14,'Interpreter','latex');
xlim([1 6]);
ylim([1 5]);

%% plot estimated parameters
figure
plot(t,x(:,17),'Color',[0 0.4470 0.7410],'linewidth',1.2);
hold on
grid on
plot(t,x(:,18), '--','Color',[0 0.4470 0.7410],'linewidth',1.2);
plot(t,x(:,19), 'Color',[0.8500 0.3250 0.0980],'linewidth',1.2);
plot(t,x(:,20), '--','Color',[0.8500 0.3250 0.0980],'linewidth',1.2);
plot(t,x(:,21), 'Color',[0.9290 0.6940 0.1250],'linewidth',1.2);
plot(t,x(:,22), '--','Color',[0.9290 0.6940 0.1250],'linewidth',1.2);
plot(t,x(:,23), 'Color',[0.4940 0.1840 0.5560],'linewidth',1.2);
plot(t,x(:,24), '--','Color',[0.4940 0.1840 0.5560],'linewidth',1.2);
legend1 = legend('$\hat{a}_{11}$','$\hat{a}_{12}$','$\hat{a}_{21}$','$\hat{a}_{22}$','$\hat{a}_{31}$','$\hat{a}_{32}$','$\hat{a}_{41}$','$\hat{a}_{42}$');
set(legend1,...
    'Position',[0.427097690121392 0.895672295648407 0.473879048371837 0.0792857159205843],...
    'Orientation','horizontal',...
    'NumColumns',4,...
    'Interpreter','latex',...
    'EdgeColor',[0.8 0.8 0.8]);
xlabel('Time (s)');
ylabel('$\hat{a}$','FontSize',14,'Interpreter','latex');
