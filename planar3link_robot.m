%% 机器人模型-平面三连杆

theta1 = 0; d1 = 0;  a1 = 1;  alpha1 = 0; 
theta2 = 0; d2 = 0;  a2 = 1;  alpha2 = 0; 
theta3 = 0; d3 = 0;  a3 = 1; alpha3 = 0; 

% 定义link1
L(1) = Link([theta1 d1 a1 alpha1  0 0]);
L(1).m = 1;
L(1).I = [0,1/12,1/12,0,0,0];
L(1).Jm = 0;
L(1).r = [1/2,0,0];
% 定义link2
L(2) = Link([theta2 d2 a2 alpha2  0 0]);
L(2).m = 1;
L(2).I = [0,1/12,1/12,0,0,0];
L(2).Jm = 0;
L(2).r = [1/2,0,0];
% 定义link3
L(3) = Link([theta3 d3 a3 alpha3  0 0]);
L(3).m = 1;
L(3).I = [0,1/12,1/12,0,0,0];
L(3).Jm = 0;
L(3).r = [1/2,0,0];

robot = SerialLink(L, 'name', '平面三连杆');
robot.gravity = [0 -9.81 0];

% 检验DH坐标系建立是否正确
robot.teach();
save('planar3link_robot.mat', 'robot');