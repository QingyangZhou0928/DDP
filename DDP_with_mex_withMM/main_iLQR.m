clear;
clc;
import casadi.*
load('planar3link_robot.mat');  % 会恢复变量 robot

%% 主函数

h = 0.05;  % RK4离散化
Tfinal = 5.0;  % 终止时间
Nt = fix(Tfinal/h)+1;  % 时间步数，即1~N

%% Cost weights

% stage cost l=1/2*(x-x_goal)'*Q(x-x_goal)+1/2*u*R*u
Q = diag([0,0,0,0.00001,0.00001,0.00001,0.1,0.1,0.1,1]);
R = diag([1,1,1]);
% terminal cost lN=1/2*(x_N-x_goal)'*QN(x_N-x_goal)
QN = diag([0,0,0,5,5,5,20,20,20,10]);

%% Initial guess
Nu = 3;  % 输入控制个数
u0 = [0.01;0.02;0.03];
utraj = kron(u0,ones(1,Nt-1));

x0 = [0; 0; 0;0;0;0;3;0;0];
xgoal = [0;0;0;0;0;0;0; 2.5; pi/2;10];

%% calculate omega
for i = 1:Nu
    link = robot.links(i);
    DH(i, :) = [link.theta, link.d, link.a, link.alpha];
end  % DH矩阵
DH =MX(DH);

x =MX.sym('x',Nu,1);
Jacob = Jacob_geo(x,DH);
Jacobgeo = Function('Jacobgeo',{x},{Jacob});

J=full(Jacobgeo(x0(1:3)));
J = J([1,2,6],:);
x0=[x0;sqrt(det(J*J'))];

Nx = length(x0);  % 状态变量个数

%%  DDP

%if isempty(gcp('nocreate'))
    %parpool(4);  % 只在没有并行池时才开启
%end


%[xtraj,utraj,J,Jhist] = iLQR(Nx, Nu, Nt, xtraj,h, utraj, xgoal, Q, QN,R,robot,dynamics_rk4_step,dfdx,dfdu,dAdx,dAdu,dBdx,dBdu);
[xtraj,utraj,J,Jhist] = DDP(Nx, Nu, Nt, xtraj,h, utraj, xgoal, Q, QN,R,robot,dynamics_rk4_step,dfdx,dfdu,dAdx,dAdu,dBdx,dBdu);
