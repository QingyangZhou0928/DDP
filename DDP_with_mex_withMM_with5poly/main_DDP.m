clear;
clc;
import casadi.*
load('planar3link_robot.mat');  % 会恢复变量 robot
load('xtraj.mat');
load('utraj.mat');
utraj = utraj(1:3,1:100);

%% 主函数

h = 0.05;  % RK4离散化
Tfinal = 5.0;  % 终止时间
Nt = fix(Tfinal/h)+1;  % 时间步数，即1~N

%% Cost weights

% stage cost l=1/2*(x-x_goal)'*Q(x-x_goal)+1/2*u*R*u
Q = diag([0,0,0,0.00001,0.00001,0.00001,10,10,10,1]);
R = diag([1,1,1]);
% terminal cost lN=1/2*(x_N-x_goal)'*QN(x_N-x_goal)
QN = diag([0,0,0,70,70,70,200,200,100,5]);

%% Initial guess
Nu = 3;  % 输入控制个数


x0 = xtraj(:,1);
xgoal = xtraj(:,Nt);

Nx = length(x0);  % 状态变量个数

%%  DDP

if isempty(gcp('nocreate'))
    parpool(4);  % 只在没有并行池时才开启
end


[xtraj,utraj,J,Jhist] = DDP(Nx, Nu, Nt, xtraj,utraj, xgoal, Q, QN,R);
