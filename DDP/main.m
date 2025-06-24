clear;
clc;
import casadi.*

%% 主函数

h = 0.05;  % RK4离散化
Tfinal = 5.0;  % 终止时间
Nt = fix(Tfinal/h)+1;  % 时间步数，即1~N

%% Cost weights

% stage cost l=1/2*(x-x_goal)'*Q(x-x_goal)+1/2*u*R*u
Q = diag([0,0,0,0.00001,0.00001,0.00001,10,10,10]);
R = diag([1,1,1]);
% terminal cost lN=1/2*(x_N-x_goal)'*QN(x_N-x_goal)
QN = diag([0,0,0,70,70,70,200,200,200]);

%% Initial guess

x0 = [0.01; 0.02; 0.03;0;0;0;3;0;0];
xgoal = [0;0;0;0;0;0;0; 2.5; pi/2];
Nx = length(x0);  % 状态变量个数
xtraj = kron(ones(1,Nt), x0);

Nu = 3;  % 输入控制个数
u0 = [0.01;0.02;0.03];
utraj = kron(u0,ones(1,Nt-1));

%%  DDP
load('planar3link_robot.mat');  % 会恢复变量 robot
%robot.teach();  % 继续使用

%if isempty(gcp('nocreate'))
    %parpool(4);  % 只在没有并行池时才开启
%end

dynamics_jiexi=dynamics(robot,Nu);
dynamics_rk4_step = rk4_step( h,dynamics_jiexi,Nu);
dfdx = dfdx(Nx,Nu,dynamics_rk4_step);
dfdu = dfdu(Nx,Nu,dynamics_rk4_step);
dAdx = dAdx(Nx,Nu,dfdx);
dAdu = dAdu(Nx,Nu,dfdx);
dBdx = dBdx(Nx,Nu,dfdu);
dBdu = dBdu(Nx,Nu,dfdu);

%[xtraj,utraj,J,Jhist] = iLQR(Nx, Nu, Nt, xtraj,h, utraj, xgoal, Q, QN,R,robot,dynamics_rk4_step,dfdx,dfdu,dAdx,dAdu,dBdx,dBdu);
[xtraj,utraj,J,Jhist] = DDP(Nx, Nu, Nt, xtraj,h, utraj, xgoal, Q, QN,R,robot,dynamics_rk4_step,dfdx,dfdu,dAdx,dAdu,dBdx,dBdu);
