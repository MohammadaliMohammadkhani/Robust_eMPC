clc
clear all
close all
warning off


%=============== This is a simulation for implimenting the eMPC control algorithm with
%disturbance and state observer ====================================================== 

%===========================  SIMULATION PARAMETERS =========================

TS=.05;   % Sampling time (sec)
tend=15; % Duration of simulation (sec)
tdist=0.5*tend; % time of disturbance activation
distAmp=0.5;

yref=0.5;

%=========================== System Configuration ==========================
% this part introduces system in state space format dx=A*x+B*u+T*theta and  y=C*x+D*u+Tprime*theta
% we define a structure format for system matrices 
% A=sysStruct.A,B=sysStruct.B,C=sysStruct.C,T=sysStruct.T,Tprime=sysStruct.Tprime
% sysStruct.xmax, sysStruct.xmin, sysStruct.umax, sysStruct.umin defines
% bounds for states and input signal
sysStruct.A=[1 TS;0 1];%[0.7326 -0.0861;0.1722 0.5909]; 
sysStruct.B=[TS^2/2;TS];%[0.1025;0.2]; 
sysStruct.C=[1 0];
sysStruct.D=0;
sysStruct.T1=[1;1];
sysStruct.T2=[0];
sysStruct.Adist=1;


% this part defines states nd input bounds
sysStruct.xmax = [5;5]; % states maximum bound
sysStruct.xmin = [-5;-5]; % states minimum bound
sysStruct.ymax = []; % output maximum bound
sysStruct.ymin = []; % output minimum bound
sysStruct.umax = [20]; % input maximum bound
sysStruct.umin = [-20]; % input minimum bound

% ========================= MPC parameters ==============================
% MPC solves an optimal problem which solves sigma(xT Q x +uT R u)
% we have to define Q, R , Nu, Nc and N
probStruct.N = 2; % prediction horizon
probStruct.Nc = 2; % Constraints Horizon
probStruct.Nu = 2; % Input Horizon
probStruct.Q=100*eye(2);   % weights on states
probStruct.P=100*eye(1);   % weights on states
probStruct.R=0.1*eye(1);   % weights on inputs
probStruct.norm=2; % kind of norms are 1,2 and inf

% ctrl=mpt_control(sysStruct,probStruct,'offline')

% ========================= Matrices for state partitioning =============
% this part generates matrices for partitioning
tStart = tic;Matrices2 = My_eMPC_Matrices_x_t(sysStruct,probStruct);t_mpc1=toc(tStart)
tStart = tic;Matrices3 = My_eMPC_Matrices_x(sysStruct,probStruct);t_mpc2=toc(tStart)
%============================= partitioning ==============================
tStart = tic;[Pn,Fi,Gi,activeConstraints,Phard,details]=mpt_mpqp(Matrices2);t_obsever=toc(tStart)
tStart = tic;[Pn0,Fi0,Gi0,activeConstraints0,Phard0,details0]=mpt_mpqp(Matrices3);t_obsever=toc(tStart)

A=sysStruct.A;B=sysStruct.B;C=sysStruct.C;D=sysStruct.D;T=sysStruct.T1;Tprime=sysStruct.T2;Adist=sysStruct.Adist;
%=================target calculator
astinv=inv([A-eye(2) B;C D]);

% ======================== observer design ==============================
state_dec=2; % 1= exact states, 2= estimated states
Aa=[0.7326 -0.0861 1;0.1722 0.9909 1;0 0 1]; % augmented matrices A
Ba=[sysStruct.B;0]; % augmented matrices B
Ca=[1 0 0]; % augmented matrices C
% Plant = ss(Aa,Ba,Ca,0,-1);
% [kalmf,L,P,M]  = kalman(Plant,1,1)

% ================ using genetic algorithm to minimization of the observer
% gains ==================================================================
% A=sysStruct.A;B=sysStruct.B;C=sysStruct.C;D=sysStruct.D;T=sysStruct.T1;Tprime=sysStruct.T2;Adist=sysStruct.Adist;
% tStart = tic;[X,FVAL,OUTPUT] =  Genetic_for_observer_matrices(A,C,T,Tprime,Adist,Matrices2,TS)
% t_genetic=toc(tStart)
% Ld=X(1);
% Lx=[X(2);X(3)];
% 
% tStart = tic;[X1,FVAL,OUTPUT] =  Genetic_for_observer_matrices1(A,C,T,Tprime,Adist,Matrices2,TS)
% t_genetic=toc(tStart)
% Ld1=X1(1);
% Lx1=[X1(2);X1(3)];
load Ld
load Ld1
load Lx;load Lx1
%======================== observer dynamics ==============================
Aobs=[A+Lx*C T;Ld*C 1];
observer_eig=eig(Aobs) %eigen valuse of the observer


% ======================= Runing Simulation ==============================
theta0=0; % initial values of disturbance
x0 = [-0.5;0.6]; % initial values of the states
x0_aug = [-0.5;0.6;0]; % initial values of the states
sim('Simulation_empc_with_observer');




