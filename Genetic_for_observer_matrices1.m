function [X2,FVAL,OUTPUT] =  Genetic_for_observer_matrices1(A1,C1,T1,Tprime1,Adist1,Matrices,Ts1)
global Matrices_Optimization A T Tprime C Adist Ts;
Matrices_Optimization=Matrices;
A=A1;
Tprime=Tprime1;
C=C1;
Adist=Adist1;
T=T1;
Ts=Ts1;

Aeq=[];
Beq=[];
Aineq=[];
Bineq=[];


% alpha=[A;A^2];beta=[B zeros(2,2);A*B B];gamma=[T;A*T+T];
% Y1=alpha'*Qhat*alpha;Y2=gamma'*Qhat*gamma;Y3=alpha'*Qhat*gamma;
% H=beta'*Qhat*beta+Rhat;
% F1=alpha'*Qhat*beta;F2=gamma'*Qhat*beta;
% F=[F1' F2']';


nvars = 3;
LB = [-50 -50 -50];
UB = [50 50 50];
%Nonlinear constraints
nonlconFunction = [];

%Start with default options
options = gaoptimset;
options = gaoptimset(options,'TolFun',1e-10,'TolCon',1e-10);
%     options = gaoptimset(options,'InitialPopulation',[0.47,.57]);

%%Modify some parameters
options = gaoptimset(options,'Display', 'diagnose');

options = gaoptimset(options,'PopInitRange',[LB;UB],'StallGenLimit',15,...
    'StallTimeLimit',400000,'Generations',50,'PopulationSize',50,'InitialPenalty',10,'PenaltyFactor',100);
nonlconFunction=@Constraint_for_observer_matrices;
[X2,FVAL,OUTPUT] = ga(@Optimization_for_observer_matrices,nvars,Aineq,Bineq,Aeq,Beq,LB,UB,nonlconFunction,options)