
% My own M-code for generating Multi-Parametric-Programming matrices:
% Dynamics:       x(k+1)=Ax(k)+Bu(k)
%                   y(k)=Cx(k)+Du(k)
% Constraints:    ymin  <=    y(k)     <= ymax
%                 umin  <=    u(k)     <= umax

% % matrices of the problem:  Matrices = struct[G,W,E,H,F,Y]
% %
% %        J=min_U  (0.5 U' H U + x(0)'FU + x(0)Yx(0))
% %                 G U <= W + E x(0)

function Matrices = My_eMPC_Matrices_x_t(sysStruct,probStruct)

A=sysStruct.A;
B=sysStruct.B;
C=sysStruct.C;
D=sysStruct.D;
T1=sysStruct.T1;
T2=sysStruct.T2;
Adist=sysStruct.Adist;

ymin=sysStruct.ymin;
ymax=sysStruct.ymax;
xmax=sysStruct.xmax;
xmin=sysStruct.xmin;
% xt=sysStruct.xt;
% ut=sysStruct.ut;
umin=sysStruct.umin;
umax=sysStruct.umax;

Q=probStruct.Q;
R=probStruct.R;
N=probStruct.N;

[n, m]=size(B);
[n1, m1]=size(T1);


%% -------------------------------------------------------------------------
% ***** Calculation of infinite invariant set Xf:
[KLQR,PLQR] = mpt_dlqr(A, B, Q, R);
KLQR=-KLQR;  % since u(t+k) = + KLQ * x(t+k)

% total_A = [KLQR; -KLQR; eye(2); -eye(2)];
% total_B = [umax; -umin; xmax; -xmin];

% X1 = polytope(total_A, total_B);
A_cl = A+ B*KLQR;              %closed loop dynamics

% [PinvSet1,tstar,fd] = mpt_infset(A_cl,X1, 100);
%% -------------------------------------------------------------------------
% load OPTION
% load Pinv
%% -------------------------------------------------------------------------
Ymin = []; Ymax = [];
Xmin = []; Xmax = [];
Umin = []; Umax = [];
Rhat = [];
Qhat = [];
Ahat = [];
Bhat = [];
Chat = [];
Ghat = [];
T1hat=[];
T2hat=[];
Xt=[];
Ut=[];
Adisthat=[];
Bt = zeros(n,N*m);
Gt = zeros(n1,m1);
P = PLQR;
for i = 1:N
%     Xt=[xt;Xt];
%     Ut=[ut;Ut];
    Rhat = blkdiag(Rhat, R);
    if i < N
        Qhat = blkdiag(Qhat, Q);
    else
        Qhat = blkdiag(Qhat, P);
    end
    Ahat = [Ahat; A^i];
    Adisthat=[Adisthat;Adist^i];
    Bt = [A^(i-1)*B  Bt(:,1:end-m)];
    Bhat = [Bhat; Bt];
    Gt = [A^(i-1)*T1+Gt];
    Ghat = [Ghat; Gt];
    Chat = blkdiag(Chat, C);
    T1hat = blkdiag(T1hat, T1);
    T2hat = blkdiag(T2hat, T2);
    Xmin = [Xmin; xmin];  Xmax = [Xmax; xmax];
    Ymin = [Ymin; ymin];  Ymax = [Ymax; ymax];
    Umin = [Umin; umin];  Umax = [Umax; umax];
end


H = (Rhat + Bhat'*Qhat*Bhat);     H = 0.5*(H+H');
F3=-Qhat*Bhat;
F1= Ahat'*Qhat*Bhat;
F2=Ghat'*Qhat*Bhat;
F4=-Rhat;

F=[F1;F2;F3;F4];
Y1 = Q + Ahat'*Qhat*Ahat;     Y1 = 0.5*(Y1+Y1');
Y2=Ghat'*Qhat*Ghat;  Y2=0.5*(Y2+Y2');
Y3=Qhat;Y3=0.5*(Y3+Y3');
Y4=Rhat;
Y=blkdiag(Y1,Y2,Y3,Y4);
BhatN = Bt;

%% -------------------------------------------------------------------------

if (length(Ymin))>=1
    G = [Chat*Bhat; -Chat*Bhat; eye(N*m); -eye(N*m)];
    W = [Ymax; -Ymin; Umax; -Umin];
E = [-Chat*Ahat -Chat*Ghat-T2hat*Adisthat zeros(N*m,(n+m)*(N)); Chat*Ahat Chat*Ghat+T2hat*Adisthat zeros(N*m,(n+m)*(N)); zeros(N*m,n+m1+(n+m)*N); zeros(N*m,n+m1+(n+m)*N)];
else
    G = [-Bhat;Bhat;eye(N*m); -eye(N*m)];
    W = [-Xmin;Xmax;Umax; -Umin];
E = [Ahat Ghat zeros(N*n,(n+m)*N);-Ahat -Ghat zeros(N*n,(n+m)*N);zeros(N*m,n+m1+(n+m)*N); zeros(N*m,n+m1+(n+m)*N)];
end

% ------- Generate output:
Matrices.G = G;
Matrices.W = W;
Matrices.E = E;
Matrices.Xt=Xt; 
Matrices.Ut=Ut; 
Matrices.H = H;
Matrices.F = F;
Matrices.F1 = F1;
Matrices.F2 = F2;
Matrices.Y = Y;

