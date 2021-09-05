
% My own M-code for generating Multi-Parametric-Programming matrices:
% Dynamics:       x(k+1)=Ax(k)+Bu(k)
%                   y(k)=Cx(k)+Du(k)
% Constraints:    ymin  <=    y(k)     <= ymax
%                 umin  <=    u(k)     <= umax

% % matrices of the problem:  Matrices = struct[G,W,E,H,F,Y] 
% % 
% %        J=min_U  (0.5 U' H U + x(0)'FU + x(0)Yx(0))
% %                 G U <= W + E x(0)

function Matrices = My_eMPC_Matrices_x(sysStruct,probStruct)

A=sysStruct.A;
B=sysStruct.B;
C=sysStruct.C;
D=sysStruct.D;

xmin=sysStruct.xmin;
xmax=sysStruct.xmax;
umin=sysStruct.umin;
umax=sysStruct.umax;

Q=probStruct.Q;
R=probStruct.R;
N=probStruct.N;

[n, m]=size(B);


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
Xmin = []; Xmax = []; 
Umin = []; Umax = []; 
Rhat = [];
Qhat = [];
Ahat = [];
Bhat = [];
Chat = [];
Bt = zeros(n,N*m);
P = PLQR;
for i = 1:N
    Rhat = blkdiag(Rhat, R);
    if i < N
        Qhat = blkdiag(Qhat, Q);
    else
        Qhat = blkdiag(Qhat, P);
    end
    
    Ahat = [Ahat; A^i];
    
    Bt = [A^(i-1)*B  Bt(:,1:end-m)];
    Bhat = [Bhat; Bt];
    
    Chat = blkdiag(Chat, C);

    Xmin = [Xmin; xmin];  Xmax = [Xmax; xmax];
    Umin = [Umin; umin];  Umax = [Umax; umax];
end


H = (Rhat + Bhat'*Qhat*Bhat);     H = 0.5*(H+H');
F = Ahat'*Qhat*Bhat;
Y = Q + Ahat'*Qhat*Ahat;     Y = 0.5*(Y+Y');

BhatN = Bt;
%% -------------------------------------------------------------------------

%% -------------------------------------------------------------------------
% [Gf, hf] = double(PinvSet1);
G = [Bhat; -Bhat; eye(N*m); -eye(N*m)];
W = [Xmax; -Xmin; Umax; -Umin];
E = [-Ahat; Ahat; zeros(N*m,n); zeros(N*m,n)];

% % Take out limits which are Inf
% aux=find(isinf(W));
% G(aux,:)=[];
% W(aux,:)=[];
% E(aux,:)=[];
% 
% aux=find(all(([G E]==zeros(size([G E])))')'); % Rows which are all 0
% G(aux,:)=[];
% W(aux,:)=[];
% E(aux,:)=[];
% 
% GEW=polytope([G -E],W);
% [GE,W]=double(GEW);
% G=GE(:,1:(end-n));
% E=-GE(:,(end-n+1):end);

%% Important Note:  the obtained matrices are not equal to the MPT ones, because some rows may exchanged. However, the result is the same.       
%% -------------------------------------------------------------------------
% adapting to the MPT:
% tt= G(6,:); tt2= G(7,:); tt3= G(8,:);  tt4= G(9,:);  G(6:9,:)=[];
% G=[ G(1,:); tt;  G(2,:); tt2;  G(3,:); tt3;  G(4,:); tt4;  G(5:end,:)];
% tt= W(6,:); tt2= W(7,:); tt3= W(8,:);  tt4= W(9,:);  W(6:9,:)=[];
% W=[ W(1,:); tt;  W(2,:); tt2;  W(3,:); tt3;  W(4,:); tt4;  W(5:end,:)];
% tt= E(6,:); tt2= E(7,:); tt3= E(8,:);  tt4= E(9,:);  E(6:9,:)=[];
% E=[ E(1,:); tt;  E(2,:); tt2;  E(3,:); tt3;  E(4,:); tt4;  E(5:end,:)];

% ------- Generate output:       
Matrices.G = G;
Matrices.W = W;
Matrices.E = E;

Matrices.H = H;
Matrices.F = F;
Matrices.Y = Y;


