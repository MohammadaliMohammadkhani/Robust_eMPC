%======================== Result =========================================
tend=10; % Duration of simulation (sec)
tdist=(1/3)*tend; % time of disturbance activation
%======================== First Simulation Run ===========================
state_dec=2; % 1= exact states, 2= estimated states 3=optimal states
sim('Simulation_empc_with_observer');
accr=1e-5;
% ======================= Squeezing signals for plot =====================
y.signals.values=squeeze(y.signals.values);
ref.signals.values=squeeze(ref.signals.values);
u.signals.values=squeeze(u.signals.values);
disthat.signals.values=squeeze(disthat.signals.values);

figure (1)
subplot(2,2,1)
stairs(ref.signals.values,'r','LineWidth',2);hold on;grid on;
stairs(y.signals.values,'--b','LineWidth',2)
subplot(2,2,2)
stairs(u.signals.values,'--b','LineWidth',2);hold on;grid on;
subplot(2,2,3)
stairs(dist.signals.values,'r','LineWidth',2);hold on;grid on;
stairs(disthat.signals.values,'--b','LineWidth',2)

%===numerical metrices
Ae=[A+Lx*C T+Lx*Tprime;Ld*C Adist+Ld*Tprime];
Fbar=[Matrices2.F1;Matrices2.F2];
Qprime=Fbar*inv(Matrices2.H)*Fbar';
SQ=sqrtm(Qprime);
z = tf('z',TS);
sizeAe=length(Ae);
TFe=SQ*inv(z*eye(sizeAe)-Ae);
optimal_norm5=norm(TFe,2)
sum_error.signals.values(end)
%======================== Second Simulation Run ==========================
state_dec=3; % 1= exact states, 2= estimated states 3=optimal states
sim('Simulation_empc_with_observer');

% ======================= Squeezing signals for plot =====================
y.signals.values=squeeze(y.signals.values);
u.signals.values=squeeze(u.signals.values);
disthat1.signals.values=squeeze(disthat1.signals.values);

figure (1)
subplot(2,2,1)
hold on;grid on;stairs(y.signals.values,':k','LineWidth',2)
xlabel('time instants');ylabel('y(t)')
legend('Ref.','MPC 1','MPC 2')
title('(a)')
subplot(2,2,2)
stairs(u.signals.values,':k','LineWidth',2);hold on;grid on;
xlabel('time instants');ylabel('u(t)')
legend('MPC 1','MPC 2')
title('(b)')
subplot(2,2,3)
hold on;grid on;stairs(disthat1.signals.values,':k','LineWidth',2)
legend('\theta(t)','Obsr 1','Obsr 2')
xlabel('time instants');ylabel('\theta(t)')
title('(c)')
subplot(2,2,4)
plot(Pn0)
title('(d)')

%===numerical metrices
Ae=[A+Lx1*C T+Lx1*Tprime;Ld1*C Adist+Ld1*Tprime];
sizeAe=length(Ae);
TFe=SQ*inv(z*eye(sizeAe)-Ae);
sum_error.signals.values(end)
optimal_norm5=norm(TFe,2)