function j = Optimization_for_observer_matrices(p)
global Matrices_Optimization A T Tprime C Adist Ts;
Ld=[p(1)];
Lx=[p(2);p(3)];
Fbar=[Matrices_Optimization.F1;Matrices_Optimization.F2];
Qprime=Fbar*inv(Matrices_Optimization.H)*Fbar';
SQ=sqrtm(Qprime);
Ae=[A+Lx*C T+Lx*Tprime;Ld*C Adist+Ld*Tprime];

z = tf('z',Ts);
sizeAe=length(Ae);
TFe=SQ*inv(z*eye(sizeAe)-Ae);
j=norm(TFe,2);

