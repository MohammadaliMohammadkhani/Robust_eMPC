function [c, ceq] = Constraint_for_observer_matrices(p)
global A C T Tprime Adist
Ld=[p(1)];
Lx=[p(2);p(3)];
Ae=[A+Lx*C T+Lx*Tprime;Ld*C Adist+Ld*Tprime];
c=[abs(eig(Ae))-ones(length(Ae),1)];ceq=[];
