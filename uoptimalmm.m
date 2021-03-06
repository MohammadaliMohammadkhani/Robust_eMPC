function U=uoptimalmm(x1)
eml.extrinsic('evalin');
eml.extrinsic('isinside');
Pn = evalin('base', 'Pn');
[inn1,k1,opt1]=isinside(Pn,x1);
Fi= evalin('base', 'Fi');
Gi= evalin('base', 'Gi');
U1 = Fi{k1}*x1 + Gi{k1};
U=U1(1,1);