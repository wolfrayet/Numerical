syms y(t)
Dy = diff(y,t,1);
DDy = diff(y,t,2);
ode = diff(y,t,3)+101*diff(y,t,2)+111*diff(y,t,1)+100*y-100 ==0;
cond1 = y(0)==0;
cond2 = Dy(0)==0;
cond3 = DDy(0)==0;
conds = [cond1,cond2,cond3];
ySol = dsolve(ode,conds);
fplot(ySol);