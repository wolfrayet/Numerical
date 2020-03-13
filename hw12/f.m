% check differential equation
syms y1(t) y2(t) y3(t)
ode1 = diff(y1,t) == -y1/(0.01*0.01)-y2/0.01+1/(0.01*0.01);
ode2 = diff(y2,t) == y1-y2-y3;
ode3 = diff(y3,t) == y2;
odes = [ode1,ode2,ode3];
cond1 = y1(0) == 0;
cond2 = y2(0) == 0;
cond3 = y3(0) == 0;
conds = [cond1,cond2,cond3];
[y1Sol,y2Sol,y3Sol] = dsolve(odes,conds);
