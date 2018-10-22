clc; clear;

pkg load symbolic;
syms p t s
Rp=[[1,0,0];[0,cos(p),sin(p)];[0,-sin(p),cos(p)]];
Rt=[[cos(t),0,-sin(t)];[0,1,0];[sin(t),0,cos(t)]];
Rs=[[cos(s),sin(s),0];[-sin(s),cos(s),0];[0,0,1]];
Tbg=Rp*Rt*Rs

Tgb=Rs'*Rt'*Rp'
