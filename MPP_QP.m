function[x,res,C_ineq,yalmiptime,solvertime]=MPP_QP(theta,A,H,E,b,d,B,f,F,C)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function MPP_QP solves the following QP:
%
%               min     1/2*x*H*x' + (C*theta+d)*x
%               s.to    A*x <= E*theta+b
%                       B*x = F*theta+f.
%
%This function returns the minimizer (x), problem with the QP (if any)
%returned by YALMIP (res), the inequality constraint satisfaction (C_ineq),
%and the times reported by YALMIP (yalmiptime) and the chosen solver (solvertime).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=sdpvar(size(A,2),1);

Obj = (1/2)*x'*H*x + d'*x + (C*theta)'*x;

C1 = A*x - E*theta-b <= 0;
C2 = B*x == F*theta+f;

constraints=[C1;C2];

ops = sdpsettings('solver','sedumi','debug',0,'verbose',0,'sedumi.eps',1e-12,'sedumi.numtol',1e-12);

sol = optimize(constraints,Obj,ops);
res = sol.problem;

yalmiptime = sol.yalmiptime;
solvertime = sol.solvertime;

x=value(x);
C_ineq=(- A*x + E*theta + b);

end

