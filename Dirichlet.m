% Imposes Dirichlet b.c. according to a list of nodes on the boundary bndr
% handles inhomogeneous boundary conditiona
% sol_ex must contain the b.c. in the boundary nodes
% but other etries can be zero
function [K,rhs] = Dirichlet(K,rhs,sol_ex,bndr)

n = size(K,1);
r = 1:n;
nbndr = setdiff(r,bndr);

if size(sol_ex,2) == 3 % for the gamma method, requires the tree time levels
    r1 = sol_ex(:,1); r1(nbndr)=0;
    ww1 = K*r1;

    r2 = sol_ex(:,2); r2(nbndr)=0;
    ww2 = K*r2;

    r3 = sol_ex(:,3); r3(nbndr)=0;
    ww3 = K*r3;
else %for compute first step, only uses one timestep
    r3 = sol_ex; r3(nbndr)=0;
    ww1 = K*r3;
    ww2 = 0;
    ww3 = 0;
end

e = ones(n,1); e(nbndr)=0;
D = spdiags(e,0,n,n);
K(bndr,:) = 0;
K(:,bndr) = 0;
K = K + D;

%computes the new right hand side
rhs = rhs - ww3 + 2*ww2 - ww1;
rhs(bndr) = r3(bndr);

return
