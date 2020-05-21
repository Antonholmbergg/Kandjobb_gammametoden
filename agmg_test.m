% Tests the scaling of agmg
T_agmg = [];

for l = 1:9
k=2^(-l-1);
h=k;
n = 1/h;
nd = n+1;
N = nd^2;
restart = 1;
maxit = n;
tol=1e-6;
I2 = speye(N,N);
gamma = 0.5;

c = 1;
e = ones(N,1);
A0 = spdiags([e,e,-4*e,e,e],[-nd,-1,0,1,nd],N,N);
rho = (c^2*k/h)^2;


for j = 1:nd-1
    A0(j*nd,j*nd+1) = 0;
    A0(j*nd+1,j*nd) = 0;    
end

A = I2 - gamma*rho*A0;

% create rhs
bcur = rand(N,1);

%constructing the solver
agmg(A,[],[],[],[],0,[],1);

%meassure the time to solve
tic;
[xi_it,fl,relres,iter,resvec] = agmg(A,bcur);
tagmg = toc;
T_agmg(l) = tagmg

end

