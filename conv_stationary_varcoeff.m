% test the convergence of the stationary problem
k=2^(-5);    % timestep
h=k;  % space discretization step
%H(l) = h;
n = 1/h;  % number of intervals per coordinate direction

nd = n+1; % number of node points in each direction
N = nd^2; % matrix size

xc = 0:h:1;
yc = 0:h:1;

[xg,yg]=meshgrid(xc,yc);
x = reshape(xg,N,1);
y = reshape(yg,N,1);
% Find the boundary points
bnd1=find(x==0);
bnd2=find(y==0);
bnd3=find(x==1);
bnd4=find(y==1);
bnda=union(bnd1,bnd2,'stable');
bndb=union(bnd3,bnd4,'stable');
bndr=union(bnda,bndb,'stable');

% exact solution, scuorce term and function for the wave speed
% q = @(x,y) x+0.2;
% frhs_stat = @(x,y) (-9*x.^2 + 0.4 + 2.8*x).*(y.^3-y.^2) -(6*y-2).*(x+0.2).*(x.^3-x.^2);
% uex = @(x,y,t) (x.^3 - x.^2).*(y.^3 - y.^2);

% q = @(x,y) 2+0*x.*y;
% frhs_stat = @(x,y,t) -2*(12*x.^2).*(y.^4-y.^1) - 2*(12*y.^2).*(x.^4-x.^1);
% uex = @(x,y,t) (x.^4-x.^1).*(y.^4-y.^1);

% A way to compute the scource term with symbolic functions, it doesn't
% work perfectly though
%
% syms q(xs,ys) uex(xs,ys,ts) frhs_stat(xs,ys)
% q(xs,ys) = atan(xs).*atan(ys)*5;
% uex(xs,ys) = (xs.^4 - xs.^1).*(ys.^4 - ys.^1); 
% frhs_stat = matlabFunction(-divergence(q(xs,ys)*gradient(uex(xs,ys)),[xs ys]));
% q = matlabFunction(q(xs,ys));
% uex = matlabFunction(uex(xs,ys));

q = @(x,y) 0.4^2 + x.*0 + y.*0;
uex = @(x,y) sin(pi.*x).*sin(pi.*y);
frhs_stat = @(x,y) 0.4^2*2.*pi.^2.*sin(pi.*x).*sin(pi.*y);

%create the matrix
variable_coeff_matrix


rhs0 = -h^2*frhs_stat(x,y);

uex0 = uex(x,y);


restart = 1;
maxit = n;
tol=1e-6;
[A,rhs]=Dirichlet(A0,rhs0,uex0,bndr);
tic;
u = A\rhs;
toc
%tic
%[xi_it,fl,relres,iter,resvec] = agmg(A0,rhs,restart,tol,maxit,1);
%toc
%u(bndr) = uex0(bndr);
err = norm(uex0-u)/norm(uex0)

