% Just some code for seeing propagation in time
k=2^(-6);    % timestep
h=k/2;  % space discretization step
%H(l) = h;
n = 1/h;  % number of intervals per coordinate direction

nd = n+1; % number of node points in each direction
N = nd^2; % matrix size
xc = 0:h:1;
yc = 0:h:1;


gamma = 0.4;

T = 2;    % time interval [0,T]
nk = ceil(T/k); % number of timesteps

% create the space discretization matrix A-discrete negative Laplacian

[xg,yg]=meshgrid(xc,yc);

%q = @(x,y) x+0.2;
q = @(x,y) ((x > 0.5)*3+1);



variable_coeff_matrix


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



% frhs = @(x,y,t) (x.^2-x).*(y.^2-y).*12.*t.^2 - (4.*x-1.4).*(y.^2-y).*(t.^4+1) - (x.^2-x).*(x+0.4).*2.*y.*(t.^4+1);
% uex = @(x,y,t) (x.^2-x).*(y.^2-y).*(t.^4+1);
% udex = @(x,y,t) (x.^2-x).*(y.^2-y).*4.*t^3;

% frhs = @(x,y,t) 2.*(x.^2-x).*(y.^2-y) - (4.*x-0.6).*(y.^2-y).*(t.^2+1) - (x.^3-x.^2-0.2.*x).*2.*(t.^2+1);
% uex = @(x,y,t) (x.^2-x).*(y.^2-y).*(t.^2+1);
% udex = @(x,y,t) (x.^2-x).*(y.^2-y).*2.*t;

uex = @(x,y,t) 0.*x;
udex = @(x,y,t) 0.*x;
frhs = @(x,y,t)  0.*x;
u0f = @(x,y) sin(pi*x).*sin(pi*y);
u0 = reshape(u0f(xg,yg),N,1);


I2 = speye(N,N);
w0 = udex(x,y,0);   % derivative at t=0
rho = (k/h)^2;
c = 1;
compute_first_step_trapets  % compute u1



err_abs = zeros(nk,1);
err_rel = zeros(nk,1);
err_abs(1,1) = norm(uex1-u1);
err_rel(1,1) = norm(uex1-u1)/norm(uex1);
cnt = 1;

u2 = zeros(N,1); % just allocation
A = I2 - gamma*rho*A0;

restart = 1;
maxit = n;
tol=1e-6;
for tcur = k:k:T-k
    uex_cur=uex(x,y,tcur+k);
    b = rho*A0*u1 + k^2*(gamma*frhs(x,y,tcur+k) + (1 - 2*gamma)*frhs(x,y,tcur) + gamma*frhs(x,y,tcur-k));
    [Acur,bcur]=Dirichlet(A,b,uex_cur,bndr);
    xi = Acur\bcur; 
    %[xi_it,fl,relres,iter,resvec] = agmg(Acur,bcur,restart,tol,maxit,1);
    
    u1(bndr) = 0; % to get the boundary rigth
    u0(bndr) = 0;
    
    u2 = xi + 2*u1 - u0;  % solution at t_cur+k
    u0 = u1;
    u1 = u2;
    
    
    pause(0.01)
    surf(reshape(u2,nd,nd), 'EdgeColor', 'none')
    %surf(reshape(uex_cur,nd,nd), 'EdgeColor', 'none')
    set(gca,'Zlim', [-1, 2])
    drawnow  
end