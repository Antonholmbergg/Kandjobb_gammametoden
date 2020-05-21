% Use one step with the Trapets method, second order accurate
% Works for all constant coefficient cases the cases

uex1 = uex(x,y,k);
udex1 = udex(x,y,k);

B = [I2, -k/2*I2; -k/(2*h^2)*c^2*A0, I2];
B0= [I2, +k/2*I2; +k/(2*h^2).*c^2*A0, I2];
b0= k/2*[zeros(N,1);(1-2*gamma)*frhs(x,y,k)+gamma*frhs(x,y,0)+gamma*frhs(x,y,2*k)];
d0 =[u0;w0];
rhs0=B0*d0+b0;
[B,rhs]=Dirichlet(B,rhs0,[uex1;udex1],[bndr;bndr+N]);
d  = B\rhs;
u1 = d(1:N,1);
