% Use one Euler backwards step, only first order accurate
% Works for the polynomial solutions in the constant case and variable case

uex1 = uex(x,y,k);

f1 = frhs(x,y,k);
Acur = I2- (k/h)^2*c^2*A0; 

bcur = u0+k*w0+k^2*f1;
[Acur,bcur] = Dirichlet(Acur,bcur,uex1,bndr);
u1   = Acur\bcur;