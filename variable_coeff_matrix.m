%creates the matrices corresponding to the different diagonals
e_u2 = q(xg+h/2,yg);
e_u1 = q(xg,yg+h/2);
e_d =  q(xg,yg-h/2) + q(xg+h/2,yg) + q(xg-h/2,yg) + q(xg,yg+h/2);
e_l1 = q(xg,yg-h/2);
e_l2 = q(xg-h/2,yg);

% second upper diagonal
e_u2 = reshape(e_u2,N,1); %reshape to vector
e_u2(end-nd+1:end) = []; % remove elements 
e_u2 = [zeros(nd,1); e_u2]; % pad with zeros so that spdiags works

% first upper diagonal
e_u1 = reshape(e_u1,N,1);
e_u1(end) = [];
e_u1 = [0; e_u1];

% main diagonal
e_d = reshape(e_d,N,1);

% first lower diagonal
e_l1 = reshape(e_l1,N,1);
e_l1(1) = [];
e_l1 = [e_l1; 0];

% second lower diagonal
e_l2 = reshape(e_l2,N,1);
e_l2(1:nd) = [];
e_l2 = [e_l2; zeros(nd,1)];

%spdiags removes the correct elements from the shorter diagonals
A0 = spdiags([e_l2, e_l1, -e_d, e_u1, e_u2],[-nd,-1,0,1,nd],N,N);

% sets elements that are outeside the domain to zero
for j = 1:nd-1
    A0(j*nd,j*nd+1) = 0;
    A0(j*nd+1,j*nd) = 0;    
end

