% N=20; h=1/(N+1);
% G=numgrid('S',N+2);
%
% % SzyldFrommer
% eta=inline('x.^2.*cos(x+y).^2','x','y'); alpha=inline('(x+y).^2.*exp(x-y)','x','y'); beta=inline('(x+y).^2.*exp(x-y)','x','y'); nu=inline('1.5 + 0*x','x','y'); mu=inline('-0.5 + 0*x','x','y');
% MyFuns_SzyldFrommer = {eta,alpha,beta,nu,mu};
% A_SzyldFrommer = ReactAdvDiff_Sqr_FD_test('SzyldFrommer',G,MyFuns_SzyldFrommer);
% % SzyldGander
% eta=inline('x.^2.*cos(x+y).^2','x','y'); a=inline('(x+y).^2.*exp(x-y)','x','y'); b1=inline('1.5 + 0*x','x','y'); b2=inline('-0.5 + 0*x','x','y');
% MyFuns_SzyldGander = {eta,a,b1,b2};
% A_SzyldGander = ReactAdvDiff_Sqr_FD_test('SzyldGander',G,MyFuns_SzyldGander);
% % look where the two differ - we adapted 'SzyldGander2012' so we assume it's correct and check against it for the case of constant advection (so that both PDEs are the same).
% err = A_SzyldFrommer - A_SzyldGander; spy(err)
% 
% 
% % SzyldFrommer - simple example for checking non-constant Advection being
% % correct
% eta=inline('x.^2.*cos(x+y).^2','x','y'); alpha=inline('(x+y).^2.*exp(x-y)','x','y'); beta=inline('(x+y).^2.*exp(x-y)','x','y'); nu=inline('1*x','x','y'); mu=inline('1*y','x','y');
% MyFuns_SzyldFrommer = {eta,alpha,beta,nu,mu};
% A_SzyldFrommer = ReactAdvDiff_Sqr_FD_test('SzyldFrommer',G,MyFuns_SzyldFrommer);

% %%% do non-symmetric AdvecDiff based SzyldFrommer eqn (15) and p.660
% G=numgrid('S',N+2);
% eta = @(x,y) 0.*x+0.*y; alpha = @(x,y) 1 + 0.*x+0.*y; beta = @(x,y) 1 + 0.*x+0.*y; nu = @(x,y) 150.*(1.*x.*(x-1).*(1-2.*y)); mu = @(x,y) 150.*(-1.*y.*(y-1).*(1-2.*x));
% A=ReactAdvDiff_Sqr_FD_test('SzyldFrommer',G,{eta,alpha,beta,nu,mu});
% rhs = ones(size(A,1),1); u_init = zeros(size(A,1),1); u_exact = A\rhs;
% u_plot = zeros(N+2); x_mesh = 0:h:1; angle1 = 224; angle2 = 35;
% % plot the approximate solution
% nmb_int_gridcols = length(x_mesh)-2; u_plot(2:nmb_int_gridcols+1,2:nmb_int_gridcols+1)=reshape(u_exact,nmb_int_gridcols,nmb_int_gridcols)';
% mesh(x_mesh,x_mesh,u_plot); xlabel('x');ylabel('y'); view(angle1,angle2); title('Exact Solution'); pause(2);




function SysMtrx = ReactAdvDiff_Sqr_FD(PDElayout, G, MyFuns)
% ReactAdvDiff_Sqr_FD assembles the centered finite difference discretization of a reaction-advection-diffusion equation on a rectangle
%   in the sparse format, namely for the operator (eta - div([alpha,beta].grad) + AdvectionTerm) on the grid G
%   The "AdvectionTerm" has two options: (i) 'SzyldGander2012' -> b_vect.grad;  (ii) 'SzyldFrommer2023' -> d/dx(nu*) + d/dy(mu*);
%   The grid G can be generated by NUMGRID (e.g., numgrid('S',N+2) )
%   Original file due to C. Moler, 91-07-16, adapted by Gander&Szyld, 08-05-13; Current adaptation by Outrata 23-12-19
if strcmp(PDElayout,'SzyldGander')
    SysMtrx = delsqnonsymmetric_SzyldGander(MyFuns,G);
elseif strcmp(PDElayout,'SzyldGander_SclrFeval')
    SysMtrx = delsqnonsymmetric_SzyldGander_SclrFeval(MyFuns,G);
elseif strcmp(PDElayout,'SzyldFrommer')
    SysMtrx = delsqnonsymmetric_SzyldFrommer(MyFuns,G);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = delsqnonsymmetric_SzyldFrommer(FunsPackeged,G)
% delsqnonsymmetric_SzyldFrommer finite difference reaction-advection-diffusion, see equation 15 in 
% Szyld, Frommer: On the convergence of randomized and greedy relaxation schemes for solving nonsingular linear systems of equations (Numerical Algoprtihms, Springer)

eta = FunsPackeged{1}; alpha = FunsPackeged{2}; beta = FunsPackeged{3}; nu = FunsPackeged{4}; mu = FunsPackeged{5};

[m,n] = size(G); h=1/(n-1); [X,Y]=meshgrid(0:h:1,0:h:(m-1)*h);
E=feval(eta,X,Y);
Am=feval(alpha,X-h/2,Y); Ap=feval(alpha,X+h/2,Y);
Bm=feval(beta,X,Y-h/2); Bp=feval(beta,X,Y+h/2);
Np=feval(nu,X+h,Y); Nm=feval(nu,X-h,Y); 
Mp=feval(mu,X,Y+h); Mm=feval(mu,X,Y-h);

% Indices of interior points
p = find(G);

% The diagonal contributions -> contribution from reaction "eta" and the analogue of the "+4" on the diagonal of the Laplacian.
i = G(p); j = G(p);
s = E(p)*h^2 + Am(p)+Ap(p) + Bm(p)+Bp(p) ;

% The off-diagonal contributions -> contribution from advection "nu,mu" and the analogue of the "-1" on the off-diagonals of the Laplacian.
% The off-diagonals for any fixed gridpoint (i,j) come from "north, east, south, west". 
% In the 1D ordering of the gridpoints we have (i,j) corresponding to some k and then: north ~ k-1; west ~ k-m; south ~ k+1; east ~ k+m;

% Possible neighbors in the north direction
ell=-1; Q = G(p+ell);
% Index of points with interior neighbors in the north direction
q = find(Q);
% Connect interior points to neighbors
i = [i; G(p(q))]; j = [j; Q(q)];
s = [s; -Bm(p(q)) - h/2*Mm(p(q))];

% Possible neighbors in the south direction
ell=1; Q = G(p+ell);
q = find(Q);
i = [i; G(p(q))]; j = [j; Q(q)];
s = [s; -Bp(p(q)) + h/2*Mp(p(q))];

% Possible neighbors in the east direction
ell=m; Q = G(p+ell);
q = find(Q);
i = [i; G(p(q))]; j = [j; Q(q)];
s = [s; -Ap(p(q)) + h/2*Np(p(q))];

% Possible neighbors in the west direction
ell=-m; Q = G(p+ell);
q = find(Q);
i = [i; G(p(q))]; j = [j; Q(q)];
s = [s; -Am(p(q)) - h/2*Nm(p(q))];

D = sparse(i,j,s)/h^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = delsqnonsymmetric_SzyldGander(FunsPackeged,G)
% delsqnonsymmetric_SzyldGander finite difference reaction-advection-diffusion, see Figure 2.1
% Gander, Loisel, Szyld: AN OPTIMAL BLOCK ITERATIVE METHOD AND PRECONDITIONER FOR BANDED MATRICES WITH APPLICATIONS TO PDES ON IRREGULAR DOMAINS (Journal on Mtrx Analysis and Appl, SIAM)
eta = FunsPackeged{1}; a = FunsPackeged{2}; b1 = FunsPackeged{3}; b2 = FunsPackeged{4};

[m,n] = size(G);
h=1/(n-1);
[X,Y]=meshgrid(0:h:1,0:h:(m-1)*h);
E=feval(eta,X,Y);
Axm=feval(a,X-h/2,Y); Axp=feval(a,X+h/2,Y);
Aym=feval(a,X,Y-h/2); Ayp=feval(a,X,Y+h/2);
B1=feval(b1,X,Y);
B2=feval(b2,X,Y);
% quiver(X,Y,B1,B2)

% Indices of interior points

p = find(G);

% Connect interior points to themselves with 4's.
i = G(p);
j = G(p);
s = (E(p)*h^2+Axm(p)+Axp(p)+Aym(p)+Ayp(p));

% for k = north, east, south, west
k=-1;
% Possible neighbors in k-th direction
Q = G(p+k);
% Index of points with interior neighbors
q = find(Q);
% Connect interior points to neighbors
i = [i; G(p(q))];
j = [j; Q(q)];
s = [s; -Aym(p(q))-h/2*B2(p(q))];
k=1;
Q = G(p+k);
q = find(Q);
i = [i; G(p(q))];
j = [j; Q(q)];
s = [s; -Ayp(p(q))+h/2*B2(p(q))];
k=m;
Q = G(p+k);
q = find(Q);
i = [i; G(p(q))];
j = [j; Q(q)];
s = [s; -Axp(p(q))+h/2*B1(p(q))];
k=-m;
Q = G(p+k);
q = find(Q);
i = [i; G(p(q))];
j = [j; Q(q)];
s = [s; -Axm(p(q))-h/2*B1(p(q))];

D = sparse(i,j,s)/h^2;
end











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = delsqnonsymmetric_SzyldGander_SclrFeval(FunsPackeged,G)
% delsqnonsymmetric_SzyldGander finite difference reaction-advection-diffusion, see Figure 2.1
% Gander, Loisel, Szyld: AN OPTIMAL BLOCK ITERATIVE METHOD AND PRECONDITIONER FOR BANDED MATRICES WITH APPLICATIONS TO PDES ON IRREGULAR DOMAINS (Journal on Mtrx Analysis and Appl, SIAM)
eta = FunsPackeged{1}; a = FunsPackeged{2}; b1 = FunsPackeged{3}; b2 = FunsPackeged{4};

[m,n] = size(G); h=1/(n-1);
x_mesh = 0:h:1; y_mesh = 0:h:(m-1)*h;
for x_ind = 1:length(x_mesh)
    for y_ind = 1:length(y_mesh)
        E(x_ind,y_ind)=feval(eta,x_mesh(x_ind),y_mesh(y_ind));
        Axm(x_ind,y_ind)=feval(a,x_mesh(x_ind)-h/2,y_mesh(y_ind)); Axp(x_ind,y_ind)=feval(a,x_mesh(x_ind)+h/2,y_mesh(y_ind));
        Aym(x_ind,y_ind)=feval(a,x_mesh(x_ind),y_mesh(y_ind)-h/2); Ayp(x_ind,y_ind)=feval(a,x_mesh(x_ind),y_mesh(y_ind)+h/2);
        B1(x_ind,y_ind)=feval(b1,x_mesh(x_ind),y_mesh(y_ind));
        B2(x_ind,y_ind)=feval(b2,x_mesh(x_ind),y_mesh(y_ind));
    end
end
E = E'; Axm = Axm'; Axp = Axp'; Aym = Aym'; Ayp = Ayp'; B1= B1'; B2 = B2';

% %%% check
% [X,Y]=meshgrid(0:h:1,0:h:(m-1)*h);
% E_check=feval(eta,X,Y);
% Axm_check=feval(a,X-h/2,Y); Axp_check=feval(a,X+h/2,Y);
% Aym_check=feval(a,X,Y-h/2); Ayp_check=feval(a,X,Y+h/2);
% B1_check=feval(b1,X,Y);
% B2_check=feval(b2,X,Y);
% disp(num2str(norm(E-E_check,'fro')))
% disp(num2str(norm(Axm-Axm_check,'fro')))
% disp(num2str(norm(Axp-Axp_check,'fro')))
% disp(num2str(norm(Aym-Aym_check,'fro')))
% disp(num2str(norm(Ayp-Ayp_check,'fro')))
% disp(num2str(norm(B1-B1_check,'fro')))
% disp(num2str(norm(B2-B2_check,'fro')))
% quiver(X,Y,B1,B2)

% Indices of interior points

p = find(G);

% Connect interior points to themselves with 4's.
i = G(p);
j = G(p);
s = (E(p)*h^2+Axm(p)+Axp(p)+Aym(p)+Ayp(p));

% for k = north, east, south, west
k=-1;
% Possible neighbors in k-th direction
Q = G(p+k);
% Index of points with interior neighbors
q = find(Q);
% Connect interior points to neighbors
i = [i; G(p(q))];
j = [j; Q(q)];
s = [s; -Aym(p(q))-h/2*B2(p(q))];
k=1;
Q = G(p+k);
q = find(Q);
i = [i; G(p(q))];
j = [j; Q(q)];
s = [s; -Ayp(p(q))+h/2*B2(p(q))];
k=m;
Q = G(p+k);
q = find(Q);
i = [i; G(p(q))];
j = [j; Q(q)];
s = [s; -Axp(p(q))+h/2*B1(p(q))];
k=-m;
Q = G(p+k);
q = find(Q);
i = [i; G(p(q))];
j = [j; Q(q)];
s = [s; -Axm(p(q))-h/2*B1(p(q))];

D = sparse(i,j,s)/h^2;
end
