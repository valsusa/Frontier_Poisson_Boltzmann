function res = PB_FH_residual(xkrylov)
global L; global dx; global NG;
global DT;

global N;
global WP; global QM; global Q; 
global rho_back e_over_kT lde;
global x0; global v0;
global E0 dphi phi0;
global xc xv t k gamma 
global Debye

% calculate the x at n+1/2 time level
x_average = x0 + xkrylov(1:N)*DT/2;
out=(x_average<0); x_average(out)=x_average(out)+L;
out=(x_average>=L);x_average(out)=x_average(out)-L;

% calculate the interpolation functions
p=1:N;p=[p p];
g1=floor(x_average/dx)+1;     
g=[g1;g1+1];
fraz1=1-abs(x_average(1:N)/dx-g1+1); 	
out=(g<1);g(out)=g(out) + NG;
out=(g>NG);g(out)=g(out)- NG;

% calculate the J
fraz=[(fraz1).*xkrylov(1:N);(1-fraz1).*xkrylov(1:N)];	
mat=sparse(p,g,fraz,N,NG);
J = full((Q/dx)*sum(mat))';
J(NG+1)=J(1);

divJ=zeros(NG,1);
divJ(1:NG)=(J(2:NG+1)-J(1:NG))/dx;

%manufactured solution

%divJ=-(k^2+1/lde^2)*sin(k*xc)*gamma*exp(-gamma*t);

res = zeros(N,1);


un=ones(NG,1);

    if(Debye)
        olde=1.0/lde^2*dx^2;
    else    
        olde=exp(e_over_kT*(phi0(1:NG)))/lde^2*dx^2;
    end 

%olde=1.0/lde^2*dx^2;
Poisson=spdiags([un -2*un-olde un],[-1 0 1],NG,NG);
Poisson(1,NG)=1;
Poisson(NG,1)=1;

dphi=Poisson\(divJ(1:NG)*dx^2*DT);%dphi=[dphi;0];

dE=zeros(NG+1,1);
dE(2:NG) = -(dphi(2:NG)-dphi(1:NG-1))/dx;
dE(1) = -(dphi(1)-dphi(NG))/dx;
dE(NG+1) = dE(1);



% residual for the average velocity
fraz=[(fraz1);1-fraz1];	
mat=sparse(p,g,fraz,N,NG);

%Boltzmann electrons
res(1:N) = xkrylov(1:N) - v0 - 0.25*mat*QM*(2*E0(1:NG)+dE(1:NG))*DT;

%inifinitely hot electrons electrons
%res(1:N) = xkrylov(1:N) - v0 - 0.25*mat*QM*(E0 +E1)*DT;



 
   
   
  


