function res = PB_JFNK_residuals(xkrylov)
global L; global dx; global NG;
global DT;

global N;
global WP; global QM; global Q; 
global rho_back e_over_kT lde;
global x0; global v0;
global E0 E1 ni Eel phi0 divJ Je;
global Nonlinear



% calculate the x at n+1/2 time level
x_average = x0 + xkrylov(1:N)*DT/2;
out=(x_average<0); x_average(out)=x_average(out)+L;
out=(x_average>=L);x_average(out)=x_average(out)-L;

% calculate the interpolation functions
p=1:N;p=[p p];
g1=floor(x_average/dx)+1;     
g=[g1;g1+1];
fraz1=1-abs(x_average(1:N)/dx-g1+1); 	
%out=(g<1);g(out)=g(out) + NG;
%out=(g>NG);g(out)=g(out)- NG;
g=mod(g-1,NG)+1;
%[max(g) min(g)]
% calculate the J
fraz=[(fraz1).*xkrylov(1:N);(1-fraz1).*xkrylov(1:N)];	
mat=sparse(p,g,fraz,N,NG);
J = full((Q/dx)*sum(mat))';
J(NG+1)=J(1);

divJ=zeros(NG,1);
divJ(1:NG)=(J(2:NG+1)-J(1:NG))/dx;

res = zeros(N+NG,1);

dphi= xkrylov(N+1:N+NG);
dE=zeros(NG+1,1);
dE(2:NG) = -(dphi(2:NG)-dphi(1:NG-1))/dx;
dE(1) = -(dphi(1)-dphi(NG))/dx;
dE(NG+1) = dE(1);


% Residual botlzmann
res(N+2:N+NG-1) = ( dphi(1:NG-2)+dphi(3:NG)-2*dphi(2:NG-1) )/dx^2 ;
res(N+1) = ( dphi(NG)+dphi(2)-2*dphi(1) )/dx^2 ;
res(N+NG) = ( dphi(NG-1)+dphi(1)-2*dphi(NG) )/dx^2 ;
%res(N+1:N+NG) = res(N+1:N+NG) + rho_back.* (ni  - exp(e_over_kT*dphi));

if(Nonlinear) 
    res(N+1:N+NG) = res(N+1:N+NG) - divJ*DT - exp(e_over_kT*(phi0+ dphi/2)) .* dphi /lde^2;
else
    res(N+1:N+NG) = res(N+1:N+NG) - divJ*DT - exp(e_over_kT*(phi0+ 0* dphi/2)) .* dphi /lde^2;
end    
%res(N+1:N+NG) = res(N+1:N+NG) - divJ*DT - exp(e_over_kT*(phi0+ dphi/2)) .* dphi /lde^2;
%res(N+1:N+NG) = res(N+1:N+NG) - divJ*DT - dphi./lde^2;
%res(N+NG) = dphi(NG);


% residual for the average velocity
fraz=[(fraz1);1-fraz1];	
mat=sparse(p,g,fraz,N,NG);

%Boltzmann electrons
res(1:N) = xkrylov(1:N) - v0 - 0.25*mat*QM*(2*E0(1:NG)+dE(1:NG))*DT;

%inifinitely hot electrons electrons
%res(1:N) = xkrylov(1:N) - v0 - 0.25*mat*QM*(E0 +E1)*DT;



 
   
   
  


