function [time_elapsed, average_energy_error, average_iterations, max_iterations] = PB_JFNK_main(resolution)

%%%%%%%%%%%%%%%%%%%%%%%
% NK PIC
%%%%%%%%%%%%%%%%%%%%%%%


close all

global L; global dx; global NG;
global DT;

global N;
global WP; global QM; global Q; global rho_back e_over_kT lde;
global x0; global v0;
global E0 E1 ni Eel phi0 divJ Je;
global Nonlinear

graphics = false;
Nonlinear = false;

% parameters
L=2*pi*10*resolution;
DT=.125*10;%/resolution;
NT=200;%*resolution;
NTOUT=25;
NG=64*resolution;
N=10000*resolution;
WP=1;
QM=1;
V0=.2;
VT=0.01;
% perturbation 
XP1=0.0*N/L; 
V1=0.1;
mode=2*resolution;
% charge and gris parameter
Q=WP^2/(QM*N/L);
rho_back=Q*N/L;
dx=L/NG;

lde=2*pi*10/1000;
e_over_kT=1/lde^2/WP^2;

E2f=[];

% 2 Stream instability
x0=linspace(0,L-L/N,N)';
v0=VT*randn(N,1); % maxwellian


pm=[1:N]';
pm=1-2*mod(pm,2);
v0=v0+pm.*V0;

% Perturbation
v0=v0+V1*sin(2*pi*x0/L*mode);
x0=x0+XP1*(L/N)*sin(2*pi*x0/L*mode);
out=(x0<0); x0(out)=x0(out)+L;
out=(x0>=L); x0(out)=x0(out)-L;

E0 = zeros(NG+1,1);
phi0=zeros(NG,1);
rhoi=ones(NG,1)*rho_back;

% grid

xv=linspace(0,L,NG+1);xv=reshape(xv,NG+1,1);
xc=linspace(dx/2,L-dx/2,NG);xc=reshape(xc,NG,1);

% absolute and relative tollerance
tol = [1E-14, 1E-14]/10;

% prepare xkrylov vector
histEnergy = [];
histEnergyP = [];
histEnergyK = [];
histSolverIteration = [];
spettro=[];
test=[];
Ee=0;
t=0;

time_start = clock();
for it=1:NT
   it;
   xkrylov = [v0; phi0]; 
   [sol, it_hist, ierr] = nsolgm(xkrylov,'PB_JFNK_residual',tol,[100 100 .9]);
   % res=residueEC_test(sol)
   v_average = sol(1:N);
   dphi=sol(N+1:N+NG);
   phi0=phi0+dphi;
   rhoi=rhoi-DT*divJ;
   Em=E0;
   E0(2:NG) = E0(2:NG)-(dphi(2:NG)-dphi(1:NG-1))/dx;
   E0(1) = E0(1) - (dphi(1)-dphi(NG))/dx;
   E0(NG+1) = E0(1);
  x_average = x0 + v_average*DT/2;   
  out=(x_average<0); x_average(out)=x_average(out)+L;
out=(x_average>=L);x_average(out)=x_average(out)-L;
   v0 = 2*v_average - v0;
   x0 = x0 + v_average*DT;
   out=(x0<0); x0(out)=x0(out)+L;
   out=(x0>=L);x0(out)=x0(out)-L;

   %E0 = E1;
   test=[test [sum(E0);  it_hist(end,2)]];

   %E2f=[E2f ; E0'];
   Ek = 0.5*abs(Q)*sum(v0.^2);
  % Ep = 0.5*sum(diff(phi0).^2);
   %Ep = -0.5*sum(phi0.*(rhoi-rho_back*exp(e_over_kT*phi0)))*dx/lde^2+0.5*sum(E0(1:NG).^2)*dx;
   
   if(Nonlinear) 
     Ep = 0.5*sum(rho_back*exp(e_over_kT*(phi0-dphi/2)).*phi0.*phi0)*dx./lde^2+0.5*sum(E0(1:NG).^2)*dx;
     divJe=-exp(e_over_kT*(phi0- dphi/2)) .* dphi /lde^2/DT;
   else
     Ep = 0.5*sum(rho_back*exp(e_over_kT*(phi0-2*dphi/2)).*phi0.*phi0)*dx./lde^2+0.5*sum(E0(1:NG).^2)*dx;
     divJe=-exp(e_over_kT*(phi0- 2*dphi/2)) .* dphi /lde^2/DT;
   end 
  

   
%      divJe=-exp(e_over_kT*(phi0- 0* dphi/2)) .* dphi /lde^2/DT;
      
%         divJe=- dphi /lde^2/DT;
un=ones(NG-1,1);
Poisson=spdiags([un -2*un un],[-1 0 1],NG-1,NG-1);
psi=Poisson\(divJe(1:NG-1)*dx^2);psi=[psi;0];
Je=zeros(NG+1,1);
Je(2:NG) = -(psi(2:NG)-psi(1:NG-1))/dx;
Je(1) = -(psi(1)-psi(NG))/dx;
Je(NG+1) = Je(1);

   Ee=Ee+sum(0.5*(E0(1:NG)+Em(1:NG)).*Je(1:NG))*dx*DT;
  Ep=Ee+0.5*sum(E0(1:NG).^2)*dx;
   
   %Ep = 0.5*sum(E0(1:NG).^2+phi0.^2./lde^2)*dx;
  
   CFL=max(abs(v_average))*DT/dx;
  
   Etot = Ek + Ep;
   histEnergy = [histEnergy Etot];
   histEnergyP = [histEnergyP Ep];
   histEnergyK = [histEnergyK Ek];
   histSolverIteration = [histSolverIteration it_hist(end,2)];
   spettro= [spettro phi0];
if(graphics & mod(it,round(NT/NTOUT))==0)
subplot(2,3,1:2)
plot(x0,v0,'.','markersize',[1])
xlabel('x')
ylabel('v')
axis tight
title(['\omega_{pe}t = ' num2str(CFL) '     CFL = ' num2str(CFL) '      \Delta x/\lambda_{De} = ' num2str(dx/VT)])
subplot(2,3,3)
plot(xc,phi0)
xlabel('x')
ylabel('\phi')
axis tight
subplot(2,3,4)
plot(1:it,histEnergy,1:it,histEnergyP,1:it,histEnergyK)
xlabel('Cycle')
legend('Total','Field','Part','location','southeast')
subplot(2,3,5)
plot((histEnergy-histEnergy(1))/histEnergy(1))
xlabel('Cycle')
axis tight
title('\delta E')
subplot(2,3,6)
plot(test(2,:))
xlabel('Cycle')
axis tight
title('Number iterations')
pause(.1)  
end

t=t+DT;
end

% timing

time_end = clock();
resolution
time_elapsed = etime(time_end,time_start)

average_energy_error=mean((histEnergy-histEnergy(1))/histEnergy(1))
average_iterations=mean(test(2,:))
max_iterations=max(test(2,:))


figure
imagesc([0 t], [0 L], spettro)
xlabel('x/d_e','fontsize',15)
xlabel('\omega_p t','fontsize',15)
set(gca, 'fontsize',15)
colorbar
print('-depsc','vanillaJFNK_spacetime.eps')

figure
pcolor(abs(fft2(spettro)))
shading interp

if(~graphics & mod(it,round(NT/NTOUT))==0)
subplot(2,3,1:2)
plot(x0,v0,'.','markersize',[1])
xlabel('x')
ylabel('v')
axis tight
title(['\omega_{pe}t = ' num2str(CFL) '     CFL = ' num2str(dx/VT) '      \Delta x/\lambda_{De} = ' num2str(dx/VT)])
subplot(2,3,3)
plot(xc,phi0)
xlabel('x')
ylabel('\phi')
axis tight
subplot(2,3,4)
plot(1:it,histEnergy,1:it,histEnergyP,1:it,histEnergyK)
xlabel('Cycle')
legend('Total','Field','Part','location','southeast')
subplot(2,3,5)
plot((histEnergy-histEnergy(1))/histEnergy(1))
xlabel('Cycle')
axis tight
title('\delta E')
subplot(2,3,6)
plot(test(2,:))
xlabel('Cycle')
axis tight
title('Number iterations')
print('-depsc','vanillaJFNK.eps')
end


save ECpic.txt histEnergy;
save PhaseSpaceEC.txt x0 v0;
save SolverIterations.txt histSolverIteration;
save TimeEC.txt time_elapsed;

end
