

clear;clc;close all
M=7

d=[ 1.24275, -0.113493, 0.0272555, -0.00731429, 0.00173544, -0.000303651, 0.0000284348];



f=  [ 1.83185, -0.350241, 0.097898, -0.0276368, 0.00670127, -0.00119173, 0.000114199];
f1= -3.11498;

h=10;
v=2500;
v=5000;
v=8500;
tau=0.001
dt=tau;
r=v*tau/h;
kh=linspace((pi)/(100),(pi),100);

for jj=1:5
    xita=(jj-1)*pi/16;
    d1kz=0;
    for m=1:M
        d1kz=d1kz+2*d(m)*sin((m-0.5)*kh*sin(xita));
    end
%     d1kz=d1kz./h;
    
    d1kx=0;
    for m=1:M
        d1kx=d1kx+2*d(m)*sin((m-0.5)*kh*cos(xita));
    end
%     d1kx=d1kx./h;
    
    d2kx=f1;
    for m=1:M
        d2kx=d2kx+2*f(m)*(cos(m*kh*cos(xita))    );
    end
%     d2kx=d2kx./h^2;
    
    d2kz=f1;
    for m=1:M
        d2kz=d2kz+2*f(m)*(cos(m*kh*sin(xita))   );
    end
%      d2kz=d2kz./h^2;
    %%%%%%%%%%%%%%%%%%%
    temp1=d2kx+d2kz;
    
    g11=-d1kx.*(d1kx+r^2*1/24*d1kx.*temp1);
    g12=-d1kz.*(d1kz+r^2*1/24*d1kz.*temp1);
%     g1=v^2*dt^2*(g11+g12);
    
    g1=r^2*(g11+g12);
    
%     g2=g1.*temp1;
    atemp=(g1)+r^2/24*(g1.*temp1);
    
    
    %%%%%%%%%%%%%%%
    
    
    a1=real(h/v*(1./(acos(1+1/2*atemp)./(kh*r))-1));
    
  if (jj==1)
        figure;plot(100*kh/(pi),a1,'m','linewidth',2)
        hold on
    elseif jj==2
        plot(100*kh/(pi),a1,'r--','linewidth',2)
    elseif jj==3
        plot(100*kh/(pi),a1,'c:','linewidth',2)
    elseif jj==4
        plot(100*kh/(pi),a1,'k-.','linewidth',2)
    else
        plot(100*kh/(pi),a1,'b','linewidth',2)
    end
    
end

axis([0 100 -0.5*10^-5 7*10^-5])
% axis([0 100 -0.5*10^-5 3*10^-5])
grid on
xlabel('percentage of kh')
legend('\theta=0', '\theta=¦Ð/16','\theta=2¦Ð/16','\theta=3¦Ð/16','\theta=4¦Ð/16')

ylabel('\delta_2_ (\theta)')