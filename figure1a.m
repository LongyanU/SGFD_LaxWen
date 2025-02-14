

clear;clc;
close all
M=8
f= [ 2.50048, -0.812038, 0.280435, -0.0938613, 0.0278821, -0.00670525, 0.00114084, -0.000101696];
d=[ 1.56764, -0.283078, 0.0822817, -0.0259518, 0.00760084, -0.00184479, 0.000321364, -0.0000296118];


h=10;
v=2500;
v=5000;
v=8500;
tau=0.001
r=v*tau/h;
kh=linspace((pi)/(100),(pi),100);

for jj=1:5
    xita=(jj-1)*pi/16;
    temp=0;
    for m=1:M
        temp=temp-2*d(m)*sin((m-0.5)*kh*sin(xita))*2.*sin(0.5*kh*sin(xita));
    end
    
    temp1=0;
    for m=1:M
        temp1=temp1-2*d(m)*sin((m-0.5)*kh*cos(xita))*2.*sin(0.5*kh*cos(xita));
    end
    
    
    temp2=0;
    for m=1:M
        temp2=temp2-2*f(m)*sin((m-0.5)*kh*sin(xita))*2.*sin(0.5*kh*sin(xita)).*( 2*cos(kh*sin(xita))-2+ 2*cos(kh*cos(xita))-2     );
    end
    
    temp3=0;
    for m=1:M
        temp3=temp3-2*f(m)*sin((m-0.5)*kh*cos(xita))*2.*sin(0.5*kh*cos(xita)) .*( 2*cos(kh*sin(xita))-2+ 2*cos(kh*cos(xita))-2     );
    end
    
    atemp=r^2*(temp1+temp)+r^4/12*(temp3+temp2);
    
%     a1=(acos(1+1/2*atemp)./(kh*r));
    
    a1=h/v*(1./(acos(1+1/2*atemp)./(kh*r))-1);
%      atemp=r^2*(temp1+temp);
%     if jj==1
%         figure;plot(atemp)
%     else
%         hold on;plot(atemp);
%     end
%     
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
grid on
xlabel('percentage of kh')
legend('\theta=0', '\theta=¦Ð/16','\theta=2¦Ð/16','\theta=3¦Ð/16','\theta=4¦Ð/16')

ylabel('\delta_1_ (\theta)')