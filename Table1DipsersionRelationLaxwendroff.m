
clear;
clc;
close all


M=7;

kh=linspace(0.01,0.59*pi,M);
% M=9;
% kh=linspace(0.01,0.65*pi,M);
AA=zeros(M,M);
b=zeros(M,1);

for ii=1:M %kx=1
    
    for kk=1:M
        AA(ii,kk)=2*(sin((kk-0.5)*kh(ii)));
    end
  
    b(ii)=(kh(ii));
end

c=AA\b;%求系数
length(c)
temp1=-2*sum(c);
digits(6)
vpa(c)'


kh=linspace((pi)/(100),(pi),100);

% a=0;
% for m=1:length(c)
%     a=a+2*c(m)*(sin(m*kh));
% end
% 
% figure;plot(a-kh,'b','LineWidth',2);grid on
% axis([0 100 -2.5*10^-4 2.5*10^-4])


a=0;
for m=1:length(c)
    a=a+2*c(m)*(sin((m-0.5)*kh));
end

figure;plot(a.^2-kh.^2,'b','LineWidth',2);grid on
% axis([0 100 -2.5*10^-4 2.5*10^-4])

%%%%%%%%%%%%%%%%%%%%%%%%%
M=8;

kh=linspace(0.01,0.62*pi,M);
% M=7;
% kh=linspace(0.01,0.57*pi,M);
AA=zeros(M,M);
b=zeros(M,1);

for ii=1:M %kx=1
    
    for kk=1:M
        AA(ii,kk)=2*(sin((kk-0.5)*kh(ii)))*2*(sin(0.5*kh(ii)));
    end
  
    b(ii)=(kh(ii))^2;
end

c=AA\b;%求系数
length(c)
% temp1=-2*sum(c);
digits(6)
vpa(c)'


kh=linspace((pi)/(100),(pi),100);

a=0;
for m=1:length(c)
    a=a+2*c(m)*(sin((m-.5)*kh))*2.*(sin(0.5*kh));
end

hold on;plot(a-kh.^2,'k','LineWidth',2);grid on
axis([0 100 -4.5*10^-4 1.5*10^-4])
legend('The previous SGFD method','The proposed SGFD method')
xlabel('The percentage of kh')
ylabel('Error')
