
clear;
clc;
close all


M=7;

kh=linspace(0.0001,0.6*pi,M+1);
% M=9;
% kh=linspace(0.01,0.65*pi,M);
AA=zeros(M+1,M+1);
b=zeros(M+1,1);

for ii=1:M+1 %kx=1
    AA(ii,1)=1;
    for kk=1:M
        AA(ii,kk+1)=2*cos(kk*kh(ii));
    end
  
    b(ii)=-(kh(ii)).^2;
end

c=AA\b;%ÇóÏµÊý
length(c)
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


a=c(1);
for m=1:length(c)-1
    a=a+2*c(m+1)*(cos((m)*kh));
end

figure;plot(-a-kh.^2,'b','LineWidth',2);grid on
axis([0 100 -2.5*10^-4 2.5*10^-4])

a=c(1);
for m=1:length(c)-1
    a=a+2*c(m+1)*(cos((m)*kh));
end

c=[ 1.24223, -0.113062, 0.0269687, -0.00716431, 0.00167719, -0.000288665, 0.0000265238];

a1=0;
for m=1:length(c)
    a1=a1+2*c(m)*(sin((m-0.5)*kh));
end

figure;plot(-a.*a1.*a1-kh.^4,'b','LineWidth',2);grid on
axis([0 100 -2.5*10^-4 5.5*10^-4])
