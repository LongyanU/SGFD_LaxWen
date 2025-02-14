clear;
 close all;
 clc

load('figure3a.mat')


plotimage((1:nt-2)*2,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')

seis_recorda=seis_record;


load('figure3bb.mat')

plotimage((1:nt-2)*2,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')
seis_recordb=seis_record;






load('figure3c.mat')

% load('figure5c3M8.mat')

plotimage((1:nt-2)*2,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')

seis_recordc=seis_record;


% load('figure5c3M8_2.mat')
% plotimage(seis_record)
% seis_recorde=seis_record;


load('figure3d.mat')
plotimage(    seis_record)
NT=nt;
B = ifft(exp(-2i*sin([0:NT-1]*pi/(2*NT))'*[0:NT-1]).*cos([0:NT-1]'*pi/(2*NT)),2*NT,'symmetric');
I = B(1:NT,1:NT)'; % <- The Inverse Time Dispersion Transform matrix

temp=I*seis_record(:,75);




figure;plot((1:nt-2)*2,seis_recorda((1:nt-2),75)./max(abs(seis_recorda((1:nt-2),75))),'b','linewidth',2);

hold on;plot((1:nt-2)*2,seis_recordb((1:nt-2),75)./max(abs(seis_recordb((1:nt-2),75))),'r','linewidth',2)
hold on;plot((1:nt-2)*2,seis_recordc((1:nt-2),75)./max(abs(seis_recordc((1:nt-2),75))),'k','linewidth',2)

hold on;plot((1:nt-2)*2,temp(1:nt-2)./max(abs(temp)),'m','linewidth',2)

% hold on;plot(seis_recorde(:,75)./max(abs(seis_recorde(:,75))),'g','linewidth',2)

% figure;
% % plot(seis_recorda(:,75)./max(abs(seis_recorda(:,75))),'b','linewidth',2);
% 
% % hold on;
% plot(1000*(seis_recordb(:,75)./max(abs(seis_recordb(:,75)))-seis_recordc(:,75)./max(abs(seis_recordc(:,75)))),'r','linewidth',2)
% % hold on;plot(seis_recordc(:,75)./max(abs(seis_recordc(:,75)))-(temp./max(abs(temp))),'k','linewidth',2)



% hold on;plot(temp-1*10^-5,'b','linewidth',1);

xlabel('time(ms)')
ylabel('Amp')
legend('FD scheme without L-W','traditional FD scheme for L-W','proposed FD scheme for L-W','Time dispersion elimination')
grid on
% axis([ 2190 2250 -8*10^-3 6*10^-3])
