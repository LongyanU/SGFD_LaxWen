clear;
 close all;
 clc



load('figure5bb.mat')

plotimage((1:nt-2)*3.5,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')

seis_recordb=seis_record;






load('figure5c.mat')

% load('figure5c3M8.mat')

plotimage((1:nt-2)*3.5,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')

seis_recordc=seis_record;






figure;

hold on;plot((1:nt-2)*3.5,seis_recordb((1:nt-2),350)./max(abs(seis_recordb((1:nt-2),350))),'r','linewidth',2)
hold on;plot((1:nt-2)*3.5,seis_recordc((1:nt-2),350)./max(abs(seis_recordc((1:nt-2),350))),'k','linewidth',2)
hold on;plot((1:nt-2)*3.5,seis_recordc((1:nt-2),350)./max(abs(seis_recordc((1:nt-2),350)))-seis_recordb((1:nt-2),350)./max(abs(seis_recordb((1:nt-2),350))),'b','linewidth',2)


xlabel('time(ms)')
ylabel('Amp')
legend('traditional FD scheme for L-W','proposed FD scheme for L-W')
grid on
axis([ 0 5150 -0.65 1.01])
box on
grid on

% load('figure5bb_2.mat')
% 
% % load('figure5c3M8.mat')
% 
% % plotimage((1:nt-2)*3.5,seis_record(1:nt-2,50:end-50))
% % xlabel('x/dx')
% % ylabel('Time(ms)')
% % title('')
% 
% seis_recordd=seis_record;
% hold on;plot((3:2:nt-2)*3.5/2,seis_recordd((3:2:nt-2),350)./max(abs(seis_recordd((3:2:nt-2),350))),'b','linewidth',2)


% hold on;plot(seis_recorde(:,350)./max(abs(seis_recorde(:,350))),'g','linewidth',2)

% figure;
% % plot(seis_recorda(:,350)./max(abs(seis_recorda(:,350))),'b','linewidth',2);
% 
% % hold on;
% plot(1000*(seis_recordb(:,350)./max(abs(seis_recordb(:,350)))-seis_recordc(:,350)./max(abs(seis_recordc(:,350)))),'r','linewidth',2)
% % hold on;plot(seis_recordc(:,350)./max(abs(seis_recordc(:,350)))-(temp./max(abs(temp))),'k','linewidth',2)



% hold on;plot(temp-1*10^-5,'b','linewidth',1);

