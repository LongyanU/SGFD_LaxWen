
% 时间已过 118.110363 秒。
clear
clc %%%%%%%
close all
% Elapsed time is 9.253353 seconds.
nt=1500;    % number of time steps
eps=.6;     % stability
isnap=60;    % snapshot sampling


load('vv.mat')
c1=flipud(c);
c=c1;

[nz,nx]=size(c1);
v=c1;


nz=nz+50;

vv=zeros(nz,nx);
for ii=1:nz-50
    for jj=1:nx
        vv(ii+50,jj)=v(ii,jj);
    end
end

for ii=1:50  %%top
    for jj=1:nx
        vv(ii,jj)=vv(51,jj);
    end
end



clear v
v=vv;



p=zeros([nz nx]); pboundarynew=p;pdan=p;
dx=20;  %calculate space increment
h=dx;
x=(0:(nx-1))*dx;
z=(0:(nz-1))*dx;  % initialize space coordinates
dt=0.0038; % calculate time step from stability criterion
dt=0.0038; % calculate time step from stability criterion

tau=dt;
r2=v.*v.*dt*dt/dx/dx;

f0=38;
t=(1:nt)*dt;
t0=4/f0;                       % initialize time axis
src=10^6*exp(-f0^2*(t-t0).*(t-t0));              % source time function
src=-diff((src))/dx^2;				% time derivative to obtain


seis_record=zeros(nt,nx);
p=zeros([nz nx]); pnew=p; pold=p;
d2px=p;
d2pz=p;

% Source location
zs=51;
xs=nx/2;

h=dx;
r=v*dt/h;


% M=8;
%
% kh=linspace(0.01,0.7*pi,M);
%
% kh=linspace(0.01,0.57*pi,M);
coeff1= [ 2.50048, -0.812038, 0.280435, -0.0938613, 0.0278821, -0.00670525, 0.00114084, -0.000101696];

for loop=1:1
    
    tic
    
    p=zeros([nz nx]); Vx=p; Vz=p;
    
    coeff=[ 1.56764, -0.283078, 0.0822817, -0.0259518, 0.00760084, -0.00184479, 0.000321364, -0.0000296118];
    
    
    for it=1:nt-2,
        
        d2px11=Vx-circshift(Vx,[0 1]);
        d2px12=(circshift(Vx,[0 -1])-circshift(Vx,[0 2]));
        d2px13=(circshift(Vx,[0 -2])-circshift(Vx,[0 3]));
        d2px14=(circshift(Vx,[0 -3])-circshift(Vx,[0 4]));
        d2px15=(circshift(Vx,[0 -4])-circshift(Vx,[0 5]));
        d2px16=(circshift(Vx,[0 -5])-circshift(Vx,[0 6]));
        d2px17=(circshift(Vx,[0 -6])-circshift(Vx,[0 7]));
        d2px18=(circshift(Vx,[0 -7])-circshift(Vx,[0 8]));
        
        
        
        d2pz11=Vz-circshift(Vz,[1 0]);
        d2pz12=(circshift(Vz,[-1 0])-circshift(Vz,[2 0]));
        d2pz13=(circshift(Vz,[-2 0])-circshift(Vz,[3 0]));
        d2pz14=(circshift(Vz,[-3 0])-circshift(Vz,[4 0]));
        d2pz15=(circshift(Vz,[-4 0])-circshift(Vz,[5 0]));
        d2pz16=(circshift(Vz,[-5 0])-circshift(Vz,[6 0]));
        d2pz17=(circshift(Vz,[-6 0])-circshift(Vz,[7 0]));
        d2pz18=(circshift(Vz,[-7 0])-circshift(Vz,[8 0]));
        
        
        d2px=coeff(1)*d2px11+coeff(2)*d2px12+coeff(3)*d2px13+coeff(4)*d2px14+coeff(5)*d2px15+coeff(6)*d2px16...
            +coeff(7)*d2px17+coeff(8)*d2px18;
        d2pz=coeff(1)*d2pz11+coeff(2)*d2pz12+coeff(3)*d2pz13+coeff(4)*d2pz14+coeff(5)*d2pz15+coeff(6)*d2pz16...
            +coeff(7)*d2pz17+coeff(8)*d2pz18;
        
        %%%%%%%%%%%
        d2pxLax=coeff1(1)*d2px11+coeff1(2)*d2px12+coeff1(3)*d2px13+coeff1(4)*d2px14+coeff1(5)*d2px15+coeff1(6)*d2px16...
            +coeff1(7)*d2px17+coeff1(8)*d2px18;
        d2pzLax=coeff1(1)*d2pz11+coeff1(2)*d2pz12+coeff1(3)*d2pz13+coeff1(4)*d2pz14+coeff1(5)*d2pz15+coeff1(6)*d2pz16...
            +coeff1(7)*d2pz17 +coeff1(8)*d2pz18;
        
        temp12=-v.^2.*(d2pzLax+d2pxLax);
        
        d2px4=(circshift(temp12,[0 -1])+circshift(temp12,[0 1]));
        d2pz4 =(circshift(temp12,[-1 0])+circshift(temp12,[1 0]));
        
        d2pxz4=v.^2.*(d2px4+d2pz4-4*temp12)*dt^3/h^3/12;
        %%%%%%%%%%%
        
        
        p=p-dt*v.^2.*(d2px+d2pz)/h +d2pxz4;
        p(zs,xs)= p(zs,xs)+src(it)*dt^2;
        
        d2px1=(circshift(p,[0 -1])-circshift(p,[0 0]));
        %         d2pz2=(circshift(p,[0 -2])-circshift(p,[0 1]));
        %         d2pz3=(circshift(p,[0 -3])-circshift(p,[0 2]));
        %         d2pz4=(circshift(p,[0 -4])-circshift(p,[0 3]));
        %         d2pz5=(circshift(p,[0 -5])-circshift(p,[0 4]));
        %         d2pz6=(circshift(p,[0 -6])-circshift(p,[0 5]));
        %         d2pz7=(circshift(p,[0 -7])-circshift(p,[0 6]));
        
        
        
        d2pz1=(circshift(p,[-1])-circshift(p,[0]));
        %         d2px2=(circshift(p,[-2])-circshift(p,[1]));
        %         d2px3=(circshift(p,[-3])-circshift(p,[2]));
        %         d2px4=(circshift(p,[-4])-circshift(p,[3]));
        %         d2px5=(circshift(p,[-5])-circshift(p,[4]));
        %         d2px6=(circshift(p,[-6])-circshift(p,[5]));
        %         d2px7=(circshift(p,[-7])-circshift(p,[6]));
        %            d2px=coeff(1)*d2px1+coeff(2)*d2px2+coeff(3)*d2px3+coeff(4)*d2px4+coeff(5)*d2px5+coeff(6)*d2px6...
        %             +coeff(7)*d2px7;
        %         d2pz=coeff(1)*d2pz1+coeff(2)*d2pz2+coeff(3)*d2pz3+coeff(4)*d2pz4+coeff(5)*d2pz5+coeff(6)*d2pz6...
        %             +coeff(7)*d2pz7;
        
        
        Vx=Vx-dt*d2px1/h;
        Vz=Vz-dt*d2pz1/h;
        
        
        
        [Vx,Vz]=spongeABC(Vx,Vz,nx,nz,50,50,0.007);
        %         Vx=taper.*Vx;
        %         Vz=Vz.*taper;
        if rem(it,isnap)== 0,
            imagesc(x,z,p), axis equal
            colormap gray
            xlabel('x'),ylabel('z')
            title(sprintf(' Time step: %i - Max ampl: %g ',it,max(max(p))))
            drawnow
        end
        
        
        
        seis_record(it,:)=p(zs,:);
        
        
    end
    
end
toc
save('figure5c_2.mat')
% figure;plot((a),'linewidth',2);hold on;plot((b),'r','linewidth',2)
% xlabel('time(ms)')
% ylabel('Amp')
% legend('receiver A','receiver B')
% grid on

