% ʱ���ѹ� 63.562920 �롣
clear
clc %%%%%%%
close all
% Elapsed time is 9.253353 seconds.
nt=1500;    % number of time steps
eps=.6;     % stability
isnap=60;    % snapshot sampling

nx=450;
nz=450;

v=ones(nz,nx)*3500;
v(1:nz/2,:)=2800;


p=zeros([nz nx]); pboundarynew=p;pdan=p;
dx=20;  %calculate space increment
h=dx;
x=(0:(nx-1))*dx;
z=(0:(nz-1))*dx;  % initialize space coordinates
dt=0.002; % calculate time step from stability criterion
tau=dt;
r2=v.*v.*dt*dt/dx/dx;

f0=21.5*pi;
t=(1:nt)*dt;
t0=4/f0;                       % initialize time axis
src=10^6*exp(-f0^2*(t-t0).*(t-t0));              % source time function
src=-diff((src))/dx^2;				% time derivative to obtain gaussian


nt=nt-1;
B = ifft(exp(-2i*sin([0:nt-1]*pi/(2*nt))'*[0:nt-1]),2*nt,'symmetric');
T = B(1:nt,1:nt) ; % <- The Forward Time Dispersion Transform matrix
src = T * src(:);    % <- Transforming the 1xN source-time vector f

seis_record=zeros(nt,nx);
p=zeros([nz nx]); pnew=p; pold=p;
d2px=p;
d2pz=p;

% Source location
zs=60;
xs=nz/2;

h=dx;
r=v*dt/h;


for loop=1:1
    tic
    
    p=zeros([nz nx]); Vx=p; Vz=p;
    
    
    
    coeff=[ 1.24275, -0.113493, 0.0272555, -0.00731429, 0.00173544, -0.000303651, 0.0000284348];
    for it=1:nt-2,
        
        d1px11=Vx-circshift(Vx,[0 1]);
        d1px12=(circshift(Vx,[0 -1])-circshift(Vx,[0 2]));
        d1px13=(circshift(Vx,[0 -2])-circshift(Vx,[0 3]));
        d1px14=(circshift(Vx,[0 -3])-circshift(Vx,[0 4]));
        d1px15=(circshift(Vx,[0 -4])-circshift(Vx,[0 5]));
        d1px16=(circshift(Vx,[0 -5])-circshift(Vx,[0 6]));
        d1px17=(circshift(Vx,[0 -6])-circshift(Vx,[0 7]));
        
        
        
        d1pz11=Vz-circshift(Vz,[1 0]);
        d1pz12=(circshift(Vz,[-1 0])-circshift(Vz,[2 0]));
        d1pz13=(circshift(Vz,[-2 0])-circshift(Vz,[3 0]));
        d1pz14=(circshift(Vz,[-3 0])-circshift(Vz,[4 0]));
        d1pz15=(circshift(Vz,[-4 0])-circshift(Vz,[5 0]));
        d1pz16=(circshift(Vz,[-5 0])-circshift(Vz,[6 0]));
        d1pz17=(circshift(Vz,[-6 0])-circshift(Vz,[7 0]));
        
        d1px=coeff(1)*d1px11+coeff(2)*d1px12+coeff(3)*d1px13+coeff(4)*d1px14+coeff(5)*d1px15+coeff(6)*d1px16...
            +coeff(7)*d1px17;
        d1pz=coeff(1)*d1pz11+coeff(2)*d1pz12+coeff(3)*d1pz13+coeff(4)*d1pz14+coeff(5)*d1pz15+coeff(6)*d1pz16...
            +coeff(7)*d1pz17;
        
        %         Q1=-v.^2.*(d1px+d1pz)/h;
        
        
        p=p-dt*v.^2.*(d1px+d1pz)/h;
        
        p(zs,xs)= p(zs,xs)+src(it)*dt^2;
        
        
        
        
        d1px1=(circshift(p,[0 -1])-circshift(p,[0 0]));
        d1px2=(circshift(p,[0 -2])-circshift(p,[0 1]));
        d1px3=(circshift(p,[0 -3])-circshift(p,[0 2]));
        d1px4=(circshift(p,[0 -4])-circshift(p,[0 3]));
        d1px5=(circshift(p,[0 -5])-circshift(p,[0 4]));
        d1px6=(circshift(p,[0 -6])-circshift(p,[0 5]));
        d1px7=(circshift(p,[0 -7])-circshift(p,[0 6]));
        
        d1pz1=(circshift(p,[-1])-circshift(p,[0]));
        d1pz2=(circshift(p,[-2])-circshift(p,[1]));
        d1pz3=(circshift(p,[-3])-circshift(p,[2]));
        d1pz4=(circshift(p,[-4])-circshift(p,[3]));
        d1pz5=(circshift(p,[-5])-circshift(p,[4]));
        d1pz6=(circshift(p,[-6])-circshift(p,[5]));
        d1pz7=(circshift(p,[-7])-circshift(p,[6]));
        
        
        
        d1px=coeff(1)*d1px1+coeff(2)*d1px2+coeff(3)*d1px3+coeff(4)*d1px4+coeff(5)*d1px5+coeff(6)*d1px6...
            +coeff(7)*d1px7;
        d1pz=coeff(1)*d1pz1+coeff(2)*d1pz2+coeff(3)*d1pz3+coeff(4)*d1pz4+coeff(5)*d1pz5+coeff(6)*d1pz6...
            +coeff(7)*d1pz7;
        
        Vx=Vx-dt*d1px/h;
        Vz=Vz-dt*d1pz/h;
        
        
        
        [Vx,Vz]=spongeABC(Vx,Vz,nx,nz,50,50,0.007);
        
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

% figure;plot((aa),'linewidth',2);hold on;plot((b),'r','linewidth',2)
% xlabel('time(ms)')
% ylabel('Amp')
% legend('receiver A','receiver B')
% grid on
save('figure3d.mat')