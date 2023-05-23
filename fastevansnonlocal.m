tsamp = [tout2(10:end-1);-flip(tout(10:end))+tout(end)+tout2(end)];


yp = yout(:,1:4);
vq = griddedInterpolant(tout(10:end),yp(10:end,1),'pchip');

yp2 = yout2(:,1:4);
vq2 = griddedInterpolant(tout2(10:end),yp2(10:end,1),'pchip');

vq3 = griddedInterpolant(tsamp,[yp2(10:end-1,1);flip(yp(10:end,1))],'pchip');
vq3inv = griddedInterpolant(flip(tsamp(end)-tsamp),flip([yp2(10:end-1,1);flip(yp(10:end,1))]),'pchip');

opts = odeset('RelTol',1e-7,'AbsTol',1e-8);
% opts = odeset('RelTol',1e-12,'AbsTol',1e-14);


 cc = 0.19686;
lam = -0.8;
lam = 1+1i;
ee = 1e-4;

% ic = [0,-1.721504100117682,+2.624889820251053];
D0 = 6.0*(-7/12)*(-0.75);
D1 = 6.0*(1.0-7.0/12.0)*(1.0-0.75);


ic = [-1.620151172036441,0,0,0,0,0];
ic2=[sqrt(D1),0,0,0,0,0];
% tout = tout;

clear anglelin Wdiff

anglelin = linspace(0,2*pi-0.01,50);


rmin = -0.8;
%  rmax = 3718.025;
% rmax = 3700.0;
rmax = 800.0;
% rmax = 4000.0;

centrept = (rmin+rmax)/2;
radi = (rmax-rmin)/2;


for i = 1:1

lam = radi*(cos(anglelin(i))+1i*sin(anglelin(i)))+centrept;



a=1
 [t,y] = ode23s(@(t,y) nonlocal(t,y,vq,cc,lam,ee),[0 tout(end)],ic,opts);
% 
a=2
% 
 [t2,y2] = ode23s(@(t2,y2) nonlocal2(t2,y2,vq2,cc,lam,ee),[0 tout2(end)],ic2,opts);

 a=3
% 
   [t3,y3] = ode23s(@(t3,y3) nonlocal2(t3,y3,vq3,cc,lam,ee),[tsamp(1) tsamp(end)],ic2,opts);
% %  
% a=4
% % 
%   [t4,y4] = ode23s(@(t4,y4) nonlocal(t4,y4,vq3inv,cc,lam,ee),[tsamp(1) tsamp(end)],ic,opts);

% ic2=[0.240557412906295,0.051113852040661,0.014889296339861];
% [t3,y3] = ode23s(@(t3,y3) nonlocal2(t3,y3,vq3,cc,lam,ee),[3.95747 3.96],ic2,opts);


 uout2 = vq(t);
 uout3 = vq2(t2);

 
%  Wdiff(i) = dot(y2(end,1:3)+1i*y2(end,4:6),y(end,1:3)+1i*y(end,4:6));

Wdiff(i) = y2(end,1)+1i*y2(end,4)-(y(end,1)+1i*y(end,4))
 
end

%  uoutall = vq3(t3);
%   uoutall2 = vq3inv(t4);

% 
%  figure(2); hold on;
% %  plot(uout2,y(:,2)./y(:,3),'b')
% %  plot(uout3,y2(:,2)./y2(:,3),'r')
% % plot(uoutall,y3(:,2)./y3(:,3),'g')
% 

% figure(2)
%  plot(uout2,y(:,1),'b')
%  plot(uout3,y2(:,1),'r')

%  plot(uoutall2,y4(:,1),'m')
%  plot(uoutall,y3(:,1),'g')

%  hold on; plot(uoutall,y3(:,2)./y3(:,3),'g--','LineWidth',2);
%  hold on; plot(uoutall2,y4(:,2)./y4(:,3),'m--','LineWidth',2);

% hold on
% plot(uout2,y(:,1),'b','LineWidth',2)
% plot([0.522331 0.522331],[-2 2],'k--','LineWidth',2)
% plot([0.811006 0.811006],[-2 2],'k--','LineWidth',2)


function dydt = nonlocal(t,y,vq,cc,lam,ee)
dydt = zeros(6,1);
u=vq(t);
Du = 6.0*(u-7.0/12.0)*(u-0.75);
Rp = -15.0*u^2+12.0*u-1.0;
einv = 1/ee;
mu = real(lam);
omega = imag(lam);

dydt(1) = einv*(Du-y(1)^2-y(3)+y(4)^2);
dydt(2) = -Rp+einv*(y(4)*y(5)-y(1)*y(2))+mu;
dydt(3) = -cc+y(2)+einv*(y(4)*y(6)-y(1)*y(3));
dydt(4) = -einv*(2*y(1)*y(4)+y(6));
dydt(5) = -einv*(y(2)*y(4)+y(1)*y(5))+omega;
dydt(6) = y(5)-einv*(y(3)*y(4)+y(1)*y(6));

dydt = -dydt;
end

function dydt = nonlocal2(t,y,vq,cc,lam,ee)
dydt = zeros(6,1);
u=vq(t);
Du = 6.0*(u-7.0/12.0)*(u-0.75);
Rp = -15.0*u^2+12.0*u-1.0;
einv = 1/ee;
mu = real(lam);
omega = imag(lam);

dydt(1) = einv*(Du-y(1)^2-y(3)+y(4)^2);
dydt(2) = -Rp+einv*(y(4)*y(5)-y(1)*y(2))+mu;
dydt(3) = -cc+y(2)+einv*(y(4)*y(6)-y(1)*y(3));
dydt(4) = -einv*(2*y(1)*y(4)+y(6));
dydt(5) = -einv*(y(2)*y(4)+y(1)*y(5))+omega;
dydt(6) = y(5)-einv*(y(3)*y(4)+y(1)*y(6));
end
