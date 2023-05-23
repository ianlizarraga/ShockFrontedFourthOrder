tsamp = [tout2(10:end-1);-flip(tout(10:end))+tout(end)+tout2(end)];


yp = yout(:,1:4);
vq = griddedInterpolant(tout(10:end),yp(10:end,1),'pchip');

yp2 = yout2(:,1:4);
vq2 = griddedInterpolant(tout2(10:end),yp2(10:end,1),'pchip');
% 
% vq3 = griddedInterpolant(tsamp,[yp2(10:end-1,1);flip(yp(10:end,1))],'pchip');
% vq3inv = griddedInterpolant(flip(tsamp(end)-tsamp),flip([yp2(10:end-1,1);flip(yp(10:end,1))]),'pchip');

opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
% opts = odeset('RelTol',1e-9,'AbsTol',1e-8);
% opts = odeset('RelTol',1e-12,'AbsTol',1e-14);


%  cc = 0.19686;

ee = 0.0001;
% cc = 0.19826;
% delta = 0.099398797595190;


cc = c0;

% cc = 0.19686;
% delta = 0.0;

%  delta = 0.07;
a = delta/c0;

clear anglelin Wdiffdet

anglelin = linspace(0,2*pi-0.01,30);




rmin = -0.82;
rmax = -0.78;

rmin = 150000;
rmax = 200000;


centrept = (rmin+rmax)/2;
radi = (rmax-rmin)/2;


RR = 1e5;

angleclock = linspace(pi/2,-pi/2,30);

lineimag = linspace(RR,1.1,500);

% lineimag = linspace(RR,1.1,10000);

% lineimag = [lineimag(1:50:986),linspace(lineimag(987),lineimag(1000),100)];
% lineimag = [lineimag(1:10:986),linspace(lineimag(987),lineimag(1042),100),lineimag(987:10:end)];
% lineimag = logspace(log(RR)/log(10),log(1.1)/log(10),300);


% lineimag = [logspace(log(RR)/log(10),log(2000)/log(10),20),linspace(1800,1000,50),logspace(log(1000)/log(10),log(1.1)/log(10),20)];
angleclock = linspace(pi/2,-pi/2,10);
angleclock2 = linspace(-pi/2,pi/2,60);

lam2 = [1i*lineimag,cos(angleclock)+1i*sin(angleclock)];
lam2 = [lam2,-1i*flip(lineimag)];
% lam2 = [lam2,1*flip(lineimag)];
lam2 = [lam2,RR*(cos(angleclock2)+1i*sin(angleclock2))];

 
% angleclock2 = linspace(0,2*pi,100);
% lam2 = (10)*(cos(angleclock2)+1i*sin(angleclock2))+20.1;
% 
%  lam2 = linspace(-0.4,1000,1000);
%  
%  lam2 = linspace(0,100,100);
%  
%  lam2 = [linspace(-0.9,0,100),linspace(0.0,95,100)];
 
%  lam2 = linspace(-0.9,95,100);

%  lam2 = 0.0;
 
for i = 1:length(lam2)
    
   lam=lam2(i)
    
    i
    
% lam = radi*(cos(anglelin(i))+1i*sin(anglelin(i)))+centrept;

% ic = [0,-1.721504100117682,+2.624889820251053];
D0 = 6.0*(-7/12)*(-0.75);
D1 = 6.0*(1.0-7.0/12.0)*(1.0-0.75);
Rp1 = -4.0;
Rp0 = -1.0;

ic = [-1.620151172036441,0,0];
ic2=[sqrt(D1),0,0];


Jac1 = [0,1/ee,0,0;(D1-a*ee*lam)/ee,-delta/ee,0,1/ee;Rp1-lam,0,0,0;cc,0,1,0];
Jac0 = [0,1/ee,0,0;(D0-a*ee*lam)/ee,-delta/ee,0,1/ee;Rp0-lam,0,0,0;cc,0,1,0];

[V1,D1] = eigs(Jac1);
[V0,D0] = eigs(Jac0);

[~,m0]=sort(real(diag(D0)),'descend');
[~,m1]=sort(real(diag(D1)),'descend');



% TT = [-1i,0,1,0;0,1i,0,1;0,0,1i,0;0,0,0,-1i];
TT = [-1i,0,1,0;0,1i,0,1;0,0,1,0;0,0,0,1];

PlaneUnst=TT*[V1(:,m1(1)) V1(:,m1(2))];

PlaneUnst0 = TT*[V0(:,m0(1)) V0(:,m0(2))];
PlaneStab0 = TT*[V0(:,m0(3)) V0(:,m0(4))];

WPlane = (PlaneUnst(3:4,1:2))*inv(PlaneUnst(1:2,1:2));

WPlane0 = (PlaneUnst0(3:4,1:2))*inv(PlaneUnst0(1:2,1:2));
SPlane0 = (PlaneStab0(3:4,1:2))*inv(PlaneStab0(1:2,1:2));

weakpart = cc+sqrt(10.0+cc^2+(5/2)*lam);

% icriccati = [0.8*weakpart*D1,-0.8*weakpart*sqrt(D1),D1,-sqrt(D1),0,0,0,0];

icriccati = [real(WPlane(1,1)),real(WPlane(1,2)),real(WPlane(2,1)),real(WPlane(2,2)),imag(WPlane(1,1)),imag(WPlane(1,2)),imag(WPlane(2,1)),imag(WPlane(2,2))];

icriccatis = [real(SPlane0(1,1)),real(SPlane0(1,2)),real(SPlane0(2,1)),real(SPlane0(2,2)),imag(SPlane0(1,1)),imag(SPlane0(1,2)),imag(SPlane0(2,1)),imag(SPlane0(2,2))];
% tout = tout;

a1=1
sflag = -1;
 [t,y] = ode23s(@(t,y) nonlocal(t,y,vq,cc,lam,ee,a,delta,sflag),[0 tout(end)],icriccatis,opts);
% 
a1=2
% 
sflag = 1;
 [t2,y2] = ode23s(@(t2,y2) nonlocal(t2,y2,vq2,cc,lam,ee,a,delta,sflag),[0 tout2(end)],icriccati,opts);
 
% 
% a=3
% 
%   [t3,y3] = ode23s(@(t3,y3) nonlocal2(t3,y3,vq3,cc,lam,ee),[tsamp(1) tsamp(end)],icriccati,opts);
%  
% a=4
% % 
%   [t4,y4] = ode23s(@(t4,y4) nonlocal2(t4,y4,vq3inv,cc,lam,ee),[tsamp(1) tsamp(end)],icriccati,opts);

% ic2=[0.240557412906295,0.051113852040661,0.014889296339861];
% [t3,y3] = ode23s(@(t3,y3) nonlocal2(t3,y3,vq3,cc,lam,ee),[3.95747 3.96],ic2,opts);


% uout2 = vq(t);
% uout3 = vq2(t2);

%  uoutall = vq3(t3);
%   uoutall2 = vq3inv(t4);
% 

Wdiff = y2(end,:)-y(end,:);
% reshape(Wdiff(1:4),2,2);
% reshape(Wdiff(1:4),2,2)'+1i*reshape(Wdiff(5:8),2,2)';
Wdiffdet(i) = det(reshape(Wdiff(1:4),2,2)'+1i*reshape(Wdiff(5:8),2,2)');

Wdiffdet(i)
end
%  figure(2); hold on;
%  plot(uout2,y(:,2)./y(:,3),'b')
%  plot(uout3,y2(:,2)./y2(:,3),'r')
% plot(uoutall,y3(:,2)./y3(:,3),'g')

% plot(uout2,y(:,1),'b')
% plot(uout3,y2(:,1),'r')

%  plot(uoutall2,y4(:,1),'m')
%  plot(uoutall,y3(:,1),'g')

%  hold on; plot(uoutall,y3(:,2)./y3(:,3),'g--','LineWidth',2);
%  hold on; plot(uoutall2,y4(:,2)./y4(:,3),'m--','LineWidth',2);

% hold on
% plot(uout2,y(:,1),'b','LineWidth',2)
% plot([0.522331 0.522331],[-2 2],'k--','LineWidth',2)
% plot([0.811006 0.811006],[-2 2],'k--','LineWidth',2)


function dydt = nonlocal(t,y,vq,cc,lam,ee,a,delta,sflag)
dydt = zeros(8,1);
u=vq(t);

y1=y(1);
y2=y(2);
y3=y(3);
y4=y(4);
y5=y(5);
y6=y(6);
y7=y(7);
y8=y(8);

Rp=-15.0*u^2+12.0*u-1.0;
Du=6.0*(u-7.0/12.0)*(u-0.75);
einv = 1.0/ee;
lamr = real(lam);
lami = imag(lam);

dydt(1) =lami+(-2).*lami.*y1+lami.*y1.^2+Du.*einv.*y2+(-1).*a.*lamr.*y2+( ...
  -1).*y1.*y2+(-1).*Du.*einv.*y1.*y2+a.*lamr.*y1.*y2+(-1).*einv.* ...
  y1.*y3+(-1).*delta.*einv.*y2.*y3+(-2).*lamr.*y5+2.*Rp.*y5+2.* ...
  lamr.*y1.*y5+(-2).*Rp.*y1.*y5+(-1).*cc.*y2.*y5+(-1).*a.*lami.*y2.* ...
  y5+(-1).*lami.*y5.^2+cc.*y6+a.*lami.*y6+(-1).*cc.*y1.*y6+(-1).*a.* ...
  lami.*y1.*y6+einv.*y3.*y6+y5.*y6+Du.*einv.*y5.*y6+(-1).*a.*lamr.* ...
  y5.*y6+einv.*y2.*y7+einv.*y5.*y7+delta.*einv.*y6.*y7;



dydt(5) =(-1).*lamr+Rp+2.*lamr.*y1+(-2).*Rp.*y1+(-1).*lamr.*y1.^2+Rp.* ...
  y1.^2+(-1).*cc.*y2+(-1).*a.*lami.*y2+cc.*y1.*y2+a.*lami.*y1.*y2+( ...
  -1).*einv.*y2.*y3+(-2).*lami.*y5+2.*lami.*y1.*y5+(-1).*y2.*y5+(-1) ...
  .*Du.*einv.*y2.*y5+a.*lamr.*y2.*y5+(-1).*einv.*y3.*y5+lamr.*y5.^2+ ...
  (-1).*Rp.*y5.^2+Du.*einv.*y6+(-1).*a.*lamr.*y6+(-1).*y1.*y6+(-1).* ...
  Du.*einv.*y1.*y6+a.*lamr.*y1.*y6+(-1).*delta.*einv.*y3.*y6+(-1).* ...
  cc.*y5.*y6+(-1).*a.*lami.*y5.*y6+(-1).*einv.*y1.*y7+(-1).*delta.* ...
  einv.*y2.*y7+einv.*y6.*y7;

dydt(2) = einv.*y1+delta.*einv.*y2+(-1).*lami.*y2+lami.*y1.*y2+(-1).*y2.^2+( ...
  -1).*Du.*einv.*y2.^2+a.*lamr.*y2.^2+(-1).*einv.*y1.*y4+(-1).* ...
  delta.*einv.*y2.*y4+lamr.*y2.*y5+(-1).*Rp.*y2.*y5+(-1).*lamr.*y6+ ...
  Rp.*y6+lamr.*y1.*y6+(-1).*Rp.*y1.*y6+(-2).*cc.*y2.*y6+(-2).*a.* ...
  lami.*y2.*y6+einv.*y4.*y6+(-1).*lami.*y5.*y6+y6.^2+Du.*einv.* ...
  y6.^2+(-1).*a.*lamr.*y6.^2+einv.*y2.*y8+einv.*y5.*y8+delta.*einv.* ...
  y6.*y8;

dydt(6) = lamr.*y2+(-1).*Rp.*y2+(-1).*lamr.*y1.*y2+Rp.*y1.*y2+cc.*y2.^2+a.* ...
  lami.*y2.^2+(-1).*einv.*y2.*y4+einv.*y5+lami.*y2.*y5+(-1).*einv.* ...
  y4.*y5+delta.*einv.*y6+(-1).*lami.*y6+lami.*y1.*y6+(-2).*y2.*y6+( ...
  -2).*Du.*einv.*y2.*y6+2.*a.*lamr.*y2.*y6+(-1).*delta.*einv.*y4.* ...
  y6+lamr.*y5.*y6+(-1).*Rp.*y5.*y6+(-1).*cc.*y6.^2+(-1).*a.*lami.* ...
  y6.^2+(-1).*einv.*y1.*y8+(-1).*delta.*einv.*y2.*y8+einv.*y6.*y8;

dydt(3) = y1+(-1).*lami.*y3+lami.*y1.*y3+(-1).*einv.*y3.^2+Du.*einv.*y4+(-1) ...
  .*a.*lamr.*y4+(-1).*y1.*y4+(-1).*Du.*einv.*y1.*y4+a.*lamr.*y1.*y4+ ...
  (-1).*delta.*einv.*y3.*y4+cc.*y5+lamr.*y3.*y5+(-1).*Rp.*y3.*y5+( ...
  -1).*cc.*y4.*y5+(-1).*a.*lami.*y4.*y5+(-1).*lamr.*y7+Rp.*y7+lamr.* ...
  y1.*y7+(-1).*Rp.*y1.*y7+einv.*y4.*y7+(-1).*lami.*y5.*y7+einv.* ...
  y7.^2+cc.*y8+a.*lami.*y8+(-1).*cc.*y1.*y8+(-1).*a.*lami.*y1.*y8+ ...
  einv.*y3.*y8+y5.*y8+Du.*einv.*y5.*y8+(-1).*a.*lamr.*y5.*y8+delta.* ...
  einv.*y7.*y8;

dydt(7) = cc+(-1).*cc.*y1+lamr.*y3+(-1).*Rp.*y3+(-1).*lamr.*y1.*y3+Rp.*y1.* ...
  y3+(-1).*cc.*y4+(-1).*a.*lami.*y4+cc.*y1.*y4+a.*lami.*y1.*y4+(-1) ...
  .*einv.*y3.*y4+y5+lami.*y3.*y5+(-1).*y4.*y5+(-1).*Du.*einv.*y4.* ...
  y5+a.*lamr.*y4.*y5+(-1).*lami.*y7+lami.*y1.*y7+(-2).*einv.*y3.*y7+ ...
  (-1).*delta.*einv.*y4.*y7+lamr.*y5.*y7+(-1).*Rp.*y5.*y7+Du.*einv.* ...
  y8+(-1).*a.*lamr.*y8+(-1).*y1.*y8+(-1).*Du.*einv.*y1.*y8+a.*lamr.* ...
  y1.*y8+(-1).*delta.*einv.*y3.*y8+(-1).*cc.*y5.*y8+(-1).*a.*lami.* ...
  y5.*y8+einv.*y7.*y8;



dydt(4) = y2+einv.*y3+lami.*y2.*y3+delta.*einv.*y4+(-1).*y2.*y4+(-1).*Du.* ...
  einv.*y2.*y4+a.*lamr.*y2.*y4+(-1).*einv.*y3.*y4+(-1).*delta.* ...
  einv.*y4.^2+cc.*y6+lamr.*y3.*y6+(-1).*Rp.*y3.*y6+(-1).*cc.*y4.*y6+ ...
  (-1).*a.*lami.*y4.*y6+lamr.*y2.*y7+(-1).*Rp.*y2.*y7+(-1).*lami.* ...
  y6.*y7+(-1).*cc.*y2.*y8+(-1).*a.*lami.*y2.*y8+2.*einv.*y4.*y8+y6.* ...
  y8+Du.*einv.*y6.*y8+(-1).*a.*lamr.*y6.*y8+einv.*y7.*y8+delta.* ...
  einv.*y8.^2;


dydt(8) = (-1).*cc.*y2+(-1).*lamr.*y2.*y3+Rp.*y2.*y3+cc.*y2.*y4+a.*lami.* ...
  y2.*y4+(-1).*einv.*y4.^2+y6+lami.*y3.*y6+(-1).*y4.*y6+(-1).*Du.* ...
  einv.*y4.*y6+a.*lamr.*y4.*y6+einv.*y7+lami.*y2.*y7+(-1).*einv.* ...
  y4.*y7+lamr.*y6.*y7+(-1).*Rp.*y6.*y7+delta.*einv.*y8+(-1).*y2.*y8+ ...
  (-1).*Du.*einv.*y2.*y8+a.*lamr.*y2.*y8+(-1).*einv.*y3.*y8+(-2).* ...
  delta.*einv.*y4.*y8+(-1).*cc.*y6.*y8+(-1).*a.*lami.*y6.*y8+einv.* ...
  y8.^2;


dydt = sflag*dydt;
end

function dydt = nonlocal2(t,y,vq,cc,lam,ee,a,delta)
dydt = zeros(8,1);
u=vq(t);

y1=y(1);
y2=y(2);
y3=y(3);
y4=y(4);
y5=y(5);
y6=y(6);
y7=y(7);
y8=y(8);

Rp=-15.0*u^2+12.0*u-1.0;
Du=6.0*(u-7.0/12.0)*(u-0.75);
einv = 1.0/ee;
lamr = real(lam);
lami = imag(lam);

dydt(1) =lami+(-2).*lami.*y1+lami.*y1.^2+Du.*einv.*y2+(-1).*a.*lamr.*y2+( ...
  -1).*y1.*y2+(-1).*Du.*einv.*y1.*y2+a.*lamr.*y1.*y2+(-1).*einv.* ...
  y1.*y3+(-1).*delta.*einv.*y2.*y3+(-2).*lamr.*y5+2.*Rp.*y5+2.* ...
  lamr.*y1.*y5+(-2).*Rp.*y1.*y5+(-1).*cc.*y2.*y5+(-1).*a.*lami.*y2.* ...
  y5+(-1).*lami.*y5.^2+cc.*y6+a.*lami.*y6+(-1).*cc.*y1.*y6+(-1).*a.* ...
  lami.*y1.*y6+einv.*y3.*y6+y5.*y6+Du.*einv.*y5.*y6+(-1).*a.*lamr.* ...
  y5.*y6+einv.*y2.*y7+einv.*y5.*y7+delta.*einv.*y6.*y7;



dydt(5) =(-1).*lamr+Rp+2.*lamr.*y1+(-2).*Rp.*y1+(-1).*lamr.*y1.^2+Rp.* ...
  y1.^2+(-1).*cc.*y2+(-1).*a.*lami.*y2+cc.*y1.*y2+a.*lami.*y1.*y2+( ...
  -1).*einv.*y2.*y3+(-2).*lami.*y5+2.*lami.*y1.*y5+(-1).*y2.*y5+(-1) ...
  .*Du.*einv.*y2.*y5+a.*lamr.*y2.*y5+(-1).*einv.*y3.*y5+lamr.*y5.^2+ ...
  (-1).*Rp.*y5.^2+Du.*einv.*y6+(-1).*a.*lamr.*y6+(-1).*y1.*y6+(-1).* ...
  Du.*einv.*y1.*y6+a.*lamr.*y1.*y6+(-1).*delta.*einv.*y3.*y6+(-1).* ...
  cc.*y5.*y6+(-1).*a.*lami.*y5.*y6+(-1).*einv.*y1.*y7+(-1).*delta.* ...
  einv.*y2.*y7+einv.*y6.*y7;

dydt(2) = einv.*y1+delta.*einv.*y2+(-1).*lami.*y2+lami.*y1.*y2+(-1).*y2.^2+( ...
  -1).*Du.*einv.*y2.^2+a.*lamr.*y2.^2+(-1).*einv.*y1.*y4+(-1).* ...
  delta.*einv.*y2.*y4+lamr.*y2.*y5+(-1).*Rp.*y2.*y5+(-1).*lamr.*y6+ ...
  Rp.*y6+lamr.*y1.*y6+(-1).*Rp.*y1.*y6+(-2).*cc.*y2.*y6+(-2).*a.* ...
  lami.*y2.*y6+einv.*y4.*y6+(-1).*lami.*y5.*y6+y6.^2+Du.*einv.* ...
  y6.^2+(-1).*a.*lamr.*y6.^2+einv.*y2.*y8+einv.*y5.*y8+delta.*einv.* ...
  y6.*y8;

dydt(6) = lamr.*y2+(-1).*Rp.*y2+(-1).*lamr.*y1.*y2+Rp.*y1.*y2+cc.*y2.^2+a.* ...
  lami.*y2.^2+(-1).*einv.*y2.*y4+einv.*y5+lami.*y2.*y5+(-1).*einv.* ...
  y4.*y5+delta.*einv.*y6+(-1).*lami.*y6+lami.*y1.*y6+(-2).*y2.*y6+( ...
  -2).*Du.*einv.*y2.*y6+2.*a.*lamr.*y2.*y6+(-1).*delta.*einv.*y4.* ...
  y6+lamr.*y5.*y6+(-1).*Rp.*y5.*y6+(-1).*cc.*y6.^2+(-1).*a.*lami.* ...
  y6.^2+(-1).*einv.*y1.*y8+(-1).*delta.*einv.*y2.*y8+einv.*y6.*y8;

dydt(3) = y1+(-1).*lami.*y3+lami.*y1.*y3+(-1).*einv.*y3.^2+Du.*einv.*y4+(-1) ...
  .*a.*lamr.*y4+(-1).*y1.*y4+(-1).*Du.*einv.*y1.*y4+a.*lamr.*y1.*y4+ ...
  (-1).*delta.*einv.*y3.*y4+cc.*y5+lamr.*y3.*y5+(-1).*Rp.*y3.*y5+( ...
  -1).*cc.*y4.*y5+(-1).*a.*lami.*y4.*y5+(-1).*lamr.*y7+Rp.*y7+lamr.* ...
  y1.*y7+(-1).*Rp.*y1.*y7+einv.*y4.*y7+(-1).*lami.*y5.*y7+einv.* ...
  y7.^2+cc.*y8+a.*lami.*y8+(-1).*cc.*y1.*y8+(-1).*a.*lami.*y1.*y8+ ...
  einv.*y3.*y8+y5.*y8+Du.*einv.*y5.*y8+(-1).*a.*lamr.*y5.*y8+delta.* ...
  einv.*y7.*y8;

dydt(7) = cc+(-1).*cc.*y1+lamr.*y3+(-1).*Rp.*y3+(-1).*lamr.*y1.*y3+Rp.*y1.* ...
  y3+(-1).*cc.*y4+(-1).*a.*lami.*y4+cc.*y1.*y4+a.*lami.*y1.*y4+(-1) ...
  .*einv.*y3.*y4+y5+lami.*y3.*y5+(-1).*y4.*y5+(-1).*Du.*einv.*y4.* ...
  y5+a.*lamr.*y4.*y5+(-1).*lami.*y7+lami.*y1.*y7+(-2).*einv.*y3.*y7+ ...
  (-1).*delta.*einv.*y4.*y7+lamr.*y5.*y7+(-1).*Rp.*y5.*y7+Du.*einv.* ...
  y8+(-1).*a.*lamr.*y8+(-1).*y1.*y8+(-1).*Du.*einv.*y1.*y8+a.*lamr.* ...
  y1.*y8+(-1).*delta.*einv.*y3.*y8+(-1).*cc.*y5.*y8+(-1).*a.*lami.* ...
  y5.*y8+einv.*y7.*y8;



dydt(4) = y2+einv.*y3+lami.*y2.*y3+delta.*einv.*y4+(-1).*y2.*y4+(-1).*Du.* ...
  einv.*y2.*y4+a.*lamr.*y2.*y4+(-1).*einv.*y3.*y4+(-1).*delta.* ...
  einv.*y4.^2+cc.*y6+lamr.*y3.*y6+(-1).*Rp.*y3.*y6+(-1).*cc.*y4.*y6+ ...
  (-1).*a.*lami.*y4.*y6+lamr.*y2.*y7+(-1).*Rp.*y2.*y7+(-1).*lami.* ...
  y6.*y7+(-1).*cc.*y2.*y8+(-1).*a.*lami.*y2.*y8+2.*einv.*y4.*y8+y6.* ...
  y8+Du.*einv.*y6.*y8+(-1).*a.*lamr.*y6.*y8+einv.*y7.*y8+delta.* ...
  einv.*y8.^2;


dydt(8) = (-1).*cc.*y2+(-1).*lamr.*y2.*y3+Rp.*y2.*y3+cc.*y2.*y4+a.*lami.* ...
  y2.*y4+(-1).*einv.*y4.^2+y6+lami.*y3.*y6+(-1).*y4.*y6+(-1).*Du.* ...
  einv.*y4.*y6+a.*lamr.*y4.*y6+einv.*y7+lami.*y2.*y7+(-1).*einv.* ...
  y4.*y7+lamr.*y6.*y7+(-1).*Rp.*y6.*y7+delta.*einv.*y8+(-1).*y2.*y8+ ...
  (-1).*Du.*einv.*y2.*y8+a.*lamr.*y2.*y8+(-1).*einv.*y3.*y8+(-2).* ...
  delta.*einv.*y4.*y8+(-1).*cc.*y6.*y8+(-1).*a.*lami.*y6.*y8+einv.* ...
  y8.^2;

dydt = dydt;
end
