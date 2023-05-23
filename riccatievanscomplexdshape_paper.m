tsamp = [tout2(10:end-1);-flip(tout(10:end))+tout(end)+tout2(end)];


yp = yout(:,1:4);
yp(:,3) = -yp(:,3);
yp(:,4) = -yp(:,4);
vq = griddedInterpolant(tout(10:end),yp(10:end,1),'pchip');

yp2 = yout2(:,1:4);
yp2(:,3) = -yp2(:,3);
yp2(:,4) = -yp2(:,4);
vq2 = griddedInterpolant(tout2(10:end),yp2(10:end,1),'pchip');

vq3 = griddedInterpolant(tsamp,[yp2(10:end-1,1);flip(yp(10:end,1))],'pchip');
vq3inv = griddedInterpolant(flip(tsamp(end)-tsamp),flip([yp2(10:end-1,1);flip(yp(10:end,1))]),'pchip');

opts = odeset('RelTol',1e-9,'AbsTol',1e-11);
opts = odeset('RelTol',1e-6,'AbsTol',1e-6);



 cc = 0.19686;

ee = 1e-4;

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

lineimag = linspace(RR,1.1,300);
lineimag = logspace(log(RR)/log(10),log(1.1)/log(10),300);


lineimag = [logspace(log(RR)/log(10),log(2000)/log(10),20),linspace(1800,1000,50),logspace(log(1000)/log(10),log(1.1)/log(10),20)];
angleclock = linspace(pi/2,-pi/2,30);

lam2 = [1i*lineimag,cos(angleclock)+1i*sin(angleclock)];

lam2 = [lam2,-1i*flip(lineimag)];
lam2 = [lam2,RR*(cos(flip(angleclock))+1i*sin(flip(angleclock)))];

        
for i = 1:length(lam2)
    
   lam=lam2(i)
    
    i
    

D0 = 6.0*(-7/12)*(-0.75);
D1 = 6.0*(1.0-7.0/12.0)*(1.0-0.75);
Rp1 = -4.0;
Rp0 = -1.0;

ic = [-1.620151172036441,0,0];
ic2=[sqrt(D1),0,0];


Jac1 = [0,1/ee,0,0;D1/ee,0,0,-1/ee;lam-Rp1,0,0,0;-cc,0,1,0];
Jac0 = [0,1/ee,0,0;D0/ee,0,0,-1/ee;lam-Rp0,0,0,0;-cc,0,1,0];

[V1,D1] = eigs(Jac1);
[V0,D0] = eigs(Jac0);


TT = [-1i,0,1,0;0,1i,0,1;0,0,1,0;0,0,0,1];

PlaneUnst=TT*[V1(:,1) V1(:,4)];

PlaneUnst0 = TT*[V0(:,1) V0(:,4)];
PlaneStab0 = TT*[V0(:,2) V0(:,3)];

WPlane = (PlaneUnst(3:4,1:2))*inv(PlaneUnst(1:2,1:2));

WPlane0 = (PlaneUnst0(3:4,1:2))*inv(PlaneUnst0(1:2,1:2));
SPlane0 = (PlaneStab0(3:4,1:2))*inv(PlaneStab0(1:2,1:2));

weakpart = cc+sqrt(10.0+cc^2+(5/2)*lam);



icriccati = [real(WPlane(1,1)),real(WPlane(1,2)),real(WPlane(2,1)),real(WPlane(2,2)),imag(WPlane(1,1)),imag(WPlane(1,2)),imag(WPlane(2,1)),imag(WPlane(2,2))];

icriccatis = [real(SPlane0(1,1)),real(SPlane0(1,2)),real(SPlane0(2,1)),real(SPlane0(2,2)),imag(SPlane0(1,1)),imag(SPlane0(1,2)),imag(SPlane0(2,1)),imag(SPlane0(2,2))];
% tout = tout;

a=1
 [t,y] = ode23s(@(t,y) nonlocal(t,y,vq,cc,lam,ee),[0 tout(end)],icriccatis,opts);
% 
a=2
% 
 [t2,y2] = ode23s(@(t2,y2) nonlocal2(t2,y2,vq2,cc,lam,ee),[0 tout2(end)],icriccati,opts);

% 

Wdiff = y2(end,:)-y(end,:);

Wdiffdet(i) = det(reshape(Wdiff(1:4),2,2)'+1i*reshape(Wdiff(5:8),2,2)');

end


figure(2); clf; hold on;
plot(real(Wdiffdet),imag(Wdiffdet))



function dydt = nonlocal(t,y,vq,cc,lam,ee)
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
Dp=6.0*(u-7.0/12.0)*(u-0.75);
einv = 1.0/ee;
lamr = real(lam);
lami = imag(lam);

dydt(1) = (-1)*lami+2*lami*y1+(-1)*lami*y1^2+Dp*einv*y2+(-1)*y1* ...
  y2+(-1)*Dp*einv*y1*y2+(-1)*einv*y1*y3+2*lamr*y5+(-2)* ...
  Rp*y5+(-2)*lamr*y1*y5+2*Rp*y1*y5+cc*y2*y5+lami*y5^2+( ...
  -1)*cc*y6+cc*y1*y6+(-1)*einv*y3*y6+y5*y6+Dp*einv*y5*y6+ ...
  (-1)*einv*y2*y7+einv*y5*y7;

dydt(5) = lamr+(-1)*Rp+(-2)*lamr*y1+2*Rp*y1+lamr*y1^2+(-1)*Rp* ...
  y1^2+cc*y2+(-1)*cc*y1*y2+einv*y2*y3+2*lami*y5+(-2)* ...
  lami*y1*y5+(-1)*y2*y5+(-1)*Dp*einv*y2*y5+(-1)*einv*y3* ...
  y5+(-1)*lamr*y5^2+Rp*y5^2+Dp*einv*y6+(-1)*y1*y6+(-1)* ...
  Dp*einv*y1*y6+cc*y5*y6+(-1)*einv*y1*y7+(-1)*einv*y6*y7; 

dydt(2) = einv*y1+lami*y2+(-1)*lami*y1*y2+(-1)*y2^2+(-1)*Dp*einv* ...
  y2^2+(-1)*einv*y1*y4+(-1)*lamr*y2*y5+Rp*y2*y5+lamr*y6+( ...
  -1)*Rp*y6+(-1)*lamr*y1*y6+Rp*y1*y6+2*cc*y2*y6+(-1)* ...
  einv*y4*y6+lami*y5*y6+y6^2+Dp*einv*y6^2+(-1)*einv*y2* ...
  y8+einv*y5*y8;

dydt(6) = (-1)*lamr*y2+Rp*y2+lamr*y1*y2+(-1)*Rp*y1*y2+(-1)*cc* ...
  y2^2+einv*y2*y4+einv*y5+(-1)*lami*y2*y5+(-1)*einv*y4*y5+ ...
  lami*y6+(-1)*lami*y1*y6+(-2)*y2*y6+(-2)*Dp*einv*y2*y6+( ...
  -1)*lamr*y5*y6+Rp*y5*y6+cc*y6^2+(-1)*einv*y1*y8+(-1)* ...
  einv*y6*y8;

dydt(3) = y1+lami*y3+(-1)*lami*y1*y3+(-1)*einv*y3^2+Dp*einv*y4+(-1) ...
  *y1*y4+(-1)*Dp*einv*y1*y4+(-1)*cc*y5+(-1)*lamr*y3*y5+ ...
  Rp*y3*y5+cc*y4*y5+lamr*y7+(-1)*Rp*y7+(-1)*lamr*y1*y7+ ...
  Rp*y1*y7+(-1)*einv*y4*y7+lami*y5*y7+einv*y7^2+(-1)*cc* ...
  y8+cc*y1*y8+(-1)*einv*y3*y8+y5*y8+Dp*einv*y5*y8;

dydt(7) = (-1)*cc+cc*y1+(-1)*lamr*y3+Rp*y3+lamr*y1*y3+(-1)*Rp*y1* ...
  y3+cc*y4+(-1)*cc*y1*y4+einv*y3*y4+y5+(-1)*lami*y3*y5+(-1) ...
  *y4*y5+(-1)*Dp*einv*y4*y5+lami*y7+(-1)*lami*y1*y7+(-2)* ...
  einv*y3*y7+(-1)*lamr*y5*y7+Rp*y5*y7+Dp*einv*y8+(-1)*y1* ...
  y8+(-1)*Dp*einv*y1*y8+cc*y5*y8+(-1)*einv*y7*y8;


dydt(4) = y2+einv*y3+(-1)*lami*y2*y3+(-1)*y2*y4+(-1)*Dp*einv*y2* ...
  y4+(-1)*einv*y3*y4+(-1)*cc*y6+(-1)*lamr*y3*y6+Rp*y3*y6+ ...
  cc*y4*y6+(-1)*lamr*y2*y7+Rp*y2*y7+lami*y6*y7+cc*y2*y8+( ...
  -2)*einv*y4*y8+y6*y8+Dp*einv*y6*y8+einv*y7*y8;

dydt(8) = cc*y2+lamr*y2*y3+(-1)*Rp*y2*y3+(-1)*cc*y2*y4+einv*y4^2+ ...
  y6+(-1)*lami*y3*y6+(-1)*y4*y6+(-1)*Dp*einv*y4*y6+einv* ...
  y7+(-1)*lami*y2*y7+(-1)*einv*y4*y7+(-1)*lamr*y6*y7+Rp* ...
  y6*y7+(-1)*y2*y8+(-1)*Dp*einv*y2*y8+(-1)*einv*y3*y8+cc* ...
  y6*y8+(-1)*einv*y8^2;
dydt = -dydt;
end

function dydt = nonlocal2(t,y,vq,cc,lam,ee)
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
Dp=6.0*(u-7.0/12.0)*(u-0.75);
einv = 1.0/ee;
lamr = real(lam);
lami = imag(lam);

dydt(1) = (-1)*lami+2*lami*y1+(-1)*lami*y1^2+Dp*einv*y2+(-1)*y1* ...
  y2+(-1)*Dp*einv*y1*y2+(-1)*einv*y1*y3+2*lamr*y5+(-2)* ...
  Rp*y5+(-2)*lamr*y1*y5+2*Rp*y1*y5+cc*y2*y5+lami*y5^2+( ...
  -1)*cc*y6+cc*y1*y6+(-1)*einv*y3*y6+y5*y6+Dp*einv*y5*y6+ ...
  (-1)*einv*y2*y7+einv*y5*y7;

dydt(5) = lamr+(-1)*Rp+(-2)*lamr*y1+2*Rp*y1+lamr*y1^2+(-1)*Rp* ...
  y1^2+cc*y2+(-1)*cc*y1*y2+einv*y2*y3+2*lami*y5+(-2)* ...
  lami*y1*y5+(-1)*y2*y5+(-1)*Dp*einv*y2*y5+(-1)*einv*y3* ...
  y5+(-1)*lamr*y5^2+Rp*y5^2+Dp*einv*y6+(-1)*y1*y6+(-1)* ...
  Dp*einv*y1*y6+cc*y5*y6+(-1)*einv*y1*y7+(-1)*einv*y6*y7; 

dydt(2) = einv*y1+lami*y2+(-1)*lami*y1*y2+(-1)*y2^2+(-1)*Dp*einv* ...
  y2^2+(-1)*einv*y1*y4+(-1)*lamr*y2*y5+Rp*y2*y5+lamr*y6+( ...
  -1)*Rp*y6+(-1)*lamr*y1*y6+Rp*y1*y6+2*cc*y2*y6+(-1)* ...
  einv*y4*y6+lami*y5*y6+y6^2+Dp*einv*y6^2+(-1)*einv*y2* ...
  y8+einv*y5*y8;

dydt(6) = (-1)*lamr*y2+Rp*y2+lamr*y1*y2+(-1)*Rp*y1*y2+(-1)*cc* ...
  y2^2+einv*y2*y4+einv*y5+(-1)*lami*y2*y5+(-1)*einv*y4*y5+ ...
  lami*y6+(-1)*lami*y1*y6+(-2)*y2*y6+(-2)*Dp*einv*y2*y6+( ...
  -1)*lamr*y5*y6+Rp*y5*y6+cc*y6^2+(-1)*einv*y1*y8+(-1)* ...
  einv*y6*y8;

dydt(3) = y1+lami*y3+(-1)*lami*y1*y3+(-1)*einv*y3^2+Dp*einv*y4+(-1) ...
  *y1*y4+(-1)*Dp*einv*y1*y4+(-1)*cc*y5+(-1)*lamr*y3*y5+ ...
  Rp*y3*y5+cc*y4*y5+lamr*y7+(-1)*Rp*y7+(-1)*lamr*y1*y7+ ...
  Rp*y1*y7+(-1)*einv*y4*y7+lami*y5*y7+einv*y7^2+(-1)*cc* ...
  y8+cc*y1*y8+(-1)*einv*y3*y8+y5*y8+Dp*einv*y5*y8;

dydt(7) = (-1)*cc+cc*y1+(-1)*lamr*y3+Rp*y3+lamr*y1*y3+(-1)*Rp*y1* ...
  y3+cc*y4+(-1)*cc*y1*y4+einv*y3*y4+y5+(-1)*lami*y3*y5+(-1) ...
  *y4*y5+(-1)*Dp*einv*y4*y5+lami*y7+(-1)*lami*y1*y7+(-2)* ...
  einv*y3*y7+(-1)*lamr*y5*y7+Rp*y5*y7+Dp*einv*y8+(-1)*y1* ...
  y8+(-1)*Dp*einv*y1*y8+cc*y5*y8+(-1)*einv*y7*y8;


dydt(4) = y2+einv*y3+(-1)*lami*y2*y3+(-1)*y2*y4+(-1)*Dp*einv*y2* ...
  y4+(-1)*einv*y3*y4+(-1)*cc*y6+(-1)*lamr*y3*y6+Rp*y3*y6+ ...
  cc*y4*y6+(-1)*lamr*y2*y7+Rp*y2*y7+lami*y6*y7+cc*y2*y8+( ...
  -2)*einv*y4*y8+y6*y8+Dp*einv*y6*y8+einv*y7*y8;

dydt(8) = cc*y2+lamr*y2*y3+(-1)*Rp*y2*y3+(-1)*cc*y2*y4+einv*y4^2+ ...
  y6+(-1)*lami*y3*y6+(-1)*y4*y6+(-1)*Dp*einv*y4*y6+einv* ...
  y7+(-1)*lami*y2*y7+(-1)*einv*y4*y7+(-1)*lamr*y6*y7+Rp* ...
  y6*y7+(-1)*y2*y8+(-1)*Dp*einv*y2*y8+(-1)*einv*y3*y8+cc* ...
  y6*y8+(-1)*einv*y8^2;
end
