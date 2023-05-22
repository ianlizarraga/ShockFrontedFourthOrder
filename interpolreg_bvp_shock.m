%% 

% This notebook approximates a heteroclinic orbit of the 4D fast-slow travelling wave system for small eps > 0,
% with initial data a singular heteroclinic orbit defined using the reduced  2D system on the critical manifold.

clear all;

load('reducedwave')

eps = 0.0001;


%% Newton/BVP combination method to compute heteroclinic for eps > 0


for j = 1:15
j
    
c0 = c0lin;

whlin = linspace(wh-1e-10,wh+1e-10,2);
cclin = linspace(c0lin-1e-10,c0lin+1e-10,2);


for i = 1:length(whlin)

yb4 = whlin(i);


xlin = linspace(0,1.0,100);


 yinit2 = [xlin(1:end-1) 1; 0*xlin(1:end-1) 1; xlin(1:end-1) 1];


solinit = bvpinit(xlin, @guess);



solinit.x = -ts1;
VV = [pts1(:,1),pts1(:,1)*0,pts1(:,2),-beta*((pts1(:,1).^3)/3-(g1+g2)*(pts1(:,1).^2)/2+g1*g2*pts1(:,1))]';


ya1 = VV(1,1);
ya3 = VV(3,1);

solinit.y = VV;
yinit2 = solinit.y; yinit2 = yinit2';

options = bvpset('RelTol',1e-5,'Nmax',150000);

solinit.yinit=[];
sol = bvp4c(@(x,y) bvpfcn(x,y,eps,delta,c0,beta,g1,g2,fmag,alpha),@(ya,yb) bcfcn(ya,yb,ya1,ya3,yb4), solinit, options);
yout = sol.y';
tout = sol.x';
% 
solinit2 = bvpinit(xlin, @guess);
solinit2.x = ts2;
VV2 = [pts2(:,1),pts2(:,1)*0,pts2(:,2),-beta*((pts2(:,1).^3)/3-(g1+g2)*(pts2(:,1).^2)/2+g1*g2*pts2(:,1))]';

ya1 = VV2(1,1);
ya3 = VV2(3,1);

solinit2.y = VV2;
 yinit3 = solinit2.y; yinit3 = yinit3';
 solinit2.yinit=[];
% 
 sol2 = bvp4c(@(x,y) bvpfcn2(x,y,eps,delta,c0,beta,g1,g2,fmag,alpha),@(ya,yb) bcfcn2(ya,yb,ya1,ya3,yb4), solinit2, options);
 yout2 = sol2.y';
 tout2 = sol2.x';

Jacwh(:,i) = yout2(end,2:3)'-yout(end,2:3)';

end



yb4 = wh;

for i = 1:length(cclin)

c0=cclin(i);


xlin = linspace(0,1.0,100);


 yinit2 = [xlin(1:end-1) 1; 0*xlin(1:end-1) 1; xlin(1:end-1) 1];


solinit = bvpinit(xlin, @guess);



solinit.x = -ts1;
VV = [pts1(:,1),pts1(:,1)*0,pts1(:,2),-beta*((pts1(:,1).^3)/3-(g1+g2)*(pts1(:,1).^2)/2+g1*g2*pts1(:,1))]';


ya1 = VV(1,1);
ya3 = VV(3,1);

solinit.y = VV;
yinit2 = solinit.y; yinit2 = yinit2';

options = bvpset('RelTol',1e-5,'Nmax',150000);

solinit.yinit=[];
sol = bvp4c(@(x,y) bvpfcn(x,y,eps,delta,c0,beta,g1,g2,fmag,alpha),@(ya,yb) bcfcn(ya,yb,ya1,ya3,yb4), solinit, options);
yout = sol.y';
tout = sol.x';
% 
solinit2 = bvpinit(xlin, @guess);
solinit2.x = ts2;
VV2 = [pts2(:,1),pts2(:,1)*0,pts2(:,2),-beta*((pts2(:,1).^3)/3-(g1+g2)*(pts2(:,1).^2)/2+g1*g2*pts2(:,1))]';

ya1 = VV2(1,1);
ya3 = VV2(3,1);

solinit2.y = VV2;
 yinit3 = solinit2.y; yinit3 = yinit3';
 solinit2.yinit=[];
% 
 sol2 = bvp4c(@(x,y) bvpfcn2(x,y,eps,delta,c0,beta,g1,g2,fmag,alpha),@(ya,yb) bcfcn2(ya,yb,ya1,ya3,yb4), solinit2, options);
 yout2 = sol2.y';
 tout2 = sol2.x';


Jaccc(:,i) = yout2(end,2:3)'-yout(end,2:3)';

end

ContJac = [(Jacwh(:,2)-Jacwh(:,1)) (Jaccc(:,2)-Jaccc(:,1))]/(2e-10);



c0 = c0lin;
yb4 = wh;

for i = 1:1



xlin = linspace(0,1.0,100);

 yinit2 = [xlin(1:end-1) 1; 0*xlin(1:end-1) 1; xlin(1:end-1) 1];


solinit = bvpinit(xlin, @guess);


solinit.x = -ts1;
VV = [pts1(:,1),pts1(:,1)*0,pts1(:,2),-beta*((pts1(:,1).^3)/3-(g1+g2)*(pts1(:,1).^2)/2+g1*g2*pts1(:,1))]';


ya1 = VV(1,1);
ya3 = VV(3,1);

solinit.y = VV;
yinit2 = solinit.y; yinit2 = yinit2';

options = bvpset('RelTol',1e-5,'Nmax',150000);

solinit.yinit=[];
sol = bvp4c(@(x,y) bvpfcn(x,y,eps,delta,c0,beta,g1,g2,fmag,alpha),@(ya,yb) bcfcn(ya,yb,ya1,ya3,yb4), solinit, options);
yout = sol.y';
tout = sol.x';
% 
solinit2 = bvpinit(xlin, @guess);
solinit2.x = ts2;
VV2 = [pts2(:,1),pts2(:,1)*0,pts2(:,2),-beta*((pts2(:,1).^3)/3-(g1+g2)*(pts2(:,1).^2)/2+g1*g2*pts2(:,1))]';

ya1 = VV2(1,1);
ya3 = VV2(3,1);

solinit2.y = VV2;
 yinit3 = solinit2.y; yinit3 = yinit3';
 solinit2.yinit=[];
% 
 sol2 = bvp4c(@(x,y) bvpfcn2(x,y,eps,delta,c0,beta,g1,g2,fmag,alpha),@(ya,yb) bcfcn2(ya,yb,ya1,ya3,yb4), solinit2, options);
 yout2 = sol2.y';
 tout2 = sol2.x';


Fout = yout2(end,2:3)'-yout(end,2:3)'

end


pnew = [wh;c0lin]-ContJac\Fout;

pdiff = pnew-[wh;c0lin];

wh = pnew(1);
c0lin = pnew(2);



end


%% Plotting


figure(1); clf; hold on;


[X,Y] = meshgrid(linspace(0.0,1,50),linspace(-0.5,0.5,50));
Z = -beta*((X.^3)/3-((g1+g2)/2)*X.^2 + g1*g2*X);
s = surf(X,Y,-Z,'FaceAlpha',0.5);



plot3(yout(:,1),-yout(:,3),-yout(:,4),'r-o')
plot3(yout2(:,1),-yout2(:,3),-yout2(:,4),'r-o')





%% BVP Functions

function dydx = bvpfcn(x,y,eps,delta,c0,beta,g1,g2,fmag,alpha)


dydx = zeros(4,1);
dydx = [(1/eps)*y(2)
       (1/eps)*(y(4)+beta*((y(1)^3)/3 - (g1 + g2)*(y(1)^2)/2+g1*g2*y(1))-delta*y(2))
       fmag*y(1)*(y(1)-alpha)*(1-y(1))
      y(3)+c0*y(1)];
    dydx = -dydx;
end

function dydx2 = bvpfcn2(x,y,eps,delta,c0,beta,g1,g2,fmag,alpha)

dydx2 = zeros(4,1);
dydx2 = [(1/eps)*y(2)
       (1/eps)*(y(4)+beta*((y(1)^3)/3 - (g1 + g2)*(y(1)^2)/2+g1*g2*y(1))-delta*y(2))
       fmag*y(1)*(y(1)-alpha)*(1-y(1))
      y(3)+c0*y(1)];
  dydx2 = dydx2;
end


function res = bcfcn(ya,yb,ya1,ya3,yb4)

res = [ya(1)-ya1
       ya(3)-ya3
       yb(1)-0.7
       yb(4)-yb4
];   
end

function res2 = bcfcn2(ya,yb,ya1,ya3,yb4)
 
res2 = [ya(1)-ya1
       ya(3)-ya3
       yb(1)-0.7
       yb(4)-yb4
]; 
end



function g = guess(x)
g = [0.3,x(2:end-1),1;
      0,0*x(2:end-1),1;
       0,x(2:end-1),1];
   
   g = [x(1:end);
   x(1:end);
   x(1:end)];

end



