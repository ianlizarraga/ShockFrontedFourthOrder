% xmesh = linspace(0,pi/2,20);

clear all;

load('reduceddata3.mat');


xlin = linspace(0,1.0,100);

% xlin = ts2Dup;

% (21.0/8.0)*x[0]-4.0*x[0]*x[0]+2.0*x[0]*x[0]*x[0]

 yinit2 = [xlin(1:end-1) 1; 0*xlin(1:end-1) 1; xlin(1:end-1) 1];


solinit = bvpinit(xlin, @guess);


% solinit.x = ts2;
% solinit.y = [pts2(:,1),pts2(:,1)*0,pts2(:,2),(21.0/8.0)*pts2(:,1)-4.0*pts2(:,1).^2+2.0*pts2(:,1).^3]';
% yinit2 = solinit.y; yinit2 = yinit2';

solinit.x = -ts;
 VV = [pts(:,1),pts(:,1)*0,pts(:,2),(21.0/8.0)*pts(:,1)-4.0*pts(:,1).^2+2.0*pts(:,1).^3,upout(:,1)*0,Duout.*pts(:,3)./pts(:,4),Duout]';
%  VV = [pts(:,1),pts(:,1)*0,pts(:,2),(21.0/8.0)*pts(:,1)-4.0*pts(:,1).^2+2.0*pts(:,1).^3,-sqrt(Duout),0*pts(:,3),0*pts(:,3)]';
solinit.y = VV;
yinit2 = solinit.y; yinit2 = yinit2';

options = bvpset('RelTol',1e-8,'Nmax',15000);

solinit.yinit=[];
sol = bvp4c(@bvpfcn, @bcfcn, solinit, options);
yout = sol.y';
tout = sol.x';

 solinit2 = bvpinit(xlin, @guess);
 solinit2.x = ts2;
 VV = [pts2(:,1),pts2(:,1)*0,pts2(:,2),(21.0/8.0)*pts2(:,1)-4.0*pts2(:,1).^2+2.0*pts2(:,1).^3,Duoutu*0,Duoutu.*pts2(:,3)./pts2(:,4),Duoutu]';
 solinit2.y = VV;
 yinit3 = solinit2.y; yinit3 = yinit3';
 solinit2.yinit=[];
% 
  sol2 = bvp4c(@bvpfcn2, @bcfcn2, solinit2, options);
 yout2 = sol2.y';
tout2 = sol2.x';

%  yout2(end,4)-yout(end,4)
 
figure(1); hold on;

% plot3(yinit2(:,1),yinit2(:,3),yinit2(:,4),'b-o')
%  plot3(yinit3(:,1),yinit3(:,3),yinit3(:,4),'b-o')
% 
%  figure(1); hold on;
% plot3(yout(:,1), yout(:,3),yout(:,4), 'r-o') 
% plot3(yout2(:,1), yout2(:,3),yout2(:,4), 'r-o') 
% 
% 
% [X,Y] = meshgrid(linspace(0.0,1,50),linspace(-0.5,0.5,50));
% Z = 2.0*X.^3-4.0*X.^2+(21.0/8.0)*X;
% s = surf(X,Y,Z,'FaceAlpha',0.5);
% 

 plot(yout(:,1), yout(:,6)./yout(:,7), 'r-')
plot(yout2(:,1), yout2(:,6)./yout2(:,7), 'r-')


% figure(1); hold on
% 
% plot3(yinit(1,:)',yinit(2,:)',yinit(3,:)','b-o')
% 
% 
% 
% plot3(yout(:,1), yout(:,2),yout(:,3), 'r-o')
% 



function dydx = bvpfcn(x,y)
eps = 0.0001;
% c0 = 0.1968109995;
c0 = 0.19686; %really good agreement for eps = 1e-4;
% lambda = -0.25;
   lambda = 100.0;
%   lambda = 0.00637;

dydx = zeros(7,1);
dydx = [(1/eps)*y(2)
       (1/eps)*((21.0/8.0)*y(1)-4.0*y(1)^2+2.0*y(1)^3-y(4))
       -5.0*y(1)*(1.0-y(1))*(y(1)-0.2)
       y(3)-c0*y(1)
       (1/eps)*(6.0*(y(1)-7.0/12.0)*(y(1)-0.75)-y(7)-y(5)^2)
       lambda-(-15.0*y(1)^2+12.0*y(1)-1.0)-(1/eps)*y(5)*y(6)
       y(6)-c0-(1/eps)*y(5)*y(7)];
   dydx = -dydx;
end

function dydx2 = bvpfcn2(x,y)
eps = 0.0001;
% c0 = 0.1968109995;
c0 = 0.19686; %really good agreement for eps = 1e-4;
% lambda = -0.25;
  lambda = 100.0;
%   lambda = 0.00637;
dydx2 = zeros(7,1);
dydx2 = [(1/eps)*y(2)
       (1/eps)*((21.0/8.0)*y(1)-4.0*y(1)^2+2.0*y(1)^3-y(4))
       -5.0*y(1)*(1.0-y(1))*(y(1)-0.2)
       y(3)-c0*y(1)
       (1/eps)*(6.0*(y(1)-7.0/12.0)*(y(1)-0.75)-y(7)-y(5)^2)
       lambda-(-15.0*y(1)^2+12.0*y(1)-1.0)-(1/eps)*y(5)*y(6)
       y(6)-c0-(1/eps)*y(5)*y(7)];
end


function res = bcfcn(ya,yb)
res = [ya(1)-0.3128739859*(1e-4)
       ya(3)+0.4770594156*(1e-4)
       yb(1)-0.6666666666666
       yb(3)-0.02243179      
       ya(6)+1.524708677764532
       ya(7)-2.624902415909677
       yb(5)-0.0002]; 
% res = [ya(1)-0.3128739859*(1e-10)
%        ya(3)+0.4770594156*(1e-10)
%        yb(1)-0.6666666666666
%        yb(3)-0.02243179  
%        ya(6)-0.0
%        ya(7)-0.0
%        yb(7)-0.01]; 
end

function res2 = bcfcn2(ya,yb)
% res2 = [ya(1)-0.999
%        ya(3)-0.1970
%        yb(1)-0.7
%        yb(3)-0.022];  
res2 = [ya(1)-0.999997028414785
       ya(3)-0.196805999500000
       yb(1)-0.6666666666666
       yb(3)-0.02243179
       ya(6)-1.682458654596142
       ya(7)-0.624946157771746
       yb(5)-0.0009]; 
end



function g = guess(x)
g = [0.3,x(2:end-1),1;
      0,0*x(2:end-1),1;
       0,x(2:end-1),1];
   
   g = [x(1:end);
   x(1:end);
   x(1:end)];
   
%    size(g)
%  p = [1;1;1];
%  g = [g,p];
end



