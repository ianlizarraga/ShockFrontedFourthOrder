options = dopset('AbsTol',1e-16,'RelTol',1e-16,'EventTol',1e-16,'Events',2,'MaxStepSize',1e-2,'MaxSteps',1e7);

timet = 1e3;

eps = 0.01;
c0 = 0.0;
% beta = 10;
% g1 = 0.3;
% g2 = 0.75;

beta = 6;
g1 = 7/12;
g2 = 3/4;
fmag = 5.0;
alpha = 0.2;


% delta = 0.3;
% % alpha = (g1+g2)/2;
% % alpha = 0.4;
% wh = -0.216562499999982;
crosssec=1.5;

 load('delwhcparams2.mat')

 load('delparams3.mat')
 
deltest = 0.24;

deltest = 0.0;

[minValue,closestIndex] = min(abs(whdely2(:,1)-deltest))


delta = whdely2(closestIndex,1);
wh = whdely2(closestIndex,2);

c0lin = linspace(0.18,0.25,10);

c0lin = [c0lin, 0.198166808404494];

c0lin = 0.198166808404494;

c0lin = whdely2(closestIndex,3);


% c0lin = 0.1;

% c0lin = [c0lin,0.006107470317511];

% c0lin = whdely(closestIndex,3);

misma =c0lin*0;

for i = 1:length(c0lin)
    
    c0 = c0lin(i);
par = [eps,c0,beta,g1,g2,delta,alpha,wh,crosssec,fmag];


Du0 = beta*g1*g2;
Du1 = beta*(1-g1)*(1-g2);
fp0 = -alpha;
fp1 = -1+alpha;

Jac0 = [-c0,-1; Du0*fp0,0];
Jac1 = [-c0,-1; Du1*fp1,0];

[V0,D0] = eigs(Jac0);
[V1,D1] = eigs(Jac1);


[~,m0]=sort(real(diag(D0)));
[~,m1]=sort(real(diag(D1)),'descend');

vstab = V0(:,m0(1));
vstab = vstab*sign(vstab(1))*(1e-8);

vunst = V1(:,m1(1));
vunst = vunst*(-sign(vunst(1)))*(1e-8);


pinit1 = [0;0]+vstab;
pinit2 = [1;-c0]+vunst;


% wh = -0.216562499999982;

Phiminus = 0.075;
Phiplus = 0.975;

 equhat = sort(roots([beta/3 -(g1+g2)*beta/2 beta*g1*g2 wh]));
 
 inflpt = (g1+g2)/(2);
 Phiinfl = beta*((inflpt^3)/3-(g1+g2)*(inflpt^2)/2+g1*g2*inflpt);
 equalrule = sort(roots([beta/3 -(g1+g2)*beta/2 beta*g1*g2 -Phiinfl]));
 
 
 
 Phiplus = beta*(g2^3/3-(g1+g2)*(g2^2)/2+g1*g2*g2);
 
par(9) = equhat(1);
% par(9) = equalrule(1);

% par(9) = g1;

[ts1,pts1,te1,ye1,ie1,stats1] = dop853('interpred2',[0 -timet],pinit1,options,par); 



par(9) = equhat(3);
% par(9) = equalrule(3);

% par(9) = g2;

[ts2,pts2,te2,ye2,ie2,stats2] = dop853('interpred2',[0 timet],pinit2,options,par);

misma(i) = ye1(2)-ye2(2);

figure(1); hold on

plot(pts1(:,1),pts1(:,2),'b')
plot(pts2(:,1),pts2(:,2),'r')

end



save('reducedwave','pts1','pts2','ts1','ts2','c0','beta','g1','g2','fmag','alpha','delta','wh','c0lin')