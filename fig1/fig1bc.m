clear;
clc;
close all;
global mu10 mu20 mu1 mu2 gamma D eta kappa;

gamma=0.99;
D=0.2;
eta=0.1;
kappa=0.005;

mu1s=0:0.005:1;
mu2s=0:0.005:1;
mu10=0.5;
mu20=0.5;
initial=[0.5 0.5];
initialHGT=[0.5 0.5 0 0];
timespan=0:0.1:80;
S1Abun=0*ones(length(mu1s),length(mu2s));
S2Abun=0*ones(length(mu1s),length(mu2s));
S1AbunHGT=0*ones(length(mu1s),length(mu2s));
S2AbunHGT=0*ones(length(mu1s),length(mu2s));

thresh=0.02;
for i=1:length(mu1s)
    i
    mu1=mu1s(i);
    for j=1:length(mu2s)
        mu2=mu2s(j);
        [t1 y1]=ode45(@TwoSpecies,timespan,initial);
        [t2 y2]=ode45(@TwoSpeciesHGT,timespan,initialHGT);
        S1Abun(i,j)=y1(end,1);
        S2Abun(i,j)=y1(end,2);
        S1AbunHGT(i,j)=y2(end,1);
        S2AbunHGT(i,j)=y2(end,2);
    end
end

tt=15;
figure(1);
XX=S1Abun;
YY=S2Abun;
for i=1:length(mu1s)
    for j=1:length(mu2s)
        if min(XX(i,j),YY(i,j))>thresh
            plot(mu1s(i),mu2s(j),'.','color',[185,221,245]/256,'markersize',tt);hold on;
        end
    end
end
set(gca,'fontsize',16);
xlabel('\mu_1','fontsize',20);
ylabel('\mu_2','fontsize',20);
set(gcf,'position',[100 100 300 300]);
saveas(gcf,'PhaseDiagram_1.fig');
saveas(gcf,'PhaseDiagram_1.pdf');


figure(2);
XX=S1AbunHGT;
YY=S2AbunHGT;
for i=1:length(mu1s)
    for j=1:length(mu2s)
        if min(XX(i,j),YY(i,j))>thresh
            plot(mu1s(i),mu2s(j),'.','color',[185,221,245]/256,'markersize',tt);hold on;
        end
    end
end
set(gca,'fontsize',16);
xlabel('\mu_1','fontsize',20);
ylabel('\mu_2','fontsize',20);
set(gcf,'position',[100 100 300 300]);
save('PhaseDiagram.mat');
saveas(gcf,'PhaseDiagram_2.fig');
saveas(gcf,'PhaseDiagram_2.pdf');

    
    
function dydt=TwoSpecies(t,y)
global mu1 mu2 gamma D ;
s1=y(1);
s2=y(2);
dydt=[mu1*s1*(1-s1-gamma*s2)-D*s1;
    mu2*s2*(1-s2-gamma*s1)-D*s2];
end

function dydt=TwoSpeciesHGT(t,y)
global mu10 mu20 mu1 mu2 gamma D eta kappa;
s1=y(1);
s2=y(2);
p1=y(3);
p2=y(4);
lambda1=mu1/mu10-1;
lambda2=mu2/mu20-1;

dydt=[mu1/s1*(s1+lambda2*p1)*s1*(1-s1-gamma*s2)-D*s1;
    mu2/s2*(s2+lambda1*p2)*s2*(1-gamma*s1-s2)-D*s2;
    mu1*(1+lambda2)*p1*(1-s1-gamma*s2)+eta*(s2+p1)*(s1-p1)-kappa*p1-D*p1;
    mu2*(1+lambda1)*p2*(1-gamma*s1-s2)+eta*(s1+p2)*(s2-p2)-kappa*p2-D*p2];
end