clear;
clc;
close all;
global mu10 mu20 lambda1 lambda2 alpha D eta kappa;
mu10=0.5;
mu20=0.5;
mu1s=0:0.1:1;
mu2s=0:0.1:1;
lambda1s=mu1s/mu10-1;
lambda2s=mu2s/mu20-1;
alpha=0.99;
kappa=0.05;
D=0.2;
initial=[0.5 0 0.5 0];
timespan=0:0.1:20;

etas=[0.05 0.2 0.5];

for jkl=1:length(etas)
    eta=etas(jkl);
    openfig('PhaseDiagram_1.fig');
    for i=1:length(mu1s)
        i
        lambda1=lambda1s(i);
        for j=1:length(mu2s)
            lambda2=lambda2s(j);
            [t,y]=ode45(@TwoSpecies,timespan,initial);
            mu1e=mu10*(1+lambda1)./y(end,1).*(y(end,1)+lambda2*y(end,2));
            mu2e=mu20*(1+lambda2)./y(end,3).*(y(end,3)+lambda1*y(end,4));
            quiver(mu1s(i),mu2s(j),mu1e-mu1s(i),mu2e-mu2s(j),'ShowArrowHead','on','MaxHeadSize',0.5);hold on;
        end
    end
    axis([0 1 0 1]);
    saveas(gcf,sprintf('ArrowPhase_%d.fig',jkl));
    saveas(gcf,sprintf('ArrowPhase_%d.pdf',jkl));
end

function dydt=TwoSpecies(t,y)
global mu10 mu20 lambda1 lambda2 alpha D eta kappa;
s1=y(1);
p1=y(2);
s2=y(3);
p2=y(4);
mu1=mu10*(1+lambda1);
mu2=mu20*(1+lambda2);

dydt=[mu1/s1*(s1+lambda2*p1)*s1*(1-s1-alpha*s2)-D*s1;
    mu1*(1+lambda2)*p1*(1-s1-alpha*s2)+eta*(s2+p1)*(s1-p1)-kappa*p1-D*p1;
    mu2/s2*(s2+lambda1*p1)*s2*(1-alpha*s1-s2)-D*s2;
    mu2*(1+lambda1)*p2*(1-alpha*s1-s2)+eta*(s1+p2)*(s2-p2)-kappa*p2-D*p2;];
end