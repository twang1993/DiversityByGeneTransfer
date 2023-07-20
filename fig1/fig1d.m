clear;
clc;
close all;
global mu10 mu20 mu1 mu2 gamma D eta kappa;
mu10=0.5;
mu20=0.5;

gammas=[0.8 0.9 0.99];
kappa=0.005;
D=0.2;
initial=[0.5 0.5 0 0];
timespan=0:0.1:200;

num=2000;
mu1s=rand(num,1);
mu2s=rand(num,1);

etas=0:0.1:0.5;
S1=0*ones(length(gammas),length(etas),num);
S2=0*ones(length(gammas),length(etas),num);

for i=1:length(gammas)
    i
    gamma=gammas(i);
    for j=1:length(etas)
        j
        eta=etas(j);
        for k=1:num
            mu1=mu1s(k);
            mu2=mu2s(k);
            [t,y]=ode45(@TwoSpecies,timespan,initial);
            S1(i,j,k)=y(end,1);
            S2(i,j,k)=y(end,2);
        end
    end
end

C=linspecer(length(gammas));
thresh=0.01;
prob=0*ones(length(gammas),length(etas));
for i=1:length(gammas)
    for j=1:length(etas)
        pin=0;
        for k=1:num
            if min(S1(i,j,k),S2(i,j,k))>thresh
                pin=pin+1;
            end
        end
        prob(i,j)=pin/num;
    end
end

for i=1:length(gammas)
    plot(etas,prob(i,:),'o-','color',C(i,:),'markersize',10,'linewidth',1.5);hold on;
end

set(gca,'fontsize',16);
xlabel('\eta','fontsize',20);
ylabel('coexistence feasibility','fontsize',20);
axis([0 max(etas) 0 0.62]);
set(gcf,'position',[100 100 270 270]);
saveas(gcf,'CoexistenceProb.fig');
saveas(gcf,'CoexistenceProb.eps');
save('CoexistenceProb.mat');

function dydt=TwoSpecies(t,y)
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