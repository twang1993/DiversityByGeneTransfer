clear;
clc;
close all;
global mu10 mu20 lambda1 lambda2 alpha D eta kappa;
mu10=0.5;
mu20=0.5;
lambda1=-0.2;
lambda2=0.2;
alpha=0.8;
kappa=0.05;%0.01;
D=0.2;%0.05;
initial=[0.5 0 0.5 0];
timespan=0:0.1:500;%800;

figure(1);
eta=0;
[t,y]=ode45(@TwoSpecies,timespan,initial);
patch('XData',[t' fliplr(t')],'YData',[0*t' fliplr(transpose(y(:,1)./(y(:,1)+y(:,3))))],'FaceColor',[217, 109, 58]/256,'FaceAlpha',1,'EdgeColor','none');hold on;
patch('XData',[t' fliplr(t')],'YData',[transpose(y(:,1)./(y(:,1)+y(:,3))) fliplr(ones(1,length(t)))],'FaceColor',[173, 201, 59]/256,'FaceAlpha',1,'EdgeColor','none');hold on;
axis([0 max(t) 0 1]);
set(gca,'fontsize',16);
xlabel('time','fontsize',20);
ylabel('fraction','fontsize',20);
set(gcf,'position',[100 100 300 200]);
saveas(gcf,'TwoSpecies_1.fig');
saveas(gcf,'TwoSpecies_1.pdf');

figure(2);
eta=0.5;
[t,y]=ode45(@TwoSpecies,timespan,initial);
patch('XData',[t' fliplr(t')],'YData',[0*t' fliplr(transpose(y(:,1)./(y(:,1)+y(:,3))))],'FaceColor',[217, 109, 58]/256,'FaceAlpha',1,'EdgeColor','none');hold on;
patch('XData',[t' fliplr(t')],'YData',[transpose(y(:,1)./(y(:,1)+y(:,3))) fliplr(ones(1,length(t)))],'FaceColor',[173, 201, 59]/256,'FaceAlpha',1,'EdgeColor','none');hold on;
set(gca,'fontsize',16);
axis([0 max(t) 0 1]);
xlabel('time','fontsize',20);
ylabel('fraction','fontsize',20);
set(gcf,'position',[100 100 300 200]);
saveas(gcf,'TwoSpecies_2.fig');
saveas(gcf,'TwoSpecies_2.pdf');



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