clear;
clc;
close all;
global NumSpecies mu lambda gamma D eta kappa;
mu0=0.5;
NumSpecies=20;
etas=[0.05 0.1 0.2];

gamma=0.9;
kappa=0.005;
D=0.2;

initial=0*ones(NumSpecies^2+NumSpecies,1);
for i=1:NumSpecies
    initial(i)=1/NumSpecies;
    initial(NumSpecies+(i-1)*NumSpecies+i)=initial(i);
end
timespan=0:0.1:200;

CC=colormap('lines');
repeat=1;
for jkl=1:length(etas)
    eta=etas(jkl);
    lambda=(-0.2+0.4*rand(NumSpecies,1));
    mu=mu0*(1+lambda);
    [t,y]=ode45(@MultiSpecies,timespan,initial);
    
    ss=0*ones(length(t),NumSpecies);
    for i=1:NumSpecies
        ss(:,i)=mu(i);
        for j=1:NumSpecies
            if i~=j
                ss(:,i)=ss(:,i).*(1+lambda(j)*y(:,NumSpecies+(i-1)*NumSpecies+j)./y(:,i));
            end
        end
    end
    yy=std(ss,0,2);
    subplot(3,1,jkl);
    patch('XData',[t' fliplr(t')],'YData',[0*yy' fliplr(yy')],'FaceColor',[46,133,198]/256,'FaceAlpha',0.4,'LineStyle','none');hold on;
    patch('XData',[t' fliplr(t')],'YData',[0*yy' fliplr(-yy')],'FaceColor',[46,133,198]/256,'FaceAlpha',0.4,'LineStyle','none');hold on;
    plot(t,0*t,'color',[46,133,198]/256);hold on;
    % plot(t,yy,'color',[46,133,198]/256);hold on;
    % plot(t,-yy,'color',[46,133,198]/256);hold on;
    set(gca,'fontsize',10);
    box on;
    H=gca;
    H.LineWidth=1;
end
subplot(3,1,3);
% xlabel('time','fontsize',16);
% ylabel('std(\mu)','fontsize',16);
set(gcf,'position',[100 100 200 500]);
% axis([0 200 0 0.12]);
saveas(gcf,'MuChange.fig');
saveas(gcf,'MuChange.pdf');
saveas(gcf,'MuChange.eps');

function dydt=MultiSpecies(t,y)
    global NumSpecies mu lambda gamma D eta kappa;
    dydt(NumSpecies*(1+NumSpecies),1)=0;
    thresh=0;
    for i=1:NumSpecies
            ss=1;
            sumy=0;
            for j=1:NumSpecies
                if i~=j
                    ss=ss*(1+lambda(j)*y(NumSpecies+(i-1)*NumSpecies+j)/y(i));
                end
                sumy=sumy+y(j);
            end
            dydt(i,1)=mu(i)*y(i)*ss*(1-(gamma*sumy-gamma*y(i)+y(i)))-D*y(i);
        for j=1:NumSpecies
                if j==i
                    dydt(NumSpecies+(i-1)*NumSpecies+j,1)=dydt(i,1);
                else
                    ss=1;
                    for k=1:NumSpecies
                        if k~=i&&k~=j
                            ss=ss*(1+lambda(k)*y(NumSpecies+(i-1)*NumSpecies+k)/y(i));
                        end
                    end
                    donor=0;
                    for k=1:NumSpecies
                        donor=donor+y(NumSpecies+(k-1)*NumSpecies+j);
                    end
                    dydt(NumSpecies+(i-1)*NumSpecies+j,1)=mu(i)*y(NumSpecies+(i-1)*NumSpecies+j)*ss*(1+lambda(j))*(1-(gamma*sumy-gamma*y(i)+y(i)))+eta*(y(i)-y(NumSpecies+(i-1)*NumSpecies+j))*donor-(kappa+D)*y(NumSpecies+(i-1)*NumSpecies+j);
                end
        end
    end

end