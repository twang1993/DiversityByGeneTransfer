clear;
clc;
close all;
global NumSpecies mu lambda gamma D eta kappa;
mu0=0.5;
Ns=[2:2:26];
etas=[0 0.1 0.2];
repeat=500;
diversity=0*ones(length(Ns),length(etas),repeat);
richness=0*ones(length(Ns),length(etas),repeat);
CoexistProb=0*ones(length(Ns),length(etas));
thresh=0.01;
for uio=1:length(Ns)
    NumSpecies=Ns(uio);
    gamma=0.9;
    kappa=0.005; 
    D=0.2;

    initial=0*ones(NumSpecies^2+NumSpecies,1);
    for i=1:NumSpecies
        initial(i)=1/NumSpecies;
        initial(NumSpecies+(i-1)*NumSpecies+i)=initial(i);
    end
    timespan=0:0.1:200;

    for i=1:length(etas)
        eta=etas(i);
        for hjk=1:repeat
            (uio-1)*length(etas)*repeat+(i-1)*repeat+hjk
            lambda=(-0.2+0.4*rand(NumSpecies,1));
            mu=mu0*(1+lambda);
            [t,y]=ode45(@MultiSpecies,timespan,initial);
            temp=y(end,1:NumSpecies);
            temp=temp/sum(temp);
            diversity(uio,i,hjk)=exp(-sum(temp.*log(temp)));
            richness(uio,i,hjk)=sum(y(end,1:NumSpecies)>thresh);
        end
        CoexistProb(uio,i)=sum(richness(uio,i,:)==NumSpecies)/repeat;
    end
end
figure(1);
C=linspecer(length(etas));
for i=1:length(etas)
plot(Ns,CoexistProb(:,i),'o-','markersize',10,'color',C(i,:),'linewidth',1.5);hold on;
end
set(gca,'fontsize',16);
xlabel('species number','fontsize',20);
ylabel('coexistence probability','fontsize',20);
set(gcf,'position',[100 100 300 300]);
axis([1 max(Ns) 0 1]);
saveas(gcf,'MultiSpecies_1.fig');
saveas(gcf,'MultiSpecies_1.pdf');

figure(2);
C=linspecer(length(Ns));
for i=1:length(Ns)
plot(etas,CoexistProb(i,:),'o-','markersize',10,'color',C(i,:),'linewidth',1.5);hold on;
end
set(gca,'fontsize',16);
xlabel('transfer rate','fontsize',20);
ylabel('coexistence probability','fontsize',20);
set(gcf,'position',[100 100 300 300]);
axis([0 max(etas) 0 1]);
saveas(gcf,'MultiSpecies_2.fig');
saveas(gcf,'MultiSpecies_2.pdf');
save('MultiSpecies.mat');


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