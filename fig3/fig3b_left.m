clear;
clc;
global NumSpecies mu alpha D;
C=linspecer(3);
NumSpecies=20;
mu0=0.5;
D=0.2;
alpha=0.9;
for repeat=1:5000
    repeat
    range=rand;
    mu=mu0*(1-range/2+range*rand(NumSpecies,1));

    initial(NumSpecies)=0;
    for i=1:NumSpecies
        initial(i)=1/NumSpecies;
    end
    timespan=0:300;
    [t,y]=ode45(@MultiSpecies,timespan,initial);
    temp=y(end,1:NumSpecies);
    temp=temp(temp>0);
    temp=temp/sum(temp);
    diversity(repeat)=exp(-sum(temp.*log(temp)));
    SS(repeat)=std(mu);
end

save('MuStdDiversity.mat');

function dydt=MultiSpecies(t,y)
    global NumSpecies mu alpha D;
    dydt(NumSpecies,1)=0;
    for i=1:NumSpecies
            sumy=0;
            for j=1:NumSpecies
                sumy=sumy+y(j);
            end
            dydt(i,1)=mu(i)*y(i)*(1-(alpha*sumy-alpha*y(i)+y(i)))-D*y(i);
    end
end