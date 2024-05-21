clc
clear all
close all

layers = 17;%set amount of layers to generate in step below

A = TriHexGen(layers);%generate graph object containing lattice of triangles which make up hexagons. neighbors() can be used to find nearest neibors for nodes

A = ConnectOppositeEdges(A);%connect opposit edges of lattics to deal with boundary issues

sigma = initialState(A);%generate initial state of all points in lattice

stepmax = 500;
stepcount = 2;

tol = 10e-6;

Ehist(1) = inf;

T_init = 2;
T = T_init;
T_final = 3;
T_step = 0.03;
candMult = 3/5;% 60% of points are candidates
Tvals = [T_init:T_step:T_final];

J = 1;
H = 0;

E_temp = zeros(size(Tvals));

    sigmasum = 0;
    sample = sigma;%randperm(length(sigma),1000);%select sample of 1000 points
for(i=1:length(sample))
    sigmasum = sigmasum + sigma(i)*sum(sigma(neighbors(A,i)));
end
E = ((-J*sigmasum)-(H*sum(sigma)));
flip_percents = zeros(1,stepmax);

for(stepcount=1:stepmax)%allow for an initial period then check for settling
    flip_count = 0;
    sigma_rand = randperm(length(sigma),floor(length(sigma)*candMult));%select random set of spins to check for possible flip
    for(i=1:length(sigma_rand))
        
        sigma_new = P_flip(A,sigma,T,sigma_rand(i));%mutate sigma according to probability equation, wait until every flip has been decided before changing the actual sigma values.
        if(sigma(sigma_rand(i))~=sigma_new(sigma_rand(i)))
            flip_count = flip_count + 1;
        end
    end
    %flip_percents(stepcount) = 100*(flip_count/length(sigma));

    sigma = sigma_new;
    
    sigmasum = 0;
    sample = sigma;%randperm(length(sigma),floor(length(sigma)/2));%select sample of 1000 points
    for(i=1:length(sample))
        sigmasum = sigmasum + sigma(i)*sum(sigma(neighbors(A,i)));
    end
    
    Ehist(stepcount)=E;
    E = ((-J*sigmasum)-(H*sum(sigma)))/length(sample);
end

for(T=Tvals)
    break
    for(iter=1:10)
        flip_count = 0;
        sigma_rand = randperm(length(sigma),floor(length(sigma)*candMult));%select random set of spins to check for possible flip
    for(i=1:length(sigma_rand))
        sigma_new = P_flip(A,sigma,T,sigma_rand(i));%mutate sigma according to probability equation, wait until every flip has been decided before changing the actual sigma values.
        if(sigma(sigma_rand(i))~=sigma_new(sigma_rand(i)))
            flip_count = flip_count + 1;
        end
    end
    flip_percents(stepcount) = 100*(flip_count/length(sigma));
    sigma = sigma_new;
    sigmasum = 0;
    sample = randperm(length(sigma),1000);%select sample of 1000 points
    for(i=1:length(sample))
        sigmasum = sigmasum + sigma(sample(i))*sum(sigma(neighbors(A,sample(i))));
    end
    

    E(iter) = ((-J*sigmasum)-(H*sum(sample)));
    M(iter)=abs(sum(sigma));
    CT = E(iter)^2;
    XT = M(iter)^2;
    
    end
    E_temp(Tvals==T)=mean(E)/A.numnodes;
    M_temp(Tvals==T)=mean(M)/A.numnodes;
    CT_temp(Tvals==T)=(mean(CT)-(mean(E)^2))/(A.numnodes*(T^2));
    XT_temp(Tvals==T)=(mean(XT)-(mean(M)^2))/(A.numnodes*(T^2));
end

T_end = stepcount;

T = T_init;%reset T for next loop

for(k=1:stepmax*0)%allow for an initial period then check for settling
    stepcount = stepcount+1;
    flip_count = 0;
    sigma_rand = randperm(length(sigma),floor(length(sigma)*candMult));%select random set of spins to check for possible flip
    for(i=1:length(sigma_rand))
        
        sigma_new = P_flip(A,sigma,T,sigma_rand(i));%mutate sigma according to probability equation, wait until every flip has been decided before changing the actual sigma values.
        if(sigma(sigma_rand(i))~=sigma_new(sigma_rand(i)))
            flip_count = flip_count + 1;
        end
    end
    flip_percents(stepcount) = 100*(flip_count/length(sigma));

    sigma = sigma_new;
    
    sigmasum = 0;
    for(i=1:length(sigma))
        sigmasum = sigmasum + sigma(i)*sum(sigma(neighbors(A,i)));
    end
    
    Ehist(stepcount) = E;
    E = ((-J*sigmasum)-(H*sum(sigma)));
end

figure
plot(2:length(Ehist),Ehist(2:end),'bo-','MarkerSize',2)
axis([0,stepcount,-0.4,0.4])

xlabel("itterations")
ylabel("average energy per node")

%for(i=2:length(Tvals)-1)
%    E_avg_temp(i-1) = mean(E_temp(i-1:i+1));%take a moving average
%end

%figure
%plot(Tvals(1:length(Tvals)),E_temp,'g-',Tvals(1:length(Tvals)),E_temp,'bo')
%xlabel("Temperature")
%ylabel("average sample energy per node")
%figure
%plot(Tvals(1:length(Tvals)),M_temp,'g-',Tvals(1:length(Tvals)),M_temp,'bo')
%xlabel("Temperature")
%ylabel("sample magnetization per node")
%figure
%plot(Tvals(1:length(Tvals)),CT_temp,'g-',Tvals(1:length(Tvals)),CT_temp,'bo')
%xlabel("Temperature")
%ylabel("sample heat capacity per node")
%figure
%plot(Tvals(1:length(Tvals)),XT_temp,'g-',Tvals(1:length(Tvals)),XT_temp,'bo')
%xlabel("Temperature")
%ylabel("sample susceptibility per node")

%axis([T_init,T_final,-0.2,0.2])

function sigma_new = P_flip(A,sigma,T,i)
sigma_new = sigma;

deltaEk = -2*sigma(i)*(sum(sigma(neighbors(A,i))));%calculates Esigmak-Esigmahat

if(deltaEk<0)
    P = 1;%flip
else
    P = exp(deltaEk/T);%maybe flip
end

if P>=rand%check probability condition
    sigma_new(i)=-sigma(i);%if condition is met then flip the sigma value at i, if not then the flip is invalid. either way return new sigma matrix after completion
end

end