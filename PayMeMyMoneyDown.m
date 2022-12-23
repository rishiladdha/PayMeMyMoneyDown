%please read accompanying PDFs to understand what's going on.

%%

%%MONTE CARLO BASICS

%%

%Magician Works Three Hours a day

% Monte Carlo Simulation of Compound Poisson process
clc
clear all
format long
 
time = 3;

lambda = 5;
N = 10000; % Number of Monte Carlo Steps
beer_price = 350;
beer = zeros(N,1);
for i=1:N
    n = poissrnd(lambda*time); %total no of people who are coming 
    if n > 0
        coins = zeros(n,1);
        for j=1:n
            U = rand(1);
            if U<=0.4
                coins(j) = 5;
            elseif U>=0.8
                coins(j) = 20;
            else
                coins(j) = 10;
            end
        end
    end
    if sum(coins)>=beer_price
        beer(i) = 1;
    else
        beer(i) = 0;
    end
end
l_hat = mean(beer)% l_hat = P(X3 >= beer_price)

%% Plotting probability of purchasing beer vs number of hours worked

% Monte Carlo Simulation of Compound Poisson process
clc
clear all
format long
 
arr=[];
time=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];

for t = 1:15
lambda = 5;
N = 10^3; % Number of Monte Carlo Steps
beer_price = 350;
beer = zeros(N,1);
for i=1:N
    n = poissrnd(lambda*t); %total no of people who are coming 
    if n > 0
        coins = zeros(n,1);
        for j=1:n
            U = rand(1);
            if U<=0.4
                coins(j) = 5;
            elseif U>=0.8
                coins(j) = 20;
            else
                coins(j) = 10;
            end
        end
    end
    if sum(coins)>=beer_price
        beer(i) = 1;
    else
        beer(i) = 0;
    end
end
l_hat = mean(beer)% l_hat = P(X3 >= beer_price)
arr = [arr l_hat]

end
% relErr_hat = std(beer) / (l_hat * sqrt(N)); % relative error of l_hat by crude

disp arr_plot;
plot(arr,time);
xlabel("No of hours");
ylabel("probability of getting beer");

%% Insurance Problem 

%%

%%Crude Monte Carlo Simulation for P(y_2)>5
lambda=3;
N = 10^6;
quartertotalclaim = 5;
claimprob = zeros(N,1);
for i = 1:N
    n=poissrnd(lambda);
    claimconcat = zeros(n,1);
    if n>=0
        for j = 1:n
            U = rand(1);
            if U < 2/3
                claimconcat(j) = 1;
            else
                claimconcat(j) = 2;
            end
        end
    end

    if sum(claimconcat)>quartertotalclaim
        claimprob(i)=1;
    else
        claimprob(i)=0;
    end

end

ans = mean(claimprob)
relErr_hat = std(claimprob) / (ans*sqrt(N))
%%

%Crude Monte Carlo Simulation for P(Y_3 > 5)

lambda=1;
N = 10^6;
quartertotalclaim = 5;
claimprob = zeros(N,1);
for i = 1:N
    n=poissrnd(lambda);
    claimconcat = zeros(n,1);
    if n>0
        
        for j = 1:n
            U = rand(1);
            if U < 2/3
                claimconcat(j) = 1;
            %elseif U >= 0.4 && U < 0.8
                %claimconcat(j) = 10;
            else
                claimconcat(j) = 2;
            end
        end
    end

    if sum(claimconcat)>quartertotalclaim
        claimprob(i)=1;
    else
        claimprob(i)=0;
    end

end

ans = mean(claimprob) 
relErr_hat = std(claimprob) / (ans*sqrt(N))

%% Risk analysis
clc
clear all
lambda=[2 3 1 3];
N = 10^4;

claimprob = [];
for i = 1:N
    total=0;
    for x=1:4
        claimconcat = [];
        temp=lambda(x);
        n=poissrnd(temp);
        if n>0
            for j = 1:n
                U = rand(1);
                if U < 2/3
                    claimconcat(j) = 1;
                else
                    claimconcat(j) = 2;
                end
            end
        end
        total=total+sum(claimconcat);
    end

    if total>10
        claimprob(i)=1;
    end

end

ans = mean(claimprob)
