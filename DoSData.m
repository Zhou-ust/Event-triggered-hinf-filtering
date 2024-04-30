% Generate the DoS attack data governed by Bernoulli distribution
% The data are generated when Matlab runs for the first time.

clear
close all
clc

alphabar=0.5;% The probability of DoS attack is alphabar
Bernoullirnd=rand(1,200);
alpha=Bernoullirnd<alphabar;% alpha is governed by Bernoulli distribution
save('DoSData.mat','alpha');