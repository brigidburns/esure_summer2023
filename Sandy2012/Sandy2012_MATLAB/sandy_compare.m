clear all
clc
close all

load('sandy2012_variables.mat'); 
load('sandy2012_processed_variables.mat');  

figure(1)
hold on
plot(sondexvec,sondeyvec,'rx')
plot(sondexvec_processed,sondeyvec_processed,'bo')
xlabel('x')
ylabel('y')
axis([-max_rad max_rad -max_rad max_rad])