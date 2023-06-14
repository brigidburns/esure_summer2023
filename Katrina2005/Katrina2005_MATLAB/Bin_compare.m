clear all
clc
close all

load('katrina2005_variables.mat'); 
load('katrina2005_processed_variables.mat'); 
load('katrina2005_dc_variables.mat'); 

figure(1)
hold on
plot(sondexvec,sondeyvec,'rx')
plot(sondexvec_processed,sondeyvec_processed,'bx')
plot(sondexvec_dc,sondeyvec_dc,'go')
xlabel('x')
ylabel('y')
axis([-max_rad max_rad -max_rad max_rad])