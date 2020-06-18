clear all; close all; clc;

data = csvread('5_original_2019-05-14-0957_2019-05-14-1006_eng.csv');

indFinal = 2000;

velocity = data(1:indFinal,26);
cof = data(1:indFinal,24);
pressure = data(1:indFinal,5);
acc = data(1:indFinal, 11);


tstep = 1/10000;
t = 0:tstep:(length(data)-1)*tstep;
t = 0:tstep:(indFinal-1)*tstep;


##plot(t, velocity)
##figure()
##plot(t, pressure)
##figure()
##plot(t, cof)

##plot(velocity, cof,'o')
##plot(pressure, cof, 'o')

x = [t', velocity, pressure, acc, cof];
plotmatrix(x)

figure()
scatter3(pressure, velocity.^2, cof)
xlabel('Pressure')
ylabel('Velocity Squared')
zlabel('COF')

data(:,90) = t';