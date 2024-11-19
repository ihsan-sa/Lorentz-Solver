clc, clearvars, clear all

addpath("C:\Users\ihsan\OneDrive\Documents\coding\ECE150\Lorentz-Solver")

p1 = readtable("p1.csv");
p2 = readtable("p2.csv");
p3 = readtable("p3.csv");
p4 = readtable("p4.csv");

figure;

%plot3(wire_data.x, wire_data.y, wire_data.z, 'Color',[.7,.2,0.7]); %plot a
%wire
%hold on
plot3(p1.x, p1.y, p1.z, '.', 'Color',[1,0,0]);
hold on
%plot3(data.x, data.y, data.z, '.', 'Color',[1,0,0]);
%plot3(p2.x, p2.y, p2.z, '.', 'Color','b');
%plot3(p3.x, p3.y, p3.z, '.', 'Color','g');
%plot3(p4.x, p4.y, p4.z, '.', 'Color','y');
%plot3(-1,0,0,'b*') %plot particle start
plot3(1,0,0,'r*') %plot particle start

%plot3(-3, 4,3.5, 'b*'); %plot a charge
xlabel('x');
ylabel('y');
zlabel('z');
grid on;



