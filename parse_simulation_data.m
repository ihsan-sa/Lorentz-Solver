clc, clearvars, clear all

addpath("C:\Users\ihsan\OneDrive\Documents\coding\ECE150\Other_CPP\Lorentz-Solver")
simulation_data = readtable("data.csv");
wire_data = readtable("wire_data.csv");

figure;

%plot3(wire_data.x, wire_data.y, wire_data.z, 'Color',[.7,.2,0.7]);
%hold on

plot3(simulation_data.x, simulation_data.y, simulation_data.z, '.', 'Color',[1,0,0]);
hold on
plot3(0,0,0,'gx')
plot3(0, -600, 100, 'b*');
plot3(0, -400, 750, 'b*');
plot3(-3000, 0, 0, 'b*');
daspect([1,1,1]);
xlabel('x');
ylabel('y');
zlabel('z');
grid on;



