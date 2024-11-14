clc, clearvars, clear all

simulation_data = readtable("data.csv");
wire_data = readtable("wire_data.csv");

figure;

plot3(wire_data.x, wire_data.y, wire_data.z, 'Color',[.7,.2,0.7]);
hold on

plot3(simulation_data.x, simulation_data.y, simulation_data.z, '.', 'Color',[1,0,0]);

daspect([1,1,1]);
xlabel('x');
ylabel('y');
zlabel('z');
grid on;



