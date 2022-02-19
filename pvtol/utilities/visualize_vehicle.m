function visualize_vehicle(vehicle)

pos = vehicle.pos;
phi = vehicle.phi;

left = vehicle.pos- [vehicle.l*cos(phi),vehicle.l*sin(phi)];
right = vehicle.pos + [vehicle.l*cos(phi),vehicle.l*sin(phi)];

plot(left,right,'Linewidth', 2);
