n = plant.n; nu = plant.nu;

% nominal trajectory
simuLen = length(times);
xnomTraj = zeros(n,simuLen);
unomTraj = zeros(nu,simuLen);
for t =1:simuLen
    xnomTraj(:,t) = x_nom_fcn(times(t));
    unomTraj(:,t) = u_nom_fcn(times(t));
end
if controller.type == CtrlDesignOpts.ccm
    addText = 'CCM';
else
    addText = 'RCCM';
end
close all;
figure(1); clf;
subplot(2,2,1); hold on;
plot(times,xnomTraj(1,:),'b--',times,xnomTraj(2,:),'r--');
plot(times,xTraj(1,:),'b-',times,xTraj(2,:),'r-');
xlabel('Time (s)')
ylabel('x & z (m)')
legend('x_{nom}', 'z_{nom}',['x: ' addText],['z: ' addText]);

subplot(2,2,2); hold on;
plot(times,xnomTraj(3,:)*180/pi,'--');
plot(times,xTraj(3,:)*180/pi);
xlabel('Time (s)')
ylabel('\phi (deg)')
legend('Nominal',addText);

subplot(2,2,3);hold on;
plot(times,unomTraj(1,:),'b--',times,unomTraj(2,:),'r--');
plot(times,uTraj(1,:),'b-',times,uTraj(2,:),'r-');
xlabel('Time (s)');
ylabel('u (N)')
legend('u_{nom,1}', 'u_{nom,2}',['u_1: ' addText],['u_2: ' addText]);

subplot(2,2,4)
plot(times, energyTraj);
xlabel('Time (s)');
ylabel('Riemann energy')
legend(addText)

figure(2);clf;
hold on;
h1 = plot(xnomTraj(1,:),xnomTraj(2,:),'b--');
h2 = plot(xTraj(1,:),xTraj(2,:),'r-');
scatter(xnomTraj(1,1),xnomTraj(1,1),'bo')
scatter(xnomTraj(1,end),xnomTraj(1,end),'b*')
scatter(xTraj(1,end),xTraj(2,end),'r*')
xlabel('z (m)')
ylabel('y (m)')
legend('Nominal', addText);


