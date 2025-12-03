function example8_main
% Oscilation of cantilever beam. Plane stress, rectangular domain
% Implicit Newmark method (example8_1.m) and explicit Newmark's method
% (example8_2.m).
% The results are compared by plotting deflection of the lower-right corner
% of the domain.

elementSize=[20, 6];
tfinal=3;
[ts1,u1,v1,a1]=example8_1(elementSize,tfinal); %#ok<*ASGLU>
[ts2,u2,v2,a2]=example8_2(elementSize,tfinal); %#ok<*ASGLU>
figure(2);clf;
hlID=elementSize(1)+1;
plot(ts1,u1(:,hlID*2),'-.',ts2,u2(:,hlID*2),'-.');
title({'Deflection of the lower-right corner of the beam'});
legend('Implicit','Explicit');
xlabel('time');
ylabel('u_y');
grid on;
end
