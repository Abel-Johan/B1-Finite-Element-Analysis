function example6_main
% Deflection of a bimetal study
% Solves steady state mechanical equilibrium and equilibrium heat equation on a
% rectangular domain (linear plane stress problem). The only aspect of
% thermomechanical coupling considered is heat expansion. See example6.m
% The location of the interface between two metals (steel and brass) within
% the domain is varied. The deflection of the tip due to equilibrium
% temperature distribution 

y0s=0:.1:2;
Tr=300;
d1s=zeros(numel(y0s),2);
d2s=zeros(numel(y0s),2);
figure(1);
for i=1:numel(y0s)
    [d1s(i,:),d2s(i,:)]=example6(y0s(i),Tr);
    disp(i);
end

figure(2);clf;
subplot(2,1,1)
plot(y0s,d1s(:,1),y0s,d2s(:,1));
title('Tip displacement x');
legend('upper corner','lower corner')
grid minor;
subplot(2,1,2)
plot(y0s,d1s(:,2),y0s,d2s(:,2));
title('Tip displacement y');
legend('upper corner','lower corner')
xlabel('y0 (interface position)')
grid minor;

end