%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation of Examples from Asymptotic Properties of Set Dynamical
% Systems
% Example:  Artstein's Circles, First simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
% Initial Set
d1=10;
x10=linspace(0.4,0.5,d1);
% n1=length(x10);
x20=linspace(0.4,0.5,d1);
[X10,X20] = meshgrid(x10,x20);
% X0 contains all the possible combinations of initial conditions given x10
% and x20 
X0 = [X10; X20];

% Number of iterations 
J=250;

% Simulation
% Xf: J x d matrix that will contains set evolution
Xf=X0;
u=0.01;

for j=1:J
    % Separate Xf recursively into 2 matrices that contains the state Xij,
    % where i=1,2 and j=1...J.
    C = mat2cell(Xf(:,:,j), [d1 d1], [d1]);
    x1=C{1,1};
    x2=C{2,1};
    % x1plus=g1(x1,x2)
    x1plus= (-(x1.^2+x2.^2).*u+x1)./(1+(x1.^2+x2.^2).*u^2-2*x1.*u);
    % x2plus=g2(x1,x2)
    x2plus= x2./(1+(x1.^2+x2.^2)*u^2-2*x1.*u);
    % check if solution is in D=[0,infty) times [0,infty) before adding it
    % to Xf
    x1sign = (x1plus<0); % compute entries with negative x1
    x2sign = (x2plus<0); % compute entries with negative x2
    inD = isempty(find(x1sign)) * isempty(find(x2sign)); % inD = 1 if both terms are 1 (meaning empty find)
    if inD == 0 % not in D
       break 
    end
    % Add G(Xj) to Xf
    Xf(:,:,j+1)=[x1plus;x2plus];
end
% Plots
figure(1), clf
%figure(2), clf

Jfinal = size(Xf);
Jfinal = Jfinal(3);

% Computation of Xf and plot
for j=1:Jfinal
XF = mat2cell(Xf(:,:,j), [d1 d1], [d1]);
x_1=XF{1,1};
x_1=x_1(:)';
x_2=XF{2,1};
x_2=x_2(:)';
figure(1)
% plot3(x_1,x_2,x_1*0+j-1,'*'); grid on; hold on; xlabel('x1j');ylabel('x2j');zlabel('j');hold on
% figure(2)
plot(x_1,x_2,'*'); grid on; hold on; xlabel('x1j');ylabel('x2j'); hold on; axis equal
end
axis([-0.1 0.65 -0.1 1.2])
%axis equal
print -depsc -tiff -r300 ArtsteinCirclesSim2
%print(fig,'-dpdf','ArtsteinCircles.pdf')
