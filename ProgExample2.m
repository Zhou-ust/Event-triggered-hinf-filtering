% Simulaton results of example2 under adaptive event-triggered scheme

clear
close all
clc

% The parameters of discrete time system
A=[0.75 -0.5 0;0.7 0 0;0 1 0];
B=[0.2;0.3;0.6];
F=eye(3);
C=[0.1 0 0];
D=0.4;
G=0;
L=[0.01 0.04 0.01];

% The given parameters
theta1=0.2;
theta2=0.8;
theta3=1;
Mtilde=[0.4;0.2;0.3];
Mtilde2=0.1;
N1tilde=[0.3 0.1 0.2];
N2tilde=0.2;
N3tilde=[0.3 0.4 0.2];
gamma=1.1;
Fbar=[0.03 0 0;0 0.01 0;0 0 0.02];
Gbar=[0 0 0;0 0 0;0 0 0];
alphabar=0.5;
a=(1-alphabar)*alphabar;
H=Fbar'*Fbar+Gbar'*Gbar;

% Solve the LMI of theorem2.
setlmis([]);
P1=lmivar(1,[2 1]);
P2=lmivar(1,[2 1]);
phi=lmivar(1,[1 1]);
Afhat=lmivar(2,[2 2]);
Bfhat=lmivar(2,[2 1]);
Cftilde=lmivar(2,[1 2]);
rho=lmivar(1,[1 1]);

P1=lmivar(1,[3 1]);
P2=lmivar(1,[3 1]);
phi=lmivar(1,[1 1]);
Afhat=lmivar(2,[3 3]);
Bfhat=lmivar(2,[3 1]);
Cftilde=lmivar(2,[1 3]);
rho=lmivar(1,[1 1]);

lmiterm([1 1 1 P1],-1,1);
lmiterm([1 1 1 0],H);
lmiterm([1 1 1 phi],theta2*C',C);

lmiterm([1 2 1 P2],-1,1);
lmiterm([1 2 2 P2],-1,1);

lmiterm([1 3 1 phi],theta2*D',C);
lmiterm([1 3 3 phi],theta2*D',D);
lmiterm([1 3 3 0],-gamma^2);

lmiterm([1 4 4 0],-1);

lmiterm([1 5 1 phi],theta2*G',C);
lmiterm([1 5 3 phi],theta2*G',D);
lmiterm([1 5 5 0],-1);
lmiterm([1 5 5 phi],theta2*G',G);

lmiterm([1 6 6 phi],-1,1);

lmiterm([1 7 1 P1],1,A);
lmiterm([1 7 1 Bfhat],(1-alphabar),C);
lmiterm([1 7 2 Afhat],1,1);
lmiterm([1 7 3 P1],1,B);
lmiterm([1 7 3 Bfhat],(1-alphabar),D);
lmiterm([1 7 4 P1],1,F);
lmiterm([1 7 4 Bfhat],(1-alphabar),G);
lmiterm([1 7 6 Bfhat],(alphabar-1),1);
lmiterm([1 7 7 P1],-1,1);

lmiterm([1 8 1 P2],1,A);
lmiterm([1 8 1 Bfhat],(1-alphabar),C);
lmiterm([1 8 2 Afhat],1,1);
lmiterm([1 8 3 P2],1,B);
lmiterm([1 8 3 Bfhat],(1-alphabar),D);
lmiterm([1 8 4 P2],1,F);
lmiterm([1 8 4 Bfhat],(1-alphabar),G);
lmiterm([1 8 6 Bfhat],(alphabar-1),1);
lmiterm([1 8 7 P2],-1,1);
lmiterm([1 8 8 P2],-1,1);

lmiterm([1 9 1 Bfhat],sqrt(a),C);
lmiterm([1 9 3 Bfhat],sqrt(a),D);
lmiterm([1 9 4 Bfhat],sqrt(a),G);
lmiterm([1 9 6 Bfhat],-sqrt(a),1);
lmiterm([1 9 9 P1],-1,1);

lmiterm([1 10 1 Bfhat],sqrt(a),C);
lmiterm([1 10 3 Bfhat],sqrt(a),D);
lmiterm([1 10 4 Bfhat],sqrt(a),G);
lmiterm([1 10 6 Bfhat],-sqrt(a),1);
lmiterm([1 10 9 P2],-1,1);
lmiterm([1 10 10 P2],-1,1);

lmiterm([1 11 1 0],L);
lmiterm([1 11 2 Cftilde],-1,1);
lmiterm([1 11 11 0],-1);

lmiterm([1 12 1 rho],1-alphabar,N2tilde*C);
lmiterm([1 12 2 rho],1,N1tilde);
lmiterm([1 12 3 rho],1-alphabar,N2tilde*D);
lmiterm([1 12 4 rho],1-alphabar,N2tilde*G);
lmiterm([1 12 6 rho],1-alphabar,N2tilde);
lmiterm([1 12 12 rho],-1,1);

lmiterm([1 13 1 rho],1,N2tilde*C);
lmiterm([1 13 3 rho],1,N2tilde*D);
lmiterm([1 13 4 rho],1,N2tilde*G);
lmiterm([1 13 6 rho],-1,N2tilde);
lmiterm([1 13 13 rho],-1,1);

lmiterm([1 14 2 rho],1,N3tilde);
lmiterm([1 14 14 rho],-1,1);

lmiterm([1 15 7 P2],Mtilde',1);
lmiterm([1 15 8 P2],Mtilde',1);
lmiterm([1 15 15 rho],-1,1);

lmiterm([1 16 9 P2],sqrt(a)*Mtilde',1);
lmiterm([1 16 10 P2],sqrt(a)*Mtilde',1);
lmiterm([1 16 16 rho],-1,1);

lmiterm([1 17 11 0],-Mtilde2');
lmiterm([1 17 17 rho],-1,1);

lmiterm([-2 1 1 P1],1,1);
lmiterm([-2 2 1 P2],1,1);
lmiterm([-2 2 2 P2],1,1)

lmiterm([-3 1 1 P1],1,1);
lmiterm([-4 1 1 P2],1,1);
lmiterm([-5 1 1 phi],1,1);
lmiterm([-6 1 1 rho],1,1);

lmisys=getlmis;
[tmin,xfeasp]=feasp(lmisys)

% Solve the matrix variables.
MP1=dec2mat(lmisys,xfeasp,P1)
MP2=dec2mat(lmisys,xfeasp,P2)
Mphi=dec2mat(lmisys,xfeasp,phi)
MAfhat=dec2mat(lmisys,xfeasp,Afhat)
MBfhat=dec2mat(lmisys,xfeasp,Bfhat)
MCftilde=dec2mat(lmisys,xfeasp,Cftilde)
Mrho=dec2mat(lmisys,xfeasp,rho)

Aftilde=inv(MP2)*MAfhat
Bftilde=inv(MP2)*MBfhat

load('DoSData.mat');% Load the DoS attack data.

x(:,1)=[0;0;0];
xf(:,1)=[0;0;0];
for k=1:length(alpha)
    if k==1
        if alpha(k)==0 %no DoS attack
            omega(k)=sin(k)*exp(-k);
            f(:,k)=[0.03*sin(x(1,k));0.01*sin(x(2,k));0.02*sin(x(3,k))];
            x(:,k+1)=A*x(:,k)+B*omega(k)+F*f(:,k);
            y(:,k)=C*x(:,k)+D*omega(k)+0;
            z(:,k)=L*x(:,k);
            % Event trigger
            sl=1;
            releaseinstant=sl;
            yf(:,k)=y(:,sl);

            delta(k)=sin(k);
            Af=Aftilde+Mtilde*delta(k)*N1tilde;
            Bf=Bftilde+Mtilde*delta(k)*N2tilde;
            Cf=MCftilde+Mtilde2*delta(k)*N3tilde;
            xf(:,k+1)=Af*xf(:,k)+Bf*yf(:,k);
            zf(:,k)=Cf*xf(:,k);

        else %there is DoS attack
            omega(k)=sin(k)*exp(-k);
            f(:,k)=[0.03*sin(x(1,k));0.01*sin(x(2,k));0.02*sin(x(3,k))];
            x(:,k+1)=A*x(:,k)+B*omega(k)+F*f(:,k);
            y(:,k)=C*x(:,k)+D*omega(k)+0;
            z(:,k)=L*x(:,k);
            % Event trigger
            sl=1;
            releaseinstant=sl;
            yf(:,k)=0;

            delta(k)=sin(k);
            Af=Aftilde+Mtilde*delta(k)*N1tilde;
            Bf=Bftilde+Mtilde*delta(k)*N2tilde;
            Cf=MCftilde+Mtilde2*delta(k)*N3tilde;
            xf(:,k+1)=Af*xf(:,k)+Bf*yf(:,k);
            zf(:,k)=Cf*xf(:,k);
        end
    end
    if k>=2 & k<=length(alpha)
        if alpha(k)==0 %no DoS attack
            omega(k)=sin(k)*exp(-k);
            f(:,k)=[0.03*sin(x(1,k));0.01*sin(x(2,k));0.02*sin(x(3,k))];
            x(:,k+1)=A*x(:,k)+B*omega(k)+F*f(:,k);
            y(:,k)=C*x(:,k)+D*omega(k)+0;
            z(:,k)=L*x(:,k);
            % Event trigger
            theta(k)=theta2+(theta1-theta2)*2/pi*atan(theta3*sqrt((y(:,k)-y(:,sl))'*(y(:,k)-y(:,sl))));
            if (y(:,k)-y(:,sl))'*Mphi*(y(:,k)-y(:,sl))>theta(k)*y(:,k)'*Mphi*y(:,k)
                sl=k;
                releaseinstant=[releaseinstant sl];
                yf(:,k)=y(:,sl);
            else
                yf(:,k)=0;
            end

            delta(k)=sin(k);
            Af=Aftilde+Mtilde*delta(k)*N1tilde;
            Bf=Bftilde+Mtilde*delta(k)*N2tilde;
            Cf=MCftilde+Mtilde2*delta(k)*N3tilde;
            xf(:,k+1)=Af*xf(:,k)+Bf*yf(:,k);
            zf(:,k)=Cf*xf(:,k);
        else % There is DoS attack.
            omega(k)=sin(k)*exp(-k);
            f(:,k)=[0.03*sin(x(1,k));0.01*sin(x(2,k));0.02*sin(x(3,k))];
            x(:,k+1)=A*x(:,k)+B*omega(k)+F*f(:,k);
            y(:,k)=C*x(:,k)+D*omega(k)+0;
            z(:,k)=L*x(:,k);
            % Event trigger
            theta(k)=theta2+(theta1-theta2)*2/pi*atan(theta3*sqrt((y(:,k)-y(:,sl))'*(y(:,k)-y(:,sl))));
            if (y(:,k)-y(:,sl))'*Mphi*(y(:,k)-y(:,sl))>theta(k)*y(:,k)'*Mphi*y(:,k)
                sl=k;
                releaseinstant=[releaseinstant sl];
                yf(:,k)=0;
            else
                yf(:,k)=0;
            end

            delta(k)=sin(k);
            Af=Aftilde+Mtilde*delta(k)*N1tilde;
            Bf=Bftilde+Mtilde*delta(k)*N2tilde;
            Cf=MCftilde+Mtilde2*delta(k)*N3tilde;
            xf(:,k+1)=Af*xf(:,k)+Bf*yf(:,k);
            zf(:,k)=Cf*xf(:,k);
        end
    end
end

length(releaseinstant)% triggering times

e=z-zf;
sqrt((e*e'))/sqrt((omega*omega'))% the steady ratio value

% Plot the instant of DoS attacks.
figure
set(gcf,'unit','centimeters','Position',[2 10 6 5]);
set(gca, 'FontName','Times New Roman','fontsize',8);
plot([1:200],alpha,'.')
axis([0,200,-0.2,1.2]);
xlabel('\fontname{Times New Roman}\fontsize{8}Time (k)')
ylabel('\fontname{Times New Roman}\fontsize{8}\alpha(k)')

% Plot the event-triggered instant and release interval.
figure
set(gcf,'unit','centimeters','Position',[10 10 6 5]);
set(gca, 'FontName','Times New Roman','fontsize',8);
releaseinstant1=releaseinstant-1;
releaseinterval=releaseinstant1-[0 releaseinstant1(1:length(releaseinstant)-1)];
stem(releaseinstant1,releaseinterval)
box on
xlabel('\fontname{Times New Roman}\fontsize{8}Time (k)')
ylabel('\fontname{Times New Roman}\fontsize{8}Release interval')

% Plot the responses of signals z(k) and zf(k).
figure
set(gcf,'unit','centimeters','Position',[18 10 6 5]);
set(gca, 'FontName','Times New Roman','fontsize',8);
plot([1:200],z(1,1:200),'-','LineWidth',2);hold on
% subplot(2,2,1)
plot([1:200],zf(1,1:200),'r--','LineWidth',2)
legend('z(k)','z_f(k)');
xlabel('\fontname{Times New Roman}\fontsize{8}Time (k)')
ylabel('\fontname{Times New Roman}\fontsize{8}z(k) and its estimation z_f(k)')

% PLot the filtering error.
figure
set(gcf,'unit','centimeters','Position',[2 2 6 5]);
set(gca, 'FontName','Times New Roman','fontsize',8);
plot([1:200],e(1,1:200),'Linewidth',2)
xlabel('\fontname{Times New Roman}\fontsize{8}Time (k)')
ylabel('\fontname{Times New Roman}\fontsize{8}Filtering error')

% PLot the evolution trajectory of adaptive threshold.
figure
set(gcf,'unit','centimeters','Position',[10 2 6 5]);
set(gca, 'FontName','Times New Roman','fontsize',8);
plot([1:200],theta(1,1:200),'-','LineWidth',2);
axis([0,200,0,0.9]);
set(gca,'YTick',[0:0.2:0.9]);
xlabel('\fontname{Times New Roman}\fontsize{8}Time (k)')
ylabel('\fontname{Times New Roman}\fontsize{8}Threshold \theta(k)')

%Plot the ratio.
figure
set(gcf,'unit','centimeters','Position',[18 2 6 5]);
set(gca, 'FontName','Times New Roman','fontsize',8);
Ratio = [];
for  i = 1:1:size(e,2)
    e_t = e(1,1:i); omega_t = omega(1,1:i);
    Ratio_t = sqrt((e_t*e_t'))/sqrt((omega_t*omega_t'));
    Ratio = [Ratio Ratio_t];
end
plot([1:200],Ratio, 'LineWidth',2)
xlabel('\fontname{Times New Roman}\fontsize{8}Time (k)')
ylabel('\fontname{Times New Roman}\fontsize{8}Ratio')

