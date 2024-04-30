% Simulaton results of example2 under static event-triggered scheme

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

Aftilde=[0.6128 -0.4507 -0.1835;0.4267 0.0294 -0.1297;-0.0533 0.6753 0.0078];
Bftilde=[-1.4031;-1.3963;-1.1148];
MCftilde=[-0.0020 -0.0266 -0.0056];
Mphi=1.6059;

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

        else % There is DoS attack.
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
            theta(k)=0.2;
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
            theta(k)=0.2;
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


