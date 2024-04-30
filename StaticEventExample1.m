% Simulaton results of example1 under static event-triggered scheme

clear
close all
clc

% The parameters of continuous time system
Ac=[-0.1 50;-1 -10];
Bc=[0;1];
Fc=[-0.5 0;0 0];
Cc=[1 0];
Dc=1;
Gc=0;
Lc=[1 0];

% Discretize the continuous time system.
T=0.01;  
A=eye(2)+T*Ac;
B=T*Bc;
F=T*Fc;
C=Cc;
D=Dc;
G=Gc;
L=Lc;

% The given parameters
theta1=0.2;
theta2=0.8;
theta3=1;
Mtilde=[0.002;0.001];
Mtilde2=-0.3;
N1tilde=[0.3 0.1];
N2tilde=0.3;
N3tilde=[0.3 0.1];
gamma=1.8;
Fbar=[0.1 0;0 0];
Gbar=[0 0;0 0];
alphabar=0.5;
a=(1-alphabar)*alphabar;
H=Fbar'*Fbar+Gbar'*Gbar;

Aftilde=[0.9116 0.3927;-0.0099 0.7645];
Bftilde=[-0.0431;-0.0047];
MCftilde=[-0.9798 -0.0942];
Mphi=0.4721;

load('DoSData.mat');% Load the DoS attack data.

x(:,1)=[0;0];
xf(:,1)=[0;0];
for k=1:length(alpha)
    if k==1
        if alpha(k)==0 % no DoS attack
            omega(k)=0.1*cos(0.3*k)*exp(-0.1*k);
            f(:,k)=[x(1,k)^3;0];
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
            omega(k)=0.1*cos(0.3*k)*exp(-0.1*k);
            f(:,k)=[x(1,k)^3;0];
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
        if alpha(k)==0 % no DoS attack
            omega(k)=0.1*cos(0.3*k)*exp(-0.1*k);
            f(:,k)=[x(1,k)^3;0];
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
            omega(k)=0.1*cos(0.3*k)*exp(-0.1*k);
            f(:,k)=[x(1,k)^3;0];
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


