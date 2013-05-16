clc
clear all
close all
format compact


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%   Constant Pile   %%%%%%%%%%%%%%%
lamda0=1310E-9;
Fm1=300;            % modulation frequence
Fm2=25000;                        
Fd= 300;           % vibration freq
N0=3.0;            % para for initialise deplacement 100000
n=1;               % in the air
Am1=0.00;%E-5;       % modulation amplitude
Am2=0.00; 
k1=1;              % on pratice : k1=coef_reflect*coef_trans
k2=0;
theta1=pi/3;
theta2=pi/6;
C=1;                % contrast of signal,


Fe=8E5;             % sampling freq
t=0:1/Fe:0.008;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%   Creat Currant Signal   %%%%%%%%%%%%%%%%%%
d=N0*lamda0*sin(2*pi*Fd.*t)/2;
D=10*lamda0/2+d;                    % total distance
dr=d/(lamda0/2);                    % reduced deplacement : dr=d/(lamda0/2)   

m1=Am1*N0*sin(2*pi*Fm1.*t);         % modulation of lamda, in order to creat a signal with a dephase pi/2
m2=Am2*N0*sin(2*pi*Fm2.*t);

phi=2*pi*n.*(dr+m1+m2);             % phi= delta phi

Km1=k1*sin(2*pi*Fm1.*t+theta1);
Km2=k2*sin(2*pi*Fm2.*t+theta2);% --> 0
K=Km1+Km2;                % la modulation de l’amplitude du signal interférométrique                       
s=C*cos(phi);      % recieved signal on the PD   s=K.*(1+C*cos(phi))    
s1=C*sin(phi);

tan=s1./s;
teta=atan(tan);
plot(teta);
figure(2)
polar(2*teta,0:1/length(teta):1-1/length(teta))
figure(3)
plot(s,s1)
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%   Figure section   %%%%%%%%%%%%%%%%%%%%%%%
figure(1);

plot(t,dr,'r');  hold on;
plot(t,phi);    hold on;    legend('Vibration','delta phi');

figure(2);
plot(t,dr/10,'r');  hold on;    
plot(t,s,'b');  hold on;
plot(t,s1,'g');  hold on;
legend('vibration','signal on PD'); Ylim([-2 2])





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%   Dephasage Function   %%%%%%%%%%%%%%%%%%%
decalage=linspace(0,10,length(s));
figure(3);
ang=acos(s)+pi;
rayon=sqrt(s.^2+s1.^2)+decalage;
polar(ang,rayon,'.');

% for count=1:length(s)
%     polar(ang,rayon,'.');
% %    plot(s1,s,'.');
% end
title('sin---cos');



% legend('FFT de signal');
% trace_fft(s,Fe,length(s));
%}

