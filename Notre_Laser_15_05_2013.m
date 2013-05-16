clc
clear all
close all
format compact

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Modélisation Laser
%---------------------------------------------
 
%*************************************************************************
%**
%**          Données et constantes
%**
%*************************************************************************
     n=1;%indice de réfraction 

     %run donnees
     c=299792458; %célérité de la lumière 
     Lambda0 = 1310E-9; %longueur d'onde dans le vide (nm)    
     Teta = 0; %angle d'incidence des rayons
     f0=c/Lambda0; %fréquence onde émise par le laser
     omega_0_E1=0; %phase à l'origine signal E1
     omega_0_E2=0; %phase à l'origine signal E2
     omega_0_E2y=0; %phase à l'origine signal E2y (quadrature)
     a0=1; %amplitude du champ

     %---> simulation des champs électriques:
        T1=96.6/100;  %end of fiber transmission
        R1=3.4/100;   %end of fiber reflection
        R2=90/100;    %reflector reflection
        T2=3.4/100;   %transmission back

    %--> numérisation HF:
        N=5e5; %nombre échantillon HF
        fe=f0*5; %fréquence d'échantillonage ( pour respecter shanon) : 100
        Te=1/fe;
        n_sample_period_HF=fe/f0;
        pas=Te;
        t=pas:pas:N*pas; %variable temps

    
     %--->  distance cible:
        L=600E-3; %distance de la cible initiale  (doit être > à l'amplitude de la vibration)
        f_deplacement=50e5; %fréquence vibration
        Amplitude_deplacement = 1965e-9;
        omega_0_BF=0; %phase à l'origine vibration
    
     %---> numérisation BF:  
        Temps_simulation = 1/f_deplacement; % (en secondes)--> 1 période de la vibration
         Fe_BF=f_deplacement*10; % Fréquence d'échantillonage pour les signeaux de basses-fréquences : *10
         Te_BF=1/Fe_BF; %Période d'échantillonage pour les signeaux de basses-fréquences
 

 pas_BF = Te_BF/Te; %pas d'échantillonage pour les BF par rapport au signal HF: nombre d'échantillon du signal à fe à compter entre chaque échantillon du signal BF
 signal_BF = zeros(1,Temps_simulation/Te_BF);
 figure(1)
 Nb_point_BF=250;
 j=1; %indice    
 
      %---->  signal quadrature
         Delta_lambda=Lambda0^2/(8*n*L); 
         Lambda1=Lambda0-Delta_lambda;
         
 
 %---> %  affichage:
 cpt=0;
 loop_aff=25; %affiche tous les x echantillons (BF)
 subplot(211)
 hold on  
 
 
 
%for i=1:1:(Temps_simulation*fe/N),

%test::::
wt=2*pi.*t*f0;
wt_deplacement=2*pi*f_deplacement.*t;

t_Ebase= pas: pas: (N + fe/f0)*pas; %(ajout d'une période en plus)
E_base=cos(2*pi.*t_Ebase*f0);

%{
test 2 :
% N_ref_cos=1024;
% ref_cos = cos(2*pi*(0:1/N_ref_cos:1-1/N_ref_cos));
% Pas_E2=0;
% pointeur=1;
% E2=zeros(1,N);
%f_sin=2;

% for i=1:1:N/8,
%         Pas_E2=f0+phi(i)/(2*pi);
%         pointeur=mod(pointeur+Pas_E2,1024);
%         pointeur=fix(pointeur);
%         testcos(i)=ref_cos(pointeur+1);
% end;
% plot(1:1:N/8,testcos)
 %}
 

%I_max=(R1*a0)^2+((T1*R2*0.05)*a0)^2+2*(R1*a0)*((T1*R2*0.05)*a0);

 test_d = zeros(1,Temps_simulation/Te_BF);
 
 
 
%**************************************************************************
%**
%**                      Boucle Simulation
%**
%**************************************************************************
for i=0:1:Nb_point_BF,

%---> distance:
    delta_d =Amplitude_deplacement*cos(wt_deplacement + omega_0_BF);  %déplacement (vibrations : cos...)
    d = (L+delta_d); %distance de la cible
    omega_0_BF=omega_0_BF+ (wt_deplacement(length(wt_deplacement))); %sauvegarde de la dernière valeur de la phase
   
%---> déphasage:
    delta_phi=4*pi*n*delta_d/Lambda0; %déphasage au réflecteur du au déplacement : vibrations
    phi=4*pi*n*d/Lambda0; %déphasage total au réflecteur
    phi1=4*pi*n*d/Lambda1; %déphasage total au réflecteur
    
%---> les champs réfléchis
   %E1=R1*a0*cos(wt + omega_0_E1);    %onde réfléchie par la fibre
  
  %test:::::::::::::::
   nbr_ech_dephasage = omega_0_E1 * n_sample_period_HF/(2*pi); %  wt%(2pi) * nbre_ech_periode / 2pi -> nombre d'échantillon correspondant au déphasage ramené sur 1 période

   E1=R1*a0*E_base(nbr_ech_dephasage+1 : (length(t)+nbr_ech_dephasage) ); %ok
   E2=(T1*R2*0.05)*a0*cos(wt + omega_0_E2 +phi); %onde réfléchie pas l'objet (T2=0.05 en réalité car pertes)
   E2y=(T1*R2*0.05)*a0*cos(wt + omega_0_E2y +phi1);
%     for ii=1:1:N, %essayer de faire rentre fe!!!
%             Pas_E2=fix((f0+phi(ii)/(2*pi))*N_ref_cos/fe);
%             pointeur=mod(pointeur+Pas_E2,N_ref_cos);
%            % pointeur=fix(pointeur);
%             E2(ii)=(T1*R2*0.05)*a0*ref_cos(pointeur+1);
%     end;
   

   
   
   
   
%---> %sauvegarde des phases finales:
 %   omega_0_E1 = omega_0_E1 + wt(length(wt));
    omega_0_E1 = mod( omega_0_E1 + wt(length(wt)) , 2*pi);
    omega_0_E2 = mod( omega_0_E2 + wt(length(wt)) , 2*pi);
    omega_0_E2y= mod( omega_0_E2y + wt(length(wt)) , 2*pi);
    
 
%    % AFFICHAGE champs
%     figure(i+4)
%     %subplot(311)
%     plot(t*fe,E1)
%     hold on
%     plot(t*fe,E2,'r')



% ---> simulation courant récepteur:
    I=(abs(E1+E2)).^2;
    Iy=(abs(E1+E2y)).^2;
    
    %I=((E1.*E3));
    I2=(R1*a0)^2+((T1*R2*0.05)*a0)^2+2*(R1*a0)*((T1*R2*0.05)*a0)*cos(phi); %avec les amplitude
    I2y=(R1*a0)^2+((T1*R2*0.05)*a0)^2+2*(R1*a0)*((T1*R2*0.05)*a0)*cos(phi1);%Vy
    
    
  %  I_norm=2*(I-I_max/2)./I_max; %courant normalisé
        
    signal_BF(j) = I(round(mod((i+pas_BF),N)+1));
    signal_Iy(j) = Iy(round(mod((i+pas_BF),N)+1));
    signal_I2(j) = I2(round(mod((i+pas_BF),N)+1));
    signal_I2y(j) = I2y(round(mod((i+pas_BF),N)+1));
    
    test_d(j)=d(round(mod((i+pas_BF),N)+1));
    
    
    %figure(2)
     %plot(I_test,'r')
     if cpt==loop_aff,
         
        figure(1)
        hold on
        %plot([j-1 j],[signal_BF(j-1) signal_BF(j)])
        plot(j-loop_aff:1:j,signal_I2(j-loop_aff:j))
        plot(j-loop_aff:1:j,signal_I2y(j-loop_aff:j),'r')
        cpt=0;
     end;


  j=j+1;
  cpt=cpt+1;
end; 




subplot(212)
I_filtre=filter(filtreIIR_But,signal_BF);
Iy_filtre=filter(filtreIIR_But,signal_Iy);

% --> Normalisation I: 
 I_norm =2*(I_filtre-max(I_filtre)/2)./max(I_filtre); %courant normalisé
 Iy_norm=2*(Iy_filtre-max(Iy_filtre)/2)./max(Iy_filtre); %courant normalisé



%plot(0:1:Nb_point_BF,I_filtre,'r')
hold on
plot(0:1:Nb_point_BF,Iy_norm,'g')
plot(0:1:Nb_point_BF,I_norm,'r')

figure(2)
plot(test_d-0.6,'r')
title ('mouvement');




% traçage du cercle:
figure(4)
polar(0,0) %pour avoir la grid
hold on
plot(I_norm,Iy_norm)

% et voila!


%   astuce, mesure temps code: 
%         tic
%         {code}
%         toc

% % AFFICHAGE courant:
%      figure(2)  
%      plot(t,I_norm)
    
% ---> FFT
%{
    Fs = 1000;                    % freq d'échantillonage
    T = 1/Fs;                     % periode d'echantillonage
    L = 1000;                     % longeur du signal
    t1 = (0:T:L-1);               % vecteur temps
    f3=50;                        % freq du cosinus

    x=cos(2*pi*f3.*t1);

    NFFT = 2^nextpow2(L);         %prochaine puiss de 2 à partir de L
    Y = fft(x,NFFT)/L;
    f = Fs/2*linspace(0,1,NFFT/2+1);

    figure(3)
    %plot de la FFT
    plot(f,2*abs(Y(1:NFFT/2+1))) 
    xlabel('Frequency (Hz)')
    ylabel('|Y(f)|')
 
%}
    

    % --- Filtre Passe-Bas --- %
             
        % Filtre RIF: définition du gabarit
%         f1 = f0/20;
%         f2 = f0/10;
%         F=[0 f1 f2 fe/2];
%         F=F*2/fe; %normalisation
%         A=[0 0 -30 -30]; %en dB
%         A=10.^(A/20); %normalisation
%         
%         Ordre_du_filtre = 5000;
%         b=fir2(Ordre_du_filtre,F,A); %retourne l'ordre minimum garentissant le gabarit   
%         I_detecteur = filter(b,A,I_test); %application du filtre
        
        % ------ affichage de la transmittance du filtre:
                %{
        [H,w]=freqz(b,1,Ordre_du_filtre); %retourne la réponse fréquentielle
       
        pas_filtre=1/Ordre_du_filtre;
        f=0:pas_filtre:1-pas_filtre;
        
        figure(9)
        semilogx(f*fe/2,20*log10(abs(H)));%trace la transmittance
        grid on;
        grid minor;
        hold on;
    
        %affichage du gabarit:
        plot([1e14 f1-1 f1+1 fe/2],...
             [0 0 -30 -30],'r')   
        plot([ f2-1 f2+1],[0 -30],'r')
    
        %affichage f0:
        plot(f0,-30,'bX')

        title('Pas-Bas RIF');
        xlabel('freq');
        ylabel('gain');
                 %}  

 
  %subplot(313)
   % plot(t,d)
 %   plot(I_detecteur,'r')%affichage I filtré
%      Axes = axis;
%      Axes(3)=-0.5;
%      Axes(4)=0.5;
%      axis(Axes);




