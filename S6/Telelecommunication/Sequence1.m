%%--------------------------------------------------------------------------
%**************************************************************************
%-----------------Etude des modulateurs de base ---------------------------
%------------------------(Séquence 1)--------------------------------------

%--------------------------------------------------------------------------
%********************Auteur : AL GHARRAS IBRAHIM***************************
%--------------------------------------------------------------------------



%%-----------------------------Paramètres----------------------------------
clear 
close all;
clc;

%nombre de bits
nb_bits = 1000;

%Fréquence d'échantillonnage
Fe = 24000;
Te = 1/Fe;

%Débit binaire
Rb = 3000;
Tb = 1/Rb;

%Ordres du modulateurs
ord_m1 = 2;
ord_m2 = 4;


% --------------------Implantation des modulateurs-------------------------
%----------------------------Modulateur1-----------------------------------
 %Debit Symbole
Rs_1 = Rb/log2(ord_m1);
Ts_1 = 1/Rs_1;
Ns_1 = floor(Ts_1/Te);

    %Generation des bits
bits_1 = randi([0,1],1,nb_bits);
  
    %Mapping
Symboles_1 = 2*bits_1-1;

    %Surechantillonnage
Suite_diracs_1 = kron(Symboles_1, [1 zeros(1,Ns_1-1)]);

    %Filtrage de mise en forme
h = ones(1,Ns_1);
 
    %Filtrage de mise en forme
x_1 = filter(h,1,Suite_diracs_1);  

    %DSP
DSP_1 = abs(fft(x_1)).^2;
f =linspace(-Fe/2,Fe/2,length(DSP_1));
DSP_1_theo = Ts_1 .* (sin(pi.*f.*Ts_1)./pi.*f.*Ts_1).^2;

   %Tracé du signal transmis et de sa DSP
figure(1);
plot((0:Te:length(x_1)*Te-Te),x_1);
%axis([0 nb_bits-1 -1.5 1.5])
xlabel("Temps en secondes");
ylabel("Signal transmis");
title("Signal 1 transmis en fonction du temps");

figure(2);
semilogy(linspace(-Fe/2,Fe/2,length(DSP_1)),fftshift(DSP_1));
xlabel("Fréquence en Hz");
ylabel("DSP1");
title("DSP 1 en fonction de la fréquence");
hold on
semilogy(linspace(-Fe/2,Fe/2,length(DSP_1)),fftshift(DSP_1_theo));


%----------------------------Modulateur2-----------------------------------

    %Debit Symbole
Rs_2 = Rb/log2(ord_m2);
Ts_2 = 1/Rs_2;
Ns_2 = floor(Ts_2/Te);

    %Generation des bits
bits_2 = randi([0,1],1,nb_bits);
  

    %Mapping
Symboles_2 = 2*(bi2de((reshape(bits_2,length(bits_2)/2,2))))-3;%%%

    %Surechantillonnage
Suite_diracs_2 = kron(Symboles_2, [1 zeros(1,Ns_2-1)]);

    %Filtrage de mise en forme
h = ones(1,Ns_2);
 
 
x_2 = filter(h,1,Suite_diracs_2);  

    %DSP
Dsp_2 = abs(fft(x_2)).^2;
DSP_2_theo = Ts_2 .* (sin(pi.*f.*Ts_2)./pi.*f.*Ts_2).^2;%%

    %Tracé du signal transmis et de sa DSP
figure(3);
plot((0:Te:length(x_2)*Te-Te),x_2);
xlabel("Temps en secondes");
ylabel("Signal transmis");
title("Signal 2 transmis en fonction du temps");

figure(4);
semilogy(linspace(-Fe/2,Fe/2,length(Dsp_2)),fftshift(Dsp_2));
xlabel("Fréquence en Hz");
ylabel("DSP_2");
title("DSP_2 en fonction de la fréquence");
hold on
semilogy(linspace(-Fe/2,Fe/2,length(DSP_2_theo)),fftshift(DSP_2_theo));



%-------------------------%Modulateur3-------------------------------------
 
%Debit Symbole
Rs_3 = Rb/log2(ord_m1);
Ts_3 = 1/Rs_3;
Ns_3 = floor(Ts_3/Te);
alpha = rand;

L = floor((nb_bits-1)/Ns_3);
%Generation des bits
bits_3 = randi([0,1],1,nb_bits);

%Mapping
Symboles_3 = 2*bits_3-1;

%Surechantillonnage
Suites_diracs_3 = kron(Symboles_3 , [1 zeros(1,Ns_3 - 1)]);

%Filtrage de mise en forme
h = rcosdesign(alpha, L, Ns_3);
x_3 = filter(h,1,Suites_diracs_3);

%La densité spectrale de puissance
Dsp_3 = abs(fft(x_3)).^2;

    %Tracé du signal transmis et de sa DSP
figure(5);
plot((0:Te:length(x_3)*Te-Te),x_3);
xlabel("Temps en secondes");
ylabel("Signal transmis");
title("Signal 3 transmis en fonction du temps");


%----------------Comparaison des deux tracés-------------------------------
               %(Dsp obtenu avec Dsp théorique)

% DSP théorique du signal transmis 1 : 
Rs1 = Rb ; 



% DSP théorique du signal transmis 3 : 
F3 = -Fe/2:Fe/(Fe-1):Fe/2;
Ts = 1/Rb ; 
DSP3_theo = zeros(1,Fe);
cond1 = (abs(F3) <= (1+alpha)/(2*Ts)) & (abs(F3) >= (1-alpha)/(2*Ts));
ind_3 = find(cond1 == 1);
cond2 = (abs(F3)<(1-alpha)/(2*Ts));
DSP3_theo(cond1==1) = (1/2).* ( ones(1,length(ind_3)) + cos((pi*Ts/alpha).*(abs(F3(ind_3)) - ((1-alpha)/(2*Ts)).*ones(1,length(ind_3)))));
DSP3_theo(cond2==1) = 1;

% Comparaison des DSP
figure(6);
hold on
Sign = linspace(-Fe/2,Fe/2,length(Dsp_3));
plot(Sign,fftshift(Dsp_3 / max(Dsp_3)));
xlabel("Fréquence en Hz");
ylabel("DSP_3");
plot(F3 , DSP3_theo / max(DSP3_theo) , 'DisplayName' , 'DSPtheorique3');
hold off
title("Comparaison DSPs du Modulateur 3");

%------------------------------FIN-----------------------------------------
%------------------------------FIN-----------------------------------------
