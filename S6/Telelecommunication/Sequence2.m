%--------------------------------------------------------------------------
%**************************************************************************
%-----------------Etude des interférences entre symbole -------------------
%------------------et du critère de Nyquist (Séquence 2)------------------

%--------------------------------------------------------------------------
%********************Auteur : AL GHARRAS IBRAHIM***************************
%--------------------------------------------------------------------------



%%-----------------------------Paramètres----------------------------------

clear ;
close all ;

%nombre de bits
n_bits = 1000 ;

%fréquence d'échantillonage
Fe = 24000 ;

%Période d'échantillonnage
Te=1/Fe ;

%Débit binaire
Rb = 3000 ;
Tb = 1/Rb ;

%Mapping binnaire
M1 = 2;

%%Débit symbole
Rs1 = Rb/log2(M1) ;
Ts1 = 1/Rs1 ;
Ns = floor(Ts1/Te) ;

%---------------------Etude sans canal de propagation----------------------
%--------------------------------------------------------------------------
%Implantation du bloc modulateur/démodulateur proposé 

%Modulateur 

%Generation des bits et des symboles
bits_1 = randi([0,1],1,n_bits) ;
symboles_1 = 2*bits_1 -1 ;

%Réponse impulsionnele du filtre de mise en forme
h = ones(1,Ns) ;

%Réponse impulsionnele du filtre de récéption
hr = h ;

%Génération de la suite de Diracs pondérés par les symboles(suréchantillonage)
diracs = kron(symboles_1, [1 zeros(1,Ns-1)]) ;

%signal de sortie du modulateur :
x1 =filter(h,1,diracs);
z_1 = filter(hr,1,x1);

%Tracé du signal en sortie du filtre de réception
figure(1);
plot((0:Te:length(z_1)*Te-Te),z_1);
xlabel("Temps en secondes");
ylabel("Signal transmis");
title("Signal 1 en sortie du filtre du récéption en fonction du temps");

%réponse impulsionnelle global de la chaine de la transmission sans canal (g) :
g = conv(h,hr) ;
figure(2) 
plot(g) 
xlabel("Temps en secondes");
ylabel("g(t)");
n0 = Ns ;
%d'après la figure de g , g a valeur maximale en Ts1 = Te*Ns , et g(Ts1 + m*Ts1) = 0 pour tout m entier (g fenetre rectangulaire de durée 2*Ts1)

%diagramme de l'oeil
figure(3)
plot(reshape(z_1(n0 +1:end),[Ns,(length(z_1(n0 +1:end))/Ns)]))
title("diagramme de l'oeil pour la première chaine de transmission")

%échantillonage de z 
z_echant = z_1(n0:Ns:end);

%Détecteur de seuil
i_bit1 = find(z_echant > 0);
i_bit0 = find(z_echant < 0);
symgenere = -ones(1,length(bits_1));
symgenere(i_bit1) = ones(1,length(i_bit1));
symgenere(i_bit0) = zeros(1,length(i_bit0));
sym_err = bits_1 - symgenere;
indicediff = find(sym_err ~= 0);
TEB = length(indicediff)/length(bits_1); 

%On reprend la même démarche en prennant n0_2 = 3
n0_2 = 3;
z_echant_2 = z_1(n0_2:Ns:end);
i_bit1_2 = find(z_echant_2 > 0);
i_bit0_2 = find(z_echant_2 < 0);
symgenere_2 = -ones(1,length(bits_1));
symgenere_2(i_bit1_2) = ones(1,length(i_bit1_2));
symgenere_2(i_bit0_2) = zeros(1,length(i_bit0_2));
sym_err_2 = bits_1 - symgenere_2;
indicediff_2 = find(sym_err_2 ~= 0);
TEB_2 = length(indicediff_2)/length(bits_1);


%------------------Etude avec canal de propagation sans bruit--------------
%--------------------------------------------------------------------------


%ordre de filtre passe de canal
N=length(x1);

%réponse en fréquence du filtre de mise en forme
H = abs(fftshift(fft(h,N)));

%réponse en fréquence du filtre de réception
Hr = H; 

%Bande de fréquence dot on va afficher les réponses en fréquence
Bdf = linspace(-(Fe/2), Fe/2, N); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pour BW = 8000Hz :
fc1 = 8000;

%Réponse impulstionnelle du filtre passe bas
hc1 = (2*fc1/Fe)*sinc(2*(fc1/Fe)*[-(N-1)/2 : (N-1)/2]);

% réponse impulsionnelle global de la chaine de la transmission sans canal (gtotal) :
gtotal1 = conv(h,conv(hc1,h));

% Traçage de la réponse impulsionnelle globale de la chaine de transmission
figure(4) 
plot(gtotal1)
xlabel("temps en secondes") ; 
ylabel("gtotal(t)") ; 
title("Réponse impulsionnelle globale de la chaine avec canal de propagation") ;

% Traçage du diagramme de l'oeil à la sortie du filtre de réception
%filtrage
z_2 = filter(hc1,1,z_1);
figure(5)
plot(reshape(z_2(n0 +1:end),[Ns,(length(z_2(n0 +1:end))/Ns)]))
title("diagramme de l'oeil à la sortie du filtre de réception") ;

%Réponse fréquentielle du filtre passe bas
Hc1 = abs(fftshift(fft(hc1,N)));
%normalisation
Hc_n = Hc1 / max(Hc1) ; 

%réponse en fréquence du filtre de canal de propagation
% représentation de |H(f)Hr(F)| et |Hc(f)| sur une même figure : 
H = abs(fftshift(fft(h , length(x1)))) ; 
Prod = H.^ 2 ; 

% normalisation  
Prod_n = Prod / max(Prod) ; 

figure(6) ; 
hold on ; 
plot( linspace(-Fe/2,Fe/2,length(x1)) , Prod_n) ; 
plot( linspace(-Fe/2,Fe/2,length(x1)) , Hc_n ) ;
legend('|Hc(f)|','|H(f) x Hr(f)|');
xlabel("f en Hz") ; 
title("représentation de |H(f)Hr(F)| et |Hc(f)|");
hold off ; 




% Détermination du TEB pour BW = 8000 Hz : 

% échantillonage : 
t0 = Ns; 
Z2_echant = z_2(t0 : Ns : end) ; 

% detecteur de seuil en 0: 
Chapeau_2 = sign(Z2_echant) ;  

% demapping : 
Sym_chap_2 = ceil((Chapeau_2 + 1)/2) ; 

% TEB :  
TEB2_8000 = sum( (Sym_chap_2 - bits_1 )~=0 )/n_bits

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Pour BW = 1000Hz :
fc2 = 1000;

%Réponse impulstionnelle du filtre passe bas
hc1 = (2*fc2/Fe)*sinc(2*(fc2/Fe)*[-(N-1)/2 : (N-1)/2]);

% réponse impulsionnelle global de la chaine de la transmission sans canal (gtotal) :
gtotal1 = conv(h,conv(hc1,h));

% Traçage de la réponse impulsionnelle globale de la chaine de transmission
figure(7) 
plot(gtotal1)
xlabel("temps en secondes") ; 
ylabel("gtotal(t)") ; 
title("Réponse impulsionnelle globale de la chaine avec canal de propagation") ;

% Traçage du diagramme de l'oeil à la sortie du filtre de réception
%filtrage
z_2 = filter(hc1,1,z_1);
figure(8)
plot(reshape(z_2(n0 +1:end),[Ns,(length(z_2(n0 +1:end))/Ns)]))
title("diagramme de l'oeil à la sortie du filtre de réception") ;

%Réponse fréquentielle du filtre passe bas
Hc1 = abs(fftshift(fft(hc1,N)));
%normalisation
Hc_n = Hc1 / max(Hc1) ; 

%réponse en fréquence du filtre de canal de propagation
% représentation de |H(f)Hr(F)| et |Hc(f)| sur une même figure : 
H = abs(fftshift(fft(h , length(x1)))) ; 
Prod = H.^ 2 ; 

% normalisation  
Prod_n = Prod / max(Prod) ; 

figure(9) ; 
hold on ; 
plot( linspace(-Fe/2,Fe/2,length(x1)) , Prod_n) ; 
plot( linspace(-Fe/2,Fe/2,length(x1)) , Hc_n ) ;
legend('|Hc(f)|','|H(f) x Hr(f)|');
xlabel("f en Hz") ; 
title("représentation de |H(f)Hr(F)| et |Hc(f)|");
hold off ; 



% Détermination du TEB pour BW = 1000 Hz : 

% échantillonage : 
t0 = Ns; 
Z2_echant = z_2(t0 : Ns : end) ; 

% detecteur de seuil en 0: 
Chapeau_2 = sign(Z2_echant) ;  

% demapping : 
Sym_chap_2 = ceil((Chapeau_2 + 1)/2) ; 

% TEB :  
TEB2_1000 = sum( (Sym_chap_2 - bits_1 )~=0 )/n_bits


%------------------------------FIN-----------------------------------------
%------------------------------FIN-----------------------------------------
