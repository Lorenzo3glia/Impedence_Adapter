clc
clear all
close all

disp('Benvenuto nel software di calcolo: Adattatore Binomiale')
disp('Realizzato da Lorenzo Treglia')
%% Fornisco i dati obbligatori
strz0 = sprintf('Inserire impedenza Z0: ');
Z0=input(strz0);

%Nel caso fornito dal professore é fissa a 120 ohm
strzc = sprintf('Inserire impedenza ZC: ');
Zc=input(strzc);

strfc= sprintf('Inserire frequenza di centrobanda: ');
fc=input(strfc);

%% Inserisco le informazioni riguardate il tipo di analisi da effettuare se sull 'ordine o sulla banda
while true
    stropt = 'Inserire il carattere B per analisi in banda oppure il carattere O per analisi su ordine: ';
    opt = input(stropt, 's'); % 's' indica che l'input è una stringa
    if strcmp(opt, 'B') || strcmp(opt, 'O')
        break; % Esce dal ciclo solo se opt è 'B' oppure 'O'
    else
        disp('Input non valido. Per favore, inserire B oppure O.');
    end
end

%% Analisi in base alla banda+ Gamma_m
if opt == 'B'
    strzbw = sprintf('Inserire banda bilatera: ');
    bw=input(strzbw);
    banda_target= bw/fc;

    strzGm = sprintf('Inserire Gamma_M ');
    Gm=input(strzGm);

    N=0; %Minimo
    non_trovato=true;

    while non_trovato
        N=N+1;
        A=2^(-(N+1))*log(Zc/Z0);
        thm=acos(0.5*(Gm/abs(A))^(1/N));
        banda(N)=2-4/pi*thm;
        non_trovato=(banda(N)<banda_target);
    end
    %Da questo while ottengo il numero di stadi
    %Prova a fare con 4 stadi e vedi se si verifica ugualmente la condizione



    Z(1)=Z0;
    for n=0:N
        gamma(n+1)=A*nchoosek(N,n); %Coefficiente binomiale
        Z(n+2)=Z(n+1)*exp(2*gamma(n+1));
    end

    Z


%% Analisi in base all'ordine + Banda
else
    strzN = sprintf('Inserire Ordine adattatore: ');
    N=input(strzN);
    strzbw = sprintf('Inserire banda monolatera: ');
    bw=input(strzbw);
    banda_target= bw/fc;
    %     strzGm = sprintf('Inserire Gamma_M ');
    %     Gm=input(strzGm);
    

    A=2^(-(N+1))*log(Zc/Z0);
    %thm=2-(pi/4)*banda_target;
    thm=(2-banda_target)*pi/4;
    %thm=acos(0.5*(Gm/abs(A))^(1/N));
    Z(1)=Z0;
    for n=0:N
        gamma(n+1)=A*nchoosek(N,n); %Coefficiente binomiale
        Z(n+2)=Z(n+1)*exp(2*gamma(n+1));
    end

    Z

end








