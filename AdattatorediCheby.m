clc
clear all
close all


disp('Benvenuto nel software di calcolo: Adattatore di Chebyshev')
disp('Realizzato da Lorenzo Treglia')


%% Fornisco i dati obbligatori
strz0 = sprintf('Inserire impedenza Z0: ');
Z0=input(strz0);


strzc = sprintf('Inserire impedenza ZC: ');
Zc=input(strzc);

strfc= sprintf('Inserire frequenza di centrobanda: ');
fc=input(strfc);

%% Inserisco le informazioni riguardate il tipo di analisi da effettuare se sull 'ordine o sulla banda
while true
    stropt = 'Inserire il carattere A per trovare il numero di elementi,inserire B per trovare il coefficente di riflessione oppure C per trovare la larghezza di banda : ';
    opt = input(stropt, 's'); % 's' indica che l'input è una stringa
    if strcmp(opt, 'A') || strcmp(opt, 'B')|| strcmp(opt, 'C')
        break; % Esce dal ciclo solo se opt è 'B' oppure 'O'
    else
        disp('Input non valido. Per favore, inserire B oppure O.');
    end
end

%% Analisi avendo la banda bilatera e il coefficente di riflessione
%In questo modo otteniamo il numero di elementi per soddisfare il vincolo
%imposto da banda e coefficiente di riflessione

if opt == 'A'
    strzbw = sprintf('Inserire banda bilatera: ');
    bw=input(strzbw);
    banda_target= bw/fc;

    %     strzGm = sprintf('Inserire Gamma_M ');
    %     Gm=input(strzGm);
    strzGm = sprintf('Inserire Gamma_M [dB] ');
    Gmdb=input(strzGm);
    Gm=10^(-Gmdb/20);

    N=0; %Minimo
    non_trovato=true;

    while non_trovato
        N=N+1;
        A=Gm;
        thm=acos(1/cosh(1/N*acosh(abs(log(Zc/Z0)/(2*Gm)))));
        banda(N)=2-4/pi*thm;
        non_trovato=(banda(N)<banda_target);
    end

    th=linspace(thm,pi-thm,N+1) %N+1 Incognite e quindi N+1 Equazioni

    for k=0:N
        ths=th(k+1);
        for n=0:N %Mi muovo lungo le colonne
            MAT(k+1,n+1)=exp(-2*1j*ths*n);
        end
        %per la riga mi serve un ciclo for
        noto(k+1)=A*exp(-1j*N*ths)*cheby(N,cos(ths)/cos(thm));
    end

    gamma=real(inv(MAT)*noto.')

    Z(1)=Z0
    for n=0:N
        Z(n+2)=Z(n+1)*exp(2*gamma(n+1));
    end

    Z

    %% Analisi avendo ordine dell'adattatore e coefficiente di riflessione massimo
    %In questo modo quello che facciamo é calcolarci la banda attraverso il
    %coefficente di riflessione e l'ordine del filtro

else if opt == 'C'
    strzN = sprintf('Inserire Ordine adattatore: ');
    N=input(strzN);
%     strzGm = sprintf('Inserire Gamma_M ');
%     Gm=input(strzGm);
    strzGm = sprintf('Inserire Gamma_M [dB] ');
    Gmdb=input(strzGm);
    Gm=10^(-Gmdb/20);

    A=Gm;
    thm=acos(1/cosh(1/N*acosh(abs(log(Zc/Z0)/(2*Gm)))));


    th=linspace(thm,pi-thm,N+1) %N+1 Incognite e quindi N+1 Equazioni
    for k=0:N
        ths=th(k+1);
        for n=0:N %Mi muovo lungo le colonne
            MAT(k+1,n+1)=exp(-2*1j*ths*n);
        end
        %per la riga mi serve un ciclo for
        noto(k+1)=A*exp(-1j*N*ths)*cheby(N,cos(ths)/cos(thm));
    end


    gamma=real(inv(MAT)*noto.')
    Z(1)=Z0
    for n=0:N
        Z(n+2)=Z(n+1)*exp(2*gamma(n+1));
    end

    Z
    %% Analisi attriaverso ordine dell'adattatore e banda bilatera
    %In questo caso invece andiamo a trovare i valori di impedenza
    %attraverso l'ordine dell'adattatore e la sua banda bilatera
else
    strzbw = sprintf('Inserire banda bilatera: ');
    bw=input(strzbw);
    banda_target= bw/fc;
    strzN = sprintf('Inserire Ordine adattatore: ');
    N=input(strzN);
  

    %thm=2-(pi/4)*banda_target
    thm=(2-banda_target)*pi/4; 
    Gm = abs(log(Zc/Z0)) / (2 * cosh(N * acosh(1/cos(thm))));


    A=Gm;

    th=linspace(thm,pi-thm,N+1) %N+1 Incognite e quindi N+1 Equazioni

    for k=0:N
        ths=th(k+1);
        for n=0:N %Mi muovo lungo le colonne
            MAT(k+1,n+1)=exp(-2*1j*ths*n);
        end
        %per la riga mi serve un ciclo for
        noto(k+1)=A*exp(-1j*N*ths)*cheby(N,cos(ths)/cos(thm));
    end

    gamma=real(inv(MAT)*noto.')

    Z(1)=Z0
    for n=0:N
        Z(n+2)=Z(n+1)*exp(2*gamma(n+1));
    end

    Z
end
end






