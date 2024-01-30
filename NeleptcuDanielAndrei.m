clc, clear all, close all;
%% Neleptcu Daniel-Andrei 332AB - Modelul 2

%Definirea parametrilor corespunzatori modelului 2
par.m = 4.1;
par.l0 = 1.7;
par.zeta = 0.9;
par.k = 0.6;
par.g = 9.8;

%Definesc bus-ul pentru a putea transmite in blocul Constant o structura
%Matlab
par_bus_info = Simulink.Bus.createObject(par);
par_bus = evalin('base', par_bus_info.busName);

%% Cerinta 1 - Implementarea modelului matematic neliniar:
%Blocul Fcn am reusit sa il scot din 'simuling_need_slupdate'
%Atat blocul Fcn cat si Matlab function se afla in acelasi fisier simulink
%Deci pentru punctele a si b:

load_system('modelul2')

%% Cerinta 2 - Crearea unui Semnal de Treaptă:
%Definire timp pentru simulare
T = 30;
set_param('modelul2', 'StopTime', num2str(T));

%Semnalul de tip treapta
t = linspace(0,T,1000);
st = double(t>=0);
usim = timeseries(st,t);

%Simularea modelului
sim('modelul2');

%Am ales un orizont de timp suficient pentru a se observa ca se ajunge
%intr-un regim stationar, graficele stabilizandu-se
figure('Name', 'Grafice verificare regim stationar');
subplot(2,1,1);
plot(simout_mat1.time, simout_mat1.data, 'r', 'LineWidth', 2);
xlabel('timp');
ylabel('x');
hold on
subplot(2,1,2);
plot(simout_mat2.time,simout_mat2.data, 'b', 'LineWidth', 2);
xlabel('timp');
ylabel('theta');

%% Cerinta 3 - Simularea Modelului:
% Simulam modelul
sim('modelul2');

%Graficele pentru blocurile corespunzatoare lui x
figure;
plot(simout_mat1.time, simout_mat1.data, 'r', 'LineWidth', 2);
hold on;
plot(simout_fcn1.time,simout_fcn1.data, 'b', 'LineWidth', 1);
title('Grafice bloc x');
xlabel('Timp (s)'), ylabel('x')

%Graficele pentru blocurile corespunzatoare lui theta
figure;
plot(simout_mat2.time, simout_mat2.data, 'r', 'LineWidth', 2);
hold on;
plot(simout_fcn2.time,simout_fcn2.data, 'b', 'LineWidth', 1);
title('Grafice bloc theta');
xlabel('Timp (s)'), ylabel('Theta (rad)')


%Salvez conform cerintei raspunsurile in variabilele simout_mat1,
%simout_mat2, respectiv simout_fcn1 si simout_fcn2
simout_fcn1 = simout_fcn1.data;
simout_fcn2 = simout_fcn2.data;

simout_mat1 = simout_mat1.data;
simout_mat2 = simout_mat2.data;

%Se poate observa ca graficele corespunzatoare blocului fcn si matlab
%function coincid, deci cel mai probabil vom avea o eroare foarte mica la
%punctul urmator

%% Cerinta 4 - Evaluarea Erorii:
err1 = norm(simout_mat1 - simout_fcn1);
err2 = norm(simout_mat2 - simout_fcn2);
disp('Eroarea fcn/matlab function pentru x');
disp(err1);
disp('Eroarea fcn/matlab function pentru theta');
disp(err2);
%Eroarea este 0, deci blocul fcn si Matlab function fac exact acelasi lucru
%in acest caz.

%% Cerinta 5 - Ilustrați caracteristica statică de funcționare a modelului matematic:

%a) 20 de intrari de tip treapta
k = linspace(0.1, 3.2, 20).';
ustar = zeros(20, 1);
y1star = zeros(20, 1);
y2star = zeros(20, 1);

%b) Simularea modelului si salvarea valorilor
for i = 1:20
    in = k(i) .* double(t>=0);
    usim = timeseries(in, t);
    
    sim('modelul2');
    ustar(i) = k(i);
    %memorarea valorilor de regim stationar pentru blocul Matlab function
    y1star(i) = simout_mat1.Data(end);
    y2star(i) = simout_mat2.Data(end);
end

%c) Calcularea polinoamelor
ustar1 = linspace(ustar(1), ustar(end), 100);

%Polinom ordin 3
p1 = polyfit(ustar, y1star, 3);
p2 = polyfit(ustar, y2star, 3);
y1starr = polyval(p1, ustar1);
y2starr = polyval(p2, ustar1);

%Polinom ordin 1
p1_lin = polyfit(ustar, y1star, 1);
p2_lin = polyfit(ustar, y2star, 1);
y1starr_lin = polyval(p1_lin, ustar1);
y2starr_lin = polyval(p2_lin, ustar1);

% Cerinta d)
figure;
subplot(2,1,1);

%punctele corespunzatoare lui b)
plot(ustar, y1star, 'xk'); 
hold on

%graficele pentru primul polinom
plot(ustar1, y1starr, 'LineWidth', 2);
plot(ustar1, y1starr_lin, 'LineWidth', 1);
xlabel('u^*'), ylabel('y^*');
legend('Puncte', 'Aproximare polinom de ordinul 3', 'Aproximare polinom de ordinul 1');
title('Figura pentru y1star - cerinta 5');

subplot(2,1,2);

%punctele corespunzatoare lui b)
plot(ustar, y2star, 'xk'); 
hold on

%graficele pentru al doilea polinom
plot(ustar1, y2starr, 'LineWidth', 2);
plot(ustar1, y2starr_lin, 'LineWidth', 1);
xlabel('u^*'), ylabel('y^*');
legend('Puncte', 'Aproximare polinom de ordinul 3', 'Aproximare polinom de ordinul 1');
title('Figura pentru y2star - cerinta 5');

%% Cerinta 6 - Aproximarea răspunsurilor:

alfa = 3.6;
beta = 4.2;
gamma = 4.9;
ustar6 = [alfa, beta, gamma];
y1star6 = polyval(p1, ustar6);
y2star6 = polyval(p2, ustar6);

figure;
subplot(2,1,1);
plot(ustar, y1star, 'xk');
hold on
plot(ustar6, y1star6, 'LineWidth', 2);
xlabel('u^*'), ylabel('y^*');
legend('Puncte', 'Aproximare polinom de ordinul 3');
title('Figura y1star - cerinta 6');

subplot(2,1,2);
plot(ustar, y2star, 'xk');
hold on
plot(ustar6, y2star6, 'LineWidth', 2);
xlabel('u^*'), ylabel('y^*');
legend('Puncte', 'Aproximare polinom de ordinul 3');
title('Figura y2star - cerinta 6');

%Putem observa grafic ca aproximarea este foarte buna, o metoda de
%verificare fiind extinderea frontului declarat la cerinta 5 pentru a
%cuprinde si alfa beta si gamma selectate la punctul 6

%% Cerinta 7 - Model simulink cu blocuri In si Out ( acelasi model )
mdlinout = 'modelul2cerinta7';

%% Cerinta 8 - Determinare PSF pentru o eroare minima
u0 = ustar(1);
[xstar8,ustar8,ystar8,~] = trim(mdlinout,[],u0,[],[],1,[]);
err_min = norm((ustar8-u0),2);
%Am considerat o eroare mai mica de 10^-3 ca fiind o eroare mica, avand in
%vedere ca eroarea in cazul meu este 0 putem concluziona deja ca este
%minima fara a mai face verificari. Exista situatii in care trim nu poate
%returna fix punctul cerut de noi, dar incearca sa returneze unul cat mai
%apropiat
if ( err_min < 1e-3) 
    disp('Eroarea este minima');
end

%% Cerinta 9 - Sistemul pe spatiul starilor
[A_lin, B_lin, C_lin, D_lin] = linmod(mdlinout, xstar8, ustar8);
sys = ss(A_lin, B_lin, C_lin, D_lin);


%% Cerinta 10 - Verificare sistem stabil si salvare spectru
% daca A are toate valorile proprii in C minus atunci sistemul este stabil
vp = eig(A_lin);
if all(vp < 0)
    disp("Sistemul este stabil");
else
    disp("Sistemul nu este stabil");
end

%% Cerinta 11 - Raspunsul modelului liniarizat

r = ustar(1);
t = linspace(0,100,10000).';
u = r.*double(t>=0);
usim = timeseries(u,t);

mdl_lin = 'modelLiniarizat';
load_system(mdl_lin);
sim(mdl_lin);

%% Cerinta 12 - Eroarea dintre sistemul liniarizat si cel neliniar

ernllin = norm(y_nl.data - y_lin.data, Inf);
disp("Eroarea dintre sistemul liniarizat si cel neliniar:")
disp(ernllin)
%Observam o eroare foarte mica ༼ つ ◕_◕ ༽つ

%% Cerinta 13 - Discretizarea sistemului
%Metodele prezentate si cunoscute de noi in cadrul laboratorului sunt zoh
%si tustin.
%Deoarece nu ni se impun cerinte legate de robustete/stabilitate sau
%pastrarea anumitor frecvente putem alege o metoda care sa ne convina mai
%mult din punct de vedere teoretic.
%zoh reprezinta o metoda pentru uz general, simpla, ce poate fi usor de
%inteles dar nu pastreaza la fel de multe informatii ca tustin
%tustin implementeaza transformata tustin, o metoda mai complicata decat
%zoh, dar permite retinerea a cat mai multe informatii.
%Mie mi se pare mai interesanta metoda tustin pentru ca pastreaza
%caracteristici importante in frecventa si ne asigura si faptul ca sistemul
%va ramane stabil ( deoarece a fost stabil la cerinta 10 ).
%deci o voi alege pe aceasta (●'◡'●)
Te = 0.1;
[numitor,numarator] = ss2tf(A_lin,B_lin,C_lin,D_lin);
SysLin = tf(numitor,numarator);
SysLinD = c2d(SysLin,Te,'tustin');

%% Cerinta 14 -
%Ecuatia cu diferente (pentru vizualizare)
equation = sprintf('y[k]=(%.10f)*((%.10f)*u[k]',1/SysLinD.den{1}(1),SysLinD.num{1}(1));
for i = 2:numel(SysLinD.num{1})
    equation = sprintf('%s + (%.10f)*u[k-%d]',equation,SysLinD.num{1}(i), i-1);
end
for i = 2:numel(SysLinD.den{1})
    equation = sprintf('%s + (%.10f)*y[k-%d]',equation,-SysLinD.den{1}(i), i-1);
end
equation = sprintf('%s)',equation);
disp('Ecuatia cu diferente:');
disp(equation);

[Adisc,Bdisc,Cdisc,Ddisc] = tf2ss(SysLinD.num{1},SysLinD.den{1});
mdl14 = 'modelul2cerinta14';
load_system(mdl14)
set_param(mdl14,'StopTime',num2str(T));
sim(mdl14);

%% Cerinta 15
figure('Name','grafice pentru sistem liniar/discret');
plot(y_discret.Time,y_discret.Data,'r');
hold on;
plot(y_neliniar.Time,y_neliniar.Data,'b');
legend('Raspuns al sistemului discretizat', 'Raspuns al sistemului neliniar');
xlabel('Timp (s)');
ylabel('Theta (rad)')
title('Diferenta sistem neliniar - discretizat');
hold off;