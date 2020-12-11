format short
clear all
clc


Ufront= @(x,y)     exp(x)*sin(y); 
Qfront= @(x,y)     exp(x)*sin(y)*x+ exp(x)*cos(y)*y; 
const=-2;

NN=16; % Número de nodos en la frontera.
NE=16; % Número de elementos en la frontera.
L=17; % Número de nodos internos.
N=NN; % 
NCI=48; % Número de celdas internas (unicamente en el cado de Uxx+Uyy=b)
NPI=7; % Número de puntos de integración para las celdas internas.
NI=4; % Número de puntos de integración para los elementos de la frontera.


% Creación de vectores y matrices para guardar datos.
X=zeros(1,N);
Y=zeros(1,N);
U=zeros(1,N);
Q=zeros(1,N);



D=zeros(1,N);
LE=zeros(1,N);
CON=zeros(N,2);
KODE=zeros(1,N);



alpha=zeros(1,N);
XY=zeros(1,N);
A=zeros(N,N);
MKJ=zeros(N,3);


DA=zeros(1,N);
CP=zeros(7,3);

% Contantes para la integración sobre elementos
POEI=[0.86113631,-0.86113631,0.33998104,-0.33998104];
FDEP=[0.34785485,0.34785485,0.65214515,0.65214515];

% Datos para la integración sobre celda
WW=[0.225,0.12592918,0.12592918,0.12592918,0.13239416,0.13239416,0.13239416];
CP=[0.33333333,0.79742699,0.10128651;0.10128651,0.05971587,0.4701426;0.47014206,0.33333333,0.10128651;0.79742699,0.10128651,0.47014206;0.05971587,0.47014206,0.33333333;0.10128651,0.10128651,0.79742699;0.47014206,0.47014206,0.05971587];

% Datos para los nodos en la frontera (LOS PUNTOS DEBEN DE COLOCARSE
% ANTIHORARIO
x=[2,1.705706,1.178800,0.597614,0,-0.597614,-1.178800,-1.705706,-2,-1.705706,-1.178800,-0.597614,0,0.597614,1.178800,1.705706];
y=[0,-0.522150,-0.807841,-0.954310,-1,-0.954310,-0.807841,-0.522150,0,0.522150,0.807841,0.954310,1,0.954310,0.807841,0.522150];


% Datos para los nodos internos
xint=[1.5,1.2,0.6,0,-0.6,-1.2,-1.5,-1.2,-0.6,0,0.6,1.2,0.9,0.3,0,-0.3,-0.9];
yint=[0,-0.35,-0.45,-0.45,-0.45,-0.35,0,0.35,0.45,0.45,0.45,0.35,0,0,0,0,0];


if const==0
    NCI=0;
%     disp('Ecuación de Laplace Uxx=Uyy')
else
    NCI=NCI;
%     disp('Ecuación de Poisson Uxx+Uyy=b')
end

for i=1:NN  % Conectando los elementos de la frontera NODO A NODO.
    
    CON(i,1)=i;
    CON(i,2)=i+1;
    CON(NN,2)=1;
    z(i)=i; % Nodos de las fronteras. 
end

for i=1:L % Agrega los nodos internos desde el nodo NN.
    zint(i)=NN+i;
end


for i=1:NE
    length(i)=sqrt((x(CON(i,2))-x(CON(i,1)))^2+(y(CON(i,2))-y(CON(i,1)))^2);
end


tabla2=[z',x',y']'; % Crea los datos de la tabla "Data for boundary nodes"
tabla1=[z',CON,length']'; % Crea la tabla "Element data"
tabla3=[zint',xint',yint']';


fprintf('%6s %12s %12s\n','Nodo','Xfront','Yfront');
fprintf('%6.0f %12.4f %12.4f\n',tabla2);
disp('--------------------------------------------------')
fprintf('%6s %11s %12s\n','Nodo',' X int','Y int');
fprintf('%6.0f %12.4f %12.4f\n',tabla3);
disp('--------------------------------------------------')
fprintf('%6s %10s %12s %12s\n','Elemento','Nodo 1','Nodo 2','Longitud');
fprintf('%6.0f %12.4f %12.4f %12.6f\n',tabla1);


for i=NN+1:NN+L
    
    KODE(i)=0;
    U(i)=0;
    Q(i)=0;
       
end




plot(x,y,'ko')
hold on 
plot(xint,yint,'bo')
grid on