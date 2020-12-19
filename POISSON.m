format short
clear all
clc


 
const=-2;

NN=16; % Número de nodos en la frontera.
NE=16; % Número de elementos en la frontera.
L=17; % Número de nodos internos.
NCI=48; % Número de celdas internas (unicamente en el cado de Uxx+Uyy=b)
NPI=7; % Número de puntos de integración para las celdas internas.
NI=4; % Número de puntos de integración para los elementos de la frontera.
VAL=2; % % Selecciona el tipo de condición de frontera VAL=1- q conocido, q desconocido; VAL=2- q conocido, u desconocido.


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


% Nodos para realizar las particiones triangulares del método

MKJ=[2,1,17;16,17,1;9,8,23;9,23,10;28,17,16;18,2,17;24,10,23;22,23,8;15,28,16;3,2,18;7,22,8;11,10,24;29,17,28;29,18,17;23,22,33;23,33,24;27,28,15;19,3,18;22,7,21;11,24,25;27,15,14;19,4,3;21,7,6;25,12,11;27,29,28;19,18,29;22,21,33;33,25,24;27,30,29;30,19,29;25,33,32;33,21,32;26,27,14;20,4,19;21,6,20;12,25,26;26,30,27;20,19,30;21,20,32;32,26,25;13,26,14;5,4,20;5,20,6;26,13,12;31,30,26;31,20,30;31,32,20;31,26,32];


if const==0
    NCI=0;
%     disp('Ecuación de Laplace Uxx=Uyy')
else
    NCI=NCI;
%     disp('Ecuación de Poisson Uxx+Uyy=b')
end

for I=1:NN  % Conectando los elementos de la frontera NODO A NODO.
    
    CON(I,1)=I;
    CON(I,2)=I+1;
    CON(NN,2)=1;
    zf(I)=I; % Nodos de las fronteras.
    
end


for I=NN+1:NN+L
    
    KODE(I)=0;
    U(I)=0;
    Q(I)=0;          
     
end


%--------------------------------------------------------------
for i=1:L % Agrega los nodos internos desde el nodo NN.
    zint(i)=NN+i;
end

% READ COORDINATES OF NODES

z=[zf';zint']; % Une los vectores en uno solo (frontera+internos)
X=[x';xint'];
Y=[y';yint'];
XYZ=[z,X,Y]'; % Vector que imprime todos los nodos y sus respectivas coordenadas.

%--------------------------------------------------------------

% Cálculo del área de las celdas internas

for J=1:NCI
    
    I1=MKJ(J,1);
    I2=MKJ(J,2);
    I3=MKJ(J,3);
    AB1=Y(I2)-Y(I3);
    AB2=Y(I3)-Y(I1);
    AB4=X(I3)-X(I2);
    AB5=X(I1)-X(I3);
    DA(J)=abs((AB1*AB5-AB2*AB4)/2);
    
end


for II=1:NN  % Selecciona el tipo de condición de frontera 1- q conocido, q desconocido; 2- q conocido, u desconocido.
    
    KODE(II)=VAL;
    
    if VAL==1
    
        U(II)=VAL;

    end
    
    if VAL==2
        
        Q(II)=VAL;

    end
end
    

% Cálculo de la longitud de cada segmento de recta que une a cada nodo.

for k=1:NE
    
    LE(k)=sqrt((x(CON(k,2))-x(CON(k,1)))^2+(y(CON(k,2))-y(CON(k,1)))^2);

end


% SUBRUTINA "ASSEM2"


for J1=1:NN+L  % Módulo para ensamblar la matriz H y G paa reducirlo a la condición AX=Y.
    
   XI=X(J1);
   YI=Y(J1);
   XY(J1)=0;

   for J=1:NN 
       
   HH(J1,J)=0;
   A(J1,J)=0;
   CC=0;
   
   for J2=1:NE
       
   LJ=LE(J2);
   N1=CON(J2,1);
   N2=CON(J2,2);
   X1=X(N1);
   X2=X(N2);
   Y1=Y(N1);
   Y2=Y(N2);
   H1=0;
   H2=0;
   G1=0;
   G2=0;
   
    for J3=1:4  % Computo de los coeficientes G y H  
          
            
        E = POEI(J3);
        w = FDEP(J3);
        XX = X1 + (1+E)*(X2-X1)/2;
        YY = Y1 + (1+E)*(Y2-Y1)/2;
        R = sqrt((XI-XX).^2 + (YI-YY).^2);
        PP = (((XI-XX).*(Y1-Y2)+(YI-YY).*(X2-X1))/(R.*R*4*pi))*w;

        % Elementos de H y G usando las formulas 2.4647 , 2.49-50
        
        H1=H1+((1-E)/2)*(PP);
        H2=H2+((1+E)/2)*(PP);
        PP=log(1/R)/(4*pi)*LJ*w;
        G1=G1+((1-E)/2)*(PP);
        G2=G2+((1+E)/2)*(PP);

         
        CC=CC-H1-H2;
        GE=LJ.*((3/2)-log(LJ))./(4.*pi); % Integración exacta para G cuando I y J son el mismo elemento.
    end
   end
   end

    
    
    if N1==J1
        % SE APLICA UNICAMENTE AL MÉTODO DRM.
        
        G1=GE;
        
    end
    
    if N2==J1
        
        G2=GE;
        
    end
        HH(J1,N1)=HH(J1,N1)+H1;
        HH(J1,N2)=HH(J1,N1)+H2;
        GG(J1,2*J2-1)=G1;
        GG(J1,2*J2)=G2;  % ###########################################################3
    
end 
% % #################################################
% % HASTA ACA ESTOY SEGURO DE QUE ESTA PERFECTO IGUAL QUE EL LIBRO
% % GUIA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     % ######################################################### DUDA
%      
%     
%    end
%    end
% end
% 
for J1=1:NN

    KK=KODE(N1);
     
    if KK==0 | KK==2
        
        XY(J1)=XY(J1)+Q(N1)*G1;
        A(J1,N1)=A(J1,N1)+H1;
    
    else
        XY(J1)=XY(J1)-U(N1)*H1;
        A(J1,N1)=A(J1,N1)-H1;
    end

  
    KK=KODE(N2);
     
    if KK==0 | KK==2
        
        XY(J1)=XY(J1)+Q(N2)*G2;
        A(J1,N2)=A(J1,N2)+H2;
    
    else
        XY(J1)=XY(J1)-U(N2)*H2;
        A(J1,N2)=A(J1,N2)-G2;
    end

    %
%     
HH(J1,J1)=CC;  % Apply boundary conditions to diagonal term
KK=KODE(J1);
% %         
%     
    if KK==0 | KK==2
        
        
        A(J1,J1)=CC;
        
    else
        
        XY(J1)=XY(J1)-U(J1)*CC;
    end
    
    end
   
% 
% % SUBRUTINA NECMOD
% 
for I=1:NN+L
    
   XI=X(I);
   YI=Y(I);
   FF=0;
   
   for J1=1:NCI
       
      AR=DA(J1);
      N1=MKJ(J1,1);
      N2=MKJ(J1,2);
      N3=MKJ(J1,3);
       
      X1=X(N1);
      X2=X(N2);
      X3=X(N3);
      Y1=Y(N1);
      Y2=Y(N2);
      Y3=Y(N3);
      
      
      for J3=1:NPI
         
         C1=CP(J3,1);
         C2=CP(J3,2);
         C3=CP(J3,3);
          
         % XP, YP COODINATES OF INTEFRATION POINT 
         
         XP=X1*C1+X2*C2+X3*C3;
         YP=Y1*C1+Y2*C2+Y3*C3;
         
         % DX, DY COMPONENTS OR R
         
         
         DX=XI-XP;
         DY=YI-YP;
         R=sqrt(DX*DX+DY*DY);
         
         
         % SUMATORIA SOBRE LAS CELDAS USANDO LA FORMULA 2.78.
         
         FF=FF-WW(J3)*AR*log(1/R)*(const/(2*pi));
         D(I)=FF;
         
         for j=1:NN
             
         XY(j)=XY(j)+FF; % ##########################################
      end
      end
   end
    
end
% 
% 
% % SUBRUTINA SOLVER
% 
for M=1:NN-1
    
   AM=(1/A(M,M));
   


   for J=M+1:NN
       
       CD=A(M,J);
       A(M,J)=CD*AM;
    
   end 
   
       for N=M+1:NN
       
       ALM=A(N,M);
       
      
       for KK=M+1:NN
           
           A(N,KK)=A(N,KK)-A(M,KK)*ALM;
           
       end
       end
end


% SUBRUTINA SOLVER

        for M=1:NN-1
            
            AM=1/A(M,M);
            XY(M)=XY(M)*AM;
            
            for N=M+1:NN
                
                ALM=A(N,M);
                XY(N)=XY(N)-XY(M)*ALM;
                XY(NN)=XY(NN)/A(NN,NN);
                NI=NN;
                NI=NI-1;
                
                
                if NI~=0
                    
                    for J=NI+1:NN
                   
                        XY(NI)=XY(NI)-A(NI,J)*XY(J);
                        
                    end 
                end
                
            end
        end

    
for I=I:NN
    
    KK=KODE(I);
    
    if KK==0 | KK==2
       U(I)=XY(I);
       
    else
        
        Q(I)=XY(I);
        
    end
end

% SUBRUTINA INTERM

for I=NN+1:NN+L
    
   U(I)=0;
   
   for K=1:NN-1
   
       U(I)=U(I)+(GG(I,2*K)+GG(I,(2*K+1)))*Q(K+1);
       
       U(I)=U(I)+(GG(I,1)+GG(I,(2*NE)))*Q(1);
       
       for K1=1:NN
          
           U(I)=U(I)-HH(I,K1)*U(K1);
           U(I)=U(I)+D(I);
           
       end
   end
    
end

tabla2=[zf',x',y']'; % Crea los datos de la tabla "Data for boundary nodes"
tabla1=[zf',CON,LE']'; % Crea la tabla "Element data"
tabla4=[MKJ,DA']';
tabla3=[zint',xint',yint']';



fprintf('%6s %12s %12s\n','Nodo','Xfront','Yfront');
fprintf('%6.0f %12.4f %12.4f\n',tabla2);
disp('--------------------------------------------------')
fprintf('%6s %11s %12s\n','Nodo',' X int','Y int');
fprintf('%6.0f %12.4f %12.4f\n',tabla3);
disp('--------------------------------------------------')
fprintf('%6s %12s %12s %12s\n','Elemento','Nodo 1','Nodo 2','Longitud');
fprintf('%6.0f %12.0f %12.0f %12.6f\n',tabla1);
% disp('--------------------------------------------------')
fprintf('%8s %8s %12s\n','Nodo','X','Y');
fprintf('%6.0f %12.4f %12.4f\n',XYZ);
disp('--------------------------------------------------')
fprintf('%8s %12s %12s %12s\n','Nodo 1','Nodo 2','Nodo 3','Área celda');
fprintf('%8.0f %12.0f %12.0f %12.4f\n',tabla4);
disp('--------------------------------------------------')
fprintf('%10s %8s %12s\n','Nodo','U','Q');
fprintf('%8.0f %12.4f %12.4f\n',[z,U',Q']');


plot(x,y,'ko')
hold on 
plot(xint,yint,'bo')
grid on
