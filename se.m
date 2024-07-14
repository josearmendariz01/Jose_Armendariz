function [c,eigenvalor,coef]=se(r,q,z)
%Calcula la función se de Mathieu de orden r, en función del parámetro q y el argumento z.
%USA METODO MATRICIAL
%Las funciones se con r=par tienen periodo pi
%Las funciones se con r=impar tienen periodo 2*pi
%Ec. de Mathieu:   y"+(a-2qcos2z)y=0
%FORMAT  see(r,q,z)
%   r=1,2,3,.. = orden integral del valor característico.
%   q=parámetro real o complejo
%   z=valor de la variable independiente en general compleja
%SALIDA [c,eigenvalor,coef]
%c = funcion de Mathieu
%eigenvalor = r-esimo eigenvalor de la funcion de Mathieu
%coef = conjunto de coeficientes de Fourier de la serie
% Programado por Julio C. Gutierrez-Vega

if r<=0
   error('Order must be real and postive')
end

if mod(r,2)==0
   N=r+50;  %Tamaño de la matriz para calcular los eigenvalores
   M = diag([2:2:2*(N+1)].^2,0) + diag(linspace(q,q,N),1) + diag(linspace(q,q,N),-1);  %Eigen matriz
   
   [A,ets]=eig(M);
   ets=diag(ets); 
   pos = r/2 ;   %Para el r'esimo orden calcula la posici´on 
   
%    [ets,index]=sort(ets);  %Ordena los eigenvalores usando la MAGNITUD
%    eigenvalor = ets(pos);  %eigenvalor para el orden r
%    coef = transpose(A(:,pos));   %coeficientes de Fourier para el orden r;  arreglados en un renglon

   [ets1,index]=sort(real(ets));  %Ordena los eigenvalores usando la parte REAL
   eigenvalor = ets(index(pos));  %eigenvalor para el orden r
   coef = transpose(A(:,index(pos)));   %coeficientes de Fourier para el orden r;  arreglados en un renglon
   
   coef = coef/sqrt(sum(coef.^2)) / sign(coef(1));  %Normalizacion y ajuste de signo
   
   %SERIE Solucion
   c=0;
   for ii=1:length(coef)
      c=c + coef(ii)*sin(2*(ii)*z);
   end
   
else
   N=r+25;  %Tamaño de la matriz para calcular los eigenvalores
   M = diag([1-q,[3:2:2*N+1].^2],0) + diag(linspace(q,q,N),1) + diag(linspace(q,q,N),-1);  %Eigen matriz
   
   [A,ets]=eig(M);
   ets=diag(ets);
   pos = (r+1)/2;   %Para el r'esimo orden calcula la posici´on 
   
%    [ets,index]=sort(ets);  %Ordena los eigenvalores usando la MAGNITUD
%    eigenvalor = ets(pos);  %eigenvalor para el orden r
%    coef = transpose(A(:,pos));   %coeficientes de Fourier para el orden r;  arreglados en un renglon

   [ets1,index]=sort(real(ets));  %Ordena los eigenvalores usando la parte REAL
   eigenvalor = ets(index(pos));  %eigenvalor para el orden r
   coef = transpose(A(:,index(pos)));   %coeficientes de Fourier para el orden r;  arreglados en un renglon
   
   coef = coef/sqrt(sum(coef.^2)) / sign(coef(1));  %Normalizacion y ajuste de signo
   
   c=0;
   for ii=0:(length(coef)-1)
      c=c + coef(ii+1)*sin((2*ii+1)*z);
   end
end
