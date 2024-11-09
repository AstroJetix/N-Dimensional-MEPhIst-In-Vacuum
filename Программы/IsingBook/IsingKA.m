function [E,Accept,M,S] = IsingKA(Nspin,J,h,NTrial,T,Sstart,Mstart,Estart)
% функция, возвращающая 
% мгновенные значения энергии системы (E), 
% число принятия (Accept), 
% мгновенный магнитный момент системы (M), 
% мгновенные конфигурации спинов системы (S)
% Nspin - число спинов системы
% J - константа обменного взаимодействия
% h - напряженность внешнего магнитного поля
% NTrial - число испытаний
% T - температура системы
% Sstart - матрица, содержащая информацию 
% об ориентации спинов в 
% момент времени t = 0
% Mstart - значение магнитного момента 
% системы в момент времени t = 0
% Estart - энергия системы в момент времени t = 0
% Copyright 2003-2009 S. V. Porshnev

Ns=Nspin.^0.5; 
M = zeros(NTrial * Nspin + 1, 1);
E = zeros(NTrial * Nspin + 1, 1);
M(1)=Mstart;
S=Sstart;
E(1)=Estart;
s=Sstart;
Accept=0;
k=2;
for i=1:NTrial
   for j=1:Nspin
      % случайный выбор номера узла 
      % для опрокидывания спина
      Ix=floor(Ns*rand(1)+1); 
      Iy=floor(Ns*rand(1)+1);
      % проверка периодических 
      % граничных условий
      if Ix==1
         Left=Ns;
      else
         Left=Ix-1;
      end
      if Ix==Ns
         Right=1;
      else
         Right=Ix+1;
      end
      if Iy==Ns
         Up=1;
      else
         Up=Iy+1;
      end
      if Iy==1
         Down=Ns;
      else
         Down=Iy-1;
      end
      % пробный переворот спина
      Temp=s(Iy,Ix);
      s(Iy,Ix)=-s(Iy,Ix);
      de=2*s(Iy,Ix)*(h+J*(s(Iy,Left)+s(Iy,Right)+s(Down,Ix)+s(Up,Ix)));
      if or(de<=0,rand(1)<=exp(-de./T)) 
         % пробный переворота спина принимается
         Accept=Accept+1;
         E(k)=E(k-1)-de;
         M(k)=M(k-1)+2*s(Iy,Ix);
      else
         % пробный переворот спина 
         % не принимается 
         s(Iy,Ix)=Temp;
         E(k)=E(k-1);
         M(k)=M(k-1);
      end
      %S=cat(3,S,s);
      k=k+1;
   end
end
Accept=Accept/(NTrial*Nspin); % среднее число 
S = s;                                       % принятия
