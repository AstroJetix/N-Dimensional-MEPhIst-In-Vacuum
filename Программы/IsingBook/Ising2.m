function [Es,Ed,SpM,A,S] = Ising2(Nspin,J,h,Esi,NTrial)
% функция, возвращающая мгновенные 
% значения: 
% энергии системы (Es), 
% энергии демона (Ed), 
% намагниченности (SpM), 
% числа принятия решений (A),
% а также мгновенные конфигурации спинов (S)
% Nspin - число спинов решетки
%  J - константа обменного взаимодействия
% h -  напряженность внешнего магнитного поля
% Esi - конечная энергия системы
% NTrial - число испытаний
% Copyright 2003-2009 S. V. Porshnev

Ns=Nspin.^0.5; % число спинов вдоль 
               % одной стороны
s=ones(Ns,Ns); % начальная конфигурация 
               % спинов
Esystem=-(J+h)*Nspin; % начальная энергия 
                      % системы
% начальная энергия демона
Edemon=4*J*floor((Esi-Esystem)/(4*J)); 

Es(1)=Esystem; 
Ed(1)=Edemon;

S=s;
k=1;
for i=1:NTrial
  Accept=0;
  for j=1:Nspin
     % случайный выбор узла сетки
     Ix=floor(Ns*rand(1)+1); 
     Iy=floor(Ns*rand(1)+1);
     % граничные условия
     if Ix==1
        Left=Ns;
     else
        Left=Ix-1;
     end
     if Ix==Ns
       Right=1;
     else;
        Right=Ix+1;
     end
     if Iy==1
        Down=Ns;
     else;  
        Down=Iy-1;
     end
     if Iy==Ns
        Up=1;
     else
        Up=Iy+1;
     end
     % пробное изменение энергии
     de=2*s(Iy,Ix)*(-h+J*(s(Iy,Left)+s(Iy,Right)+s(Down,Ix)+s(Up,Ix)));
     if de<=Edemon % принятие пробного 
                   % изменения энергии
        s(Iy,Ix)=-s(Iy,Ix);
        Accept=Accept+1;
        Edemon=Edemon-de;
        Esystem=Esystem+de;
     end
     k=k+1;
     Es(k)=Esystem;
     Ed(k)=Edemon;
     A(k-1)=Accept;
     s1=sum(s);
     SpM(k)=sum(s1);
     S=cat(3,S,s);
  end
end
A=A/NTrial;
