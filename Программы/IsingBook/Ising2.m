function [Es,Ed,SpM,A,S] = Ising2(Nspin,J,h,Esi,NTrial)
% �������, ������������ ���������� 
% ��������: 
% ������� ������� (Es), 
% ������� ������ (Ed), 
% ��������������� (SpM), 
% ����� �������� ������� (A),
% � ����� ���������� ������������ ������ (S)
% Nspin - ����� ������ �������
%  J - ��������� ��������� ��������������
% h -  ������������� �������� ���������� ����
% Esi - �������� ������� �������
% NTrial - ����� ���������
% Copyright 2003-2009 S. V. Porshnev

Ns=Nspin.^0.5; % ����� ������ ����� 
               % ����� �������
s=ones(Ns,Ns); % ��������� ������������ 
               % ������
Esystem=-(J+h)*Nspin; % ��������� ������� 
                      % �������
% ��������� ������� ������
Edemon=4*J*floor((Esi-Esystem)/(4*J)); 

Es(1)=Esystem; 
Ed(1)=Edemon;

S=s;
k=1;
for i=1:NTrial
  Accept=0;
  for j=1:Nspin
     % ��������� ����� ���� �����
     Ix=floor(Ns*rand(1)+1); 
     Iy=floor(Ns*rand(1)+1);
     % ��������� �������
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
     % ������� ��������� �������
     de=2*s(Iy,Ix)*(-h+J*(s(Iy,Left)+s(Iy,Right)+s(Down,Ix)+s(Up,Ix)));
     if de<=Edemon % �������� �������� 
                   % ��������� �������
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
