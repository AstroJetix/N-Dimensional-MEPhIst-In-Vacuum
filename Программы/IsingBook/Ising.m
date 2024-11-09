function [Es,Ed,SpM,Accept,S]=Ising(Nspin,J,h,NTrial,T)
% �������, ������������ ��������:
% ������� ������� (Es),
% ������� ������ (Ed),
% ���������� ������� ������� (SpM),
% �������� ������� (Accept),
% ���������� ������������ ������ (S)
% Nspin -����� ������
% J - ��������� ��������� ��������������
% h - ������������� �������� ���������� ����
% NTrial - ����� ���������
% T - ����������� �������

beta=1/T;
s=ones(Nspin);
S(1)=s;
M=Nspin;
Esystem=-(J+h)*Nspin;
Edemon=2*J*ceil((Esi-Esystem)/(2*J));
Es(1)=Esystem;
Ed(1)=Edemon;
SpM(1)=M;
Accept=0;
k=1;

for i=1:NTrial
  for j=1:Nspin
    Ispin=floor(Nspin*rand(1)+1); % ���������
                                  % ����� �����
    if Ispin==1
      Left=s(Nspin);
    else
      Left=s(Ispin-1);
    end  
    if Ispin==Nspin
      Right=s(1);
    else
      Right=s(Ispin+1);
    end
    de=2*s(Ispin)*(-h+J*(Left+Right));% ������� 
                                      % ��������� 
                                      % �������
    
    if (de<=0) | (exp(-beta*de)>rand(1))
      % ��������� �����������          
      s(Ispin)=-s(Ispin);
      Accept=Accept+1;
      Esystem=Esystem+de;
    end
    k=k+1;
    S=cat(3,S,s);
    Es(k)=Esystem/Nspin;
    Ed(k)=Edemon;
    SpM(k)=sum(s)/Nspin;
  end
end
Accept=Accept/(NTrial*Nspin);
