function [V,Es,Accept] = IdealGas(Vel,T,dVmax,NTrial)
% �������, ������������ 
% ���������� �������� �������� (V),
% ������� (E), 
% ������� ����� �������� (Accept)
% Vel - ��������� ��������
% T - ����������� ����
% dVmax - ������������ ��������� ��������
% NTrial - ����� ���������
% Copyright 2003-2009 S. V. Porshnev

beta=1/T;
E=Vel.^2/2; % ��������� �������
Accept=0;

for i=1:NTrial
   dV=(2*rand(1)-1)*dVmax;
   Vtrial=Vel+dV; % ��������� ��������� 
                  % ��������
   de=0.5*(Vtrial.^2-Vel.^2); % ������� 
                              % ��������� �������
   if de>0 
     if exp(-beta*de)>=rand(1)
        % ��� �����������
       Vel=Vtrial;
       Accept=Accept+1;
       E=E+de;
     end
   else
      % ��� �� �����������
      Vel=Vtrial; 
      Accept=Accept+1;
      E=E+de;
   end
   Es(i)=E;
   V(i)=Vel;
end
Accept=Accept/NTrial;