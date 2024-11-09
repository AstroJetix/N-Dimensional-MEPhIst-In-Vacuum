function [V,Es,Accept] = IdealGas(Vel,T,dVmax,NTrial)
% функция, возвращающая 
% мгновенные значения скорости (V),
% энергии (E), 
% среднее число принятия (Accept)
% Vel - начальная скорость
% T - температура газа
% dVmax - максимальное изменение скорости
% NTrial - число испытаний
% Copyright 2003-2009 S. V. Porshnev

beta=1/T;
E=Vel.^2/2; % начальная энергия
Accept=0;

for i=1:NTrial
   dV=(2*rand(1)-1)*dVmax;
   Vtrial=Vel+dV; % случайное изменение 
                  % скорости
   de=0.5*(Vtrial.^2-Vel.^2); % пробное 
                              % изменение энергии
   if de>0 
     if exp(-beta*de)>=rand(1)
        % шаг принимается
       Vel=Vtrial;
       Accept=Accept+1;
       E=E+de;
     end
   else
      % шаг не принимается
      Vel=Vtrial; 
      Accept=Accept+1;
      E=E+de;
   end
   Es(i)=E;
   V(i)=Vel;
end
Accept=Accept/NTrial;