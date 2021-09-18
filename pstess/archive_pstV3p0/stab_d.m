% m file for stabilizer design
% requires system data from svm_smse
function [tw,t1,t2,t3,t4]=d_stab(a,b,c,tstate)
[f,yphase]=stabf(a,b,c,tstate);
st_des=input('Do you wish to do a stabilizer design?Y/N[Y]','s')
if isempty(st_des);st_des='Y';end;
if st_des=='y';st_des='Y';end
if st_des=='Y';
   stab_name=input('What is the name of the generator?','s')
   point=length(f)
   ideal_phase=-yphase;
   st_type=input('Has the stabilizer power or speed input,P/S[S]:','s');
   if isempty(st_type);st_type='S';end
   if st_type=='s';st_type='S';end
   if st_type~='S'
      new_stab='Y';
      while(new_stab=='Y')
         %power input stabilizer
         [phase,tw,t1,t2,t3,t4]=pss_phse(f);
         phase=phase+90.*ones(1,point);
         subplot(111),plot(f,phase,f,ideal_phase);
         grid
         tit=['Ideal and ',stab_name,' pss phase'];
         title(tit)
         new_stab=input('Do you wish to try another stabilizer design,Y/N[N]:','s');
         if isempty(new_stab);new_stab='N';end
         if new_stab=='y';new_stab='Y';end
      end
   else
      new_stab='Y';
      while(new_stab=='Y')
         %speed input stabilizer
         [phase,tw,t1,t2,t3,t4]=pss_phse(f);
         subplot(111), plot(f,phase,f,ideal_phase);
         grid
         tit=['Ideal and ',stab_name,' phase'];
         title(tit)
         new_stab=input('Do you wish to try another stabilizer design,Y/N[N]:','s');
         if isempty(new_stab);new_stab='N';end
         if new_stab=='y';new_stab='Y';end
      end
   end
end
return