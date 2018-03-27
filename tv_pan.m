function x = tv_pan(f,iter,a,d,eps)
%TV
 [s1,s2]=size(f);
%  eps=1e-8;
%  a=0.2;
%  d=norm(f(:)-lastf(:));
 for ntv=1:iter                  
     %extend
     fe=[zeros(1,s2+2);[zeros(s1,1) f zeros(s1,1)];zeros(1,s2+2)];
     part1=((f-fe(1:s1,2:s2+1))+(f-fe(2:s1+1,1:s2)))...
         ./sqrt((f-fe(1:s1,2:s2+1)).^2+(f-fe(2:s1+1,1:s2)).^2+eps);
     part2=(-f+fe(3:s1+2,2:s2+1))...
         ./sqrt((-f+fe(3:s1+2,2:s2+1)).^2+(fe(3:s1+2,2:s2+1)-fe(3:s1+2,1:s2)).^2+eps);
     part3=(-f+fe(2:s1+1,3:s2+2))...
         ./sqrt((-f+fe(2:s1+1,3:s2+2)).^2+(fe(2:s1+1,3:s2+2)-fe(1:s1,3:s2+2)).^2+eps);
     v=part1-part2-part3;
     v=v/norm(v(:));
     x=f-a*d*v;
 end;
