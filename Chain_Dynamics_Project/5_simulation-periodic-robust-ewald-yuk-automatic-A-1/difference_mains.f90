5c5
<   integer :: nosteps=47400!99900!1024!2048!4096!8192!16384!32768
---
>   integer :: nosteps=50000!99900!1024!2048!4096!8192!16384!32768
239,240c239
< 
<   !the underlying calls are dirty, make them robust
---
>   
248c247
<   V2=0.d0!!!! 
---
> !!!V2=0.d0!!!! 
255c254
<   !v_s = 0.d0
---
>   v_s = 0.d0
425c424
<   real*8 :: r_yuk,fac,f0
---
>   real*8 :: r_yuk,fac
428,429d426
< 
<   f0=2.d0 !correction for vaue of quadrupole at 1
431c428
<   fac = f0*((-ka/(r_yuk**m_yuk))-(m_yuk/(r_yuk**(m_yuk+1))))
---
>   fac = ((-ka/(r_yuk**m_yuk))-(m_yuk/(r_yuk**(m_yuk+1))))
566c563
<   write(*,*)'From generate_af:Lx=',Lx,'Ly=',Ly
---
>   write(*,*)'Lx=',Lx,'Ly=',Ly
