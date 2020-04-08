_IF1() 
_IF1()  $Author: jmht $
_IF1()  $Date: 2008-10-06 23:36:09 +0200 (Mon, 06 Oct 2008) $
_IF1()  $Locker:  $
_IF1()  $Revision: 5733 $
_IF1()  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/assem.m,v $
_IF1()  $State: Exp $
_IF1()  
_IF(cray)
_IFN(c90)
_IF(unicos)
          ident revpri
_ENDIF
_IFN(unicos)
          ident openv
_ENDIF
*
*   coded      v. r. saunders    july 1982
*
  common disc
incr bss 1
npb bss 1
npbsz bss 1
irep bss 1
npbm1 bss 1
ipos bss 40
nam bss 40
  common discne
i40 bss 40
length bss 40
iostat bss 40
iobuff bss 40
itab bss 1600
istabf bss 10
ipribf bss 9
junibf bss 9
jpbbf bss 9
junjpb bss 9
ibfbas bss 9
igetm bss 9
iputm bss 9
mbuff bss 1
kunkpb bss 1
lun bss 1
l40 bss 1
latm bss 1
lbl bss 1
lbfbas bss 1
lbuff bss 1
  block *
  block *
_IFN(unicos)
*
*   call openv(ftable,dsname,nblock)
*   ftable=dsp(24 words) + odn(2 words) + ddl (6 words)
*
* entry openv
openv enter np=3
  a5 le@dsp+le@odn+le@ddl
  argadd a1,1                 a1=<dsp>
  vl a5
  s7 >2
  v7 0
  argadd a5,2
  a0 a1
  s5 0,a5
  ,a0,1 v7
  a2 le@dsp+le@odn
  a3 le@dsp
  a2 a2+a1                     a2=<ddl>
  a3 a3+a1                     a3=<odn>
  0,a1 s5                      dsname to dsp
  0,a2 s5                      dsname to ddl
  0,a3 s5                      dsname to odn
  2,a2 s7                      random/unblocked flags to ddl
  1,a3 a1                      <dsp> to odn
  s1 a2
  s0 f$dnt                     request dnt entry
  a7 512
  s7 0
  ex
  close a3                     to lose an existing dsp
  put,s7 s5&s6,dperr,a1        clear dperr flags
  4,a1 a7                      to pacify open
  a6 irep
  open a3,io
  s7 1,a1
  s7 s7<2
  s0 s7>53
  s7 s7>53
  1,a1 s0
  0,a6 s7                      set reply in irep
  jsn ret
  argadd a1,3
  s1 a3
  s0 f$lbn
  s7 1
  ex
  s7 s7+s2
  0,a1 s7                   store file size in nblock
ret exit
*
*
*
* entry vclose
*
*   call vclose(ftable)
*
vclose enter np=1
  argadd a1,1
  s1 a1
  s7 le@dsp
  s0 f$cls
  s1 s1+s7
  ex
  exit
*
*
*
* entry fillv
*
*  call fillv(ftable,buffer,mblock,iblock)
*
fillv enter np=4
  argadd a4,4                 a4=<iblock>
  argadd a3,3                 a3=<mblock>
  argadd a1,1                 a1=<dsp>
  argadd a2,2                 a2=<intermediate>
  s2 a2                       s2=<buffer>
  s6 >40
  s7 0
  s4 0,a4              s4=iblock
  s3 0,a3              s3=mblock
  s1 1
  s5 1,a1
  3,a1 s2              set  dpout
  s4 s4-s1
  s3 s3<9             s3=buflen
  s4 s4<24
  s5 s5&s6
  s6 <51
  s3 s3+s2                 s3=dplmt
  7,a1 s7                  set dpbpd for read
  s5 s5&s6
  s1 a1
  s0 f$rdc
  4,a1 s3                  set dplmt
  s4 s4!s2                 s4=dpibn and dpin
  s5 s5!s2
  2,a1 s4                  set dpibn and dpin
  1,a1 s5                  set dplmt
  ex
  exit
*
*
*
* entry flushv
*
*  call flushv(ftable,buffer,mblock,iblock)
*
flushv enter np=4
  argadd a4,4                    a4=<iblock>
  argadd a3,3                    a3=<mblock>
  argadd a1,1                    a1=<dsp>
  argadd a2,2
  s2 a2                          s2=<buffer>
  s6 >40
  s7 >1
  s4 0,a4                      s4=iblock
  s3 0,a3                      s3=mblock
  s1 1
  s5 1,a1
  2,a1 s2                      set dpin
  s7 s7>4
  s4 s4-s1
  s3 s3<9                     s3=buflen
  s4 s4<24
  s5 s5&s6
  s6 <51
  s3 s3+s2                     s3=dplmt
  7,a1 s7                      set dpbpd to write
  s5 s5&s6
  s1 a1
  s0 f$wdc
  4,a1 s3                      set dplmt
  s4 s4!s2               s4=dpobn+dpout
  s5 s5!s2
  3,a1 s4                set dpobn and dpout
  1,a1 s5                set dpfrst
  ex
  exit
*
*
*
* entry chekrd
*
*  call chekrd(ftable,iblock)
*
chekrd enter np=2
  argadd a2,2                     a2=<iblock>
  argadd a1,1                     a1=<dsp>
  a7 irep
  s7 0,a2                         s7=iblock
  s6 1,a1                         s6=busy+err+dpfrst
  a5 2,a1
  s7 s7<9                         s7=iblock*512
  a6 s6
  s3 s6
  s5 a5                           s5=dpin
  s4 a6                           s4=dpfrst
  s3 s3<2
  s4 s4+s7                        s4=limit
  s3 s3>53
  s5 s5-s4
  s0 s5&s6
  s1 a1
  jsp exitr
loopr s0 f$rcl
  ex
  a5 2,a1
  s6 1,a1
  s5 a5
  s5 s5-s4
  s0 s5&s6
  jsm loopr
exitr 0,a7 s3                     set reply word
  exit
*
*
*
* entry chekwr
*
*  call chekwr(ftable)
*
chekwr enter np=1
  argadd a1,1
  a7 irep
  s0 1,a1
  s1 a1
  jsp exitw
loopw s0 f$rcl
  ex
  s0 1,a1
  jsm loopw
exitw s0 s0<2
  s0 s0>53
  0,a7 s0
  exit
*
_ENDIF
*    fortran equivalent for revpri
*
*       ix=ipribf(mbuff)
*       do 1 loop=1,npb
*       if(ix.gt.ipribf(loop))ipribf(loop)=ipribf(loop)+1
* 1     continue
*       ipribf(mbuff)=0
* entry revpri
revpri enter np=0
  a6 npb,0
  a7 mbuff,0
  a0 ipribf
  s1 1
  vl a6
  a5 a7-1
  v7 ,a0,1
  v6 s1+v7
  s2 v7,a5
  v5 s2-v6
  vm v5,p
  v4 v6!v7&vm
  v4,a5 0
  ,a0,1 v4
  exit
*
*
*
  bss 2
* entry locat
*
*   m=locat(list,test)
*
locat enter np=2
  argadd a2,2
  a1 npb,0
  argadd a6,1
  a0 a6
  s2 0,a2
  vl a1
  v7 ,a0,1
  v6 s2-v7
  vm v6,z
  s1 vm
  s0 vm
  a7 zs1
  jsz retl
  a7 a7+1
  s1 a7
retl exit
  end
          ident upackh
* entry upackh
upackh enter np=3
  argadd a1,1
  argadd a2,2
  argadd a3,3
  s2 <24
  s1 0,a1
  s3 s1&s2
  s1 s1>32
  0,a3 s3
  0,a2 s1
  exit
  end
          ident locate
* entry locate
*
*   l = locate(list,nlist,ispace,pattern)
*
*  to find member of   list   =   pattern
*  or return zero if no match found
*
locate enter np=4
  argadd a2,2
  argadd a3,3
  argadd a4,4
  argadd a1,1
  a7 64
  s1 0
  a2 0,a2
  a3 0,a3
  s4 0,a4
  a6 0
  vl a7
  a0 a7-a2
  a5 a3*a7
  jap vvv
  a0 a1
  j xxx
yyy s0 vm
  a0 a1
  jsn zzz
xxx v7 ,a0,a3
  a2 a2-a7
  a1 a1+a5
  v6 s4-v7
  a0 a7-a2
  vm v6,z
  a6 a6+a7
  jam yyy
  s0 vm
  a0 a2
  jsz www
zzz s6 vm
  a6 a6-a7
aaa a7 zs6
  a6 a6+1
  a6 a6+a7
  s1 a6
out exit
vvv a0 a2
www jaz out
  a0 a1
  vl a2
  v7 ,a0,a3
  v6 s4-v7
  vm v6,z
  s0 vm
  s6 vm
  jsn aaa
  exit
  end
          ident square
*  entry square
*  call square(r,a,mrowr,n)
square enter np=4
  argadd a7,4
  argadd a6,3
  argadd a5,2
  argadd a4,1
  a3 64
  a7 0,a7
  a6 0,a6
  a0 a7
  a1 1
  a2 a4
  b70 a7
  jaz done
  a7 a6*a3
loop1 a0 a1-a3
  b71 a1
  b72 a4
  b73 a2
  vl a3
loop2 jam ex1
  a0 a5
  v7 ,a0,1
  a0 a4
  a1 a1-a3
  a5 a5+a3
  ,a0,1 v7
  a0 a2
  a4 a4+a3
  ,a0,a6 v7
  a0 a1-a3
  a2 a2+a7
  j loop2
ex1 a0 a1
  vl a1
  jaz term
  a0 a5
  v7 ,a0,1
  a0 a4
  a5 a5+a1
  ,a0,1 v7
  a0 a2
  ,a0,a6 v7
term a1 b71
  a2 b70
  a4 b72
  a0 a1-a2
  a2 b73
  a1 a1+1
  a4 a4+a6
  a2 a2+1
  jan loop1
done exit
  end
          ident sqtrip
*  entry sqtrip
*  call sqtrip (r , a , n)
sqtrip enter np=3
  argadd a7,3
  argadd a5,1
  a1 64
  a7 0,a7
  a4 a5+1
  a0 a1-a7
  vl a7
  a2 a7+1
  a3 a5+a7
  jam ex1
  v7 0
  j ex4
ex1 vl a1
  v7 0
  b70 a7
ex2 a0 a5
  ,a0,a2 v7
  a7 a7-a1
  a6 a2*a1
  a0 a7-a1
  a5 a5+a6
  jap ex2
  a0 a7
  vl a7
  a7 b70
  jaz ex3
ex4 a0 a5
  ,a0,a2 v7
ex3 argadd a6,2
  a2 a7*a1
  a5 1
loop1 a0 a5-a1
  b70 a5
  b71 a4
  b72 a3
  vl a1
  jap loop2
  vl a5
  j ex6
loop2 a0 a6
  v2 ,a0,1
  a5 a5-a1
  a6 a6+a1
  a0 a3
  v1 -fv2
  ,a0,1 v2
  a0 a4
  a3 a3+a1
  ,a0,a7 v1
  a0 a5-a1
  a4 a4+a2
  jap loop2
  a0 a5
  vl a5
  jaz ex7
ex6 a0 a6
  v2 ,a0,1
  a6 a6+a5
  a0 a3
  v1 -fv2
  ,a0,1 v2
  a0 a4
  ,a0,a7 v1
ex7 a5 b70
  a4 b71
  a3 b72
  a5 a5+1
  a4 a4+1
  a3 a3+a7
  a0 a5-a7
  jan loop1
  exit
  end
  ident sgmata
*  entry sgmata
* call sgmata(fock,dens)
  common blkin
gin bss 340
ijin bss 170
mword bss 1
  common vectem
ikyk bss 64
ijkl bss 64
l bss 64
k bss 64
  block *
  block *
sgmata enter np=2
  s1 mword,0
  a1 gin
  argadd a5,2
  a2 ijin
  argadd a6,1
  b71 a1
  b70 a2
  s6 s1
  s1 s1>1
  s3 1
  s7 s6-s1
  s5 a5
  a6 a6-a5
  s4 s7-s3
  t70 s5
  s5 s7\s1
  s4 s4>6
  t77 s5
  a1 s4
  s4 s4<6
  a0 a1
  a3 64
  s2 s7-s4
  b74 a3
  vl a3
  a2 s2
  s0 0
  b75 a2
aaaa jan www
  a3 b75
  s0 t77
  vl a3
www b73 a1
  b72 a3
* vl set
* a3=vl
* s0=svl-lvl
* even cases
  a0 b70
  a4 8
* clock 1
  v7 ,a0,1
  s5 <8
  v6 v7>a4
  v5 s5&v6
  a1 24
  a0 l
* clock 2
  v4 s5&v7
  v3 v5<a1
  v2 v3*fv3
* clock 3
  ,a0,1 v4
  v1 v2-v5
  v0 v1>1
  a0 k
* clock 4
  ,a0,1 v5
  v7 v6>a4
  v3 s5&v7
  v2 v4+v0
  a7 b71
  a2 2
  a0 ikyk
* clock 5
  v6 v7>a4
  v1 s5&v6
  ,a0,1 v0
  a0 a7+1
* clock 6
  v7 v1<a1
  v5 v7*fv7
  v0 ,a0,a2
* clock 7
  v4 v5-v1
  v7 v4>1
  s7 4.
  a7 a3-1
  a0 a3-1
  s6 t70
* clock 8
  v4 s7*rv0
  v5 v7+v3
  v0 v6>a4
* v0=ijkl v1=i v2=kl v3=j
* v4=g4 v5=ij v6=scratch v7=ikyi
* k,l,ikyk stored
  jsz pppf
  jaz ttt
  a3 a7
  a7 a7-1
* s0=svl-lvl
* a0=vl-1  a6=<a-p> a7=vl-1
* t70=<p>
* t77=lvl-svl(for last time round)
* b70=<ijklijkl> b71=<g> b72=cvl
* b73=nvord-1 b74=64 b75=lrem
* coulomb part(even cases)
pppf s0 a3
  b77 a3
* clock 9
  v6 s6+v2
  s0 s0<63
* clock 10
  v2 s6+v5
* v2=<pij> v4=g4 v6=<pkl>
  s1 v6,a7
  s2 v2,a7
  jsp ppph
  s3 v4,a7
  a1 s1
  a2 s2
  s1 -1,a1
  s2 -1,a2
  a2 a2+a6
  a1 a1+a6
  s4 -1,a2
  s6 -1,a1
  s1 s1*rs3
  s2 s2*rs3
  a0 a7-1
  a7 a7-1
  s1 s1+fs4
  s2 s2+fs6
  -1,a2 s1
  -1,a1 s2
  jam pppg
* above handles loop folded cases (odd no. of terms)
ppp s1 v6,a7
  s2 v2,a7
ppph a5 a7-1
  s3 v4,a7
  s6 v6,a5
  a1 s1
  a2 s2
  s1 -1,a1
  s2 -1,a2
  a2 a2+a6
  a1 a1+a6
  s7 v2,a5
  s4 -1,a2
  s5 -1,a1
  a3 s6
  s1 s1*rs3
  a4 s7
  s2 s2*rs3
  s6 -1,a3
  s7 -1,a4
  s3 v4,a5
  s1 s1+fs4
  s2 s2+fs5
  a4 a4+a6
  a3 a3+a6
  s6 s6*rs3
  s7 s7*rs3
  -1,a2 s1
  -1,a1 s2
  s4 -1,a4
  s5 -1,a3
  a0 a5-1
  a7 a5-1
  s6 s6+fs4
  s7 s7+fs5
  -1,a4 s6
  -1,a3 s7
  jap ppp
* above ppp loop costs about 70 clocks per two elements
*
* exchange part (even cases)
* code based on sequence
* kik=gik*rjl + kik
* kil=gil*rjk + kil
* kjk=gil*ril + kjk
* kjl=gik*rik + kjl
* clock 1
pppg a0 ijkl
  a1 24
  ,a0,1 v0
  v4 v3<a1
  v5 v4*fv4
  a0 k
* clock 2
  v2 ,a0,1
  v0 v1\v2
  v6 v5-v3
  v4 v6>1
  a0 l
* clock 3
  vm v0,z
  v6 ,a0,1
  v1 v3-v6
* clock 4
  s1 vm
  vm v1,z
  v5 v7+v2
  v0 v6<a1
* clock 5
  s2 vm
  vm v1,m
  ,a0,1 v5
  s1 s1!s2
  a0 k
* clock 6
  v1 v7+v6
  v5 v0*fv0
* clock 7
  v7 v4+v6
  ,a0,1 v1
  a0 ikyk
  a2 2
* clock 8
  v0 v5-v6
  v1 v0>1
  a7 b71
* clock 9
  v5 ,a0,1
  v6 v1+v3
  v0 v6!v7&vm
  a0 a7+1
* clock 10
  v7 ,a0,a2
  v6 v4+v2
  a0 l
* clock 11
  v1 v3-v2
  vm v1,m
* clock 12
  v2 v5+v3
  v4 v2!v6&vm
* clock 13
  v5 ,a0,1
  vm v1,z
  v2 v7+fv7
  a7 b77
  a0 k
* clock 14
  v6 ,a0,1
  v3 v2!v7&vm
  a7 a7-1
  vm s1
  s6 t70
* clock 15
  v1 v2!v7&vm
* in the above il stored in k ik stored in l
* clock 16
  v2 s6+v5
* clock 17
  v7 s6+v6
* clock 18
  v5 s6+v4
* clock 19
  v6 s6+v0
* v1=gik v2=<pik> v3=gil v5=<pjk> v6=<pjl> v7=<pil>
qqq s1 v2,a7
  s2 v7,a7
  s3 v5,a7
  s4 v6,a7
  s5 v1,a7
  s6 v3,a7
  a1 s1
  a2 s2
  s1 -1,a1
  s2 -1,a2
  a3 s3
  a4 s4
  s3 -1,a3
  s4 -1,a4
  a4 a4+a6
  a3 a3+a6
  a2 a2+a6
  s7 -1,a4
  s1 s1*rs5
  s2 s2*rs6
  s3 s3*rs6
  s6 -1,a3
  a1 a1+a6
  s4 s4*rs5
  s5 -1,a2
  s1 s7-fs1
  s7 -1,a1
  a0 a7-1
  a7 a7-1
  s2 s6-fs2
  s3 s5-fs3
  s4 s7-fs4
  -1,a4 s1
  -1,a3 s2
  -1,a2 s3
  -1,a1 s4
  jap qqq
* above costs about 105 clocks
* odd cases
  a0 ijkl
  a4 8
* clock 1
  v0 ,a0,1
  s5 <8
ttt v6 v0>a4
  v5 s5&v6
  a7 b70
  a3 b72
  a1 24
  a0 l
  a7 a7+a3
* clock 2
  v4 s5&v0
  v3 v5<a1
  v2 v3*fv3
  b70 a7
* clock 3
  ,a0,1 v4
  v1 v2-v5
  v7 v1>1
  a0 k
* clock 4
  ,a0,1 v5
  v0 v6>a4
  v3 s5&v0
  v2 v4+v7
  a0 b71
  a2 2
  s7 4.
* clock 5
  v1 ,a0,a2
  v5 s7*rv1
  v6 v0>a4
* clock 6
  v1 v7+v3
  v4 v6<a1
  v0 v4*fv4
  a7 a3-1
  a0 ikyk
* clock 7
  ,a0,1 v1
  v7 v0-v6
  v4 v7>1
  s6 t70
* clock 8
  v0 v4+v3
  s0 a3
* l,k stored
* ikykj stored in ikyk
* v1=ikykj v2=kl v3=j v4=ikyi v5=g4
* v6=i v0=ij
* coulomb part(odd cases)
* clock 9
  v7 s6+v2
  s0 s0<63
* clock 10
  v2 s6+v0
  s1 v7,a7
  s2 v2,a7
  jsp pph
  s3 v5,a7
* v2=<pij> v5=g4 v7=<pkl>
  a1 s1
  a2 s2
  s1 -1,a1
  s2 -1,a2
  a2 a2+a6
  a1 a1+a6
  s4 -1,a2
  s6 -1,a1
  s1 s1*rs3
  s2 s2*rs3
  a0 a7-1
  a7 a7-1
  s1 s1+fs4
  s2 s2+fs6
  -1,a2 s1
  -1,a1 s2
  jam ppg
pp s1 v7,a7
  s2 v2,a7
pph a5 a7-1
  s3 v5,a7
  s6 v7,a5
  a1 s1
  a2 s2
  s1 -1,a1
  s2 -1,a2
  a2 a2+a6
  a1 a1+a6
  s7 v2,a5
  s4 -1,a2
  s5 -1,a1
  a3 s6
  s1 s1*rs3
  a4 s7
  s2 s2*rs3
  s6 -1,a3
  s7 -1,a4
  s3 v5,a5
  s1 s1+fs4
  s2 s2+fs5
  a4 a4+a6
  a3 a3+a6
  s6 s6*rs3
  s7 s7*rs3
  -1,a2 s1
  -1,a1 s2
  s4 -1,a4
  s5 -1,a3
  a0 a5-1
  a7 a5-1
  s6 s6+fs4
  s7 s7+fs5
  -1,a4 s6
  -1,a3 s7
  jap pp
* v1=ikykj v3=j v4=ikyi v6=i
ppg a0 l
  a1 24
* clock 1
  v0 ,a0,1
  v5 v4+v0
  v2 v3<a1
  v7 v2*fv2
  a7 b71
  a3 b72
* clock 2
  ,a0,1 v5
  v1 v7-v3
  v2 v1>1
  a0 k
  a7 a7+a3
* clock 3
  v5 ,a0,1
  v1 v4+v5
  v7 v0<a1
  a7 a7+a3
* clock 4
  ,a0,1 v1
  v4 v6-v5
  vm v4,z
  a0 b71
  a2 2
* clock 5
  v6 v7*fv7
  v1 v6-v0
  v4 v1>1
  s7 vm
  b71 a7
* clock 6
  v7 ,a0,a2
  v1 v4+v3
* clock 7
  v6 v3-v0
  vm v6,z
  a7 a3-1
* clock 8
  v4 v2+v0
  s6 vm
  vm v6,m
  s7 s7!s6
* clock 9
  v6 v3-v5
  v0 v1!v4&vm
  a0 ikyk
* clock 10
  v3 v7+fv7
  vm v6,m
  v1 v2+v5
  v4 ,a0,1
  a0 k
* clock 11
  v5 ,a0,1
  v2 v4!v1&vm
  a0 l
* clock 12
  vm v6,z
  v4 ,a0,1
* clock 13
  v1 v3!v7&vm
* clock 14
  vm s7
  s6 t70
* clock 15
  v6 v3!v7&vm
* in the above ik stored in k il stored in l
* v1=gil v2=jk v4=il v5=ik v6=gik v0=jl
* clock 16
  v7 s6+v5
* clock 17
  v3 s6+v4
* clock 18
  v5 s6+v2
* clock 19
  v4 s6+v0
* v1=gil v3=<pil> v4=<pjl> v5=<pjk> v6=gik v7=<pik>
qq s1 v7,a7
  s2 v3,a7
  s3 v5,a7
  s4 v4,a7
  s5 v6,a7
  s6 v1,a7
  a1 s1
  a2 s2
  s1 -1,a1
  s2 -1,a2
  a3 s3
  a4 s4
  s3 -1,a3
  s4 -1,a4
  a4 a4+a6
  a3 a3+a6
  a2 a2+a6
  s7 -1,a4
  s1 s1*rs5
  s2 s2*rs6
  s3 s3*rs6
  s6 -1,a3
  a1 a1+a6
  s4 s4*rs5
  s5 -1,a2
  s1 s7-fs1
  s7 -1,a1
  a0 a7-1
  a7 a7-1
  s2 s6-fs2
  s3 s5-fs3
  s4 s7-fs4
  -1,a4 s1
  -1,a3 s2
  -1,a2 s3
  -1,a1 s4
  jap qq
  a1 b73
  a3 b74
  a0 a1-1
  a1 a1-1
  s0 0
  jap aaaa
  exit
  end
          ident proc2
* entry proc2
* call proc2(fock,dens,k,q)
  common blkin
gin bss 340
ijin bss 170
mword bss 1
  common vectem
ikyk bss 64
ijkl bss 64
l bss 64
k bss 64
  block *
  block *
count set 3
bsave bss count
proc2 enter np=4
  a1 count
  a0 bsave
  a2 gin
  ,a0 b1,a1
  s1 mword,0
  argadd a3,4
  argadd a4,3
  argadd a6,1
  argadd a5,2
  b71 a2
  s6 s1
  s1 s1>1
  a2 ijin
  s3 1
  s7 s6-s1
  a4 a4-a3
  a3 a3-a6
  a6 a6-a5
  s4 s7-s3
  s5 a5
  b70 a2
  t70 s5
  s5 s7\s1
  s4 s4>6
  t77 s5
  b76 a3
  a1 s4
  s4 s4<6
  a0 a1
  a3 64
  b77 a4
  s2 s7-s4
  b1 a6
  b74 a3
  vl a3
  a2 s2
  s0 0
  b75 a2
aaaa jan www
  a3 b75
  s0 t77
  vl a3
www b73 a1
  b72 a3
* vl set
* a3=vl
* s0=svl-lvl
* even cases
  a0 b70
  a4 8
* clock 1
  v7 ,a0,1
  s5 <8
  v6 v7>a4
  v5 s5&v6
  a1 24
  a0 l
* clock 2
  v4 s5&v7
  v3 v5<a1
  v2 v3*fv3
* clock 3
  ,a0,1 v4
  v1 v2-v5
  v0 v1>1
  a0 k
* clock 4
  ,a0,1 v5
  v7 v6>a4
  v3 s5&v7
  v2 v4+v0
  a7 b71
  a2 2
  a0 ikyk
* clock 5
  v6 v7>a4
  v1 s5&v6
  ,a0,1 v0
  a0 a7+1
* clock 6
  v7 v1<a1
  v5 v7*fv7
  v0 ,a0,a2
* clock 7
  v4 v5-v1
  v7 v4>1
  s7 4.
  a7 a3-1
  a0 a3-1
  s6 t70
* clock 8
  v4 s7*rv0
  v5 v7+v3
  v0 v6>a4
* v0=ijkl v1=i v2=kl v3=j
* v4=g4 v5=ij v6=scratch v7=ikyi
* k,l,ikyk stored
  jsz pppf
  jaz ttt
  a3 a7
  a7 a7-1
* s0=svl-lvl
* a0=vl-1  a6=<f-p> a7=vl-1
* t70=<p>
* t77=lvl-svl(for last time round)
* b70=<ijklijkl> b71=<g> b72=cvl
* b73=nvord-1 b74=64 b75=lrem
* b76=<q-f> b77=<k-q>
* b1=<f-p>
* coulomb part(even cases)
pppf s0 a3
  b3 a3
* clock 9
  v6 s6+v2
  s0 s0<63
* clock 10
  v2 s6+v5
* v2=<pij> v4=g4 v6=<pkl>
  s1 v6,a7
  s2 v2,a7
  jsp ppph
  s3 v4,a7
  a1 s1
  a2 s2
  s1 -1,a1
  s2 -1,a2
  a2 a2+a6
  a1 a1+a6
  s4 -1,a2
  s6 -1,a1
  s1 s1*rs3
  s2 s2*rs3
  a0 a7-1
  a7 a7-1
  s1 s1+fs4
  s2 s2+fs6
  -1,a2 s1
  -1,a1 s2
  jam pppg
* above handles loop folded cases (odd no. of terms)
ppp s1 v6,a7
  s2 v2,a7
ppph a5 a7-1
  s3 v4,a7
  s6 v6,a5
  a1 s1
  a2 s2
  s1 -1,a1
  s2 -1,a2
  a2 a2+a6
  a1 a1+a6
  s7 v2,a5
  s4 -1,a2
  s5 -1,a1
  a3 s6
  s1 s1*rs3
  a4 s7
  s2 s2*rs3
  s6 -1,a3
  s7 -1,a4
  s3 v4,a5
  s1 s1+fs4
  s2 s2+fs5
  a4 a4+a6
  a3 a3+a6
  s6 s6*rs3
  s7 s7*rs3
  -1,a2 s1
  -1,a1 s2
  s4 -1,a4
  s5 -1,a3
  a0 a5-1
  a7 a5-1
  s6 s6+fs4
  s7 s7+fs5
  -1,a4 s6
  -1,a3 s7
  jap ppp
* above ppp loop costs about 70 clocks per two elements
*
* exchange part (even cases)
* code based on sequence
* kik=gik*rjl + kik
* kil=gil*rjk + kil
* kjk=gil*ril + kjk
* kjl=gik*rik + kjl
* clock 1
pppg a0 ijkl
  a1 24
  ,a0,1 v0
  v4 v3<a1
  v5 v4*fv4
  a0 k
* clock 2
  v2 ,a0,1
  v0 v1\v2
  v6 v5-v3
  v4 v6>1
  a0 l
* clock 3
  vm v0,z
  v6 ,a0,1
  v1 v3-v6
* clock 4
  s1 vm
  vm v1,z
  v5 v7+v2
  v0 v6<a1
* clock 5
  s2 vm
  vm v1,m
  ,a0,1 v5
  s1 s1!s2
  a0 k
* clock 6
  v1 v7+v6
  v5 v0*fv0
* clock 7
  v7 v4+v6
  ,a0,1 v1
  a0 ikyk
  a2 2
* clock 8
  v0 v5-v6
  v1 v0>1
  a7 b71
* clock 9
  v5 ,a0,1
  v6 v1+v3
  v0 v6!v7&vm
  a0 a7+1
* clock 10
  v7 ,a0,a2
  v6 v4+v2
  a0 l
* clock 11
  v1 v3-v2
  vm v1,m
* clock 12
  v2 v5+v3
  v4 v2!v6&vm
* clock 13
  v5 ,a0,1
  vm v1,z
  v2 v7+fv7
  a7 b3
  a0 k
* clock 14
  v6 ,a0,1
  v3 v2!v7&vm
  a7 a7-1
  vm s1
  s6 t70
* clock 15
  v1 v2!v7&vm
* in the above il stored in k ik stored in l
* clock 16
  v2 s6+v5
* clock 17
  v7 s6+v6
* clock 18
  v5 s6+v4
* clock 19
  v6 s6+v0
* v1=gik v2=<pik> v3=gil v5=<pjk> v6=<pjl> v7=<pil>
qqq s1 v2,a7
  s2 v7,a7
  s3 v5,a7
  s4 v6,a7
  s5 v1,a7
  s6 v3,a7
  a1 s1
  a2 s2
  s1 -1,a1
  s2 -1,a2
  a3 s3
  a4 s4
  s3 -1,a3
  s4 -1,a4
  a4 a4+a6
  a3 a3+a6
  a2 a2+a6
  s7 -1,a4
  s1 s1*rs5
  s2 s2*rs6
  s3 s3*rs6
  t72 s6
  s6 -1,a3
  a1 a1+a6
  s4 s4*rs5
  t71 s5
  s5 -1,a2
  a6 b76
  s1 s7-fs1
  s7 -1,a1
  a5 a1+a6
  b2 a7
  a7 a2+a6
  s2 s6-fs2
  s6 -1,a5
  s3 s5-fs3
  s5 -1,a7
  s4 s7-fs4
  s7 t71
  -1,a4 s1
  -1,a3 s2
  a3 a3+a6
  a6 a6+a4
  s1 -1,a3
  s2 -1,a6
  -1,a2 s3
  s3 t72
  a4 b77
  s6 s7*rs6
  a2 a6+a4
  -1,a1 s4
  s4 -1,a2
  s5 s5*rs3
  s1 s1*rs3
  a3 a3+a4
  a1 a7+a4
  a5 a5+a4
  s3 -1,a3
  s2 s2*rs7
  s7 -1,a1
  s4 s4+fs6
  s6 -1,a5
  a7 b2
  s3 s3+fs5
  a0 a7-1
  a7 a7-1
  -1,a2 s4
  s7 s7+fs1
  s6 s6+fs2
  a6 b1
  -1,a3 s3
  -1,a1 s7
  -1,a5 s6
  jap qqq
* above costs about 105 clocks
* odd cases
  a0 ijkl
  a4 8
* clock 1
  v0 ,a0,1
  s5 <8
ttt v6 v0>a4
  v5 s5&v6
  a7 b70
  a3 b72
  a1 24
  a0 l
  a7 a7+a3
* clock 2
  v4 s5&v0
  v3 v5<a1
  v2 v3*fv3
  b70 a7
* clock 3
  ,a0,1 v4
  v1 v2-v5
  v7 v1>1
  a0 k
* clock 4
  ,a0,1 v5
  v0 v6>a4
  v3 s5&v0
  v2 v4+v7
  a0 b71
  a2 2
  s7 4.
* clock 5
  v1 ,a0,a2
  v5 s7*rv1
  v6 v0>a4
* clock 6
  v1 v7+v3
  v4 v6<a1
  v0 v4*fv4
  a7 a3-1
  a0 ikyk
* clock 7
  ,a0,1 v1
  v7 v0-v6
  v4 v7>1
  s6 t70
* clock 8
  v0 v4+v3
  s0 a3
* l,k stored
* ikykj stored in ikyk
* v1=ikykj v2=kl v3=j v4=ikyi v5=g4
* v6=i v0=ij
* coulomb part(odd cases)
* clock 9
  v7 s6+v2
  s0 s0<63
* clock 10
  v2 s6+v0
  s1 v7,a7
  s2 v2,a7
  jsp pph
  s3 v5,a7
* v2=<pij> v5=g4 v7=<pkl>
  a1 s1
  a2 s2
  s1 -1,a1
  s2 -1,a2
  a2 a2+a6
  a1 a1+a6
  s4 -1,a2
  s6 -1,a1
  s1 s1*rs3
  s2 s2*rs3
  a0 a7-1
  a7 a7-1
  s1 s1+fs4
  s2 s2+fs6
  -1,a2 s1
  -1,a1 s2
  jam ppg
pp s1 v7,a7
  s2 v2,a7
pph a5 a7-1
  s3 v5,a7
  s6 v7,a5
  a1 s1
  a2 s2
  s1 -1,a1
  s2 -1,a2
  a2 a2+a6
  a1 a1+a6
  s7 v2,a5
  s4 -1,a2
  s5 -1,a1
  a3 s6
  s1 s1*rs3
  a4 s7
  s2 s2*rs3
  s6 -1,a3
  s7 -1,a4
  s3 v5,a5
  s1 s1+fs4
  s2 s2+fs5
  a4 a4+a6
  a3 a3+a6
  s6 s6*rs3
  s7 s7*rs3
  -1,a2 s1
  -1,a1 s2
  s4 -1,a4
  s5 -1,a3
  a0 a5-1
  a7 a5-1
  s6 s6+fs4
  s7 s7+fs5
  -1,a4 s6
  -1,a3 s7
  jap pp
* v1=ikykj v3=j v4=ikyi v6=i
ppg a0 l
  a1 24
* clock 1
  v0 ,a0,1
  v5 v4+v0
  v2 v3<a1
  v7 v2*fv2
  a7 b71
  a3 b72
* clock 2
  ,a0,1 v5
  v1 v7-v3
  v2 v1>1
  a0 k
  a7 a7+a3
* clock 3
  v5 ,a0,1
  v1 v4+v5
  v7 v0<a1
  a7 a7+a3
* clock 4
  ,a0,1 v1
  v4 v6-v5
  vm v4,z
  a0 b71
  a2 2
* clock 5
  v6 v7*fv7
  v1 v6-v0
  v4 v1>1
  s7 vm
  b71 a7
* clock 6
  v7 ,a0,a2
  v1 v4+v3
* clock 7
  v6 v3-v0
  vm v6,z
  a7 a3-1
* clock 8
  v4 v2+v0
  s6 vm
  vm v6,m
  s7 s7!s6
* clock 9
  v6 v3-v5
  v0 v1!v4&vm
  a0 ikyk
* clock 10
  v3 v7+fv7
  vm v6,m
  v1 v2+v5
  v4 ,a0,1
  a0 k
* clock 11
  v5 ,a0,1
  v2 v4!v1&vm
  a0 l
* clock 12
  vm v6,z
  v4 ,a0,1
* clock 13
  v1 v3!v7&vm
* clock 14
  vm s7
  s6 t70
* clock 15
  v6 v3!v7&vm
* in the above ik stored in k il stored in l
* v1=gil v2=jk v4=il v5=ik v6=gik v0=jl
* clock 16
  v7 s6+v5
* clock 17
  v3 s6+v4
* clock 18
  v5 s6+v2
* clock 19
  v4 s6+v0
* v1=gil v3=<pil> v4=<pjl> v5=<pjk> v6=gik v7=<pik>
qq s1 v7,a7
  s2 v3,a7
  s3 v5,a7
  s4 v4,a7
  s5 v6,a7
  s6 v1,a7
  a1 s1
  a2 s2
  s1 -1,a1
  s2 -1,a2
  a3 s3
  a4 s4
  s3 -1,a3
  s4 -1,a4
  a4 a4+a6
  a3 a3+a6
  a2 a2+a6
  s7 -1,a4
  s1 s1*rs5
  s2 s2*rs6
  s3 s3*rs6
  t72 s6
  s6 -1,a3
  a1 a1+a6
  s4 s4*rs5
  t71 s5
  s5 -1,a2
  a6 b76
  s1 s7-fs1
  s7 -1,a1
  a5 a1+a6
  b2 a7
  a7 a2+a6
  s2 s6-fs2
  s6 -1,a5
  s3 s5-fs3
  s5 -1,a7
  s4 s7-fs4
  s7 t71
  -1,a4 s1
  -1,a3 s2
  a3 a3+a6
  a6 a6+a4
  s1 -1,a3
  s2 -1,a6
  -1,a2 s3
  s3 t72
  a4 b77
  s6 s7*rs6
  a2 a6+a4
  -1,a1 s4
  s4 -1,a2
  s5 s5*rs3
  s1 s1*rs3
  a3 a3+a4
  a1 a7+a4
  a5 a5+a4
  s3 -1,a3
  s2 s2*rs7
  s7 -1,a1
  s4 s4+fs6
  s6 -1,a5
  a7 b2
  s3 s3+fs5
  a0 a7-1
  a7 a7-1
  -1,a2 s4
  s7 s7+fs1
  s6 s6+fs2
  a6 b1
  -1,a3 s3
  -1,a1 s7
  -1,a5 s6
  jap qq
  a1 b73
  a3 b74
  a0 a1-1
  a1 a1-1
  s0 0
  jap aaaa
  a0 bsave
  a1 count
  b1,a1 ,a0
  exit
  end
          ident mxmbn
*  entry mxmbn
* call mxmbn(a,mcola,mrowa,b,mcolb,mrowb,r,mcolr,mrowr,
* ncol,nlink,nrow)
*
* r(ncol:mcolr,nrow:mrowr)=r(ncol:mcolr,nrow:mrowr)-
* a(ncol:mcola,nlink:mrowa)*b(nlink:mcolb,nrow:mrowb)
* matrix mult. where sparsity of b is used
*     *****r must be pre-set*****
count set 20
btemp bss count
*
* b00=internal link    b01=b(ib)
* b02=mrowb   b03=r(ir)   b04=mrowr
* b05=mcola*64   b06=mcolr*64   b07=r(irrr)
* b10=a(iaaa)   b11=nlast b12=mcolb  b13=mrowa*64
* b14=64 b15=<step1> b16=<step3> b17=<p1c>
* b20=<ste22> b21=<step2>
* b22=<p22x>  b23=<ste10>
*
* b70=a(iaa)   b71=address for delayed store of r(ir)
* b72=vector length for delayed store of r(ir)
* b73=mcolr  b74=b(1)  b75=nvast
* b76=b(ibb)  b77=mcolb*64
*
* t70=nn t71=nrow t72=nlseg
*
mxmbn enter np=12
  a0 btemp
  a1 count
  s6 1
  a6 64
  ,a0 b00,a1
  a0 ste10
  a1 p22x
  a2 ste22
  a3 step2
  a5 p1c
  a7 step1
  a4 step3
  b23 a0
  b22 a1
  b20 a2
  b21 a3
  b14 a6
  b17 a5
  b15 a7
  b16 a4
  argadd a1,10
  argadd a2,11
  argadd a3,8
  argadd a5,9
  argadd a7,12
  argadd a4,6
  a0 0
  s1 0,a1
  s2 0,a2
  a3 0,a3
  a5 0,a5
  s7 0,a7
  a4 0,a4
  argadd a1,5
  argadd a2,2
  argadd a7,3
  s4 s1-s6
  s5 s2-s6
  b73 a3
  a3 a3*a6
  b04 a5
  argadd a5,7
  t71 s7
  s4 s4>6
  s5 s5>6
  b02 a4
  a4 0,a1
  a1 0,a2
  a2 0,a7
  t72 s5
  s0 s4
  b06 a3
  argadd a3,4
  a7 a1*a6
  b12 a4
  a4 a4*a6
  b71 a0
  s4 s4<6
  s5 s5<6
  a0 a2*a6
  s1 s1-s4
  s2 s2-s5
  b05 a7
  argadd a7,1
  a6 s1
  b77 a4
  a4 s2
  b74  a3
  b13 a0
  b75 a6
  b11 a4
p9999 a4 b14
  t70 s0
  jsn p9998
  a4 b75
*
  j p9998
  align
*
p9998 b10 a7
  b07 a5
  a7 b74
  s5 t71
p1 a6 b15
  b01 a7
  b03 a5
  a3 b10
  s0 t72
  s4 t72
  s1 0
  s2 0
  b00 a6
p2 a6 b14
  jsn p2a
  a6 b11
p2a a5 b12
  a0 a7
  vl a6
  b76 a7
  v0 ,a0,a5
  a5 0
  b70 a3
  vm v0,n
  vl a4
  s3 vm
p22 a7 zs3
  s0 s3
  s3 s3<1
p22y a6 a2*a7
p22x jsz p2x
p22z a5 a5+a7
  s3 s3<a7
p22w a0 a3+a6
  a3 a3+a6
*  a1=mcola  a2=mrowa  a3=a(ia) a4=cvl for central loop
*  a5=posn in b vector(v0)
*  s1=ipar1 s2=ipar2 s3=vector mask
*  s4=kk s5=j s6=1 s7=fac
  j b00
*
* 1st possibility
*
step1 a0 b03
  a6 b73
  a7 b03
  v2 ,a0,a6
  a0 b71
  b71 a7
  a7 b72
  jaz sta
  vl a7
  ,a0,a6 v7
  vl a4
sta a0 a3
  b72 a4
  a6 b21
  j stb
*
* 2nd possibility
*
step2 v7 ,a0,a1
  a7 zs3
  s7 v0,a5
  s0 s3
  s2 -1
  a7 a7+1
  v3 s7*rv7
  r p22y
*
* 3rd possibility
*
step3 v6 ,a0,a1
  a7 zs3
  s7 v0,a5
  s0 s3
  s1 1
  a7 a7+1
  v7 s7*rv6
  a6 a2*a7
  v2 v1-fv7
  r p22x
*
* 4th possibility
*
  v5 ,a0,a1
  a7 zs3
  s7 v0,a5
  s0 s3
  s2 1
  a7 a7+1
  v6 s7*rv5
  a6 a2*a7
  v4 v3+fv6
  r p22x
*
* 5th possibility
*
  v7 ,a0,a1
  a7 zs3
  s7 v0,a5
  s0 s3
  s1 -1
  a7 a7+1
  v5 s7*rv7
  a6 a2*a7
  v1 v2-fv5
  r p22x
*
* 6th possibility
*
  v6 ,a0,a1
  a7 zs3
  s7 v0,a5
  s0 s3
  s2 -1
  a7 a7+1
  v7 s7*rv6
  a6 a2*a7
  v3 v4+fv7
  r p22x
*
* 7th possibility
*
  v5 ,a0,a1
  a7 zs3
  s7 v0,a5
  s0 s3
  s1 1
  a7 a7+1
  v6 s7*rv5
  a6 a2*a7
  v2 v1-fv6
  r p22x
*
* 8th possibility
*
  v7 ,a0,a1
  a7 zs3
  s7 v0,a5
  s0 s3
  s2 1
  a7 a7+1
  v5 s7*rv7
  a6 a2*a7
  v4 v3+fv5
  r p22x
*
* 9th possibility
*
  v6 ,a0,a1
  a7 zs3
  s7 v0,a5
  s0 s3
  s1 -1
  a0 b23
  a7 a7+1
  v7 s7*rv6
  b00 a0
  a6 a2*a7
  a5 a5+a7
  s3 s3<a7
  v1 v2-fv7
  jsn p22w
  j p2x
*
* 10th possibility
*
ste10 v5 ,a0,a1
  a7 zs3
  s7 v0,a5
  s0 s3
  s2 -1
  a7 a7+1
  v6 s7*rv5
  a6 a2*a7
  v3 v4+fv6
  r p22x
*
* 11th possibility
*
  v7 ,a0,a1
  a7 zs3
  s7 v0,a5
  s0 s3
  s1 1
  a7 a7+1
  v5 s7*rv7
  a6 a2*a7
  v2 v1-fv5
  r p22x
*
* 12th possibility
*
  v6 ,a0,a1
  a7 zs3
  s7 v0,a5
  s0 s3
  s2 1
  a7 a7+1
  v7 s7*rv6
  a6 a2*a7
  v4 v3+fv7
  r p22x
*
* 1st possibility 2nd time round
*
  a6 b20
stb v5 ,a0,a1
  a7 zs3
  s7 v0,a5
  s0 s3
  s1 -1
  a7 a7+1
  v6 s7*rv5
  b00 a6
  a6 a2*a7
  v1 v2-fv6
  j b22
*
* 2nd possibilty 2nd time round
*
ste22 v7 ,a0,a1
  a7 zs3
  s7 v0,a5
  s0 s3
  s2 -1
  a7 a7+1
  v5 s7*rv7
  a0 b16
  a6 a2*a7
  b00 a0
  v3 v4+fv5
  jsn p22z
p2x s0 s4-s6
  a7 b76
  a6 b77
  a3 b70
  a5 b13
  s4 s4-s6
  a7 a7+a6
  a3 a3+a5
  jsp p2
  s0 s1
  a7 b01
  a6 b02
  jsm p1b
  jsz p1c
*
* first part of result in v2
*
  s0 s2
  jsm p1d
* first part in v2 second in v4
  v7 v2-fv4
  j b17
* first part in v2 second in v3
p1d v7 v2-fv3
  j b17
*
* first part of result in v1
*
p1b s0 s2
  jsm p1e
  jsz p1f
* first part in v1 second in v4
  v7 v1-fv4
  j b17
* first part in v1 second in v3
p1e v7 v1-fv3
p1c a7 a7+a6
  s0 s5-s6
  a5 b03
  a6 b04
  s5 s5-s6
  a5 a5+a6
  jsn p1
  s4 t70
  a5 b07
  s0 s4-s6
  a4 b06
  a7 b10
  a6 b05
  a5 a5+a4
  a7 a7+a6
  jsp p9999
  a0 b71
  a2 b73
  jaz exit
  ,a0,a2 v7
exit a0 btemp
  a7 count
  b00,a7 ,a0
  exit
* first part in v1 no second part
p1f v7 +fv1
  j b17
  end
          ident locat1
*
*   l = locat1(list,nlist,pattern)
*
*  to find member of   list   =   pattern
*  or return zero if no match found
*
*
locat1 enter np=3
  argadd a2,2
  argadd a4,3
  argadd a1,1
  a7 64
  s1 0
  a2 0,a2
  s4 0,a4
  a6 0
  vl a7
  a0 a7-a2
  jap vvv
  a0 a1
  j xxx
yyy s0 vm
  a0 a1
  jsn zzz
xxx v7 ,a0,1
  a2 a2-a7
  a1 a1+a7
  v6 s4-v7
  a0 a7-a2
  vm v6,z
  a6 a6+a7
  jam yyy
  s0 vm
  a0 a2
  jsz www
zzz s6 vm
  a6 a6-a7
aaa a7 zs6
  a6 a6+1
  a6 a6+a7
  s1 a6
out exit
vvv a0 a2
www jaz out
  a0 a1
  vl a2
  v7 ,a0,1
  v6 s4-v7
  vm v6,z
  s0 vm
  s6 vm
  jsn aaa
  exit
  end
          ident mxmb
*
*  coded v.r. saunders   ---- june 87
*  compatible with eam to 16 mw
*  x-mp version --- not cray 1
*
*  call mxmb(a,mcola,mrowa,b,mcolb,mrowb,r,mcolr,mrowr,ncol,link,nrow)
*     r(matrix) = r(matrix) + a(matrix) * b(matrix)
*
  align
  entry mxmb
mxmb a4 12,a6
  a3 10,a6
  a1 11,a6
  a2 6,a6
  a5 9,a6
  a7 8,a6
  a4 0,a4
  a3 0,a3
  a1 0,a1
  a2 0,a2
  s5 0,a5
  a5 5,a6
  a7 0,a7
  s1 1,a6
  s2 7,a6
  a0 a3*a4
  b73 a3
  s0 a1
  b75 a1
  a1 2,a6
  b74 a2
  a2 3,a6
  a5 0,a5
  jaz ret
  jsz ret
  s3 4,a6
  t71 s5
  b77 a7
  t77 s1
  a1 0,a1
  a2 0,a2
  b76 a5
*
*  a1=mcola a2=mrowa
*  s2=<r> destroyed  s3=<b> updated
*  b77=mcolr b76=mcolb b75=link b73=ncol
*  t77=<a>
*  sets
*  t76=<a>(updated in myv128) t75=<a>(updated in bfetch)
*  t74=<b> t72=flag  t70=<r>
*  b72=link(updated in bfetch)   t73=address of last bfetch
*
top s7 0
  b71 a4
  a4 zs0
  t70 s2
  t74 s3
tmyv a5 b75
  s0 s7\s3
  vl a5
  a7 b76
  a6 vl
  a0 s3
  a5 a5-a6
  jsz exitt
  v0 ,a0,a7
  a7 a7*a6
  a6 a2*a6
  t73 s3
  s6 a7
  s5 a6
  vm v0,n
  s4 s1+s5
  s3 s3+s6
exitt a0 a4-a3
  s0 0
  a6 a4+a4
  t76 s1
  t72 s0
  jam tyv128
rx64 s7 vm
  s0 vm
  a7 zs7
  a4 ps7
  vl a3
  b72 a5
  a6 a2*a7
  s7 s7<a7
  jsz ex64
  s7 #sb&s7
  s6 a6
  s7 s7>a7
  s6 s1+s6
  s5 v0,a7
  a5 zs7
  s0 t72
  a7 2
  a0 s6
  a4 a4-a7
  a6 a2*a5
  jsn ix64
  v3 ,a0,a1
  a0 s2
  a3 b77
  s7 s7<a5
  v4 ,a0,a3
  a4 a4-1
  s0 s7
  s7 #sb&s7
  s6 a6
  v7 s5*rv3
  s6 s1+s6
  t72 s1
  s7 s7>a5
  s5 v0,a5
  a0 s6
  a5 zs7
  v1 v4+fv7
  jsz ex64
  a6 a2*a5
kx64 v5 ,a0,a1
  s7 s7<a5
  a4 a4-1
  s0 sb
  s7 #sb&s7
  s6 a6
  s7 s7>a5
  t72 s0
  s0 +a4
  v2 s5*rv5
  s5 v0,a5
  a5 zs7
  s6 s1+s6
  a6 a2*a5
  a0 s6
wx64 jsm xx64
tx64 v3 ,a0,a1
  s6 a6
  s7 s7<a5
  s6 s1+s6
  s7 #sb&s7
  a0 s6
  s7 s7>a5
  v4 ,a0,a1
  v5 s5*rv3
  a3 zs7
  s5 v0,a5
  a4 a4-a7
  a6 a2*a3
  s7 s7<a3
  v1 v1+fv5
  s7 #sb&s7
  s6 a6
  s7 s7>a3
  v6 s5*rv4
  s6 s1+s6
  a5 zs7
  s5 v0,a3
  s0 +a4
  a0 s6
  a6 a2*a5
  v2 v2+fv6
  jsm xx64
  v3 ,a0,a1
  s6 a6
  s7 s7<a5
  s6 s1+s6
  s7 #sb&s7
  a0 s6
  s7 s7>a5
  v5 ,a0,a1
  v4 s5*rv3
  a3 zs7
  s5 v0,a5
  a4 a4-a7
  a6 a2*a3
  s7 s7<a3
  v1 v1+fv4
  s7 #sb&s7
  s6 a6
  s7 s7>a3
  v6 s5*rv5
  s6 s1+s6
  a5 zs7
  s5 v0,a3
  s0 +a4
  a0 s6
  a6 a2*a5
  v2 v2+fv6
  jsp tx64
xx64 s0 s0<63
  jsz ex64
  v3 ,a0,a1
  v5 s5*rv3
  v1 v1+fv5
ex64 a0 b72
  s0 t72
  a5 b72
  a7 b76
  a3 vl
  jan sx64
  a7 b77
  a0 s2
  jsp hx64
  a4 b71
  a5 b74
  a4 a4-1
  s2 t70
  s3 t74
  v6 v1+fv2
  s0 a4
  s4 t71
  s5 a5
  s1 t77
  a3 b73
  s2 s2+s4
  s3 s3+s5
  ,a0,a7 v6
  jsn top
ret j b00
sx64 vl a5
  a0 s3
  a6 vl
  v0 ,a0,a7
  a7 a7*a6
  a4 a2*a6
  t73 s3
  s1 s4
  a5 a5-a6
  s6 a7
  s5 a4
  vm v0,n
  s4 s1+s5
  s3 s3+s6
  j rx64
ix64 jsp kx64
  s0 +a4
  j wx64
hx64 jsz return
  ,a0,a7 v1
  j return
tyv128 a0 a6-a3
  s6 a3
  s7 a3
  s6 s6>1
  t75 s4
  a6 a6+a4
  s5 s7-s6
  jam uyv128
  a4 s6
  a3 s5
rx127 s7 vm
  s0 vm
  a7 zs7
  b72 a5
  a5 a3*a1
  a6 a2*a7
  jsz ex127
  s6 a6
  s4 a5
  s6 s1+s6
  s0 t72
  s5 v0,a7
  vl a3
  a0 s6
  s6 s6+s4
  v3 ,a0,a1
  jsn ux127
  a0 s2
  a3 b77
  s7 s7<a7
  v4 ,a0,a3
  v7 s5*rv3
  s7 #sb&s7
  a0 s6
  s7 s7>a7
  a3 vl
  v1 v4+fv7
  a7 zs7
  vl a4
  s0 s7
  t72 s1
  a6 a2*a7
  v5 ,a0,a1
  s6 a6
  s6 s1+s6
  a0 s6
  s6 s6+s4
  v2 s5*rv5
  vl a3
  jsz ex127
tx127 v3 ,a0,a1
  s5 v0,a7
ux127 v4 s5*rv3
  s7 s7<a7
  a0 s6
  s7 #sb&s7
  s7 s7>a7
  v1 v1+fv4
  vl a4
  a7 zs7
  v5 ,a0,a1
  a6 a2*a7
  s0 s7
  s6 a6
  v6 s5*rv5
  s6 s1+s6
  a0 s6
  s6 s6+s4
  v2 v2+fv6
  vl a3
  jsn tx127
ex127 a0 b72
  s0 t72
  a5 b72
  a7 b76
  vl a5
  jaz sx127
  a0 s3
  a6 vl
  v0 ,a0,a7
  a7 a7*a6
  a5 a5-a6
  a6 a2*a6
  t73 s3
  s1 t75
  s6 a7
  s5 a6
  vm v0,n
  s5 s5+s1
  s3 s3+s6
  t75 s5
  j rx127
sx127 a7 b77
  a0 s2
  vl a3
  jsz return
  ,a0,a7 v1
  a5 a3*a7
  s7 a5
  vl a4
  s7 s7+s2
  a0 s7
  v5 ,a0,a7
  v6 v5+fv2
  ,a0,a7 v6
return a5 b74
  a4 b71
  s2 t70
  s4 t71
  a0 a4-1
  s3 t74
  s5 a5
  a4 a4-1
  s1 t77
  a3 b73
  s2 s2+s4
  s3 s3+s5
  jan top
  j b00
uyv128 a0 a3-a6
  a7 s5
  a6 32
  jap rx128
  a4 a7-a6
rx128 s7 vm
  s0 vm
  a7 zs7
  vl a4
  b72 a5
  a5 a4*a1
  a6 a2*a7
  jsz ex128
  s6 a6
  s4 a5
  s6 s1+s6
  s0 t72
  s5 v0,a7
  a0 s6
  s6 s6+s4
  v3 ,a0,a1
  jsn ux128
  a0 s2
  a4 b77
  s7 s7<a7
  v4 ,a0,a4
  v7 s5*rv3
  s7 #sb&s7
  a0 s6
  s7 s7>a7
  a4 vl
  v1 v4+fv7
  a7 zs7
  s0 s7
  t72 s1
  a6 a2*a7
  v5 ,a0,a1
  s6 a6
  s6 s1+s6
  a0 s6
  s6 s6+s4
  v2 s5*rv5
  jsz ex128
tx128 v3 ,a0,a1
  s5 v0,a7
ux128 a0 s6
  v4 s5*rv3
  v5 ,a0,a1
  s7 s7<a7
  v1 v1+fv4
  s7 #sb&s7
  s7 s7>a7
  a7 zs7
  s0 s7
  a6 a2*a7
  s6 a6
  v6 s5*rv5
  s6 s1+s6
  a0 s6
  s6 s6+s4
  v2 v2+fv6
  jsn tx128
ex128 a0 b72
  s0 t72
  a5 b72
  a7 b76
  jaz sx128
  vl a5
  a0 s3
  a6 vl
  v0 ,a0,a7
  a7 a7*a6
  a5 a5-a6
  a6 a2*a6
  t73 s3
  s1 t75
  s6 a7
  s5 a6
  vm v0,n
  s5 s5+s1
  s3 s3+s6
  t75 s5
  j rx128
sx128 a0 s2
  a7 b77
  jsz return
  ,a0,a7 v1
  a5 a4*a7
  a6 a4+a4
  s7 a5
  s7 s7+s2
  a3 a3-a6
  a5 a1*a6
  a6 a7*a6
  a0 s7
  s1 t76
  v5 ,a0,a7
  a4 zs0
  s6 a5
  s5 a6
  s1 s1+s6
  s2 s2+s5
  v6 v5+fv2
  s7 t73
  s3 t74
  ,a0,a7 v6
  j tmyv
  end
  ident dbuild
* sgmata vrs
* dbuild mods rjh 30/1/87
* differs from sgmata in that indices are packed two to a word
* in /craypk/ with offset of 12 bits (*4096). this eliminates
* need to handle odd and even cases. odd ones deleted.
* additional confusion from loop unrolling etc
*  entry dbuild
* call dbuild(fock,dens)
  common shlt
tol bss 1
cutoff bss 1
icount bss 1
ic4 bss 1
  common blkin
gin bss 340
  common craypk
ind bss 680
  common vectem
ikyk bss 64
ijkl bss 64
l bss 64
k bss 64
  block *
  block *
  block *
  block *
dbuild enter np=2
  s1 icount,0              icount into s1
  a1 gin                   <g> into a1
  s3 1                     s3 = 1
  argadd a5,2              a5 = <p>
  a2 ind                   a2 = <ind>
  argadd a6,1              a6 = <f>
  s1 s1-s3                 s1 = mword = no. of integrals = icount-1
  b71 a1                   b71 = <g>
  b70 a2                   b70 = <ind>
  s4 s1-s3                 s4 = mword-1
  s5 a5                    s5 = <p>
  a6 a6-a5                 a6 = <f>-<p>
  t70 s5                   t70 = <p>
  s4 s4>6                  s4 = no. of vector loads - 1
  a1 s4                    a1 =    ditto
  s4 s4<6
  a0 a1                    a0 =    ditto
  a3 64                    a3 = 64
  s2 s1-s4                 s2 = final vector length
  b74 a3                   b74 = 64
  vl a3                    vl = 64
  a2 s2                    a2 = final vector length
  b75 a2                   b75 =     ditto
aaaa jan www               if last vector set vl accordingly
  a3 b75
  vl a3
www b73 a1                 b73 = no. of vector loads - 1
  b72 a3                   b72 = current vl
* vl set
* a3=vl
* get k,l and start to form iky(k) etc
  a2 b70
  a0 a2+1                  a0 = address of first kl
  a2 2                     a2 = 2 = increment between kl
  a4 12                    a4 = 12 (4096 = 2**12)
* clock 1
  v7 ,a0,a2                v7 = kl
  s5 <12                   s5 = 12 bit mask
  v5 v7>a4                 v5 = k
  a1 24                    a1 = 24
  a0 l                     a0 = <l>
* clock 2
  v4 s5&v7                 v4 = l
  v3 v5<a1
*
* j dbug
*
  v2 v3*fv3                v2 = k*k
* clock 3
  ,a0,1 v4                 store l at <l>
  v1 v2-v5
  v0 v1>1                  v0 = iky(k)
  a0 k                     a0 = <k>
* clock 4
  ,a0,1 v5                 store k at <k>
* load ij into v7
  a0 b70
  v7 ,a0,a2                v7 = ij
  v3 s5&v7                 v3 = j
  v2 v4+v0                 v2 = iky(k)+l
  a0 ikyk                  a0 = <ikyk>
* clock 5
  v1 v7>a4                 v1 = i
  ,a0,1 v0                 store iky(k) at <ikyk>
  a0 b71                   a0 = <g>
* clock 6
  v6 v1<a1
  v5 v6*fv6                v5 = i*i
  v0 ,a0,1                 v0 = g
* clock 7
  v4 v5-v1
  v7 v4>1                  v7 = ikyk(i)
  s7 4.                    s7 = 4.0
  a7 a3-1                  a7 = vl-1
  a0 a3-1                  a0 = vl-1
  s6 t70                   s6 = <p>
* clock 8
  v4 s7*rv0                v4 = 4.0*g
  v5 v7+v3                 v5 = ikyk(i) + j
* v1=i v2=kl v3=j
* v4=g4 v5=ij v6=scratch v7=ikyi
* k,l,ikyk stored
* a0=vl-1  a6=<a-p> a7=vl-1
* t70=<p>
* b70=<ijklijkl> b71=<g> b72=cvl
* b73=nvord-1 b74=64 b75=lrem
* coulomb part(even cases)
pppf s0 a3                 a0 = vl
  b77 a3                   b77 = vl
* clock 9
  v6 s6+v2                 v6 = <p(kl)>
  s0 s0<63                 s0 has sign bit set if vl is odd
* clock 10
  v2 s6+v5                 v2 = <p(ij)>
* v2=<pij> v4=g4 v6=<pkl>
  s1 v6,a7                 s1 = <p(kl)>,vl
  s2 v2,a7                 s2 = <p(ij)>,vl
  jsp ppph                 jump if vl is even
  s3 v4,a7                 s3 = 4.*g,vl
  a1 s1                    a1 = <p(kl)>,vl
  a2 s2                    a2 = <p(ij)>,vl
  s1 -1,a1                 s1 = p(kl),vl
  s2 -1,a2                 s2 = p(ij),vl
  a2 a2+a6                 a2 = <f(ij)>,vl
  a1 a1+a6                 a1 = <f(kl)>,vl
  s4 -1,a2                 s4 = f(ij),vl
  s6 -1,a1                 s6 = f(kl),vl
  s1 s1*rs3                s1 = p(kl),vl * 4.*g,vl
  s2 s2*rs3                s2 = p(ij),vl * 4.*g,vl
  a0 a7-1                  a0 = vl-2
  a7 a7-1                  a7 =  ditto
  s1 s1+fs4                s1 = (f(ij)+p(kl)*4.g),vl
  s2 s2+fs6                s2 =    kl    ij
  -1,a2 s1                 store  f(ij),vl
  -1,a1 s2                 store    kl
  jam pppg                 jump if vl = 1
* above handles loop folded cases (odd no. of terms)
ppp s1 v6,a7
  s2 v2,a7
ppph a5 a7-1               a5 = counter -1 = offset for second element
  s3 v4,a7                 s3 = 4g,1
  s6 v6,a5
  a1 s1
  a2 s2
  s1 -1,a1                 s1 = p(kl),1
  s2 -1,a2                 s2 = p(ij),1
  a2 a2+a6
  a1 a1+a6
  s7 v2,a5
  s4 -1,a2                 s4 = f(ij),1
  s5 -1,a1                 s5 = f(kl),1
  a3 s6
  s1 s1*rs3
  a4 s7
  s2 s2*rs3
  s6 -1,a3                 s6 = p(kl),2
  s7 -1,a4                 s7 = p(ij),2
  s3 v4,a5
  s1 s1+fs4
  s2 s2+fs5
  a4 a4+a6
  a3 a3+a6
  s6 s6*rs3
  s7 s7*rs3
  -1,a2 s1                 store f(ij),1
  -1,a1 s2                 store f(kl),1
  s4 -1,a4                 s4 = f(ij),2
  s5 -1,a3                 s5 = f(kl),2
  a0 a5-1                  a0 = no. remaining-1
  a7 a5-1                  a7 = updated
  s6 s6+fs4
  s7 s7+fs5
  -1,a4 s6                 store f(ij),2
  -1,a3 s7                 store f(kl),2
  jap ppp                  jump if more to do.
* above ppp loop costs about 70 clocks per two elements
*
* exchange part (even cases)
* code based on sequence
* kik=gik*rjl + kik
* kil=gil*rjk + kil
* kjk=gil*ril + kjk
* kjl=gik*rik + kjl
* clock 1
pppg a1 24
  v4 v3<a1
  v5 v4*fv4
  a0 k
* clock 2
  v2 ,a0,1
  v0 v1\v2
  v6 v5-v3
  v4 v6>1
  a0 l
* clock 3
  vm v0,z
  v6 ,a0,1
  v1 v3-v6
* clock 4
  s1 vm
  vm v1,z
  v5 v7+v2
  v0 v6<a1
* clock 5
  s2 vm
  vm v1,m
  ,a0,1 v5
  s1 s1!s2
  a0 k
* clock 6
  v1 v7+v6
  v5 v0*fv0
* clock 7
  v7 v4+v6
  ,a0,1 v1
  a0 ikyk
  a2 2
* clock 8
  v0 v5-v6
  v1 v0>1
* clock 9
  v5 ,a0,1
  v6 v1+v3
  v0 v6!v7&vm
  a0 b71
* clock 10
  v7 ,a0,1
  v6 v4+v2
  a0 l
* clock 11
  v1 v3-v2
  vm v1,m
* clock 12
  v2 v5+v3
  v4 v2!v6&vm
* clock 13
  v5 ,a0,1
  vm v1,z
  v2 v7+fv7
  a7 b77
  a0 k
* clock 14
  v6 ,a0,1
  v3 v2!v7&vm
  a7 a7-1
  vm s1
  s6 t70
* clock 15
  v1 v2!v7&vm
* in the above il stored in k ik stored in l
* clock 16
  v2 s6+v5
* clock 17
  v7 s6+v6
* clock 18
  v5 s6+v4
* clock 19
  v6 s6+v0
* v1=gik v2=<pik> v3=gil v5=<pjk> v6=<pjl> v7=<pil>
qqq s1 v2,a7
  s2 v7,a7
  s3 v5,a7
  s4 v6,a7
  s5 v1,a7
  s6 v3,a7
  a1 s1
  a2 s2
  s1 -1,a1
  s2 -1,a2
  a3 s3
  a4 s4
  s3 -1,a3
  s4 -1,a4
  a4 a4+a6
  a3 a3+a6
  a2 a2+a6
  s7 -1,a4
  s1 s1*rs5
  s2 s2*rs6
  s3 s3*rs6
  s6 -1,a3
  a1 a1+a6
  s4 s4*rs5
  s5 -1,a2
  s1 s7-fs1
  s7 -1,a1
  a0 a7-1
  a7 a7-1
  s2 s6-fs2
  s3 s5-fs3
  s4 s7-fs4
  -1,a4 s1
  -1,a3 s2
  -1,a2 s3
  -1,a1 s4
  jap qqq
* above costs about 105 clocks
  a3 b77
  a7 b71
  a7 a7+a3
  b71 a7
  a7 b70
  a7 a7+a3
  a7 a7+a3
  b70 a7
  a1 b73
  a3 b74
  a0 a1-1
  a1 a1-1
  s0 0
  jap aaaa
dbug  s1 1
  icount,0 s1
  ic4,0 s1
  exit
  end
_IF(1s)
          ident szero
* entry szero
*  call szero(a,n)
szero enter np=2
  argadd a7,2
  argadd a1,1
  a3 64
  a7 0,a7
  vl a3
  a0 a7-a3
  jam clean
  v7 0
step a0 a1
  ,a0,1 v7
  a7 a7-a3
  a1 a1+a3
  a0 a7-a3
  jap step
  a0 a7
  vl a7
  jaz out
  a0 a1
  ,a0,1 v7
out exit
clean a0 a7
  vl a7
  jaz out
  v7 0
  a0 a1
  ,a0,1 v7
  exit
  end
          ident fmove
* entry fmove
*   call fmove(a,b,n)
*   b  =  a  for  n  elements
fmove enter np=3
  argadd a7,3
  argadd a1,2
  argadd a2,1
  a3 64
  a7 0,a7
  vl a3
  a0 a3-a7
  jap clean
* n gt 64
  a0 a2
  v7 ,a0,1
  a7 a7-a3
  a2 a2+a3
  a0 a7-a3
step jam ex1
  a0 a2
  v4 ,a0,1
  a0 a1
  a7 a7-a3
  a2 a2+a3
  ,a0,1 v7
  a0 a7-a3
  a1 a1+a3
  jam ex2
  a0 a2
  v7 ,a0,1
  a0 a1
  a7 a7-a3
  a2 a2+a3
  ,a0,1 v4
  a0 a7-a3
  a1 a1+a3
  j step
ex1 a0 a1
  ,a0,1 v7
  a1 a1+a3
  j clean
ex2 a0 a1
  ,a0,1 v4
  a1 a1+a3
clean a0 a7
  vl a7
  jaz out
  a0 a2
  v6 ,a0,1
  a0 a1
  ,a0,1 v6
out exit
  end
          ident setsto
*
*  call setsto(n,scalar,a)
*
*
setsto enter np=3
  argadd a7,1
  argadd a5,2
  argadd a1,3
  a3 64
  s7 >64
  a7 0,a7
  s5 0,a5
  vl a3
  vm s7
  a0 a7-a3
  jam clean
  v7 s5!v6&vm
step a0 a1
  ,a0,1 v7
  a7 a7-a3
  a1 a1+a3
  a0 a7-a3
  jap step
  a0 a7
  vl a7
  jaz out
  a0 a1
  ,a0,1 v7
out exit
clean a0 a7
  vl a7
  jaz out
  v7 s5!v6&vm
  a0 a1
  ,a0,1 v7
  exit
  end
_ENDIF
_IFN(1s)
          ident upak8z
*
*  call upak8z(n,in,iout,jout,kout,lout)
*  called by list,punch/serv sgmat,jkgen/util proc2/scf
*  not protected for n=0 or -ve
*  optimized for xmp not cray-1
*  compatible with eam to 16 mwords
*
*  coded v.r.saunders     march 1987
*
  align
  entry upak8z
upak8z a7 1,a6
  s2 2,a6
  s6 6,a6
  s5 5,a6
  s4 4,a6
  s3 3,a6
  s7 0,a7
  s1 1
  a2 2
  a1 8
  a3 16
  s7 s7+s1
  s6 s6+s1
  s5 s5+s1
  s4 s4+s1
  s7 s7>1
  s3 s3+s1
  a7 s7
top a0 s2
  vl a7
  s7 <8
  v0 ,a0,1  v0=in
  a0 s6
  s6 s6-s1
  v1 s7&v0  v1=lout(odd)
  ,a0,a2 v1  store lout(odd)
  a0 s5
  s5 s5-s1
  a5 vl
  v2 v0>a1  v2=in>8
  v3 s7&v2  v3=kout(odd)
  ,a0,a2 v3  store kout(odd)
  a0 s4
  s4 s4-s1
  a7 a7-a5
  v4 v0>a3  v4=in>16
  v5 s7&v4  v5=jout(odd)
  ,a0,a2 v5 store jout(odd)
  a0 s3
  s3 s3-s1
  s0 a7
  v6 v2>a3  v6=in>24
  v7 s7&v6  v7=iout(odd)
  ,a0,a2 v7  store iout(odd)
  a0 s6
  v0 v4>a3  v0=in>32
  v1 s7&v0  v1=lout(even)
  ,a0,a2 v1 store lout(even)
  a0 s5
  v2 v6>a3  v2=in>40
  v3 s7&v2  v3=kout(even)
  ,a0,a2 v3  store kout(even)
  a0 s4
  v4 v0>a3  v4=in>48
  v5 s7&v4  v5=jout(even)
  ,a0,a2 v5  save jout(even)
  s7 a5
  a0 s3
  s2 s2+s7
  s7 s7<1
  s7 s7+s1
  s6 s6+s7
  s5 s5+s7
  s4 s4+s7
  s3 s3+s7
  v6 v2>a3  v6=iout(even)
  ,a0,a2 v6  store iout(even)
  jsn top
  j b00
  end
_ENDIF
_IF(j90)
_IFN(sv1)
          ident gather
*
*  call gather(n,r,a,mapa)
*  r = a(mapa) for n elements
*  protected for n=0 but not for -ve n
*  optimized for xmp not cray-1
*  compatible with eam to 16 mwords
*
*  coded v.r.saunders     march 1987
*
  align
  entry gather
gather a7 1,a6
  s3 3,a6
  s4 4,a6
  s2 2,a6
  s6 1
  a1 zs0
  a7 0,a7
  s3 s3-s6
  vl a7
  s0 a7
  a6 vl
  a0 s4
  s1 a6
  a7 a7-a6
  jsz return
top v0 ,a0,1
  a0 s3
  s0 a7
  v1 ,a0,v0
  a0 s2
  s4 s4+s1
  s2 s2+s1
  ,a0,1 v1
  jsz return
  a0 s4
  vl a1
  s1 a1
  v2 ,a0,1
  a0 s3
  a7 a7-a1
  v3 ,a0,v2
  a0 s2
  s0 a7
  s4 s4+s1
  s2 s2+s1
  a7 a7-a1
  ,a0,1 v3
  a0 s4
  jsn top
return j b00
  end
          ident scatter
*
*  call scatter(n,r,mapr,a)
*  r(mapr) = a for n elements
*  protected for n=0 but not for -ve n
*  optimized for xmp not cray-1
*  compatible with eam to 16 mwords
*
*  coded v.r.saunders     march 1987
*
  align
  entry scatter
scatter a7 1,a6
  s2 2,a6
  s4 4,a6
  s3 3,a6
  s6 1
  a1 zs0
  a7 0,a7
  s2 s2-s6
  vl a7
  s0 a7
  a6 vl
  a0 s4
  s1 a6
  a7 a7-a6
  jsz return
top v0 ,a0,1
  a0 s3
  s0 a7
  v1 ,a0,1
  a0 s2
  s4 s4+s1
  s3 s3+s1
  ,a0,v1 v0
  jsz return
  a0 s4
  vl a1
  s1 a1
  v2 ,a0,1
  a0 s3
  a7 a7-a1
  v3 ,a0,1
  a0 s2
  s0 a7
  s4 s4+s1
  s3 s3+s1
  a7 a7-a1
  ,a0,v3 v2
  a0 s4
  jsn top
return j b00
  end
_ENDIF
_ENDIF
_IFN(1s)
          ident setsto
*
*  call setsto(n,scalar,r)
*  r = scalar for n elements
*  protected for n=0 but not for -ve n
*  optimized for xmp not cray-1
*  compatible with eam to 16 mwords
*
*  coded v.r.saunders     march 1987
*
  align
  entry setsto
setsto a4 2,a6
  a7 1,a6
  s2 3,a6
  s6 -1
  a1 0
  vm s6
  vl a1
  s5 64
  s4 0,a4
  a7 0,a7
  a0 s2
  v3 s4!v4&vm
  vl a7
  s0 a7
  a3 vl
  s7 a7
  s3 a3
  jsz return
  ,a0,1 v3
  s0 s7\s3
  s2 s2+s3
  s7 s7-s3
  vl a1
  jsz return
top a0 s2
  s0 s7\s5
  ,a0,1 v3
  s2 s2+s5
  s7 s7-s5
  jsn top
return j b00
  end
          ident szero
*
*  call szero(r,n)
*  r = 0.0 for n elements
*  protected for n=0 but not for -ve n
*  optimized for xmp not cray-1
*  compatible with eam to 16 mwords
*
*  coded v.r.saunders     march 1987
*
  align
  entry szero
szero a7 2,a6
  s2 1,a6
  a1 0
  s5 64
  vl a1
  v3 v0\v0
  a7 0,a7
  a0 s2
  vl a7
  s0 a7
  a4 vl
  s7 a7
  s4 a4
  jsz return
  ,a0,1 v3
  s0 s7\s4
  s2 s2+s4
  s7 s7-s4
  vl a1
  jsz return
top a0 s2
  s0 s7\s5
  ,a0,1 v3
  s2 s2+s5
  s7 s7-s5
  jsn top
return j b00
  end
          ident fmove
*
*  call fmove(a,r,n)
*  r = a for n elements
*  protected for n=0 but not for -ve n
*  optimized for xmp not cray-1
*  compatible with eam to 16 mwords
*
*  coded v.r.saunders     march 1987
*
  align
  entry fmove
fmove a7 3,a6
  s3 1,a6
  s2 2,a6
  a6 zs0
  a7 0,a7
  s0 s2\s3
  jsz return
  vl a7
  s0 a7
  a5 vl
  s7 a7
  s5 a5
  a0 s3
  jsz return
top v0 ,a0,1
  a0 s2
  s0 s7\s5
  s3 s3+s5
  s7 s7-s5
  s2 s2+s5
  ,a0,1 v0
  jsz return
  a0 s3
  vl a6
  s5 a6
  v4 ,a0,1
  a0 s2
  s0 s7\s5
  s3 s3+s5
  s2 s2+s5
  s7 s7-s5
  ,a0,1 v4
  a0 s3
  jsn top
return j b00
  end
          ident ujau
*
*  call ujau(n,iscalar,ir,ia)   ---  ir=ia+iscalar
*
*  protected for n=0 but not n -ve
*  optimized for xmp not cray-1
*  compatible with eam to 16 mwords
*
*  coded v.r.saunders     may 1987
*
  align
  entry ujau
ujau a1 1,a6
  a2 2,a6
  s3 4,a6
  s6 3,a6
  a7 zs0
  a1 0,a1
  s2 0,a2
  vl a1
  s0 a1
  a5 vl
  s1 a1
  s5 a5
  a0 s3
  jsz return
top v6 ,a0,1
  s0 s1\s5
  a0 s6
  s3 s3+s5
  s6 s6+s5
  s1 s1-s5
  v5 s2+v6
  ,a0,1 v5
  jsz return
  vl a7
  a0 s3
  s5 a7
  v3 ,a0,1
  a0 s6
  s0 s1\s5
  s1 s1-s5
  s3 s3+s5
  s6 s6+s5
  v2 s2+v3
  ,a0,1 v2
  a0 s3
  jsn top
return j b00
  end
          ident scaler
*
*  call scaler(n,scalar,r,a)   ---  r=a*scalar
*
*  protected for n=0 but not n -ve
*  optimized for xmp not cray-1
*  compatible with eam to 16 mwords
*
*  coded v.r.saunders     may 1987
*
  align
  entry scaler
scaler a1 1,a6
  a2 2,a6
  s3 4,a6
  s6 3,a6
  a7 zs0
  a1 0,a1
  s2 0,a2
  vl a1
  s0 a1
  a5 vl
  s1 a1
  s5 a5
  a0 s3
  jsz return
top v6 ,a0,1
  s0 s1\s5
  a0 s6
  s3 s3+s5
  s6 s6+s5
  s1 s1-s5
  v5 s2*rv6
  ,a0,1 v5
  jsz return
  vl a7
  a0 s3
  s5 a7
  v3 ,a0,1
  a0 s6
  s0 s1\s5
  s1 s1-s5
  s3 s3+s5
  s6 s6+s5
  v2 s2*rv3
  ,a0,1 v2
  a0 s3
  jsn top
return j b00
  end
          ident swapv
*
*  call swapv(n,a,b)   --- swaps the two vectors a and b
*
*  protected for n=0 but not n -ve
*  optimized for xmp not cray-1
*  compatible with eam to 16 mwords
*
*  coded v.r.saunders     april 1987
*
  align
  entry swapv
swapv a1 1,a6
  s2 2,a6
  s3 3,a6
  a7 zs0
  a1 0,a1
  s7 a7
  s0 a1
  vl a1
  a0 s2
  s4 s2
  a6 vl
  jsz return
  v0 ,a0,1
  a0 s3
  s6 a6
  v1 ,a0,1
  a0 a1-a6
  s5 s3
  s2 s2+s6
  s3 s3+s6
  a1 a1-a6
  jaz exit1
top a0 s2
  vl a7
  v2 ,a0,1
  a0 s3
  a1 a1-a7
  v3 ,a0,1
  a0 s5
  vl a6
  s5 s3
  ,a0,1 v0
  a0 s4
  s0 a1
  s4 s2
  s2 s2+s7
  s3 s3+s7
  a6 a7
  ,a0,1 v1
  jsz exit2
  a0 s2
  vl a7
  v4 ,a0,1
  a0 s3
  a1 a1-a7
  v5 ,a0,1
  a0 s5
  s5 s3
  ,a0,1 v2
  a0 s4
  s0 a1
  s4 s2
  s2 s2+s7
  s3 s3+s7
  ,a0,1 v3
  jsz exit3
  a0 s2
  a1 a1-a7
  v0 ,a0,1
  a0 s3
  s0 a1
  v1 ,a0,1
  a0 s5
  s5 s3
  ,a0,1 v4
  a0 s4
  s4 s2
  ,a0,1 v5
  jsn top
exit1 a0 s5
  vl a6
  cmr
  ,a0,1 v0
  a0 s4
  ,a0,1 v1
return j b00
exit2 a0 s5
  vl a7
  cmr
  ,a0,1 v2
  a0 s4
  ,a0,1 v3
  j b00
exit3 a0 s5
  cmr
  ,a0,1 v4
  a0 s4
  ,a0,1 v5
  j b00
  end
_ENDIF
_IF(j90)
          ident jkadd
*
***      subroutine jkadd(msmall,mbig,f,indf,ibias,p)
***      dimension f(1),indf(1),ibias(1),p(1)
***      ms=msmall
***      mb=mbig
***      ib=0
***1     m=min0(ms,mb)
***cdir$ ivdep
***      do 2 l=1,m
***2     f(indf(ib+l)+ibias(l))=f(indf(ib+l)+ibias(l))+p(ib+l)
***      mb=mb-m
***      ib=ib+ms
***      if(mb.ne.0)goto 1
***      return
***      end
*
*  msmall must be .le. 64
*  optimized for xmp not cray-1
*  compatible with eam to 16 mwords
*
*  coded v.r.saunders     march 1987
*
  align
  entry jkadd
jkadd a1 1,a6
  a2 2,a6
  s5 5,a6
  s3 3,a6
  s4 4,a6
  s6 6,a6
  s1 1
  a1 0,a1
  a2 0,a2
  a0 s5
  s3 s3-s1
  vl a1
  v0 ,a0,1  v0=ibias
  a0 a2-a1
top jap topa
  vl a2
topa a0 s4
  a7 vl
  v7 ,a0,1  v7=indf
  a0 s6
  a2 a2-a7
  s7 a7
  v6 v0+v7  v6=indf+ibias
  v5 ,a0,1  v5=p
  a0 s3
  s0 a2
  v4 ,a0,v6  v4=f
  s4 s4+s7
  s6 s6+s7
  v3 v4+fv5  v3=f+p
  ,a0,v6 v3  scatter f+p
  a0 a2-a1
  jsz return
  jap topb
  vl a2
topb a0 s4
  a7 vl
  v7 ,a0,1  v7=indf
  a0 s6
  a2 a2-a7
  s7 a7
  v2 v0+v7  v2=indf+ibias
  v1 ,a0,1  v1=p
  a0 s3
  s0 a2
  v4 ,a0,v2  v4=f
  s4 s4+s7
  s6 s6+s7
  v3 v4+fv1  v3=f+p
  ,a0,v2 v3  scatter f+p
  a0 a2-a1
  jsn top
return j b00
  end
          ident scatt
*
*  call scatt(n,scalar,v,mapv)
*  v(mapv) = scalar for n elements
*  protected for n=0 but not for n -ve
*  optimized for xmp not cray-1
*  compatible with eam to 16 mwords
*
*  coded v.r.saunders     march 1987
*
  align
  entry scatt
scatt a4 2,a6
  a7 1,a6
  s1 3,a6
  s2 4,a6
  s5 -1
  a6 zs0
  s4 0,a4
  a7 0,a7
  vm s5
  vl a6
  s1 s1+s5
  v3 s4!v2&vm
  vl a7
  s0 a7
  a5 vl
  a0 s2
  a7 a7-a5
  s5 a5
  jsz return
top v5 ,a0,1
  a0 s1
  s0 a7
  s2 s2+s5
  ,a0,v5 v3
  jsz return
  a0 s2
  vl a6
  s5 a6
  v2 ,a0,1
  a7 a7-a6
  s2 s2+s5
  a0 s1
  s0 a7
  a7 a7-a6
  ,a0,v2 v3
  a0 s2
  jsn top
return j b00
  end
          ident srotg
*
*  call srotg(n,a,b,mapb,cos,sin)
*  a = a*cos + b*sin  and  b = b*cos - a*sin for n elements
*  where b is gathered/scattered under control of mapb
*  protected for n=0 but not for n -ve
*  optimized for xmp not cray-1
*  compatible with eam to 16 mwords
*
*  coded v.r.saunders     march 1987
*
  align
  entry srotg
srotg a7 1,a6      a7=<n>
  a2 6,a6      a2=<sin>
  a1 5,a6      a1=<cos>
  s4 3,a6      s4=<b>
  s3 2,a6      s3=<a>
  s5 4,a6      s5=<mapb>
  s6 1
  a7 0,a7      a7=n
  s2 0,a2      s2=sin
  s1 0,a1      s1=cos
  s4 s4-s6     s4=<b>-1
  s0 a7
  vl a7
  a0 s3
  a5 vl
  jsz return
top v7 ,a0,1     v7=a
  a0 s5        a0=<mapb>
  a7 a7-a5
  v6 ,a0,1     v6=mapb
  v5 s2*rv7    v5=a*sin
  a0 s4        a0=<b>-1
  s6 a5
  v4 ,a0,v6    v4=b(mapb)
  v3 s1*rv4    v3=b(mapb)*cos
  s0 a7
  v2 v3-fv5    v2=updated b(mapb)
  s5 s5+s6     update <mapb>
  v1 s1*rv7    v1=a*cos
  ,a0,v6 v2    store updated b(mapb)
  a0 s3        a0=<a>
  s3 s3+s6     update <a>
  v0 s2*rv4    v0=b(mapb)*sin
  v3 v0+fv1    v3=updated a
  ,a0,1 v3     store updated a
  vl a7
  a0 s3
  a5 vl
  jsn top
return j b00
  end
          ident srotv
*
*  call srotv(n,a,b,cos,sin)
*  a = a*cos + b*sin  and  b = b*cos - a*sin for n elements
*  not protected for n=0 or -ve
*  optimized for xmp not cray-1
*  compatible with eam to 16 mwords
*
*  coded v.r.saunders     march 1987
*
  align
  entry srotv
srotv a7 1,a6
  a2 5,a6
  a1 4,a6
  s4 3,a6
  s3 2,a6
  a7 0,a7
  s2 0,a2
  s1 0,a1
top a0 s4
  vl a7
  v7 ,a0,1
  a0 s3
  a6 vl
  v6 ,a0,1
  a7 a7-a6
  s6 a6
  v5 s2*rv7
  v4 s1*rv6
  s0 a7
  v3 v4+fv5
  ,a0,1 v3
  a0 s4
  s3 s3+s6
  s4 s4+s6
  v2 s1*rv7
  v1 s2*rv6
  v0 v2-fv1
  ,a0,1 v0
  jsn top
  j b00
  end
_ENDIF
_IFN(1s)
          ident triad
*
*  call triad(n,scalar,r,a)
*  r = a*scalar + r for n elements
*  protected for n=0 but not for -ve n
*  optimized for xmp not cray-1
*  compatible with eam to 16 mwords
*
*  coded v.r.saunders     march 1987
*
  align
  entry triad
triad a7 1,a6
  a5 2,a6
  s3 4,a6
  s2 3,a6
  a6 zs0
  a7 0,a7
  s5 0,a5
  vl a7
  a0 a7
  a1 vl
  s0 s5
  s1 a1
  a7 a7-a1
  jaz return
  a0 s3
  jsz return
top v3 ,a0,1
  a0 s2
  s0 a7
  v2 ,a0,1
  v0 s5*rv3
  v1 v2+fv0
  s3 s3+s1
  s2 s2+s1
  ,a0,1 v1
  jsz return
  a0 s3
  vl a6
  s1 a6
  v7 ,a0,1
  a0 s2
  a7 a7-a6
  s3 s3+s1
  v6 ,a0,1
  v4 s5*rv7
  s0 a7
  s2 s2+s1
  v5 v6+fv4
  a7 a7-a6
  ,a0,1 v5
  a0 s3
  jsn top
return j b00
  end
          ident gtriad
*
*  call gtriad(n,scalar,r,b,a)
*  r = a*scalar + b for n elements
*  protected for n=0 but not for -ve n
*  optimized for xmp not cray-1
*  compatible with eam to 16 mwords
*
*  coded v.r.saunders     march 1987
*
  align
  entry gtriad
gtriad a7 1,a6
  a5 2,a6
  s3 5,a6
  s4 4,a6
  s2 3,a6
  a6 zs0
  a7 0,a7
  s5 0,a5
  vl a7
  s0 a7
  a1 vl
  a0 s3
  a7 a7-a1
  s1 a1
  jsz return
top v3 ,a0,1
  a0 s4
  s3 s3+s1
  v2 ,a0,1
  a0 s2
  v0 s5*rv3
  s0 a7
  v1 v0+fv2
  s4 s4+s1
  s2 s2+s1
  ,a0,1 v1
  jsz return
  a0 s3
  vl a6
  a7 a7-a6
  v7 ,a0,1
  a0 s4
  s1 a6
  v6 ,a0,1
  v4 s5*rv7
  s0 a7
  a0 s2
  a7 a7-a6
  s3 s3+s1
  s4 s4+s1
  s2 s2+s1
  v5 v4+fv6
  ,a0,1 v5
  a0 s3
  jsn top
return j b00
  end
          ident gtrian
*
*  call gtrian(n,scalar,r,b,a)
*  r = a*scalar - b for n elements
*  protected for n=0 but not for -ve n
*  optimized for xmp not cray-1
*  compatible with eam to 16 mwords
*
*  coded v.r.saunders     march 1987
*
  align
  entry gtrian
gtrian a7 1,a6
  a5 2,a6
  s3 5,a6
  s4 4,a6
  s2 3,a6
  a6 zs0
  a7 0,a7
  s5 0,a5
  vl a7
  s0 a7
  a1 vl
  a0 s3
  a7 a7-a1
  s1 a1
  jsz return
top v3 ,a0,1
  a0 s4
  s3 s3+s1
  v2 ,a0,1
  v0 s5*rv3
  a0 s2
  s0 a7
  v1 v0-fv2
  s4 s4+s1
  s2 s2+s1
  ,a0,1 v1
  jsz return
  a0 s3
  vl a6
  a7 a7-a6
  v7 ,a0,1
  a0 s4
  s1 a6
  v6 ,a0,1
  v4 s5*rv7
  s0 a7
  a0 s2
  a7 a7-a6
  s3 s3+s1
  s4 s4+s1
  s2 s2+s1
  v5 v4-fv6
  ,a0,1 v5
  a0 s3
  jsn top
return j b00
  end
          ident vvtv
*
*  call vvtv(n,r,a,b)
*  r = a * b for n elements
*  protected for n=0 but not for -ve n
*  optimized for xmp not cray-1
*  compatible with eam to 16 mwords
*
*  coded v.r.saunders     march 1987
*
  align
  entry vvtv
vvtv a7 1,a6
  s3 3,a6
  s4 4,a6
  s2 2,a6
  a6 zs0
  a7 0,a7
  vl a7
  s0 a7
  a5 vl
  s7 a7
  a0 s3
  s5 a5
  jsz return
top v3 ,a0,1
  a0 s4
  s0 s7\s5
  v2 ,a0,1
  a0 s2
  v1 v2*rv3
  s3 s3+s5
  s7 s7-s5
  s4 s4+s5
  s2 s2+s5
  ,a0,1 v1
  jsz return
  a0 s3
  vl a6
  s5 a6
  v7 ,a0,1
  a0 s4
  s0 s7\s5
  v6 ,a0,1
  a0 s2
  s3 s3+s5
  v5 v6*rv7
  s7 s7-s5
  s4 s4+s5
  s2 s2+s5
  ,a0,1 v5
  a0 s3
  jsn top
return j b00
  end
_IFN(charmm)
          ident addvec
*
*  call addvec(r,a,b,n)
*  r = a + b for n elements
*  protected for n=0 but not for -ve n
*  optimized for xmp not cray-1
*  compatible with eam to 16 mwords
*
*  coded v.r.saunders     march 1987
*
  align
  entry addvec
addvec a7 4,a6
  s3 2,a6
  s4 3,a6
  s2 1,a6
  a6 zs0
  a7 0,a7
  vl a7
  s0 a7
  a5 vl
  s7 a7
  a0 s3
  s5 a5
  jsz return
top v3 ,a0,1
  a0 s4
  s0 s7\s5
  v2 ,a0,1
  a0 s2
  v1 v2+fv3
  s3 s3+s5
  s7 s7-s5
  s4 s4+s5
  s2 s2+s5
  ,a0,1 v1
  jsz return
  a0 s3
  vl a6
  s5 a6
  v7 ,a0,1
  a0 s4
  s0 s7\s5
  v6 ,a0,1
  a0 s2
  s3 s3+s5
  v5 v6+fv7
  s7 s7-s5
  s4 s4+s5
  s2 s2+s5
  ,a0,1 v5
  a0 s3
  jsn top
return j b00
  end
_ENDIF
          ident subvec
*
*  call subvec(r,a,b,n)
*  r = a + b for n elements
*  protected for n=0 but not for -ve n
*  optimized for xmp not cray-1
*  compatible with eam to 16 mwords
*
*  coded v.r.saunders     march 1987
*
  align
  entry subvec
subvec a7 4,a6
  s3 2,a6
  s4 3,a6
  s2 1,a6
  a6 zs0
  a7 0,a7
  vl a7
  s0 a7
  a5 vl
  s7 a7
  a0 s3
  s5 a5
  jsz return
top v3 ,a0,1
  a0 s4
  s0 s7\s5
  v2 ,a0,1
  a0 s2
  v1 v3-fv2
  s3 s3+s5
  s7 s7-s5
  s4 s4+s5
  s2 s2+s5
  ,a0,1 v1
  jsz return
  a0 s3
  vl a6
  s5 a6
  v7 ,a0,1
  a0 s4
  s0 s7\s5
  v6 ,a0,1
  a0 s2
  s3 s3+s5
  v5 v7-fv6
  s7 s7-s5
  s4 s4+s5
  s2 s2+s5
  ,a0,1 v5
  a0 s3
  jsn top
return j b00
  end
_ENDIF
_IF(j90)
          ident jkuno
*
*  call jkuno(n,p,g,i,j,k,l,indf,pg)
*
*  to generate triangle matrix indices ik,jl,il,jk,ij,kl
*  and p*g contributions for k<p> j<p>
*
*  not protected for n=0
*  optimized for xmp not cray-1
*  compatible with eam to 16 mwords
*  but note that indf,pg must lie within first 8mword
*
*  coded v.r.saunders     april 1987
*
  common mapper
iky bss 1
  block *
  align
  entry jkuno
jkuno a1 1,a6
  s7 2,a6
  s1 3,a6
  s2 4,a6
  s3 5,a6
  s4 6,a6
  s5 7,a6
  a1 0,a1
  a3 8,a6
  a2 9,a6
  s6 1
  a4 6
  s7 s7-s6
  a5 iky-1
  a7 s7
  b77 a7
*
*  a1=n a2=<pg> a3=<indf> a4=6 a5=<iky>-1 b77=<p>-1
*  s1=<g> s2=<i> s3=<j> s4=<k> s5=<l>
*
top vl a1
  a0 s2
  a6 4
  v0 ,a0,1                 v0=i
  a0 s4
  a7 5
  v1 ,a0,1                 v1=k
  a0 s3
  v2 v0-v1                 v2=i-k
  vm v2,z                  vm=i.eq.k
* v0=i v1-k v2=i-k vm=ieqk
  v3 ,a0,1                 v3=j
  a0 a5
  v4 ,a0,v0                v4=ikyi
  v5 v3-v1                 v5=j-k
  a0 s5
  s6 vm                    s6=i.eq.k
  vm v5,z                  vm=j.eq.k
* v0=i v1=k v2=i-k v3=j v4=ikyi v5=j-k vm=jeqk
  v6 ,a0,1                 v6=l
  a0 a5
  v7 ,a0,v1                v7=ikyk
  v2 v4+v3                 v2=ij
  a0 a3+a6
  s7 vm                    s7=j.eq.k
  vm v5,p                  vm=j.ge.k
  ,a0,a4 v2                store ij
* v0=i v1=k v2=ij v3=j v4=ikyi v5=j-k v6=l v7=ikyk vm=jgek
  a0 s1
  v0 ,a0,1                 v0=g
  a0 b77
  v2 ,a0,v2                v2=pij
  v5 v7+v6                 v5=kl
  a0 a3+a7
  ,a0,a4 v5                store kl
* v0=g v1=k v2=pij v3=j v4=ikyi v5=kl v6=l v7=ikyk vm=jgek
  a0 b77
  v6 v6+v4                 v6=il
  v5 ,a0,v5                v5=pkl
  v1 v0*rv2                v1=jkl/2
  a0 a2+a7
  a7 3
  ,a0,a4 v1                store jkl/2
* v0=g v1=jkl v2=pij v3=j v4=ikyi v5=pkl v6=il v7=ikyk vm=jgek
  a0 a5
  v2 ,a0,v3                v2=ikyj
  v1 v0*rv5                v1=jij/2
  a0 a2+a6
  a6 2
  ,a0,a4 v1                store jij/2
* v0=g v1=jij v2=ikyj v3=j v4=ikyi v5=pkl v6=il v7=ikyk vm=jgek
  a0 s4
  v5 ,a0,1                 v5=k
  a0 a3+a6
  v7 v7+v3                 v7=xkj
  ,a0,a4 v6                store il
* v0=g v1=jij v2=ikyj v3=j v4=ikyi v5=k v6=il v7=xkj vm=jgek
  a0 b77
  v1 v2+v5                 v1=xjk
  v3 v1!v7&vm              v3=jk
  v6 ,a0,v6                v6=pil
  a0 a3+a7
  ,a0,a4 v3                store jk
* v0=g v1=xjk v2=ikyj v3=jk v4=ikyi v5=k v6=pil v7=xkj vm=jgek
  a0 b77
  vm s7                    vm=j.eq.k
  v4 v4+v5                 v4=ik
  v7 v0+fv0                v7=g2
  v3 ,a0,v3                v3=pjk
* v0=g v1=xjk v2=ikyj v3=pjk v4=ik v5=k v6=pil v7=g2 vm=jeqk
  a0 a2+a7
  a7 vl
  v1 v7!v0&vm              v1=gil
  v5 v1*rv6                v5=kjk
  ,a0,a4 v5                store kjk
* v0=g v1=gil v2=ikyj v3=pjk v4=ik v5=kjk v6=pil v7=g2 vm=jeqk
  a0 s5
  v7 ,a0,1                 v7=l
  a0 s3
  v6 ,a0,1                 v6=j
  a0 a3
  v3 v1*rv3                v3=kil
  v5 v6-v7                 v5=j-l
  vm v5,z                  vm=j.eq.l
  ,a0,a4 v4                 store ik
* v0=g v1=gil v2=ikyj v3=kil v4=ik v5=j-l v6=j v7=l vm=jeql
  a0 a5
  a1 a1-a7
  v1 ,a0,v7                v1=ikyl
  a0 a2+a6
  a6 a7*a4
  s0 a1
  v6 v1+v6                 v6=xlj
  s7 vm                    s7=j.eq.l
  vm v5,p                  vm=j.ge.l
  ,a0,a4 v3                store kil
* v0=g v1=ikyl v2=ikyj v3=kil v4=ik v5=j-l v6=xlj v7=l vm=jgel
  a0 b77
  s6 s6!s7                 s6=i.eq.k.or.j.eq.l
  v4 ,a0,v4                v4=pik
  v5 v0+fv0                v5=g2
  v1 v2+v7                 v1=xjl
  a0 a3+1
  v3 v1!v6&vm              v3=jl
  vm s6                    vm=i.eq.k.or.j.eq.l
  s7 a7
  ,a0,a4 v3                store jl
* v0=g v1=xjl v2=ikyj v3=jl v4=pik v5=g2 v6=xlj v7=l vm=ieqkorjeql
  a0 b77
  v7 ,a0,v3                v7=pjl
  v6 v5!v0&vm              v6=gik
  v2 v6*rv4                v2=kjl
  a0 a2+1
  ,a0,a4 v2                store kjl
* v0=g v1=xjl v2=kjl v3=jl v4=pik v5=g2 v6=gik v7=pjl
  a0 a2
  s1 s1+s7
  s2 s2+s7
  s3 s3+s7
  s4 s4+s7
  s5 s5+s7
  v5 v6*rv7                v5=kik
  a3 a3+a6
  a2 a2+a6
  ,a0,a4 v5                store kik
  jsn top
  j b00
  end
_ENDIF
_ELSE
*
*   opmerkingen voor c90
*
*
*ident   entry   opmerkingen
*------  ------  -------------------------------------
*setsto  setsto  aangepast voor c90
*dbuild  dbuild  aangepast, c90 max vectorlengte=64
*sgmata  sgmata  aangepast, c90 max vectorlengte=64
*proc2   proc2   aangepast, c90 max vectorlengte=64
*square  square  ok
*sqtrip  sqtrip  ok
*upak8z  upak8z  ok
*revpri  revpri  ok, voor npb<=64, c90 en ymp
*        locat   ok, voor npb<=64, c90 en ymp
*upackh  upackh  ok
*locate  locate  ok
*szero    szero   aangepast, c90 max vectorlengte=64
*fmove   fmove   ok
*ujau    ujau    ok
*scalar  scalar  ok
*swapv   swapv   ok
*jkadd   jkadd   aangepast, c90 max vectorlengte=64
*scatt   scatt   aangepast, c90 max vectorlengte=64
*srotg   srotg   ok
*srotv   srotv   ok
*triad   triad   ok
*gtriad  gtriad  ok
*gtrian  gtrian  ok
*vvtv    vvtv    ok
*addvec adddvec  ok
*subvec  subvec  ok
*jkuno   jkuno   aangepast, c90 max vectorlente=64
*gather  gather  was uitgecomment, lijkt ok
*scatter scatter was uitgecomment, lijkt ok
*
*mxmb/mxmbn van vic maart 1994
*
          ident revpri
*
*   coded      v. r. saunders    july 1982
*
  common disc
incr bss 1
npb bss 1
npbsz bss 1
irep bss 1
npbm1 bss 1
ipos bss 40
nam bss 40
  common discne
i40 bss 40
length bss 40
iostat bss 40
iobuff bss 40
itab bss 1600
istabf bss 10
ipribf bss 9
junibf bss 9
jpbbf bss 9
junjpb bss 9
ibfbas bss 9
igetm bss 9
iputm bss 9
mbuff bss 1
kunkpb bss 1
lun bss 1
l40 bss 1
latm bss 1
lbl bss 1
lbfbas bss 1
lbuff bss 1
  block *
  block *
*    fortran equivalent for revpri
*
*       ix=ipribf(mbuff)
*       do 1 loop=1,npb
*       if(ix.gt.ipribf(loop))ipribf(loop)=ipribf(loop)+1
* 1     continue
*       ipribf(mbuff)=0
* entry revpri
revpri enter np=0
  a6 npb,0
  a7 mbuff,0
  a0 ipribf
  s1 1
  vl a6
  a5 a7-1
  v7 ,a0,1
  v6 s1+v7
  s2 v7,a5
  v5 s2-v6
  vm v5,p
  v4 v6!v7&vm
  v4,a5 0
  ,a0,1 v4
  exit
*
*
*
  bss 2
* entry locat
*
*   m=locat(list,test)
*
locat enter np=2
  argadd a2,2
  a1 npb,0
  argadd a6,1
  a0 a6
  s2 0,a2
  vl a1
  v7 ,a0,1
  v6 s2-v7
  vm v6,z
  s1 vm
  s0 vm
  a7 zs1
  jsz retl
  a7 a7+1
  s1 a7
retl exit
  end
          ident upackh
* entry upackh
upackh enter np=3
  argadd a1,1
  argadd a2,2
  argadd a3,3
  s2 <24
  s1 0,a1
  s3 s1&s2
  s1 s1>32
  0,a3 s3
  0,a2 s1
  exit
  end
          ident locate
* entry locate
*
*   l = locate(list,nlist,ispace,pattern)
*
*  to find member of   list   =   pattern
*  or return zero if no match found
*
locate enter np=4
  argadd a2,2
  argadd a3,3
  argadd a4,4
  argadd a1,1
  a7 64
  s1 0
  a2 0,a2
  a3 0,a3
  s4 0,a4
  a6 0
  vl a7
  a0 a7-a2
  a5 a3*a7
  jap vvv
  a0 a1
  j xxx
yyy s0 vm
  a0 a1
  jsn zzz
xxx v7 ,a0,a3
  a2 a2-a7
  a1 a1+a5
  v6 s4-v7
  a0 a7-a2
  vm v6,z
  a6 a6+a7
  jam yyy
  s0 vm
  a0 a2
  jsz www
zzz s6 vm
  a6 a6-a7
aaa a7 zs6
  a6 a6+1
  a6 a6+a7
  s1 a6
out exit
vvv a0 a2
www jaz out
  a0 a1
  vl a2
  v7 ,a0,a3
  v6 s4-v7
  vm v6,z
  s0 vm
  s6 vm
  jsn aaa
  exit
  end
          ident square
*  entry square
*  call square(r,a,mrowr,n)
square enter np=4
  argadd a7,4
  argadd a6,3
  argadd a5,2
  argadd a4,1
  a3 64
  a7 0,a7
  a6 0,a6
  a0 a7
  a1 1
  a2 a4
  b70 a7
  jaz done
  a7 a6*a3
loop1 a0 a1-a3
  b71 a1
  b72 a4
  b73 a2
  vl a3
loop2 jam ex1
  a0 a5
  v7 ,a0,1
  a0 a4
  a1 a1-a3
  a5 a5+a3
  ,a0,1 v7
  a0 a2
  a4 a4+a3
  ,a0,a6 v7
  a0 a1-a3
  a2 a2+a7
  j loop2
ex1 a0 a1
  vl a1
  jaz term
  a0 a5
  v7 ,a0,1
  a0 a4
  a5 a5+a1
  ,a0,1 v7
  a0 a2
  ,a0,a6 v7
term a1 b71
  a2 b70
  a4 b72
  a0 a1-a2
  a2 b73
  a1 a1+1
  a4 a4+a6
  a2 a2+1
  jan loop1
done exit
  end
          ident sqtrip
*  entry sqtrip
*  call sqtrip (r , a , n)
sqtrip enter np=3
  argadd a7,3
  argadd a5,1
  a1 64
  a7 0,a7
  a4 a5+1
  a0 a1-a7
  vl a7
  a2 a7+1
  a3 a5+a7
  jam ex1
  v7 0
  j ex4
ex1 vl a1
  v7 0
  b70 a7
ex2 a0 a5
  ,a0,a2 v7
  a7 a7-a1
  a6 a2*a1
  a0 a7-a1
  a5 a5+a6
  jap ex2
  a0 a7
  vl a7
  a7 b70
  jaz ex3
ex4 a0 a5
  ,a0,a2 v7
ex3 argadd a6,2
  a2 a7*a1
  a5 1
loop1 a0 a5-a1
  b70 a5
  b71 a4
  b72 a3
  vl a1
  jap loop2
  vl a5
  j ex6
loop2 a0 a6
  v2 ,a0,1
  a5 a5-a1
  a6 a6+a1
  a0 a3
  v1 -fv2
  ,a0,1 v2
  a0 a4
  a3 a3+a1
  ,a0,a7 v1
  a0 a5-a1
  a4 a4+a2
  jap loop2
  a0 a5
  vl a5
  jaz ex7
ex6 a0 a6
  v2 ,a0,1
  a6 a6+a5
  a0 a3
  v1 -fv2
  ,a0,1 v2
  a0 a4
  ,a0,a7 v1
ex7 a5 b70
  a4 b71
  a3 b72
  a5 a5+1
  a4 a4+1
  a3 a3+a7
  a0 a5-a7
  jan loop1
  exit
  end
**    sgmata c90 version (hans nelemans (cray), march 1994)
  ident sgmata
*  entry sgmata
* call sgmata(fock,dens)
  common blkin
gin bss 340
ijin bss 170
mword bss 1
  common vectem
ikyk bss 64
ijkl bss 64
l bss 64
k bss 64
  block *
  block *
sgmata enter np=2
  s1 mword,0
  a1 gin
  argadd a5,2
  a2 ijin
  argadd a6,1
  b71 a1
  b70 a2
  s6 s1
  s1 s1>1
  s3 1
  s7 s6-s1
  s5 a5
  a6 a6-a5
  s4 s7-s3
  t70 s5
  s5 s7\s1
  s4 s4>6
  t77 s5
  a1 s4
  s4 s4<6
  a3 64
  s2 s7-s4
  a0 s2
  jan nfull
  s2 64
nfull a0 a1
  b74 a3
  vl a3
  a2 s2
  s0 0
  b75 a2
aaaa jan www
  a3 b75
  s0 t77
  vl a3
www b73 a1
  b72 a3
* vl set
* a3=vl
* s0=svl-lvl
* even cases
  a0 b70
  a4 8
* clock 1
  v7 ,a0,1
  s5 <8
  v6 v7>a4
  v5 s5&v6
  a1 24
  a0 l
* clock 2
  v4 s5&v7
  v3 v5<a1
  v2 v3*fv3
* clock 3
  ,a0,1 v4
  v1 v2-v5
  v0 v1>1
  a0 k
* clock 4
  ,a0,1 v5
  v7 v6>a4
  v3 s5&v7
  v2 v4+v0
  a7 b71
  a2 2
  a0 ikyk
* clock 5
  v6 v7>a4
  v1 s5&v6
  ,a0,1 v0
  a0 a7+1
* clock 6
  v7 v1<a1
  v5 v7*fv7
  v0 ,a0,a2
* clock 7
  v4 v5-v1
  v7 v4>1
  s7 4.
  a7 a3-1
  a0 a3-1
  s6 t70
* clock 8
  v4 s7*rv0
  v5 v7+v3
  v0 v6>a4
* v0=ijkl v1=i v2=kl v3=j
* v4=g4 v5=ij v6=scratch v7=ikyi
* k,l,ikyk stored
  jsz pppf
  jaz ttt
  a3 a7
  a7 a7-1
* s0=svl-lvl
* a0=vl-1  a6=<a-p> a7=vl-1
* t70=<p>
* t77=lvl-svl(for last time round)
* b70=<ijklijkl> b71=<g> b72=cvl
* b73=nvord-1 b74=64 b75=lrem
* coulomb part(even cases)
pppf s0 a3
  b77 a3
* clock 9
  v6 s6+v2
  s0 s0<63
* clock 10
  v2 s6+v5
* v2=<pij> v4=g4 v6=<pkl>
  s1 v6,a7
  s2 v2,a7
  jsp ppph
  s3 v4,a7
  a1 s1
  a2 s2
  s1 -1,a1
  s2 -1,a2
  a2 a2+a6
  a1 a1+a6
  s4 -1,a2
  s6 -1,a1
  s1 s1*rs3
  s2 s2*rs3
  a0 a7-1
  a7 a7-1
  s1 s1+fs4
  s2 s2+fs6
  -1,a2 s1
  -1,a1 s2
  jam pppg
* above handles loop folded cases (odd no. of terms)
ppp s1 v6,a7
  s2 v2,a7
ppph a5 a7-1
  s3 v4,a7
  s6 v6,a5
  a1 s1
  a2 s2
  s1 -1,a1
  s2 -1,a2
  a2 a2+a6
  a1 a1+a6
  s7 v2,a5
  s4 -1,a2
  s5 -1,a1
  a3 s6
  s1 s1*rs3
  a4 s7
  s2 s2*rs3
  s6 -1,a3
  s7 -1,a4
  s3 v4,a5
  s1 s1+fs4
  s2 s2+fs5
  a4 a4+a6
  a3 a3+a6
  s6 s6*rs3
  s7 s7*rs3
  -1,a2 s1
  -1,a1 s2
  s4 -1,a4
  s5 -1,a3
  a0 a5-1
  a7 a5-1
  s6 s6+fs4
  s7 s7+fs5
  -1,a4 s6
  -1,a3 s7
  jap ppp
* above ppp loop costs about 70 clocks per two elements
*
* exchange part (even cases)
* code based on sequence
* kik=gik*rjl + kik
* kil=gil*rjk + kil
* kjk=gil*ril + kjk
* kjl=gik*rik + kjl
* clock 1
pppg a0 ijkl
  a1 24
  ,a0,1 v0
  v4 v3<a1
  v5 v4*fv4
  a0 k
* clock 2
  v2 ,a0,1
  v0 v1\v2
  v6 v5-v3
  v4 v6>1
  a0 l
* clock 3
  vm v0,z
  v6 ,a0,1
  v1 v3-v6
* clock 4
  s1 vm
  vm v1,z
  v5 v7+v2
  v0 v6<a1
* clock 5
  s2 vm
  vm v1,m
  ,a0,1 v5
  s1 s1!s2
  a0 k
* clock 6
  v1 v7+v6
  v5 v0*fv0
* clock 7
  v7 v4+v6
  ,a0,1 v1
  a0 ikyk
  a2 2
* clock 8
  v0 v5-v6
  v1 v0>1
  a7 b71
* clock 9
  v5 ,a0,1
  v6 v1+v3
  v0 v6!v7&vm
  a0 a7+1
* clock 10
  v7 ,a0,a2
  v6 v4+v2
  a0 l
* clock 11
  v1 v3-v2
  vm v1,m
* clock 12
  v2 v5+v3
  v4 v2!v6&vm
* clock 13
  v5 ,a0,1
  vm v1,z
  v2 v7+fv7
  a7 b77
  a0 k
* clock 14
  v6 ,a0,1
  v3 v2!v7&vm
  a7 a7-1
  vm s1
  s6 t70
* clock 15
  v1 v2!v7&vm
* in the above il stored in k ik stored in l
* clock 16
  v2 s6+v5
* clock 17
  v7 s6+v6
* clock 18
  v5 s6+v4
* clock 19
  v6 s6+v0
* v1=gik v2=<pik> v3=gil v5=<pjk> v6=<pjl> v7=<pil>
qqq s1 v2,a7
  s2 v7,a7
  s3 v5,a7
  s4 v6,a7
  s5 v1,a7
  s6 v3,a7
  a1 s1
  a2 s2
  s1 -1,a1
  s2 -1,a2
  a3 s3
  a4 s4
  s3 -1,a3
  s4 -1,a4
  a4 a4+a6
  a3 a3+a6
  a2 a2+a6
  s7 -1,a4
  s1 s1*rs5
  s2 s2*rs6
  s3 s3*rs6
  s6 -1,a3
  a1 a1+a6
  s4 s4*rs5
  s5 -1,a2
  s1 s7-fs1
  s7 -1,a1
  a0 a7-1
  a7 a7-1
  s2 s6-fs2
  s3 s5-fs3
  s4 s7-fs4
  -1,a4 s1
  -1,a3 s2
  -1,a2 s3
  -1,a1 s4
  jap qqq
* above costs about 105 clocks
* odd cases
  a0 ijkl
  a4 8
* clock 1
  v0 ,a0,1
  s5 <8
ttt v6 v0>a4
  v5 s5&v6
  a7 b70
  a3 b72
  a1 24
  a0 l
  a7 a7+a3
* clock 2
  v4 s5&v0
  v3 v5<a1
  v2 v3*fv3
  b70 a7
* clock 3
  ,a0,1 v4
  v1 v2-v5
  v7 v1>1
  a0 k
* clock 4
  ,a0,1 v5
  v0 v6>a4
  v3 s5&v0
  v2 v4+v7
  a0 b71
  a2 2
  s7 4.
* clock 5
  v1 ,a0,a2
  v5 s7*rv1
  v6 v0>a4
* clock 6
  v1 v7+v3
  v4 v6<a1
  v0 v4*fv4
  a7 a3-1
  a0 ikyk
* clock 7
  ,a0,1 v1
  v7 v0-v6
  v4 v7>1
  s6 t70
* clock 8
  v0 v4+v3
  s0 a3
* l,k stored
* ikykj stored in ikyk
* v1=ikykj v2=kl v3=j v4=ikyi v5=g4
* v6=i v0=ij
* coulomb part(odd cases)
* clock 9
  v7 s6+v2
  s0 s0<63
* clock 10
  v2 s6+v0
  s1 v7,a7
  s2 v2,a7
  jsp pph
  s3 v5,a7
* v2=<pij> v5=g4 v7=<pkl>
  a1 s1
  a2 s2
  s1 -1,a1
  s2 -1,a2
  a2 a2+a6
  a1 a1+a6
  s4 -1,a2
  s6 -1,a1
  s1 s1*rs3
  s2 s2*rs3
  a0 a7-1
  a7 a7-1
  s1 s1+fs4
  s2 s2+fs6
  -1,a2 s1
  -1,a1 s2
  jam ppg
pp s1 v7,a7
  s2 v2,a7
pph a5 a7-1
  s3 v5,a7
  s6 v7,a5
  a1 s1
  a2 s2
  s1 -1,a1
  s2 -1,a2
  a2 a2+a6
  a1 a1+a6
  s7 v2,a5
  s4 -1,a2
  s5 -1,a1
  a3 s6
  s1 s1*rs3
  a4 s7
  s2 s2*rs3
  s6 -1,a3
  s7 -1,a4
  s3 v5,a5
  s1 s1+fs4
  s2 s2+fs5
  a4 a4+a6
  a3 a3+a6
  s6 s6*rs3
  s7 s7*rs3
  -1,a2 s1
  -1,a1 s2
  s4 -1,a4
  s5 -1,a3
  a0 a5-1
  a7 a5-1
  s6 s6+fs4
  s7 s7+fs5
  -1,a4 s6
  -1,a3 s7
  jap pp
* v1=ikykj v3=j v4=ikyi v6=i
ppg a0 l
  a1 24
* clock 1
  v0 ,a0,1
  v5 v4+v0
  v2 v3<a1
  v7 v2*fv2
  a7 b71
  a3 b72
* clock 2
  ,a0,1 v5
  v1 v7-v3
  v2 v1>1
  a0 k
  a7 a7+a3
* clock 3
  v5 ,a0,1
  v1 v4+v5
  v7 v0<a1
  a7 a7+a3
* clock 4
  ,a0,1 v1
  v4 v6-v5
  vm v4,z
  a0 b71
  a2 2
* clock 5
  v6 v7*fv7
  v1 v6-v0
  v4 v1>1
  s7 vm
  b71 a7
* clock 6
  v7 ,a0,a2
  v1 v4+v3
* clock 7
  v6 v3-v0
  vm v6,z
  a7 a3-1
* clock 8
  v4 v2+v0
  s6 vm
  vm v6,m
  s7 s7!s6
* clock 9
  v6 v3-v5
  v0 v1!v4&vm
  a0 ikyk
* clock 10
  v3 v7+fv7
  vm v6,m
  v1 v2+v5
  v4 ,a0,1
  a0 k
* clock 11
  v5 ,a0,1
  v2 v4!v1&vm
  a0 l
* clock 12
  vm v6,z
  v4 ,a0,1
* clock 13
  v1 v3!v7&vm
* clock 14
  vm s7
  s6 t70
* clock 15
  v6 v3!v7&vm
* in the above ik stored in k il stored in l
* v1=gil v2=jk v4=il v5=ik v6=gik v0=jl
* clock 16
  v7 s6+v5
* clock 17
  v3 s6+v4
* clock 18
  v5 s6+v2
* clock 19
  v4 s6+v0
* v1=gil v3=<pil> v4=<pjl> v5=<pjk> v6=gik v7=<pik>
qq s1 v7,a7
  s2 v3,a7
  s3 v5,a7
  s4 v4,a7
  s5 v6,a7
  s6 v1,a7
  a1 s1
  a2 s2
  s1 -1,a1
  s2 -1,a2
  a3 s3
  a4 s4
  s3 -1,a3
  s4 -1,a4
  a4 a4+a6
  a3 a3+a6
  a2 a2+a6
  s7 -1,a4
  s1 s1*rs5
  s2 s2*rs6
  s3 s3*rs6
  s6 -1,a3
  a1 a1+a6
  s4 s4*rs5
  s5 -1,a2
  s1 s7-fs1
  s7 -1,a1
  a0 a7-1
  a7 a7-1
  s2 s6-fs2
  s3 s5-fs3
  s4 s7-fs4
  -1,a4 s1
  -1,a3 s2
  -1,a2 s3
  -1,a1 s4
  jap qq
  a1 b73
  a3 b74
  a0 a1-1
  a1 a1-1
  s0 0
  jap aaaa
  exit
  end
**    proc2 c90 version (hans nelemans (cray), march 1994)
          ident proc2
* entry proc2
* call proc2(fock,dens,k,q)
  common blkin
gin bss 340
ijin bss 170
mword bss 1
  common vectem
ikyk bss 64
ijkl bss 64
l bss 64
k bss 64
  block *
  block *
count set 3
bsave bss count
proc2 enter np=4
  a1 count
  a0 bsave
  a2 gin
  ,a0 b1,a1
  s1 mword,0
  argadd a3,4
  argadd a4,3
  argadd a6,1
  argadd a5,2
  b71 a2
  s6 s1
  s1 s1>1
  a2 ijin
  s3 1
  s7 s6-s1
  a4 a4-a3
  a3 a3-a6
  a6 a6-a5
  s4 s7-s3
  s5 a5
  b70 a2
  t70 s5
  s5 s7\s1
  s4 s4>6
  t77 s5
  b76 a3
  a1 s4
  s4 s4<6
  a3 64
  b77 a4
  s2 s7-s4
  a0 s2
  jan nfull
  s2 64
nfull a0 a1
  b1 a6
  b74 a3
  vl a3
  a2 s2
  s0 0
  b75 a2
aaaa jan www
  a3 b75
  s0 t77
  vl a3
www b73 a1
  b72 a3
* vl set
* a3=vl
* s0=svl-lvl
* even cases
  a0 b70
  a4 8
* clock 1
  v7 ,a0,1
  s5 <8
  v6 v7>a4
  v5 s5&v6
  a1 24
  a0 l
* clock 2
  v4 s5&v7
  v3 v5<a1
  v2 v3*fv3
* clock 3
  ,a0,1 v4
  v1 v2-v5
  v0 v1>1
  a0 k
* clock 4
  ,a0,1 v5
  v7 v6>a4
  v3 s5&v7
  v2 v4+v0
  a7 b71
  a2 2
  a0 ikyk
* clock 5
  v6 v7>a4
  v1 s5&v6
  ,a0,1 v0
  a0 a7+1
* clock 6
  v7 v1<a1
  v5 v7*fv7
  v0 ,a0,a2
* clock 7
  v4 v5-v1
  v7 v4>1
  s7 4.
  a7 a3-1
  a0 a3-1
  s6 t70
* clock 8
  v4 s7*rv0
  v5 v7+v3
  v0 v6>a4
* v0=ijkl v1=i v2=kl v3=j
* v4=g4 v5=ij v6=scratch v7=ikyi
* k,l,ikyk stored
  jsz pppf
  jaz ttt
  a3 a7
  a7 a7-1
* s0=svl-lvl
* a0=vl-1  a6=<f-p> a7=vl-1
* t70=<p>
* t77=lvl-svl(for last time round)
* b70=<ijklijkl> b71=<g> b72=cvl
* b73=nvord-1 b74=64 b75=lrem
* b76=<q-f> b77=<k-q>
* b1=<f-p>
* coulomb part(even cases)
pppf s0 a3
  b3 a3
* clock 9
  v6 s6+v2
  s0 s0<63
* clock 10
  v2 s6+v5
* v2=<pij> v4=g4 v6=<pkl>
  s1 v6,a7
  s2 v2,a7
  jsp ppph
  s3 v4,a7
  a1 s1
  a2 s2
  s1 -1,a1
  s2 -1,a2
  a2 a2+a6
  a1 a1+a6
  s4 -1,a2
  s6 -1,a1
  s1 s1*rs3
  s2 s2*rs3
  a0 a7-1
  a7 a7-1
  s1 s1+fs4
  s2 s2+fs6
  -1,a2 s1
  -1,a1 s2
  jam pppg
* above handles loop folded cases (odd no. of terms)
ppp s1 v6,a7
  s2 v2,a7
ppph a5 a7-1
  s3 v4,a7
  s6 v6,a5
  a1 s1
  a2 s2
  s1 -1,a1
  s2 -1,a2
  a2 a2+a6
  a1 a1+a6
  s7 v2,a5
  s4 -1,a2
  s5 -1,a1
  a3 s6
  s1 s1*rs3
  a4 s7
  s2 s2*rs3
  s6 -1,a3
  s7 -1,a4
  s3 v4,a5
  s1 s1+fs4
  s2 s2+fs5
  a4 a4+a6
  a3 a3+a6
  s6 s6*rs3
  s7 s7*rs3
  -1,a2 s1
  -1,a1 s2
  s4 -1,a4
  s5 -1,a3
  a0 a5-1
  a7 a5-1
  s6 s6+fs4
  s7 s7+fs5
  -1,a4 s6
  -1,a3 s7
  jap ppp
* above ppp loop costs about 70 clocks per two elements
*
* exchange part (even cases)
* code based on sequence
* kik=gik*rjl + kik
* kil=gil*rjk + kil
* kjk=gil*ril + kjk
* kjl=gik*rik + kjl
* clock 1
pppg a0 ijkl
  a1 24
  ,a0,1 v0
  v4 v3<a1
  v5 v4*fv4
  a0 k
* clock 2
  v2 ,a0,1
  v0 v1\v2
  v6 v5-v3
  v4 v6>1
  a0 l
* clock 3
  vm v0,z
  v6 ,a0,1
  v1 v3-v6
* clock 4
  s1 vm
  vm v1,z
  v5 v7+v2
  v0 v6<a1
* clock 5
  s2 vm
  vm v1,m
  ,a0,1 v5
  s1 s1!s2
  a0 k
* clock 6
  v1 v7+v6
  v5 v0*fv0
* clock 7
  v7 v4+v6
  ,a0,1 v1
  a0 ikyk
  a2 2
* clock 8
  v0 v5-v6
  v1 v0>1
  a7 b71
* clock 9
  v5 ,a0,1
  v6 v1+v3
  v0 v6!v7&vm
  a0 a7+1
* clock 10
  v7 ,a0,a2
  v6 v4+v2
  a0 l
* clock 11
  v1 v3-v2
  vm v1,m
* clock 12
  v2 v5+v3
  v4 v2!v6&vm
* clock 13
  v5 ,a0,1
  vm v1,z
  v2 v7+fv7
  a7 b3
  a0 k
* clock 14
  v6 ,a0,1
  v3 v2!v7&vm
  a7 a7-1
  vm s1
  s6 t70
* clock 15
  v1 v2!v7&vm
* in the above il stored in k ik stored in l
* clock 16
  v2 s6+v5
* clock 17
  v7 s6+v6
* clock 18
  v5 s6+v4
* clock 19
  v6 s6+v0
* v1=gik v2=<pik> v3=gil v5=<pjk> v6=<pjl> v7=<pil>
qqq s1 v2,a7
  s2 v7,a7
  s3 v5,a7
  s4 v6,a7
  s5 v1,a7
  s6 v3,a7
  a1 s1
  a2 s2
  s1 -1,a1
  s2 -1,a2
  a3 s3
  a4 s4
  s3 -1,a3
  s4 -1,a4
  a4 a4+a6
  a3 a3+a6
  a2 a2+a6
  s7 -1,a4
  s1 s1*rs5
  s2 s2*rs6
  s3 s3*rs6
  t72 s6
  s6 -1,a3
  a1 a1+a6
  s4 s4*rs5
  t71 s5
  s5 -1,a2
  a6 b76
  s1 s7-fs1
  s7 -1,a1
  a5 a1+a6
  b2 a7
  a7 a2+a6
  s2 s6-fs2
  s6 -1,a5
  s3 s5-fs3
  s5 -1,a7
  s4 s7-fs4
  s7 t71
  -1,a4 s1
  -1,a3 s2
  a3 a3+a6
  a6 a6+a4
  s1 -1,a3
  s2 -1,a6
  -1,a2 s3
  s3 t72
  a4 b77
  s6 s7*rs6
  a2 a6+a4
  -1,a1 s4
  s4 -1,a2
  s5 s5*rs3
  s1 s1*rs3
  a3 a3+a4
  a1 a7+a4
  a5 a5+a4
  s3 -1,a3
  s2 s2*rs7
  s7 -1,a1
  s4 s4+fs6
  s6 -1,a5
  a7 b2
  s3 s3+fs5
  a0 a7-1
  a7 a7-1
  -1,a2 s4
  s7 s7+fs1
  s6 s6+fs2
  a6 b1
  -1,a3 s3
  -1,a1 s7
  -1,a5 s6
  jap qqq
* above costs about 105 clocks
* odd cases
  a0 ijkl
  a4 8
* clock 1
  v0 ,a0,1
  s5 <8
ttt v6 v0>a4
  v5 s5&v6
  a7 b70
  a3 b72
  a1 24
  a0 l
  a7 a7+a3
* clock 2
  v4 s5&v0
  v3 v5<a1
  v2 v3*fv3
  b70 a7
* clock 3
  ,a0,1 v4
  v1 v2-v5
  v7 v1>1
  a0 k
* clock 4
  ,a0,1 v5
  v0 v6>a4
  v3 s5&v0
  v2 v4+v7
  a0 b71
  a2 2
  s7 4.
* clock 5
  v1 ,a0,a2
  v5 s7*rv1
  v6 v0>a4
* clock 6
  v1 v7+v3
  v4 v6<a1
  v0 v4*fv4
  a7 a3-1
  a0 ikyk
* clock 7
  ,a0,1 v1
  v7 v0-v6
  v4 v7>1
  s6 t70
* clock 8
  v0 v4+v3
  s0 a3
* l,k stored
* ikykj stored in ikyk
* v1=ikykj v2=kl v3=j v4=ikyi v5=g4
* v6=i v0=ij
* coulomb part(odd cases)
* clock 9
  v7 s6+v2
  s0 s0<63
* clock 10
  v2 s6+v0
  s1 v7,a7
  s2 v2,a7
  jsp pph
  s3 v5,a7
* v2=<pij> v5=g4 v7=<pkl>
  a1 s1
  a2 s2
  s1 -1,a1
  s2 -1,a2
  a2 a2+a6
  a1 a1+a6
  s4 -1,a2
  s6 -1,a1
  s1 s1*rs3
  s2 s2*rs3
  a0 a7-1
  a7 a7-1
  s1 s1+fs4
  s2 s2+fs6
  -1,a2 s1
  -1,a1 s2
  jam ppg
pp s1 v7,a7
  s2 v2,a7
pph a5 a7-1
  s3 v5,a7
  s6 v7,a5
  a1 s1
  a2 s2
  s1 -1,a1
  s2 -1,a2
  a2 a2+a6
  a1 a1+a6
  s7 v2,a5
  s4 -1,a2
  s5 -1,a1
  a3 s6
  s1 s1*rs3
  a4 s7
  s2 s2*rs3
  s6 -1,a3
  s7 -1,a4
  s3 v5,a5
  s1 s1+fs4
  s2 s2+fs5
  a4 a4+a6
  a3 a3+a6
  s6 s6*rs3
  s7 s7*rs3
  -1,a2 s1
  -1,a1 s2
  s4 -1,a4
  s5 -1,a3
  a0 a5-1
  a7 a5-1
  s6 s6+fs4
  s7 s7+fs5
  -1,a4 s6
  -1,a3 s7
  jap pp
* v1=ikykj v3=j v4=ikyi v6=i
ppg a0 l
  a1 24
* clock 1
  v0 ,a0,1
  v5 v4+v0
  v2 v3<a1
  v7 v2*fv2
  a7 b71
  a3 b72
* clock 2
  ,a0,1 v5
  v1 v7-v3
  v2 v1>1
  a0 k
  a7 a7+a3
* clock 3
  v5 ,a0,1
  v1 v4+v5
  v7 v0<a1
  a7 a7+a3
* clock 4
  ,a0,1 v1
  v4 v6-v5
  vm v4,z
  a0 b71
  a2 2
* clock 5
  v6 v7*fv7
  v1 v6-v0
  v4 v1>1
  s7 vm
  b71 a7
* clock 6
  v7 ,a0,a2
  v1 v4+v3
* clock 7
  v6 v3-v0
  vm v6,z
  a7 a3-1
* clock 8
  v4 v2+v0
  s6 vm
  vm v6,m
  s7 s7!s6
* clock 9
  v6 v3-v5
  v0 v1!v4&vm
  a0 ikyk
* clock 10
  v3 v7+fv7
  vm v6,m
  v1 v2+v5
  v4 ,a0,1
  a0 k
* clock 11
  v5 ,a0,1
  v2 v4!v1&vm
  a0 l
* clock 12
  vm v6,z
  v4 ,a0,1
* clock 13
  v1 v3!v7&vm
* clock 14
  vm s7
  s6 t70
* clock 15
  v6 v3!v7&vm
* in the above ik stored in k il stored in l
* v1=gil v2=jk v4=il v5=ik v6=gik v0=jl
* clock 16
  v7 s6+v5
* clock 17
  v3 s6+v4
* clock 18
  v5 s6+v2
* clock 19
  v4 s6+v0
* v1=gil v3=<pil> v4=<pjl> v5=<pjk> v6=gik v7=<pik>
qq s1 v7,a7
  s2 v3,a7
  s3 v5,a7
  s4 v4,a7
  s5 v6,a7
  s6 v1,a7
  a1 s1
  a2 s2
  s1 -1,a1
  s2 -1,a2
  a3 s3
  a4 s4
  s3 -1,a3
  s4 -1,a4
  a4 a4+a6
  a3 a3+a6
  a2 a2+a6
  s7 -1,a4
  s1 s1*rs5
  s2 s2*rs6
  s3 s3*rs6
  t72 s6
  s6 -1,a3
  a1 a1+a6
  s4 s4*rs5
  t71 s5
  s5 -1,a2
  a6 b76
  s1 s7-fs1
  s7 -1,a1
  a5 a1+a6
  b2 a7
  a7 a2+a6
  s2 s6-fs2
  s6 -1,a5
  s3 s5-fs3
  s5 -1,a7
  s4 s7-fs4
  s7 t71
  -1,a4 s1
  -1,a3 s2
  a3 a3+a6
  a6 a6+a4
  s1 -1,a3
  s2 -1,a6
  -1,a2 s3
  s3 t72
  a4 b77
  s6 s7*rs6
  a2 a6+a4
  -1,a1 s4
  s4 -1,a2
  s5 s5*rs3
  s1 s1*rs3
  a3 a3+a4
  a1 a7+a4
  a5 a5+a4
  s3 -1,a3
  s2 s2*rs7
  s7 -1,a1
  s4 s4+fs6
  s6 -1,a5
  a7 b2
  s3 s3+fs5
  a0 a7-1
  a7 a7-1
  -1,a2 s4
  s7 s7+fs1
  s6 s6+fs2
  a6 b1
  -1,a3 s3
  -1,a1 s7
  -1,a5 s6
  jap qq
  a1 b73
  a3 b74
  a0 a1-1
  a1 a1-1
  s0 0
  jap aaaa
  a0 bsave
  a1 count
  b1,a1 ,a0
  exit
  end
******* 8/3/94  van vic *******
* date: fri, 4 mar 1994 17:09:36 gmt
* from: "v.r.saunders" <v.r.saunders@dl.ac.uk>
          ident mxmb
*  coded by v.r.saunders february 1994 for cray-c90
*  email v.r.saunders@dl.ac.uk   (internet)
*  code is re-entrant - may be multitasked
*  matrix multiply routine
*  r(ncol,nrow) = a(ncol,nlink)*b(nlink,nrow) + r(ncol,nrow)
*  the r-matrix   ***must***  be initialized
*  sparcity of the b-matrix is exploited
*  column and row strides are
*  mcola,mrowa   for the a-matrix
*  mcolb,mrowb   for the b-matrix
*  mcolr,mrowr   for the r-matrix
*
*  fortran equivalent (logic as close as possible)
*
*      subroutine mxmb
*     *(a,mcola,mrowa,
*     * b,mcolb,mrowb,
*     * r,mcolr,mrowr,
*     * ncol,nlink,nrow)
*      dimension a(*),b(*),r(*)
*      jr=1
*      ja=1
*11    ncoll=mod(ncol,128)
*      jb=1
*      ir=jr
*      do 1 i=1,nrow
*c fetch r-vector with vl=ncoll, base=ir, stride=mcolr
*      mlink=nlink
*      ib=ja
*      ia=ja
*22    link=mod(mlink,128)
*      mlink=mlink-link
*c fetch b-vector with vl=link, base=ib, stride=mcolb
*      do 2 j=1,link
*      fac=b(ib)
*      ib=ib+mcolb
*      if(fac)222,2,222
*222   ka=ia
*      kr=ir
*      do 3 k=1,ncoll
*      r(kr)=a(ka)*fac+r(kr)
*      ka=ka+mcola
*3     kr=kr+mcolr
*2     ia=ia+mrowa
*      if(mlink.ne.0)goto 22
*c store r-vector with vl=ncoll, base=ir, stride=mcolr
*      jb=jb+mrowb
*1     ir=ir+mrowr
*      ncol=ncol-ncoll
*      jr=mcolr*ncoll+jr
*      ja=mcola*ncoll+ja
*      if(ncol.ne.0)goto 11
*      return
*      end
*
     align
     entry mxmb
*  get arguments
mxmb a7 10,a6     a7=<ncol>
     a4 11,a6     a4=<nlink>
     a5 12,a6     a5=<nrow>
     a1 9,a6      a1=<mrowr)
     a2 8,a6      a2=<mcolr>
     a3 5,a6      a3=<mcolb>
     s4 4,a6      s4=<b>
     s3 7,a6      s3=jr
     s2 1,a6      s2=ja
*
     a7 0,a7      a7=ncol
     a4 0,a4      a4=nlink
     s6 0,a5      s6=nrow
     a5 6,a6      a5=<mrowb>
     s1 0,a1      s1=mrowr
     a1 2,a6      a1=<mcola>
     a6 3,a6      a6=<mrowa>
     a2 0,a2      a2=mcolr
     a3 0,a3      a3=mcolb
     t71 s4       t71=<b>
     s7 1         s7=1
*
     a0 a4*a7     a0=ncol*nlink
     s0 s6        s0=nrow
     b74 a4       b74=nlink
     t70 s6       t70=nrow
     t73 s1       t73=mrowr
     b73 a2       b73=mcolr
*
     jaz exit     exit if ncol or nlink = 0
     jsz exit     exit if nrow = 0
*
     s5 0,a5      s5=mrowb
     a1 0,a1      a1=mcola
     a2 0,a6      a2=mrowa
*  top section of loop p11
p11  vl a7        vl=ncoll
     t72 s3       t72=jr
     a6 vl        a6=ncoll
     b72 a7       b72=ncol
     b71 a6       b71=ncoll
*  top section of loop p1
p1   a0 s3        a0=ir
     a6 b73       a6=mcolr
     a7 b74       a7=mlink
     v1 ,a0,a6    v1=r(1)
     a5 s4        a5=ib
     a4 s2        a4=ia
*  start of loops p22 and p2
*
*  a1   mcola
*  a2   mrowa
*  a3   mcolb
*  a4   ia
*  a5   ib
*  a6   link-1
*  a7   j
*
*  s1   fac
*  s2   ja
*  s3   ir
*  s4   jb
*  s5   mrowb
*  s6   i
*  s7   1
*
*  b70  mlink
*  b71  ncoll
*  b72  ncol (via a7)
*  b73  mcolr
*  b74  nlink
*
*  t70  nrow (via s6)
*  t71  <b>  (via s4)
*  t72  jr   (via s3)
*  t73  mrowr
*
     vl a7        vl=link
     a0 a5        a0=ib
     a6 vl        a6=link
     v0 ,a0,a3    v0=b-vector
     a7 a7-a6     a7=mlink-link=mlink
     a6 a6-1      a6=link-1
     b70 a7       b70=mlink
     a7 b71       a7=ncoll
     vl a7        vl=ncoll
     a7 0         a7=j
*  look for first finite element of b-vector
af1  a0 a6-a7     a0=link-j
     s0 v0,a7     s0=fac
     s1 v0,a7     s1=fac
     a4 a4+a2     a4=ia+mrowa=ianew
     a5 a5+a3     a5=ib+mcolb=ib
     a7 a7+1      a7=j+1=j
     jaz af2
     a0 a4-a2     a0=ianew-mrowa=ia
     jsz af1
     v5 ,a0,a1    v5=a-vector
     v4 s1*rv5    v4=fac*a=r(2)
*  update r(1)
bf1  a0 a6-a7     a0=link-j
     s0 v0,a7     s0=fac
     s1 v0,a7     s1=fac
     a4 a4+a2     a4=ia+mrowa=ianew
     a5 a5+a3     a5=ib+mcolb=ib
     a7 a7+1      a7=j+1=j
     jaz bf2
     a0 a4-a2     a0=ianew-mrowa=ia
     jsz bf1
     v3 ,a0,a1    v3=a-vector
     v2 s1*rv3    v2=fac*a
     v1 v1+fv2    v1=r(1)+fac*a=r(1)
*  update r(2)
cf1  a0 a6-a7     a0=link-j
     s0 v0,a7     s0=fac
     s1 v0,a7     s1=fac
     a4 a4+a2     a4=ia+mrowa=ianew
     a5 a5+a3     a5=ib+mcolb=ib
     a7 a7+1      a7=j+1=j
     jaz cf2
     a0 a4-a2     a0=ianew-mrowa=ia
     jsz cf1
     v6 ,a0,a1    v6=a-vector
     v5 s1*rv6    v5=fac*a
     v4 v4+fv5    v4=r(2)+fac*a=r(2)
     j bf1
af2  a7 b70       a7=mlink
     jsz af3
     a0 a4-a2     a0=ianew-mrowa=ja
     v5 ,a0,a1    v5=a-vector
     a0 b70       a0=mlink
     v4 s1*rv5    v4=fac*a=r(1)
     jaz stor
     vl a7        vl=link
     a0 a5        a0=ib
     a6 vl        a6=link
     v0 ,a0,a3    v0=b-vector
     a7 a7-a6     a7=mlink-link=mlink
     a6 a6-1      a6=link-1
     b70 a7       b70=mlink
     a7 b71       a7=ncoll
     vl a7        vl=ncoll
     a7 0         a7=j
     j bf1
af3  a0 b70       a0=mlink
     jaz nsto
     vl a7        vl=link
     a0 a5        a0=ib
     a6 vl        a6=link
     v0 ,a0,a3    v0=b-vector
     a7 a7-a6     a7=mlink-link=mlink
     a6 a6-1      a6=link-1
     b70 a7       b70=mlink
     a7 b71       a7=ncoll
     vl a7        vl=ncoll
     a7 0         a7=j
     j af1
bf2  a7 b70       a7=mlink
     jsz bf3
     a0 a4-a2     a0=ianew-mrowa=ja
     v3 ,a0,a1    v3=a-vector
     a0 b70       a0=mlink
     v2 s1*rv3    v2=fac*a
     v1 v1+fv2    v1=r(1)+fac*a=r(1)
     jaz stor
     vl a7        vl=link
     a0 a5        a0=ib
     a6 vl        a6=link
     v0 ,a0,a3    v0=b-vector
     a7 a7-a6     a7=mlink-link=mlink
     a6 a6-1      a6=link-1
     b70 a7       b70=mlink
     a7 b71       a7=ncoll
     vl a7        vl=ncoll
     a7 0         a7=j
     j cf1
bf3  a0 b70       a0=mlink
     jaz stor
     vl a7        vl=link
     a0 a5        a0=ib
     a6 vl        a6=link
     v0 ,a0,a3    v0=b-vector
     a7 a7-a6     a7=mlink-link=mlink
     a6 a6-1      a6=link-1
     b70 a7       b70=mlink
     a7 b71       a7=ncoll
     vl a7        vl=ncoll
     a7 0         a7=j
     j bf1
cf2  a7 b70       a7=mlink
     jsz cf3
     a0 a4-a2     a0=ianew-mrowa=ja
     v6 ,a0,a1    v6=a-vector
     a0 b70       a0=mlink
     v5 s1*rv6    v5=fac*a
     v4 v4+fv5    v4=r(2)+fac*a=r(2)
     jaz stor
     vl a7        vl=link
     a0 a5        a0=ib
     a6 vl        a6=link
     v0 ,a0,a3    v0=b-vector
     a7 a7-a6     a7=mlink-link=mlink
     a6 a6-1      a6=link-1
     b70 a7       b70=mlink
     a7 b71       a7=ncoll
     vl a7        vl=ncoll
     a7 0         a7=j
     j bf1
cf3  a0 b70       a0=mlink
     jaz stor
     vl a7        vl=link
     a0 a5        a0=ib
     a6 vl        a6=link
     v0 ,a0,a3    v0=b-vector
     a7 a7-a6     a7=mlink-link=mlink
     a6 a6-1      a6=link-1
     b70 a7       b70=mlink
     a7 b71       a7=ncoll
     vl a7        vl=ncoll
     a7 0         a7=j
     j cf1
*  store r-vector
stor v7 v1+fv4    r=r(1)+r(2)
     a0 s3        a0=ir
     a6 b73       a6=mcolr
     ,a0,a6 v7    store  r-vector
*  bottom section of loop p1
nsto s0 s6-s7     s0=i-1=i
     s1 t73       s1=mrowr 
     s4 s4+s5     s4=jb+mrowb=jb
     s6 s6-s7     s6=i-1=i
     s3 s3+s1     s3=ir+mrowr=ir
     jsn p1
*  bottom section of loop p11
     a6 b71       a6=ncoll
     a7 b72       a7=ncol
     a4 b73       a4=mcolr
     a0 a7-a6     a0=ncol-ncoll=ncol
     a5 a1*a6     a5=mcola*ncoll
     a4 a4*a6     a4=mcolr*ncoll
     s3 t72       s3=jr
     s6 t70       s6=i
     a7 a7-a6     a7=ncol-ncoll=ncol
     s1 +a5       s1=mcola*ncoll
     s4 +a4       s4=mcolr*ncoll
     s2 s2+s1     s2=ja+mcola*ncoll=ja
     s3 s3+s4     s3=jr+mcolr*ncoll=jr
     s4 t71       s4=jb
     jan p11
*  exit from routine
exit j b00
     end
          ident mxmbn
*  coded by v.r.saunders february 1994 for cray-c90
*  email v.r.saunders@dl.ac.uk   (internet)
*  code is re-entrant - may be multitasked
*  matrix multiply routine
*  r(ncol,nrow) = r(ncol,nrow) - a(ncol,nlink)*b(nlink,nrow)
*  the r-matrix   ***must***  be initialized
*  sparcity of the b-matrix is exploited
*  column and row strides are
*  mcola,mrowa   for the a-matrix
*  mcolb,mrowb   for the b-matrix
*  mcolr,mrowr   for the r-matrix
*
*  fortran equivalent (logic as close as possible)
*
*      subroutine mxmbn
*     *(a,mcola,mrowa,
*     * b,mcolb,mrowb,
*     * r,mcolr,mrowr,
*     * ncol,nlink,nrow)
*      dimension a(*),b(*),r(*)
*      jr=1
*      ja=1
*11    ncoll=mod(ncol,128)
*      jb=1
*      ir=jr
*      do 1 i=1,nrow
*c fetch r-vector with vl=ncoll, base=ir, stride=mcolr
*      mlink=nlink
*      ib=ja
*      ia=ja
*22    link=mod(mlink,128)
*      mlink=mlink-link
*c fetch b-vector with vl=link, base=ib, stride=mcolb
*      do 2 j=1,link
*      fac=b(ib)
*      ib=ib+mcolb
*      if(fac)222,2,222
*222   ka=ia
*      kr=ir
*      do 3 k=1,ncoll
*      r(kr)=r(kr)-a(ka)*fac
*      ka=ka+mcola
*3     kr=kr+mcolr
*2     ia=ia+mrowa
*      if(mlink.ne.0)goto 22
*c store r-vector with vl=ncoll, base=ir, stride=mcolr
*      jb=jb+mrowb
*1     ir=ir+mrowr
*      ncol=ncol-ncoll
*      jr=mcolr*ncoll+jr
*      ja=mcola*ncoll+ja
*      if(ncol.ne.0)goto 11
*      return
*      end
*
     align
     entry mxmbn
*  get arguments
mxmbn a7 10,a6     a7=<ncol>
     a4 11,a6     a4=<nlink>
     a5 12,a6     a5=<nrow>
     a1 9,a6      a1=<mrowr)
     a2 8,a6      a2=<mcolr>
     a3 5,a6      a3=<mcolb>
     s4 4,a6      s4=<b>
     s3 7,a6      s3=jr
     s2 1,a6      s2=ja
*
     a7 0,a7      a7=ncol
     a4 0,a4      a4=nlink
     s6 0,a5      s6=nrow
     a5 6,a6      a5=<mrowb>
     s1 0,a1      s1=mrowr
     a1 2,a6      a1=<mcola>
     a6 3,a6      a6=<mrowa>
     a2 0,a2      a2=mcolr
     a3 0,a3      a3=mcolb
     t71 s4       t71=<b>
     s7 1         s7=1
*
     a0 a4*a7     a0=ncol*nlink
     s0 s6        s0=nrow
     b74 a4       b74=nlink
     t70 s6       t70=nrow
     t73 s1       t73=mrowr
     b73 a2       b73=mcolr
*
     jaz exit     exit if ncol or nlink = 0
     jsz exit     exit if nrow = 0
*
     s5 0,a5      s5=mrowb
     a1 0,a1      a1=mcola
     a2 0,a6      a2=mrowa
*  top section of loop p11
p11  vl a7        vl=ncoll
     t72 s3       t72=jr
     a6 vl        a6=ncoll
     b72 a7       b72=ncol
     b71 a6       b71=ncoll
*  top section of loop p1
p1   a0 s3        a0=ir
     a6 b73       a6=mcolr
     a7 b74       a7=mlink
     v1 ,a0,a6    v1=r(1)
     a5 s4        a5=ib
     a4 s2        a4=ia
*  start of loops p22 and p2
*
*  a1   mcola
*  a2   mrowa
*  a3   mcolb
*  a4   ia
*  a5   ib
*  a6   link-1
*  a7   j
*
*  s1   fac
*  s2   ja
*  s3   ir
*  s4   jb
*  s5   mrowb
*  s6   i
*  s7   1
*
*  b70  mlink
*  b71  ncoll
*  b72  ncol (via a7)
*  b73  mcolr
*  b74  nlink
*
*  t70  nrow (via s6)
*  t71  <b>  (via s4)
*  t72  jr   (via s3)
*  t73  mrowr
*
     vl a7        vl=link
     a0 a5        a0=ib
     a6 vl        a6=link
     v0 ,a0,a3    v0=b-vector
     a7 a7-a6     a7=mlink-link=mlink
     a6 a6-1      a6=link-1
     b70 a7       b70=mlink
     a7 b71       a7=ncoll
     vl a7        vl=ncoll
     a7 0         a7=j
*  look for first finite element of b-vector
af1  a0 a6-a7     a0=link-j
     s0 v0,a7     s0=fac
     s1 v0,a7     s1=fac
     a4 a4+a2     a4=ia+mrowa=ianew
     a5 a5+a3     a5=ib+mcolb=ib
     a7 a7+1      a7=j+1=j
     jaz af2
     a0 a4-a2     a0=ianew-mrowa=ia
     jsz af1
     v5 ,a0,a1    v5=a-vector
     v4 s1*rv5    v4=fac*a=r(2)
*  update r(1)
bf1  a0 a6-a7     a0=link-j
     s0 v0,a7     s0=fac
     s1 v0,a7     s1=fac
     a4 a4+a2     a4=ia+mrowa=ianew
     a5 a5+a3     a5=ib+mcolb=ib
     a7 a7+1      a7=j+1=j
     jaz bf2
     a0 a4-a2     a0=ianew-mrowa=ia
     jsz bf1
     v3 ,a0,a1    v3=a-vector
     v2 s1*rv3    v2=fac*a
     v1 v1-fv2    v1=r(1)-fac*a=r(1)
*  update r(2)
cf1  a0 a6-a7     a0=link-j
     s0 v0,a7     s0=fac
     s1 v0,a7     s1=fac
     a4 a4+a2     a4=ia+mrowa=ianew
     a5 a5+a3     a5=ib+mcolb=ib
     a7 a7+1      a7=j+1=j
     jaz cf2
     a0 a4-a2     a0=ianew-mrowa=ia
     jsz cf1
     v6 ,a0,a1    v6=a-vector
     v5 s1*rv6    v5=fac*a
     v4 v4+fv5    v4=r(2)+fac*a=r(2)
     j bf1
af2  a7 b70       a7=mlink
     jsz af3
     a0 a4-a2     a0=ianew-mrowa=ja
     v5 ,a0,a1    v5=a-vector
     a0 b70       a0=mlink
     v4 s1*rv5    v4=fac*a=r(1)
     jaz stor
     vl a7        vl=link
     a0 a5        a0=ib
     a6 vl        a6=link
     v0 ,a0,a3    v0=b-vector
     a7 a7-a6     a7=mlink-link=mlink
     a6 a6-1      a6=link-1
     b70 a7       b70=mlink
     a7 b71       a7=ncoll
     vl a7        vl=ncoll
     a7 0         a7=j
     j bf1
af3  a0 b70       a0=mlink
     jaz nsto
     vl a7        vl=link
     a0 a5        a0=ib
     a6 vl        a6=link
     v0 ,a0,a3    v0=b-vector
     a7 a7-a6     a7=mlink-link=mlink
     a6 a6-1      a6=link-1
     b70 a7       b70=mlink
     a7 b71       a7=ncoll
     vl a7        vl=ncoll
     a7 0         a7=j
     j af1
bf2  a7 b70       a7=mlink
     jsz bf3
     a0 a4-a2     a0=ianew-mrowa=ja
     v3 ,a0,a1    v3=a-vector
     a0 b70       a0=mlink
     v2 s1*rv3    v2=fac*a
     v1 v1-fv2    v1=r(1)-fac*a=r(1)
     jaz stor
     vl a7        vl=link
     a0 a5        a0=ib
     a6 vl        a6=link
     v0 ,a0,a3    v0=b-vector
     a7 a7-a6     a7=mlink-link=mlink
     a6 a6-1      a6=link-1
     b70 a7       b70=mlink
     a7 b71       a7=ncoll
     vl a7        vl=ncoll
     a7 0         a7=j
     j cf1
bf3  a0 b70       a0=mlink
     jaz stor
     vl a7        vl=link
     a0 a5        a0=ib
     a6 vl        a6=link
     v0 ,a0,a3    v0=b-vector
     a7 a7-a6     a7=mlink-link=mlink
     a6 a6-1      a6=link-1
     b70 a7       b70=mlink
     a7 b71       a7=ncoll
     vl a7        vl=ncoll
     a7 0         a7=j
     j bf1
cf2  a7 b70       a7=mlink
     jsz cf3
     a0 a4-a2     a0=ianew-mrowa=ja
     v6 ,a0,a1    v6=a-vector
     a0 b70       a0=mlink
     v5 s1*rv6    v5=fac*a
     v4 v4+fv5    v4=r(2)+fac*a=r(2)
     jaz stor
     vl a7        vl=link
     a0 a5        a0=ib
     a6 vl        a6=link
     v0 ,a0,a3    v0=b-vector
     a7 a7-a6     a7=mlink-link=mlink
     a6 a6-1      a6=link-1
     b70 a7       b70=mlink
     a7 b71       a7=ncoll
     vl a7        vl=ncoll
     a7 0         a7=j
     j bf1
cf3  a0 b70       a0=mlink
     jaz stor
     vl a7        vl=link
     a0 a5        a0=ib
     a6 vl        a6=link
     v0 ,a0,a3    v0=b-vector
     a7 a7-a6     a7=mlink-link=mlink
     a6 a6-1      a6=link-1
     b70 a7       b70=mlink
     a7 b71       a7=ncoll
     vl a7        vl=ncoll
     a7 0         a7=j
     j cf1
*  store r-vector
stor v7 v1-fv4    r=r(1)-r(2)
     a0 s3        a0=ir
     a6 b73       a6=mcolr
     ,a0,a6 v7    store  r-vector
*  bottom section of loop p1
nsto s0 s6-s7     s0=i-1=i
     s1 t73       s1=mrowr 
     s4 s4+s5     s4=jb+mrowb=jb
     s6 s6-s7     s6=i-1=i
     s3 s3+s1     s3=ir+mrowr=ir
     jsn p1
*  bottom section of loop p11
     a6 b71       a6=ncoll
     a7 b72       a7=ncol
     a4 b73       a4=mcolr
     a0 a7-a6     a0=ncol-ncoll=ncol
     a5 a1*a6     a5=mcola*ncoll
     a4 a4*a6     a4=mcolr*ncoll
     s3 t72       s3=jr
     s6 t70       s6=i
     a7 a7-a6     a7=ncol-ncoll=ncol
     s1 +a5       s1=mcola*ncoll
     s4 +a4       s4=mcolr*ncoll
     s2 s2+s1     s2=ja+mcola*ncoll=ja
     s3 s3+s4     s3=jr+mcolr*ncoll=jr
     s4 t71       s4=jb
     jan p11
*  exit from routine
exit j b00
     end
*          ident locat1
**
**   l = locat1(list,nlist,pattern)
**
**  to find member of   list   =   pattern
**  or return zero if no match found
**
**
*locat1 enter np=3
*  argadd a2,2
*  argadd a4,3
*  argadd a1,1
*  a7 64
*  s1 0
*  a2 0,a2
*  s4 0,a4
*  a6 0
*  vl a7
*  a0 a7-a2
*  jap vvv
*  a0 a1
*  j xxx
*yyy s0 vm
*  a0 a1
*  jsn zzz
*xxx v7 ,a0,1
*  a2 a2-a7
*  a1 a1+a7
*  v6 s4-v7
*  a0 a7-a2
*  vm v6,z
*  a6 a6+a7
*  jam yyy
*  s0 vm
*  a0 a2
*  jsz www
*zzz s6 vm
*  a6 a6-a7
*aaa a7 zs6
*  a6 a6+1
*  a6 a6+a7
*  s1 a6
*out exit
*vvv a0 a2
*www jaz out
*  a0 a1
*  vl a2
*  v7 ,a0,1
*  v6 s4-v7
*  vm v6,z
*  s0 vm
*  s6 vm
*  jsn aaa
*  exit
*  end
**    dbuild c90 version (hans nelemans (cray), march 1994)
  ident dbuild
* sgmata vrs
* dbuild mods rjh 30/1/87
* differs from sgmata in that indices are packed two to a word
* in /craypk/ with offset of 12 bits (*4096). this eliminates
* need to handle odd and even cases. odd ones deleted.
* additional confusion from loop unrolling etc
*  entry dbuild
* call dbuild(fock,dens)
  common shlt
tol bss 1
cutoff bss 1
icount bss 1
ic4 bss 1
  common blkin
gin bss 340
  common craypk
ind bss 680
  common vectem
ikyk bss 64
ijkl bss 64
l bss 64
k bss 64
  block *
  block *
  block *
  block *
dbuild enter np=2
  s1 icount,0              icount into s1
  a1 gin                   <g> into a1
  s3 1                     s3 = 1
  argadd a5,2              a5 = <p>
  a2 ind                   a2 = <ind>
  argadd a6,1              a6 = <f>
  s1 s1-s3                 s1 = mword = no. of integrals = icount-1
  b71 a1                   b71 = <g>
  b70 a2                   b70 = <ind>
  s4 s1-s3                 s4 = mword-1
  s5 a5                    s5 = <p>
  a6 a6-a5                 a6 = <f>-<p>
  t70 s5                   t70 = <p>
  s4 s4>6                  s4 = no. of vector loads - 1
  a1 s4                    a1 =    ditto
  s4 s4<6
  a3 64                    a3 = 64
  s2 s1-s4                 s2 = final vector length
  a0 s2
  jan nfull                jif s2 not full vector
  s2 64                    s2 = 64 (full vector length)
nfull a0 a1                a0 = no. of vector loads - 1
  b74 a3                   b74 = 64
  vl a3                    vl = 64
  a2 s2                    a2 = final vector length
  b75 a2                   b75 =     ditto
aaaa jan www               if last vector set vl accordingly
  a3 b75
  vl a3
www b73 a1                 b73 = no. of vector loads - 1
  b72 a3                   b72 = current vl
* vl set
* a3=vl
* get k,l and start to form iky(k) etc
  a2 b70
  a0 a2+1                  a0 = address of first kl
  a2 2                     a2 = 2 = increment between kl
  a4 12                    a4 = 12 (4096 = 2**12)
* clock 1
  v7 ,a0,a2                v7 = kl
  s5 <12                   s5 = 12 bit mask
  v5 v7>a4                 v5 = k
  a1 24                    a1 = 24
  a0 l                     a0 = <l>
* clock 2
  v4 s5&v7                 v4 = l
  v3 v5<a1
*
* j dbug
*
  v2 v3*fv3                v2 = k*k
* clock 3
  ,a0,1 v4                 store l at <l>
  v1 v2-v5
  v0 v1>1                  v0 = iky(k)
  a0 k                     a0 = <k>
* clock 4
  ,a0,1 v5                 store k at <k>
* load ij into v7
  a0 b70
  v7 ,a0,a2                v7 = ij
  v3 s5&v7                 v3 = j
  v2 v4+v0                 v2 = iky(k)+l
  a0 ikyk                  a0 = <ikyk>
* clock 5
  v1 v7>a4                 v1 = i
  ,a0,1 v0                 store iky(k) at <ikyk>
  a0 b71                   a0 = <g>
* clock 6
  v6 v1<a1
  v5 v6*fv6                v5 = i*i
  v0 ,a0,1                 v0 = g
* clock 7
  v4 v5-v1
  v7 v4>1                  v7 = ikyk(i)
  s7 4.                    s7 = 4.0
  a7 a3-1                  a7 = vl-1
  a0 a3-1                  a0 = vl-1
  s6 t70                   s6 = <p>
* clock 8
  v4 s7*rv0                v4 = 4.0*g
  v5 v7+v3                 v5 = ikyk(i) + j
* v1=i v2=kl v3=j
* v4=g4 v5=ij v6=scratch v7=ikyi
* k,l,ikyk stored
* a0=vl-1  a6=<a-p> a7=vl-1
* t70=<p>
* b70=<ijklijkl> b71=<g> b72=cvl
* b73=nvord-1 b74=64 b75=lrem
* coulomb part(even cases)
pppf s0 a3                 a0 = vl
  b77 a3                   b77 = vl
* clock 9
  v6 s6+v2                 v6 = <p(kl)>
  s0 s0<63                 s0 has sign bit set if vl is odd
* clock 10
  v2 s6+v5                 v2 = <p(ij)>
* v2=<pij> v4=g4 v6=<pkl>
  s1 v6,a7                 s1 = <p(kl)>,vl
  s2 v2,a7                 s2 = <p(ij)>,vl
  jsp ppph                 jump if vl is even
  s3 v4,a7                 s3 = 4.*g,vl
  a1 s1                    a1 = <p(kl)>,vl
  a2 s2                    a2 = <p(ij)>,vl
  s1 -1,a1                 s1 = p(kl),vl
  s2 -1,a2                 s2 = p(ij),vl
  a2 a2+a6                 a2 = <f(ij)>,vl
  a1 a1+a6                 a1 = <f(kl)>,vl
  s4 -1,a2                 s4 = f(ij),vl
  s6 -1,a1                 s6 = f(kl),vl
  s1 s1*rs3                s1 = p(kl),vl * 4.*g,vl
  s2 s2*rs3                s2 = p(ij),vl * 4.*g,vl
  a0 a7-1                  a0 = vl-2
  a7 a7-1                  a7 =  ditto
  s1 s1+fs4                s1 = (f(ij)+p(kl)*4.g),vl
  s2 s2+fs6                s2 =    kl    ij
  -1,a2 s1                 store  f(ij),vl
  -1,a1 s2                 store    kl
  jam pppg                 jump if vl = 1
* above handles loop folded cases (odd no. of terms)
ppp s1 v6,a7
  s2 v2,a7
ppph a5 a7-1               a5 = counter -1 = offset for second element
  s3 v4,a7                 s3 = 4g,1
  s6 v6,a5
  a1 s1
  a2 s2
  s1 -1,a1                 s1 = p(kl),1
  s2 -1,a2                 s2 = p(ij),1
  a2 a2+a6
  a1 a1+a6
  s7 v2,a5
  s4 -1,a2                 s4 = f(ij),1
  s5 -1,a1                 s5 = f(kl),1
  a3 s6
  s1 s1*rs3
  a4 s7
  s2 s2*rs3
  s6 -1,a3                 s6 = p(kl),2
  s7 -1,a4                 s7 = p(ij),2
  s3 v4,a5
  s1 s1+fs4
  s2 s2+fs5
  a4 a4+a6
  a3 a3+a6
  s6 s6*rs3
  s7 s7*rs3
  -1,a2 s1                 store f(ij),1
  -1,a1 s2                 store f(kl),1
  s4 -1,a4                 s4 = f(ij),2
  s5 -1,a3                 s5 = f(kl),2
  a0 a5-1                  a0 = no. remaining-1
  a7 a5-1                  a7 = updated
  s6 s6+fs4
  s7 s7+fs5
  -1,a4 s6                 store f(ij),2
  -1,a3 s7                 store f(kl),2
  jap ppp                  jump if more to do.
* above ppp loop costs about 70 clocks per two elements
*
* exchange part (even cases)
* code based on sequence
* kik=gik*rjl + kik
* kil=gil*rjk + kil
* kjk=gil*ril + kjk
* kjl=gik*rik + kjl
* clock 1
pppg a1 24
  v4 v3<a1
  v5 v4*fv4
  a0 k
* clock 2
  v2 ,a0,1
  v0 v1\v2
  v6 v5-v3
  v4 v6>1
  a0 l
* clock 3
  vm v0,z
  v6 ,a0,1
  v1 v3-v6
* clock 4
  s1 vm
  vm v1,z
  v5 v7+v2
  v0 v6<a1
* clock 5
  s2 vm
  vm v1,m
  ,a0,1 v5
  s1 s1!s2
  a0 k
* clock 6
  v1 v7+v6
  v5 v0*fv0
* clock 7
  v7 v4+v6
  ,a0,1 v1
  a0 ikyk
  a2 2
* clock 8
  v0 v5-v6
  v1 v0>1
* clock 9
  v5 ,a0,1
  v6 v1+v3
  v0 v6!v7&vm
  a0 b71
* clock 10
  v7 ,a0,1
  v6 v4+v2
  a0 l
* clock 11
  v1 v3-v2
  vm v1,m
* clock 12
  v2 v5+v3
  v4 v2!v6&vm
* clock 13
  v5 ,a0,1
  vm v1,z
  v2 v7+fv7
  a7 b77
  a0 k
* clock 14
  v6 ,a0,1
  v3 v2!v7&vm
  a7 a7-1
  vm s1
  s6 t70
* clock 15
  v1 v2!v7&vm
* in the above il stored in k ik stored in l
* clock 16
  v2 s6+v5
* clock 17
  v7 s6+v6
* clock 18
  v5 s6+v4
* clock 19
  v6 s6+v0
* v1=gik v2=<pik> v3=gil v5=<pjk> v6=<pjl> v7=<pil>
qqq s1 v2,a7
  s2 v7,a7
  s3 v5,a7
  s4 v6,a7
  s5 v1,a7
  s6 v3,a7
  a1 s1
  a2 s2
  s1 -1,a1
  s2 -1,a2
  a3 s3
  a4 s4
  s3 -1,a3
  s4 -1,a4
  a4 a4+a6
  a3 a3+a6
  a2 a2+a6
  s7 -1,a4
  s1 s1*rs5
  s2 s2*rs6
  s3 s3*rs6
  s6 -1,a3
  a1 a1+a6
  s4 s4*rs5
  s5 -1,a2
  s1 s7-fs1
  s7 -1,a1
  a0 a7-1
  a7 a7-1
  s2 s6-fs2
  s3 s5-fs3
  s4 s7-fs4
  -1,a4 s1
  -1,a3 s2
  -1,a2 s3
  -1,a1 s4
  jap qqq
* above costs about 105 clocks
  a3 b77
  a7 b71
  a7 a7+a3
  b71 a7
  a7 b70
  a7 a7+a3
  a7 a7+a3
  b70 a7
  a1 b73
  a3 b74
  a0 a1-1
  a1 a1-1
  s0 0
  jap aaaa
dbug  s1 1
  icount,0 s1
  ic4,0 s1
  exit
  end
          ident upak8z
*
*  call upak8z(n,in,iout,jout,kout,lout)
*  called by list,punch/serv sgmat,jkgen/util proc2/scf
*  not protected for n=0 or -ve
*  optimized for xmp not cray-1
*  compatible with eam to 16 mwords
*
*  coded v.r.saunders     march 1987
*
  align
  entry upak8z
upak8z a7 1,a6
  s2 2,a6
  s6 6,a6
  s5 5,a6
  s4 4,a6
  s3 3,a6
  s7 0,a7
  s1 1
  a2 2
  a1 8
  a3 16
  s7 s7+s1
  s6 s6+s1
  s5 s5+s1
  s4 s4+s1
  s7 s7>1
  s3 s3+s1
  a7 s7
top a0 s2
  vl a7
  s7 <8
  v0 ,a0,1  v0=in
  a0 s6
  s6 s6-s1
  v1 s7&v0  v1=lout(odd)
  ,a0,a2 v1  store lout(odd)
  a0 s5
  s5 s5-s1
  a5 vl
  v2 v0>a1  v2=in>8
  v3 s7&v2  v3=kout(odd)
  ,a0,a2 v3  store kout(odd)
  a0 s4
  s4 s4-s1
  a7 a7-a5
  v4 v0>a3  v4=in>16
  v5 s7&v4  v5=jout(odd)
  ,a0,a2 v5 store jout(odd)
  a0 s3
  s3 s3-s1
  s0 a7
  v6 v2>a3  v6=in>24
  v7 s7&v6  v7=iout(odd)
  ,a0,a2 v7  store iout(odd)
  a0 s6
  v0 v4>a3  v0=in>32
  v1 s7&v0  v1=lout(even)
  ,a0,a2 v1 store lout(even)
  a0 s5
  v2 v6>a3  v2=in>40
  v3 s7&v2  v3=kout(even)
  ,a0,a2 v3  store kout(even)
  a0 s4
  v4 v0>a3  v4=in>48
  v5 s7&v4  v5=jout(even)
  ,a0,a2 v5  save jout(even)
  s7 a5
  a0 s3
  s2 s2+s7
  s7 s7<1
  s7 s7+s1
  s6 s6+s7
  s5 s5+s7
  s4 s4+s7
  s3 s3+s7
  v6 v2>a3  v6=iout(even)
  ,a0,a2 v6  store iout(even)
  jsn top
  j b00
  end
*          ident gather
**
**  call gather(n,r,a,mapa)
**  r = a(mapa) for n elements
**  protected for n=0 but not for -ve n
**  optimized for xmp not cray-1
**  compatible with eam to 16 mwords
**
**  coded v.r.saunders     march 1987
**
*  align
*  entry gather
*gather a7 1,a6
*  s3 3,a6
*  s4 4,a6
*  s2 2,a6
*  s6 1
*  a1 zs0
*  a7 0,a7
*  s3 s3-s6
*  vl a7
*  s0 a7
*  a6 vl
*  a0 s4
*  s1 a6
*  a7 a7-a6
*  jsz return
*top v0 ,a0,1
*  a0 s3
*  s0 a7
*  v1 ,a0,v0
*  a0 s2
*  s4 s4+s1
*  s2 s2+s1
*  ,a0,1 v1
*  jsz return
*  a0 s4
*  vl a1
*  s1 a1
*  v2 ,a0,1
*  a0 s3
*  a7 a7-a1
*  v3 ,a0,v2
*  a0 s2
*  s0 a7
*  s4 s4+s1
*  s2 s2+s1
*  a7 a7-a1
*  ,a0,1 v3
*  a0 s4
*  jsn top
*return j b00
*  end
*          ident scatter
**
**  call scatter(n,r,mapr,a)
**  r(mapr) = a for n elements
**  protected for n=0 but not for -ve n
**  optimized for xmp not cray-1
**  compatible with eam to 16 mwords
**
**  coded v.r.saunders     march 1987
**
*  align
*  entry scatter
*scatter a7 1,a6
*  s2 2,a6
*  s4 4,a6
*  s3 3,a6
*  s6 1
*  a1 zs0
*  a7 0,a7
*  s2 s2-s6
*  vl a7
*  s0 a7
*  a6 vl
*  a0 s4
*  s1 a6
*  a7 a7-a6
*  jsz return
*top v0 ,a0,1
*  a0 s3
*  s0 a7
*  v1 ,a0,1
*  a0 s2
*  s4 s4+s1
*  s3 s3+s1
*  ,a0,v1 v0
*  jsz return
*  a0 s4
*  vl a1
*  s1 a1
*  v2 ,a0,1
*  a0 s3
*  a7 a7-a1
*  v3 ,a0,1
*  a0 s2
*  s0 a7
*  s4 s4+s1
*  s3 s3+s1
*  a7 a7-a1
*  ,a0,v3 v2
*  a0 s4
*  jsn top
*return j b00
*  end
          ident setsto
*
*  call setsto(n,scalar,r)
*  r = scalar for n elements
*  protected for n=0 but not for -ve n
*  optimized for xmp not cray-1
*  compatible with eam to 16 mwords
*
*  coded v.r.saunders     march 1987
*  modified hans nelemans cray research b.v. march 1994
*	this version is now cray c90 compatible
*
  align
  entry setsto
setsto a4 2,a6
  a7 1,a6
  s2 3,a6
  s6 -1
  a1 0
  vl a1                 full c90 vl
  a1 vl
  s5 a1                 save vectorlength
  vm s6                 set v-mask bit 0..63
  vm1 s6                     and bits 64..127
  s4 0,a4
  a7 0,a7
  a0 s2
  v3 s4!v4&vm
  vl a7
  s0 a7
  a3 vl
  s7 a7
  s3 a3
  jsz return
  ,a0,1 v3
  s0 s7\s3
  s2 s2+s3
  s7 s7-s3
  vl a1
  jsz return
top a0 s2
  s0 s7\s5
  ,a0,1 v3
  s2 s2+s5
  s7 s7-s5
  jsn top
return j b00
  end
*         ident szero
*
*  call szero(r,n)
*  r = 0.0 for n elements
*  protected for n=0 but not for -ve n
*  optimized for xmp not cray-1
*  compatible with eam to 16 mwords
*
*  coded v.r.saunders     march 1987
*   modified for c90 march 1994
*
*  align
*  entry szero
*szero a7 2,a6
*  s2 1,a6
*  a1 64
*  s5 a1
*  vl a1
*  v3 v0\v0
*  a7 0,a7
*  a0 s2
*  vl a7
*  s0 a7               s0=n 
*  s0 s0>6
*  jsz vlok            jif n<64
*  vl a1               set max vl
*vlok s0 a7
*  a4 vl
*  s7 a7
*  s4 a4
*  jsz return
*  ,a0,1 v3
*  s0 s7\s4
*  s2 s2+s4
*  s7 s7-s4
*  vl a1
*  jsz return
*top a0 s2
*  s0 s7\s5
*  ,a0,1 v3
*  s2 s2+s5
*  s7 s7-s5
*  jsn top
*return j b00
*  end
          ident fmove
*
*  call fmove(a,r,n)
*  r = a for n elements
*  protected for n=0 but not for -ve n
*  optimized for xmp not cray-1
*  compatible with eam to 16 mwords
*
*  coded v.r.saunders     march 1987
*
  align
  entry fmove
fmove a7 3,a6
  s3 1,a6
  s2 2,a6
  a6 zs0
  a7 0,a7
  s0 s2\s3
  jsz return
  vl a7
  s0 a7
  a5 vl
  s7 a7
  s5 a5
  a0 s3
  jsz return
top v0 ,a0,1
  a0 s2
  s0 s7\s5
  s3 s3+s5
  s7 s7-s5
  s2 s2+s5
  ,a0,1 v0
  jsz return
  a0 s3
  vl a6
  s5 a6
  v4 ,a0,1
  a0 s2
  s0 s7\s5
  s3 s3+s5
  s2 s2+s5
  s7 s7-s5
  ,a0,1 v4
  a0 s3
  jsn top
return j b00
  end
          ident ujau
*
*  call ujau(n,iscalar,ir,ia)   ---  ir=ia+iscalar
*
*  protected for n=0 but not n -ve
*  optimized for xmp not cray-1
*  compatible with eam to 16 mwords
*
*  coded v.r.saunders     may 1987
*
  align
  entry ujau
ujau a1 1,a6
  a2 2,a6
  s3 4,a6
  s6 3,a6
  a7 zs0
  a1 0,a1
  s2 0,a2
  vl a1
  s0 a1
  a5 vl
  s1 a1
  s5 a5
  a0 s3
  jsz return
top v6 ,a0,1
  s0 s1\s5
  a0 s6
  s3 s3+s5
  s6 s6+s5
  s1 s1-s5
  v5 s2+v6
  ,a0,1 v5
  jsz return
  vl a7
  a0 s3
  s5 a7
  v3 ,a0,1
  a0 s6
  s0 s1\s5
  s1 s1-s5
  s3 s3+s5
  s6 s6+s5
  v2 s2+v3
  ,a0,1 v2
  a0 s3
  jsn top
return j b00
  end
          ident scaler
*
*  call scaler(n,scalar,r,a)   ---  r=a*scalar
*
*  protected for n=0 but not n -ve
*  optimized for xmp not cray-1
*  compatible with eam to 16 mwords
*
*  coded v.r.saunders     may 1987
*
  align
  entry scaler
scaler a1 1,a6
  a2 2,a6
  s3 4,a6
  s6 3,a6
  a7 zs0
  a1 0,a1
  s2 0,a2
  vl a1
  s0 a1
  a5 vl
  s1 a1
  s5 a5
  a0 s3
  jsz return
top v6 ,a0,1
  s0 s1\s5
  a0 s6
  s3 s3+s5
  s6 s6+s5
  s1 s1-s5
  v5 s2*rv6
  ,a0,1 v5
  jsz return
  vl a7
  a0 s3
  s5 a7
  v3 ,a0,1
  a0 s6
  s0 s1\s5
  s1 s1-s5
  s3 s3+s5
  s6 s6+s5
  v2 s2*rv3
  ,a0,1 v2
  a0 s3
  jsn top
return j b00
  end
          ident swapv
*
*  call swapv(n,a,b)   --- swaps the two vectors a and b
*
*  protected for n=0 but not n -ve
*  optimized for xmp not cray-1
*  compatible with eam to 16 mwords
*
*  coded v.r.saunders     april 1987
*
  align
  entry swapv
swapv a1 1,a6
  s2 2,a6
  s3 3,a6
  a7 zs0
  a1 0,a1
  s7 a7
  s0 a1
  vl a1
  a0 s2
  s4 s2
  a6 vl
  jsz return
  v0 ,a0,1
  a0 s3
  s6 a6
  v1 ,a0,1
  a0 a1-a6
  s5 s3
  s2 s2+s6
  s3 s3+s6
  a1 a1-a6
  jaz exit1
top a0 s2
  vl a7
  v2 ,a0,1
  a0 s3
  a1 a1-a7
  v3 ,a0,1
  a0 s5
  vl a6
  s5 s3
  ,a0,1 v0
  a0 s4
  s0 a1
  s4 s2
  s2 s2+s7
  s3 s3+s7
  a6 a7
  ,a0,1 v1
  jsz exit2
  a0 s2
  vl a7
  v4 ,a0,1
  a0 s3
  a1 a1-a7
  v5 ,a0,1
  a0 s5
  s5 s3
  ,a0,1 v2
  a0 s4
  s0 a1
  s4 s2
  s2 s2+s7
  s3 s3+s7
  ,a0,1 v3
  jsz exit3
  a0 s2
  a1 a1-a7
  v0 ,a0,1
  a0 s3
  s0 a1
  v1 ,a0,1
  a0 s5
  s5 s3
  ,a0,1 v4
  a0 s4
  s4 s2
  ,a0,1 v5
  jsn top
exit1 a0 s5
  vl a6
  cmr
  ,a0,1 v0
  a0 s4
  ,a0,1 v1
return j b00
exit2 a0 s5
  vl a7
  cmr
  ,a0,1 v2
  a0 s4
  ,a0,1 v3
  j b00
exit3 a0 s5
  cmr
  ,a0,1 v4
  a0 s4
  ,a0,1 v5
  j b00
  end
          ident jkadd
*
***      subroutine jkadd(msmall,mbig,f,indf,ibias,p)
***      dimension f(1),indf(1),ibias(1),p(1)
***      ms=msmall
***      mb=mbig
***      ib=0
***1     m=min0(ms,mb)
***cdir$ ivdep
***      do 2 l=1,m
***2     f(indf(ib+l)+ibias(l))=f(indf(ib+l)+ibias(l))+p(ib+l)
***      mb=mb-m
***      ib=ib+ms
***      if(mb.ne.0)goto 1
***      return
***      end
*
*  msmall must be .le. 64
*  optimized for xmp not cray-1
*  compatible with eam to 16 mwords
*
*  coded v.r.saunders     march 1987
*   modified for c90 march 1994
*
  align
  entry jkadd
jkadd a1 1,a6
  a2 2,a6
  s5 5,a6
  s3 3,a6
  s4 4,a6
  s6 6,a6
  s1 1
  a1 0,a1
  a2 0,a2
  a4 zs0
  a0 s5
  s3 s3-s1
  vl a1
  v0 ,a0,1  v0=ibias
  a0 a2-a1
top jap topa
  vl a2
  s0 a2
  s0 s0>6   limit vl to 64
  jsz topa  jif vl<64
  vl a4
topa a0 s4
  a7 vl
  v7 ,a0,1  v7=indf
  a0 s6
  a2 a2-a7
  s7 a7
  v6 v0+v7  v6=indf+ibias
  v5 ,a0,1  v5=p
  a0 s3
  s0 a2
  v4 ,a0,v6  v4=f
  s4 s4+s7
  s6 s6+s7
  v3 v4+fv5  v3=f+p
  ,a0,v6 v3  scatter f+p
  a0 a2-a1
  jsz return
  jap topb
  vl a2
  s0 a2
  s0 s0>6    limit vl to 64
  jsz topb   jif vl<64
  vl a4
topb a0 s4
  a7 vl
  v7 ,a0,1  v7=indf
  a0 s6
  a2 a2-a7
  s7 a7
  v2 v0+v7  v2=indf+ibias
  v1 ,a0,1  v1=p
  a0 s3
  s0 a2
  v4 ,a0,v2  v4=f
  s4 s4+s7
  s6 s6+s7
  v3 v4+fv1  v3=f+p
  ,a0,v2 v3  scatter f+p
  a0 a2-a1
  jsn top
return j b00
  end
*          ident scatt
*
*  call scatt(n,scalar,v,mapv)
*  v(mapv) = scalar for n elements
*  protected for n=0 but not for n -ve
*  optimized for xmp not cray-1
*  compatible with eam to 16 mwords
*
*  coded v.r.saunders     march 1987
*   modified for c90 march 1994
*
*  align
*  entry scatt
*scatt a4 2,a6
*  a7 1,a6
*  s1 3,a6
*  s2 4,a6
*  s5 -1
*  a6 zs0
*  s4 0,a4
*  a7 0,a7
*  vm s5
*  vl a6
*  s1 s1+s5
*  v3 s4!v2&vm
*  vl a7
*  s0 a7
*  s0 s0>6
*  jsz vlok     jif vl<64
*  vl a6        set max vl to use     
*vlok s0 a7
*  a5 vl
*  a0 s2
*  a7 a7-a5
*  s5 a5
*  jsz return
*top v5 ,a0,1
*  a0 s1
*  s0 a7
*  s2 s2+s5
*  ,a0,v5 v3
*  jsz return
*  a0 s2
*  vl a6
*  s5 a6
*  v2 ,a0,1
*  a7 a7-a6
*  s2 s2+s5
*  a0 s1
*  s0 a7
*  a7 a7-a6
*  ,a0,v2 v3
*  a0 s2
*  jsn top
*return j b00
*  end
          ident srotg
*
*  call srotg(n,a,b,mapb,cos,sin)
*  a = a*cos + b*sin  and  b = b*cos - a*sin for n elements
*  where b is gathered/scattered under control of mapb
*  protected for n=0 but not for n -ve
*  optimized for xmp not cray-1
*  compatible with eam to 16 mwords
*
*  coded v.r.saunders     march 1987
*
  align
  entry srotg
srotg a7 1,a6      a7=<n>
  a2 6,a6      a2=<sin>
  a1 5,a6      a1=<cos>
  s4 3,a6      s4=<b>
  s3 2,a6      s3=<a>
  s5 4,a6      s5=<mapb>
  s6 1
  a7 0,a7      a7=n
  s2 0,a2      s2=sin
  s1 0,a1      s1=cos
  s4 s4-s6     s4=<b>-1
  s0 a7
  vl a7
  a0 s3
  a5 vl
  jsz return
top v7 ,a0,1     v7=a
  a0 s5        a0=<mapb>
  a7 a7-a5
  v6 ,a0,1     v6=mapb
  v5 s2*rv7    v5=a*sin
  a0 s4        a0=<b>-1
  s6 a5
  v4 ,a0,v6    v4=b(mapb)
  v3 s1*rv4    v3=b(mapb)*cos
  s0 a7
  v2 v3-fv5    v2=updated b(mapb)
  s5 s5+s6     update <mapb>
  v1 s1*rv7    v1=a*cos
  ,a0,v6 v2    store updated b(mapb)
  a0 s3        a0=<a>
  s3 s3+s6     update <a>
  v0 s2*rv4    v0=b(mapb)*sin
  v3 v0+fv1    v3=updated a
  ,a0,1 v3     store updated a
  vl a7
  a0 s3
  a5 vl
  jsn top
return j b00
  end
          ident srotv
*
*  call srotv(n,a,b,cos,sin)
*  a = a*cos + b*sin  and  b = b*cos - a*sin for n elements
*  not protected for n=0 or -ve
*  optimized for xmp not cray-1
*  compatible with eam to 16 mwords
*
*  coded v.r.saunders     march 1987
*
  align
  entry srotv
srotv a7 1,a6
  a2 5,a6
  a1 4,a6
  s4 3,a6
  s3 2,a6
  a7 0,a7
  s2 0,a2
  s1 0,a1
top a0 s4
  vl a7
  v7 ,a0,1
  a0 s3
  a6 vl
  v6 ,a0,1
  a7 a7-a6
  s6 a6
  v5 s2*rv7
  v4 s1*rv6
  s0 a7
  v3 v4+fv5
  ,a0,1 v3
  a0 s4
  s3 s3+s6
  s4 s4+s6
  v2 s1*rv7
  v1 s2*rv6
  v0 v2-fv1
  ,a0,1 v0
  jsn top
  j b00
  end
          ident triad
*
*  call triad(n,scalar,r,a)
*  r = a*scalar + r for n elements
*  protected for n=0 but not for -ve n
*  optimized for xmp not cray-1
*  compatible with eam to 16 mwords
*
*  coded v.r.saunders     march 1987
*
  align
  entry triad
triad a7 1,a6
  a5 2,a6
  s3 4,a6
  s2 3,a6
  a6 zs0
  a7 0,a7
  s5 0,a5
  vl a7
  a0 a7
  a1 vl
  s0 s5
  s1 a1
  a7 a7-a1
  jaz return
  a0 s3
  jsz return
top v3 ,a0,1
  a0 s2
  s0 a7
  v2 ,a0,1
  v0 s5*rv3
  v1 v2+fv0
  s3 s3+s1
  s2 s2+s1
  ,a0,1 v1
  jsz return
  a0 s3
  vl a6
  s1 a6
  v7 ,a0,1
  a0 s2
  a7 a7-a6
  s3 s3+s1
  v6 ,a0,1
  v4 s5*rv7
  s0 a7
  s2 s2+s1
  v5 v6+fv4
  a7 a7-a6
  ,a0,1 v5
  a0 s3
  jsn top
return j b00
  end
          ident gtriad
*
*  call gtriad(n,scalar,r,b,a)
*  r = a*scalar + b for n elements
*  protected for n=0 but not for -ve n
*  optimized for xmp not cray-1
*  compatible with eam to 16 mwords
*
*  coded v.r.saunders     march 1987
*
  align
  entry gtriad
gtriad a7 1,a6
  a5 2,a6
  s3 5,a6
  s4 4,a6
  s2 3,a6
  a6 zs0
  a7 0,a7
  s5 0,a5
  vl a7
  s0 a7
  a1 vl
  a0 s3
  a7 a7-a1
  s1 a1
  jsz return
top v3 ,a0,1
  a0 s4
  s3 s3+s1
  v2 ,a0,1
  a0 s2
  v0 s5*rv3
  s0 a7
  v1 v0+fv2
  s4 s4+s1
  s2 s2+s1
  ,a0,1 v1
  jsz return
  a0 s3
  vl a6
  a7 a7-a6
  v7 ,a0,1
  a0 s4
  s1 a6
  v6 ,a0,1
  v4 s5*rv7
  s0 a7
  a0 s2
  a7 a7-a6
  s3 s3+s1
  s4 s4+s1
  s2 s2+s1
  v5 v4+fv6
  ,a0,1 v5
  a0 s3
  jsn top
return j b00
  end
          ident gtrian
*
*  call gtrian(n,scalar,r,b,a)
*  r = a*scalar - b for n elements
*  protected for n=0 but not for -ve n
*  optimized for xmp not cray-1
*  compatible with eam to 16 mwords
*
*  coded v.r.saunders     march 1987
*
  align
  entry gtrian
gtrian a7 1,a6
  a5 2,a6
  s3 5,a6
  s4 4,a6
  s2 3,a6
  a6 zs0
  a7 0,a7
  s5 0,a5
  vl a7
  s0 a7
  a1 vl
  a0 s3
  a7 a7-a1
  s1 a1
  jsz return
top v3 ,a0,1
  a0 s4
  s3 s3+s1
  v2 ,a0,1
  v0 s5*rv3
  a0 s2
  s0 a7
  v1 v0-fv2
  s4 s4+s1
  s2 s2+s1
  ,a0,1 v1
  jsz return
  a0 s3
  vl a6
  a7 a7-a6
  v7 ,a0,1
  a0 s4
  s1 a6
  v6 ,a0,1
  v4 s5*rv7
  s0 a7
  a0 s2
  a7 a7-a6
  s3 s3+s1
  s4 s4+s1
  s2 s2+s1
  v5 v4-fv6
  ,a0,1 v5
  a0 s3
  jsn top
return j b00
  end
          ident vvtv
*
*  call vvtv(n,r,a,b)
*  r = a * b for n elements
*  protected for n=0 but not for -ve n
*  optimized for xmp not cray-1
*  compatible with eam to 16 mwords
*
*  coded v.r.saunders     march 1987
*
  align
  entry vvtv
vvtv a7 1,a6
  s3 3,a6
  s4 4,a6
  s2 2,a6
  a6 zs0
  a7 0,a7
  vl a7
  s0 a7
  a5 vl
  s7 a7
  a0 s3
  s5 a5
  jsz return
top v3 ,a0,1
  a0 s4
  s0 s7\s5
  v2 ,a0,1
  a0 s2
  v1 v2*rv3
  s3 s3+s5
  s7 s7-s5
  s4 s4+s5
  s2 s2+s5
  ,a0,1 v1
  jsz return
  a0 s3
  vl a6
  s5 a6
  v7 ,a0,1
  a0 s4
  s0 s7\s5
  v6 ,a0,1
  a0 s2
  s3 s3+s5
  v5 v6*rv7
  s7 s7-s5
  s4 s4+s5
  s2 s2+s5
  ,a0,1 v5
  a0 s3
  jsn top
return j b00
  end
_IFN(charmm)
          ident addvec
*
*  call addvec(r,a,b,n)
*  r = a + b for n elements
*  protected for n=0 but not for -ve n
*  optimized for xmp not cray-1
*  compatible with eam to 16 mwords
*
*  coded v.r.saunders     march 1987
*
  align
  entry addvec
addvec a7 4,a6
  s3 2,a6
  s4 3,a6
  s2 1,a6
  a6 zs0
  a7 0,a7
  vl a7
  s0 a7
  a5 vl
  s7 a7
  a0 s3
  s5 a5
  jsz return
top v3 ,a0,1
  a0 s4
  s0 s7\s5
  v2 ,a0,1
  a0 s2
  v1 v2+fv3
  s3 s3+s5
  s7 s7-s5
  s4 s4+s5
  s2 s2+s5
  ,a0,1 v1
  jsz return
  a0 s3
  vl a6
  s5 a6
  v7 ,a0,1
  a0 s4
  s0 s7\s5
  v6 ,a0,1
  a0 s2
  s3 s3+s5
  v5 v6+fv7
  s7 s7-s5
  s4 s4+s5
  s2 s2+s5
  ,a0,1 v5
  a0 s3
  jsn top
return j b00
  end
_ENDIF
          ident subvec
*
*  call subvec(r,a,b,n)
*  r = a + b for n elements
*  protected for n=0 but not for -ve n
*  optimized for xmp not cray-1
*  compatible with eam to 16 mwords
*
*  coded v.r.saunders     march 1987
*
  align
  entry subvec
subvec a7 4,a6
  s3 2,a6
  s4 3,a6
  s2 1,a6
  a6 zs0
  a7 0,a7
  vl a7
  s0 a7
  a5 vl
  s7 a7
  a0 s3
  s5 a5
  jsz return
top v3 ,a0,1
  a0 s4
  s0 s7\s5
  v2 ,a0,1
  a0 s2
  v1 v3-fv2
  s3 s3+s5
  s7 s7-s5
  s4 s4+s5
  s2 s2+s5
  ,a0,1 v1
  jsz return
  a0 s3
  vl a6
  s5 a6
  v7 ,a0,1
  a0 s4
  s0 s7\s5
  v6 ,a0,1
  a0 s2
  s3 s3+s5
  v5 v7-fv6
  s7 s7-s5
  s4 s4+s5
  s2 s2+s5
  ,a0,1 v5
  a0 s3
  jsn top
return j b00
  end
          ident jkuno
*
*  call jkuno(n,p,g,i,j,k,l,indf,pg)
*
*  to generate triangle matrix indices ik,jl,il,jk,ij,kl
*  and p*g contributions for k<p> j<p>
*
*  not protected for n=0
*  optimized for xmp not cray-1
*  compatible with eam to 16 mwords
*  but note that indf,pg must lie within first 8mword
*
*  coded v.r.saunders     april 1987
*   modified for c90 march 1994 (hans nelemans cray)
*
  common mapper
iky bss 1
  block *
  align
  entry jkuno
jkuno a1 1,a6
  s7 2,a6
  s1 3,a6
  s2 4,a6
  s3 5,a6
  s4 6,a6
  s5 7,a6
  a1 0,a1
  a3 8,a6
  a2 9,a6
  s6 1
  a4 6
  s7 s7-s6
  a5 iky-1
  a7 s7
  b77 a7
*
*  a1=n a2=<pg> a3=<indf> a4=6 a5=<iky>-1 b77=<p>-1
*  s1=<g> s2=<i> s3=<j> s4=<k> s5=<l>
*
top vl a1
  s0 a1
  s0 s0>6
  jsz vlok                 jif vl.lt.64
  a6 zs0
  vl a6                    set max vl
vlok a0 s2
  a6 4
  v0 ,a0,1                 v0=i
  a0 s4
  a7 5
  v1 ,a0,1                 v1=k
  a0 s3
  v2 v0-v1                 v2=i-k
  vm v2,z                  vm=i.eq.k
* v0=i v1-k v2=i-k vm=ieqk
  v3 ,a0,1                 v3=j
  a0 a5
  v4 ,a0,v0                v4=ikyi
  v5 v3-v1                 v5=j-k
  a0 s5
  s6 vm                    s6=i.eq.k
  vm v5,z                  vm=j.eq.k
* v0=i v1=k v2=i-k v3=j v4=ikyi v5=j-k vm=jeqk
  v6 ,a0,1                 v6=l
  a0 a5
  v7 ,a0,v1                v7=ikyk
  v2 v4+v3                 v2=ij
  a0 a3+a6
  s7 vm                    s7=j.eq.k
  vm v5,p                  vm=j.ge.k
  ,a0,a4 v2                store ij
* v0=i v1=k v2=ij v3=j v4=ikyi v5=j-k v6=l v7=ikyk vm=jgek
  a0 s1
  v0 ,a0,1                 v0=g
  a0 b77
  v2 ,a0,v2                v2=pij
  v5 v7+v6                 v5=kl
  a0 a3+a7
  ,a0,a4 v5                store kl
* v0=g v1=k v2=pij v3=j v4=ikyi v5=kl v6=l v7=ikyk vm=jgek
  a0 b77
  v6 v6+v4                 v6=il
  v5 ,a0,v5                v5=pkl
  v1 v0*rv2                v1=jkl/2
  a0 a2+a7
  a7 3
  ,a0,a4 v1                store jkl/2
* v0=g v1=jkl v2=pij v3=j v4=ikyi v5=pkl v6=il v7=ikyk vm=jgek
  a0 a5
  v2 ,a0,v3                v2=ikyj
  v1 v0*rv5                v1=jij/2
  a0 a2+a6
  a6 2
  ,a0,a4 v1                store jij/2
* v0=g v1=jij v2=ikyj v3=j v4=ikyi v5=pkl v6=il v7=ikyk vm=jgek
  a0 s4
  v5 ,a0,1                 v5=k
  a0 a3+a6
  v7 v7+v3                 v7=xkj
  ,a0,a4 v6                store il
* v0=g v1=jij v2=ikyj v3=j v4=ikyi v5=k v6=il v7=xkj vm=jgek
  a0 b77
  v1 v2+v5                 v1=xjk
  v3 v1!v7&vm              v3=jk
  v6 ,a0,v6                v6=pil
  a0 a3+a7
  ,a0,a4 v3                store jk
* v0=g v1=xjk v2=ikyj v3=jk v4=ikyi v5=k v6=pil v7=xkj vm=jgek
  a0 b77
  vm s7                    vm=j.eq.k
  v4 v4+v5                 v4=ik
  v7 v0+fv0                v7=g2
  v3 ,a0,v3                v3=pjk
* v0=g v1=xjk v2=ikyj v3=pjk v4=ik v5=k v6=pil v7=g2 vm=jeqk
  a0 a2+a7
  a7 vl
  v1 v7!v0&vm              v1=gil
  v5 v1*rv6                v5=kjk
  ,a0,a4 v5                store kjk
* v0=g v1=gil v2=ikyj v3=pjk v4=ik v5=kjk v6=pil v7=g2 vm=jeqk
  a0 s5
  v7 ,a0,1                 v7=l
  a0 s3
  v6 ,a0,1                 v6=j
  a0 a3
  v3 v1*rv3                v3=kil
  v5 v6-v7                 v5=j-l
  vm v5,z                  vm=j.eq.l
  ,a0,a4 v4                 store ik
* v0=g v1=gil v2=ikyj v3=kil v4=ik v5=j-l v6=j v7=l vm=jeql
  a0 a5
  a1 a1-a7
  v1 ,a0,v7                v1=ikyl
  a0 a2+a6
  a6 a7*a4
  s0 a1
  v6 v1+v6                 v6=xlj
  s7 vm                    s7=j.eq.l
  vm v5,p                  vm=j.ge.l
  ,a0,a4 v3                store kil
* v0=g v1=ikyl v2=ikyj v3=kil v4=ik v5=j-l v6=xlj v7=l vm=jgel
  a0 b77
  s6 s6!s7                 s6=i.eq.k.or.j.eq.l
  v4 ,a0,v4                v4=pik
  v5 v0+fv0                v5=g2
  v1 v2+v7                 v1=xjl
  a0 a3+1
  v3 v1!v6&vm              v3=jl
  vm s6                    vm=i.eq.k.or.j.eq.l
  s7 a7
  ,a0,a4 v3                store jl
* v0=g v1=xjl v2=ikyj v3=jl v4=pik v5=g2 v6=xlj v7=l vm=ieqkorjeql
  a0 b77
  v7 ,a0,v3                v7=pjl
  v6 v5!v0&vm              v6=gik
  v2 v6*rv4                v2=kjl
  a0 a2+1
  ,a0,a4 v2                store kjl
* v0=g v1=xjl v2=kjl v3=jl v4=pik v5=g2 v6=gik v7=pjl
  a0 a2
  s1 s1+s7
  s2 s2+s7
  s3 s3+s7
  s4 s4+s7
  s5 s5+s7
  v5 v6*rv7                v5=kik
  a3 a3+a6
  a2 a2+a6
  ,a0,a4 v5                store kik
  jsn top
  j b00
  end
_ENDIF
_ENDIF
_IF(cyber205)
  output
ind ident
*
*    coded    v.r. saunders ---- aug 85.
*
mapper msec 3
viky res,128 64
  msec 2
i equ 3*64
j equ 4*64
iky equ 5*64
ii equ 6*64
jj equ 7*64
r16 equ #16*64
r17 equ #17*64
r18 equ #18*64
r1a equ #1a*64
  entry ind
ind lod [r17],i
  lod [r17,r16],j
  ex iky,viky-64
  lod [i],i
  lod [j],j
  lod [iky,i],ii
  lod [iky,j],jj
  ibxgt,rel ,j,p1,i,
  addx ii,j,r18
  barb,br ,r1a
p1 addx i,jj,r18
  barb,br ,r1a
  end
ihint ident
*
*    coded    v.r. saunders ---- aug 85.
*
*      function ihint(ii,jj)
*      common/helpr/isym(200),iky(200),niky(8),m1e,m2e,m1ia,m2coul,
*     *m2exch,m2fock,m2psup,m2qsup,
*     *mbase(288),nbase(8),npairi(36),norbi(8),mbasi(8),mnbasi(8)
*     *,iap(62),mbfock(62),mbexch(36),mbcoul(36),mbpsup(36),mbqsup(36)
*     *,norbas(200)
*c...
*c... to compute index of 1-electron integral
*c...
*      if(ii.ge.jj)goto 111
*      i=jj
*      j=ii
*      goto 1
*111   i=ii
*      j=jj
*1     ihint=iky(norbas(i))+j+nbase(isym(i))
*      return
*      end
helpr msec 3
visym res,128 200*64
viky res,64 200*64
vniky res,64 8*64
vm1e res,64 64
vm2e res,64 64
vm1ia res,64 64
vm2coul res,64 64
vm2exch res,64 64
vm2fock res,64 64
vm2psup res,64 64
vm2qsup res,64 64
vmbase res,64 288*64
vnbase res,64 8*64
vnpairi res,64 36*64
vnorbi res,64 8*64
vmbasi res,64 8*64
vmnbasi res,64 8*64
viap res,64 62*64
vmbfock res,64 62*64
vmbexch res,64 36*64
vmbcoul res,64 36*64
vmbpsup res,64 36*64
vmbqsup res,64 36*64
vnorbas res,64 200*64
  msec 2
i equ 3*64
j equ 4*64
ii equ 5*64
jj equ 6*64
isym equ 7*64
norbas equ 8*64
nbase equ 9*64
iky equ #a*64
r16 equ #16*64
r17 equ #17*64
r18 equ #18*64
r1a equ #1a*64
  entry ihint
ihint lod [r17],ii
  lod [r17,r16],jj
  ex isym,visym-64
  ex norbas,vnorbas-64
  ex nbase,vnbase-64
  ex iky,viky-64
  lod [ii],ii
  lod [jj],jj
  rtor ii,i
  rtor jj,j
  ibxge,rel ,ii,ihint1,jj,
  rtor jj,i
  rtor ii,j
ihint1 lod [isym,i],isym
  lod [norbas,i],norbas
  lod [nbase,isym],nbase
  lod [iky,norbas],iky
  addx nbase,j,r18
  addx iky,r18,r18
  barb,br ,r1a
  end
identi ident
*
*    coded    v.r. saunders ---- feb 85.
*
*               fortran equivalent
*
*      subroutine identi(nrow,ncol,q,ilifq)
*      dimension q(1),ilifq(1)
*      do 1 i=1,ncol
*      call szero(q(ilifq(i)+1),nrow)
*1     q(ilifq(i)+i)=1.0
*      return
*      end
*
scr205 msec 3
vtemp res,128 64
  msec 2
q equ 3*64
ilifq equ 4*64
nrow equ 5*64
ncol equ 6*64
temp equ 7*64
one equ 8*64
six equ one
loop equ one
ilif equ 9*64
*
r16 equ #16*64
r17 equ #17*64
r1a equ #1a*64
*
  entry identi
identi lod [r17,r16],ncol
  lod [r17],nrow
  es ilifq,3
  lod [r17,ilifq],ilifq
  es q,2
  lod [r17,q],q
  ex temp,vtemp
  lod [ncol],ncol
  lod [nrow],nrow
  es six,6
  pack ncol,temp,temp
  pack ncol,ilifq,ilifq
  pack nrow,q,q
  linkv,rb
  shiftv,b [ilifq],six,temp,
  addxv,a q,[temp],temp,
  es loop,0
vrs lod [temp,loop],ilif
  vtov,a ,ilif,
  ibxlt,rel r16,loop,vrs,ncol,loop
  interval ,r16,temp,
  addxv [ilifq],[temp],temp,
  addn ,r16,one
  vtovx,b [temp],one,q
  barb,br ,r1a
  end
diag ident
*
*    coded    v.r. saunders ---- feb 85.
*
*               fortran equivalent
*
*      subroutine diag(a,iky,newbas,q,ilifq,nrow,e,thresh,
*     1 small,tei)
*      dimension a(2),iky(2),q(2),ilifq(2),e(2)
*      common/scr205/c1(512),c2(512),y(512)
*      nm1=newbas-1
*      teisq=tei*tei
*      nvrs=iky(newbas+1)-2
*      call scatt(nm1,0.0,a,iky(2))
*20500 te=absmax(nvrs,0.0,a(2))
*      if(te.lt.thresh)return
*      temp=tei
*      if(te.lt.teisq)temp=dsqrt(te)
*      te=te*temp
*      do 2 i=1,nm1
*      i1=i-1
*      ip1=i+1
*      ikyi=iky(i)
*      itest=newbas-i
*      ilifi=ilifq(i)
*      aii=e(i)
*      call fmove(a(ikyi+1),y,i)
*      call gather(itest,y(ip1),a(ip1),iky(ip1))
*      do 22 j=ip1,newbas
*      vij=y(j)
*      avij=dabs(vij)
*      if(avij.lt.te) go to 22
*      vjj=e(j)
*      j1=j-1
*      jp1=j+1
*      temp=aii-vjj
*      if(avij.ge.(small*dabs(temp)))goto 97531
*c... determine rotation parameters (3rd order pert. theory)
*      tem=vij/temp
*      temp=tem*tem
*      sint=tem-(1.5*tem*temp)
*      cost=1.0-(temp*0.5)
*      goto 86420
*c... determine rotation parameters (rutishauser algorithm)
*97531 temp=temp*0.5
*      tem=vij/(dsqrt(temp*temp+vij*vij)+dabs(temp))
*      cost=dsqrt(1.0/(tem*tem+1.0))
*      if(temp)12345,54321,54321
*12345 tem=-tem
*54321 sint=cost*tem
*86420 call drot(j1,y,a(iky(j)+1),cost,sint)
*      tem=tem*vij
*      aii=aii+tem
*      e(j)=vjj-tem
*      y(j)=0.0
*      if(jp1.gt.newbas)goto 7
*      do 8 k=jp1,newbas
*      jj=j+iky(k)
*      vij=y(k)
*      y(k)=vij*cost+a(jj)*sint
*8     a(jj)=a(jj)*cost-vij*sint
*7     call drot(nrow,q(ilifi+1),q(ilifq(j)+1),cost,sint)
*22    continue
*      call fmove(y,a(ikyi+1),i1)
*      call scatter(itest,a(ip1),iky(ip1),y(ip1))
*2     e(i)=aii
*      goto 20500
*      end
*
scr205 msec 3
c1b res,128 512*64
c2b res,64 512*64
yb res,64 512*64
  msec 2
r15 equ #15*64
r16 equ #16*64
r17 equ #17*64
r1a equ #1a*64
r1c equ #1c*64
*
one equ #13*64
tei equ #12*64
half equ #11*64
te equ #10*64
e equ #f*64
nrow equ #e*64
onep5 equ nrow
ilifq equ #d*64
q equ #c*64
newbas equ #b*64
iky equ #a*64
a equ 9*64
c1 equ 8*64
c2 equ 7*64
y equ 6*64
thresh equ 5*64
nvrs equ 4*64
yi equ nvrs
small equ 3*64
*
i1 equ #1a*64
nm1 equ #1b*64
ikyi equ #1d*64
aii equ #1e*64
ilifi equ #1f*64
i equ #20*64
j1 equ i
ikyip1 equ #21*64
aip1 equ #22*64
ilifj equ #23*64
ikyj equ #24*64
temp equ #25*64
jtemp equ temp
tem equ #26*64
yj equ tem
vij equ #27*64
cost equ #28*64
sint equ #29*64
avij equ #2a*64
atemp equ avij
aj equ avij
temph equ #2b*64
j equ temph
yip1 equ #2c*64
vjj equ #2d*64
q22e equ #2e*64
q97531 equ #2f*64
graham equ #30*64
nrowc1 equ #31*64
nrowc2 equ #32*64
alan equ #33*64
teisq equ #34*64
*
  entry diag
diag es newbas,2
  lod [r17,newbas],newbas
  es nrow,5
  lod [r17,nrow],nrow
  es thresh,7
  lod [r17,thresh],thresh
  es small,8
  lod [r17,small],small
  es tei,9
  lod [r17,tei],tei
  lod [r17,r16],iky
  lod [r17],a
  es ilifq,4
  lod [r17,ilifq],ilifq
  es q,3
  lod [r17,q],q
  es e,6
  lod [r17,e],e
  lod [newbas],newbas
  lod [nrow],nrow
  lod [thresh],thresh
  lod [small],small
  lod [tei],tei
  shifti r16,46,half
  addn ,r16,one
  ex c1,c1b
  ex c2,c2b
  ex y,yb
  elen half,#ffd1   generate  f.p.  0.5
  lod [iky,newbas],nvrs
  swap ,r15,r1c
  es alan,64
  subx a,alan,atemp
  subx newbas,r16,nm1
  pack nm1,iky,ikyi
  addx ikyi,alan,ikyi
  vtovx,b [ikyi,r16],,atemp
  subx nvrs,r16,nvrs
  pack nvrs,a,a
  max,ma [a,r16],vij,te,
  mpys tei,tei,teisq
  pack nrow,q,q
  pack nrow,c1,nrowc1
  pack nrow,c2,nrowc2
  ex alan,gary
  ex q22e,p22e
  ex q97531,p97531
  ex graham,neil
  addn one,half,onep5    generate  f.p.  1.5
p20500 abs te,te
  es i1,0
  sqrt te,temp
  bge te,thresh,alan
  swap r1c,r15,
  barb,br ,r1a
gary blt te,teisq,graham
  mpys te,tei,te
  bim ,p2
neil mpys te,temp,te
p2 lod [iky,i1],ikyi
  lod [ilifq,i1],ilifi
  lod [e,i1],aii
  addx i1,r16,i
  pack i,y,yi
  shifti ikyi,6,ikyi
  shifti ilifi,6,ilifi
  addx a,ikyi,ikyi
  vtov [ikyi],yi,
  shifti i,6,tem
  subx newbas,i,temp
  addx q,ilifi,ilifi
  addx iky,tem,ikyip1
  addx a,tem,aip1
  is aip1,-64
  addx y,tem,yip1
  pack temp,ikyip1,ikyip1
  vxtov [ikyip1],aip1,yip1
p22 lod [y,j1],vij
  lod [e,j1],vjj
  lod [iky,j1],ikyj
  lod [ilifq,j1],ilifj
  pack j1,y,y
  pack j1,c1,c1
  abs vij,avij
  subn aii,vjj,temp
  shifti ikyj,6,ikyj
  shifti ilifj,6,ilifj
  blt avij,te,q22e
  mpys small,temp,sint
  mpys half,temp,temph
  mpys vij,vij,cost
  addx a,ikyj,ikyj
  addx q,ilifj,ilifj
  blt sint,avij,q97531
*   determine rotation parameters (3rd order pert. theory)
  divs vij,temp,tem
  mpys tem,tem,temp
  mpys tem,onep5,sint
  mpys half,temp,temph
  mpys sint,temp,sint
  subn one,temph,cost
  subn tem,sint,sint
  bim ,p86420
*   determine rotation parameters (rutishauser algorithm)
p97531 mpys temph,temph,tem
  abs temph,atemp
  addn tem,cost,tem
  sqrt tem,tem
  addn tem,atemp,tem
  divs vij,tem,tem
  mpys tem,tem,cost
  addn one,cost,cost
  divs one,cost,cost
  sqrt cost,cost
  ibxge,rel ,temph,p54321,,
  subn ,tem,tem
p54321 mpys cost,tem,sint
p86420 mpysv,a cost,[y],c1,
  pack j1,c2,c2
  mpysv,a sint,[y],c2,
  mpys tem,vij,tem
  pack j1,ikyj,ikyj
  addx r16,j1,j
  linkv,ra
  mpysv,b [ikyj],sint,y,
  addnv [y],[c1],y,
  addn aii,tem,aii
  subn vjj,tem,vjj
  subx newbas,j,jtemp
  linkv,ra
  mpysv,b [ikyj],cost,ikyj,
  subnv [ikyj],[c2],ikyj,
  ibxeq,rel ,jtemp,p7,,
  shifti j,6,yj
  addx iky,yj,ikyj
  pack jtemp,ikyj,ikyj
  addx a,yj,aj
  is aj,-64
  vxtov [ikyj],aj,c1
  pack jtemp,c1,c1
  pack jtemp,c2,c2
  mpysv,a cost,[c1],c2,
  mpysv,a sint,[c1],c1,
  addx y,yj,yj
  pack jtemp,yj,yj
  linkv,rb
  mpysv,a sint,[yj],c2,
  subnv [c2],[c2],c2,
  linkv,rb
  mpysv,a cost,[yj],yj,
  addnv [c1],[yj],yj,
  vtovx [ikyj],c2,aj
p7 sto [e,j1],vjj
  sto [y,j1],
  mpysv,a cost,[ilifi],nrowc1,
  mpysv,a sint,[ilifi],nrowc2
  linkv,ra
  mpysv,b [ilifj],sint,ilifi,
  addnv [ilifi],[nrowc1],ilifi,
  linkv,ra
  mpysv,b [ilifj],cost,ilifj,
  subnv [ilifj],[nrowc2],ilifj,
p22e ibxlt,rel r16,j1,p22,newbas,j1
  sto [e,i1],aii
  vtovx [ikyip1],yip1,aip1
  ibxeq,rel ,i1,vrs,,
  pack i1,ikyi,ikyi
  vtov [yi],ikyi,
vrs ibxlt,rel r16,i1,p2,nm1,i1
  max,ma [a,r16],vij,te,
  bim ,p20500
  end
mvges ident
*
*     nbig=mvges(nvec,range,r,a,masker)
*
*    coded     v.r.saunders ---- oct 84 for vamp-205.
*
  msec 2
nvec equ 3*64
range equ 4*64
r equ 5*64
a equ 6*64
masker equ 7*64
r16 equ #16*64
r17 equ #17*64
r18 equ #18*64
r1a equ #1a*64
  entry mvges
mvges lod [r17],nvec
  lod [r17,r16],range
  es r,2
  lod [r17,r],r
  es a,3
  lod [r17,a],a
  es masker,4
  lod [r17,masker],masker
  lod [nvec],nvec
  lod [range],range
  pack nvec,a,a
  arithcps,b [a],range,r,masker
  ltor r,r18
  shifti r18,6,range
  addx r,range,r
  cpsv,z a,r,masker
  barb,br ,r1a
  end
mrgle ident
scr205 msec 3
temp res,128 64
maskvrs msec 3
maskx res,128 2048*64
  msec 2
n equ 3*64
scalar equ 4*64
scalneg equ 5*64
r equ 6*64
masker equ 7*64
ffff equ 8*64
tem equ 9*64
temlen equ #a*64
lenr equ #b*64
lenr64 equ #c*64
r16 equ #16*64
r17 equ #17*64
r1a equ #1a*64
  entry mrgle
*    call mrgle(n,scalar,r)
*    r = dsign(scalar,r) if dabs(r) .lt. scalar
*
*    coded     v.r.saunders ---- nov 84.
*
mrgle lod [r17,r16],scalar
  lod [r17],n
  es r,2
  lod [r17,r],r
  ex masker,maskx
  ex tem,temp
  lod [scalar],scalar
  lod [n],n
  elen r,#ffff
  ltor r,ffff
  subn ,scalar,scalneg
alan ibxge,rel ,n,gary,ffff,
  pack n,r,r
gary arithcps,b,n [r],scalneg,tem,masker
  ltor r,lenr
  shifti lenr,6,lenr64
  subx n,lenr,n
  ltor tem,temlen
  ibxeq,rel ,temlen,graham,,
  shifti temlen,6,temlen
  addx tem,temlen,temlen
  cmpge,b [tem],,temlen
  maskv,a,b scalar,scalneg,tem,temlen
  mrgv,sb tem,r,r,masker
graham ibxle ,n,[r1a],,
  addx r,lenr64,r
  bim ,alan
  end
mrgge ident
scr205 msec 3
temp res,128 64
maskvrs msec 3
maskx res,128 2048*64
  msec 2
n equ 3*64
scalar equ 4*64
scalneg equ 5*64
r equ 6*64
masker equ 7*64
ffff equ 8*64
tem equ 9*64
temlen equ #a*64
lenr equ #b*64
lenr64 equ #c*64
r16 equ #16*64
r17 equ #17*64
r1a equ #1a*64
  entry mrgge
*    call mrgge(n,scalar,r)
*    r = dsign(scalar,r) if dabs(r) .ge. scalar
*
*    coded     v.r.saunders ---- nov 84.
*
mrgge lod [r17,r16],scalar
  lod [r17],n
  es r,2
  lod [r17,r],r
  ex masker,maskx
  ex tem,temp
  lod [scalar],scalar
  lod [n],n
  elen r,#ffff
  ltor r,ffff
  subn ,scalar,scalneg
alan ibxge,rel ,n,gary,ffff,
  pack n,r,r
gary arithcps,b,ma [r],scalar,tem,masker
  ltor r,lenr
  shifti lenr,6,lenr64
  subx n,lenr,n
  ltor tem,temlen
  ibxeq,rel ,temlen,graham,,
  shifti temlen,6,temlen
  addx tem,temlen,temlen
  cmpge,b [tem],,temlen
  maskv,a,b scalar,scalneg,tem,temlen
  mrgv,sb tem,r,r,masker
graham ibxle ,n,[r1a],,
  addx r,lenr64,r
  bim ,alan
  end
myml ident
maskvrs msec 3
vmask res,128 2048*64
  msec 2
*
*    long-vector matrix multiply routine  (r = a * b)
*    coded    v.r. saunders ---- aug 85.
*
*      subroutine myml(a,b,r,mcolr,mrowr,ncol,link,nrow)
*      dimension a(1),b(1),r(1)
*      dimension index(1)
*      common/scr205/x(1)
*      equivalence (index(1),x(1))
*
*      common/scr205/ is actually in the dynamic stack
*
*      r(ncol,nrow) has col,row spacings mcolr,mrowr
*      a is a ncol*link matrix stored row-wise and tight packed
*      b is a link*nrow matrix stored col-wise and tight packed
*      ncol-1 extra columns have been appended to the b matrix
*      ncol  ***must*** be .ge. nrow
*
r16 equ #16*64
r17 equ #17*64
r18 equ #18*64
xx equ r18
r19 equ #19*64
yy equ r19
r1a equ #1a*64
r1b equ #1b*64
*
a equ 3*64
ncnr64 equ a
mcmr equ a
b equ 4*64
stride equ b
mink equ b
inc equ b
r equ 5*64
mcolr equ 6*64
rr equ mcolr
mrowr equ 7*64
ncol equ 8*64
link equ 9*64
xxlenh equ link
n equ link
nrow equ #a*64
ncli equ #b*64
lenh equ ncli
m equ ncli
ncli64 equ #c*64
mask equ ncli64
link64 equ #d*64
nink equ link64
limit equ link64
x equ #e*64
lenz equ x
ss equ x
loop equ #f*64
lenx equ #10*64
tt equ lenx
btemp equ lenx
y equ #11*64
uu equ y
nrem equ y
ncnr equ #13*64
ngot equ ncnr
  entry myml
myml es link,6
  lod [r17,link],link
  es ncol,5
  lod [r17,ncol],ncol
  es nrow,7
  lod [r17,nrow],nrow
  es mrowr,4
  lod [r17,mrowr],mrowr
  es mcolr,3
  lod [r17,mcolr],mcolr
  es r,2
  lod [r17,r],r
  lod [r17,r16],b
  lod [r17],a
  lod [link],link
  lod [ncol],ncol
  lod [nrow],nrow
  lod [mrowr],mrowr
  lod [mcolr],mcolr
*
*      ncli=ncol*link
*      ngot=nrow*link
*      nrem=ncli-link
*101   if(nrem.eq.0)goto 100
*      lenz=min0(ngot,nrem)
*      call fmove(b,b(ngot+1),lenz)
*      nrem=nrem-lenz
*      ngot=ngot+lenz
*      goto 101
*
  mpyx ncol,link,ncli
  ibxne,rel ,ncol,p100,r16,
  pack link,a,a
  pack link,b,b
  dotv [a],[b],r18,
  addn ,r18,r18
  addn r18,r19,r18
  sto [r],r18
  barb,br ,r1a
p100 mpyx nrow,link,ngot
  subx ncli,link,nrem
p102 rtor nrem,lenz
  ibxle,rel ,nrem,p101,ngot,
  rtor ngot,lenz
p101 shifti ngot,6,btemp
  pack lenz,b,b
  addx b,btemp,btemp
  vtov [b],btemp,
  subx nrem,lenz,nrem
  addx ngot,lenz,ngot
  ibxgt,rel ,nrem,p102,,
*
*            long vector multiply generates x-matrix
*
*100   ix=0
*      ib=0
*      do 1 loop=1,nrow
*      call vvtv(ncli,x(ix+1),a,b(ib+1))
*      ix=ix+ncli
*1     ib=ib+link
*
  shifti ncli,6,ncli64
  shifti link,6,link64
  rtor r1b,x
  es loop,1
  rtor r1b,xx
  mpyx ncli64,nrow,lenx
  pack ncli,a,a
  pack ncli,b,b
  pack ncli,x,x
  addx r1b,lenx,y
  addx r1b,lenx,yy
p1 mpysv [a],[b],x,
  addx x,ncli64,x
  addx b,link64,b
  ibxle,rel r16,loop,p1,nrow,loop
*
*            transpose x-matrix to y-matrix
*
*      lenx=ncli*nrow
*      ncnr=ncol*nrow
*      iy=lenx
*      ix=0
*      if(ncnr.ge.link)goto 2222
*      do 3 loop=1,ncnr
*      call scattr(link,x(iy+1),ncnr,x(ix+1))
*      ix=ix+link
*3     iy=iy+1
*      goto 4444
*2222  do 2 loop=1,link
*      call gathr(ncnr,x(iy+1),x(ix+1),link)
*      ix=ix+1
*2     iy=iy+ncnr
*
  mpyx ncol,nrow,ncnr
  rtor xx,x
  es loop,1
  shifti ncnr,6,ncnr64
  ibxge,rel ,ncnr,p2222,link,
  pack link,ncnr,stride
p3 vtovx,fia [stride],x,y
  addx x,link64,x
  is y,64
  ibxle,rel r16,loop,p3,ncnr,loop
  bim ,p4444
p2222 pack ncnr,link,stride
p2 vxtov,fia [stride],x,y
  addx y,ncnr64,y
  is x,64
  ibxle,rel r16,loop,p2,link,loop
*
*            condense y-matrix by bi-section
*
*4444  mink=link
*4     if(mink.eq.1)goto 5555
*      nink=mink/2
*      lenh=ncnr*nink
*      mink=mink-nink
*      lenz=mink*ncnr
*      call addvec(x(lenx+1),x(lenx+1),x(lenx+lenz+1),lenh)
*      goto 4
*
p4444 rtor yy,y
  rtor link,mink
p4 ibxeq,rel ,mink,p5555,r16,
  shifti mink,-1,nink
  mpyx ncnr,nink,lenh
  subx mink,nink,mink
  mpyx mink,ncnr64,lenz
  pack lenh,y,y
  addx y,lenz,lenz
  addnv [y],[lenz],y,
  bim ,p4
*
*            scatter y-matrix to r-matrix
*
*5555  if(nrow.ne.1)goto 6666
*      if(mcolr.ne.1)goto 555
*      call fmove(x(lenx+1),r,ncol)
*      goto 9999
*555   call scattr(ncol,r,mcolr,x(lenx+1))
*      goto 9999
*
p5555 ibxne,rel ,nrow,p6666,r16,
 ibxne,rel ,mcolr,p555,r16,
  pack ncol,yy,yy
  pack ncol,r,r
  vtov [yy],r,
  barb,br ,r1a
p555 pack ncol,mcolr,mcolr
  vtovx,fia [mcolr],yy,r
  barb,br ,r1a
*
*            diagonals of index vector
*
*6666  mcmr=mcolr+mrowr
*      call ibasgn(nrow,1,mcmr,index)
*      inc=mcolr*nrow
*      lenh=nrow
*6     loop=ncol-lenh
*      if(loop.eq.0)goto 66
*      loop=min0(loop,lenh)
*      call ujau(loop,inc,index(lenh+1),index)
*      inc=inc+inc
*      lenh=lenh+loop
*      goto 6
*
p6666 addx mcolr,mrowr,mcmr
  es loop,0
  pack nrow,xx,xx
  interval loop,mcmr,xx,
  mpyx mcolr,nrow,inc
  rtor nrow,lenh
p6 ibxeq,rel ,lenh,p66,ncol,
  subx ncol,lenh,loop
  ibxlt,rel ,loop,p6a,lenh,
  rtor lenh,loop
p6a pack loop,xx,xx
  shifti lenh,6,xxlenh
  addx xx,xxlenh,xxlenh
  addxv,b [xx],inc,xxlenh,
  addx lenh,loop,lenh
  addx inc,inc,inc
  bim ,p6
*
*            uncorrected off-diagonals of index vector
*
*66    inc=mrowr
*666   loop=ncnr-lenh
*      if(loop.eq.0)goto 7777
*      loop=min0(loop,lenh)
*      call ujau(loop,inc,index(lenh+1),index)
*      inc=inc+inc
*      lenh=lenh+loop
*      goto 666
*
p66 rtor mrowr,inc
p666 ibxeq,rel ,lenh,p7777,ncnr,
  subx ncnr,lenh,loop
  ibxlt,rel ,loop,p6b,lenh,
  rtor lenh,loop
p6b pack loop,xx,xx
  shifti lenh,6,xxlenh
  addx xx,xxlenh,xxlenh
  addxv,b [xx],inc,xxlenh,
  addx lenh,loop,lenh
  addx inc,inc,inc
  bim ,p666
*
*            correct index vector
*
*7777  limit=ncol*mcolr
*      inc=nrow*mrowr
*      m=ncol+nrow-1
*      n=ncnr-m
*      if(limit.le.mrowr)goto 8888
*
p7777 mpyx ncol,mcolr,limit
  mpyx nrow,mrowr,inc
  addx ncol,nrow,m
  is m,-1
  subx ncnr,m,n
  shifti m,6,m
  ex mask,vmask
  addx m,xx,m
  pack n,m,m
  ibxle,rel ,limit,p8888,mrowr,
*
*            awkward index correction case (ncol*mcolr .gt. mrowr)
*
  pack ncol,mask,tt
  pack nrow,,ss
  subx nrow,r16,loop
  es uu,-1
p7 pack loop,,rr
  maskz rr,ss,tt
  addx tt,ncol,tt
  ibxne,rel uu,loop,p7,,loop
  addx mask,nrow,mask
  is mask,-1
  bim ,p88
*
*            simple index correction case (ncol*mcolr .le. mrowr)
*
*8888  limit=limit+inc-mcmr+1
*      do 8 loop=1,n
*      if(index(m+loop).gt.limit)index(m+loop)=index(m+loop)-inc
*8     continue
*88    call scatter(ncnr,r,index,x(lenx+1))
*
p8888 addx limit,inc,limit
  subx limit,mcmr,limit
  is limit,1
  pack n,mask,mask
  cmpge,b [m],limit,mask
p88 subxv,b [m],inc,m,mask
  pack ncnr,xx,xx
  vtovx [xx],yy,r
*
*9999  return
*      end
*
  barb,br ,r1a
  end
mymi ident
*
*    matrix multiply by inner product algorithm --- r = a * b
*
*      subroutine mymi(a,mcola,b,mrowb,r,mcolr,mrowr,
*     1 ncol,link,nrow)
*      dimension a(1),b(1),r(1)
*
*      a(ncol,link) has col.,row spacing of mcola,1
*      b(link,nrow) has col.,row spacing of 1,mrowb
*      r(ncol,nrow) has col.,row spacing of mcolr,mrowr
*
*    coded    v.r. saunders ---- sep 1985.
*
  msec 2
a equ 3*64
mcola equ 4*64
b equ 5*64
mrowb equ 6*64
r equ 7*64
mcolr equ 8*64
mrowr equ 9*64
inc equ mrowr
ncol equ #a*64
link equ #b*64
i equ link
nrow equ #c*64
save equ #d*64
dotter equ #e*64
ir equ #f*64
q1 equ #10*64
ia equ #11*64
q50 equ #12*64
qot2 equ #13*64
*
r14 equ #14*64
r16 equ #16*64
r17 equ #17*64
qstore equ r17
r18 equ #18*64
r19 equ #19*64
r1a equ #1a*64
r1b equ #1b*64
*
r20 equ #20*64
r21 equ #21*64
bstore equ #22*64
it equ #23*64
*
plus14 equ 14
nsave equ plus14+plus14+4
*
s1 equ #24*64
s2 equ s1+64
s3 equ s2+64
s4 equ s3+64
s5 equ s4+64
s6 equ s5+64
s7 equ s6+64
s8 equ s7+64
s9 equ s8+64
s10 equ s9+64
s11 equ s10+64
s12 equ s11+64
s13 equ s12+64
scalar equ s13+64
*
m1 equ scalar+64
m2 equ m1+64
m3 equ m2+64
m4 equ m3+64
m5 equ m4+64
m6 equ m5+64
m7 equ m6+64
m8 equ m7+64
m9 equ m8+64
m10 equ m9+64
m11 equ m10+64
m12 equ m11+64
m13 equ m12+64
irr equ m13+64
*
*           fetch arguments
*
  entry mymi
mymi es link,8
  lod [r17,link],link
  es ncol,7
  lod [r17,ncol],ncol
  es mcolr,5
  lod [r17,mcolr],mcolr
  lod [r17,r16],mcola
  es mrowb,3
  lod [r17,mrowb],mrowb
  es mrowr,6
  lod [r17,mrowr],mrowr
  es nrow,9
  lod [r17,nrow],nrow
  lod [r17],a
  es b,2
  lod [r17,b],b
  es r,4
  lod [r17,r],r
  lod [link],link
  lod [ncol],ncol
  lod [mcolr],mcolr
  lod [mcola],mcola
  lod [mrowb],mrowb
  lod [mrowr],mrowr
  lod [nrow],nrow
  rtor r1b,save
  elen save,nsave
  swap ,r14,save
*
*      inc=mrowr-mcolr*ncol
*      ir=0
*      ib=0
*      do 1 j=1,nrow
*      ia=0
*      do 2 i=1,ncol
*      r(ir+1)=vecsum(a(ia+1),b(ib+1),link)
*      ir=ir+mcolr
*2     ia=ia+mcola
*      ib=ib+mrowb
*1     ir=ir+inc
*
  pack link,a,a
  pack link,b,b
  dotv [a],[b],r18,
  mpyx mcolr,ncol,ia
  ex qstore,astore
  ex qot2,dot2
  ex q50,p50
  ex q1,p1
  ex dotter,dot2
  es ir,0
  es it,plus14
  subx mrowr,ia,inc
  rtor qstore,bstore
  rtor a,ia
  rtor ncol,i
  shifti mcola,6,mcola
  shifti mrowb,6,mrowb
  barb,br ,q50
p1 rtor a,ia
  rtor ncol,i
  barb,br ,dotter
dot2 dotv [ia],[b],r20,
  addn ,r18,scalar
  bsave dotter,[bstore]
  dotv [ia],[b],r18,
  addn ,r20,scalar
  rtor qot2,dotter
  barb,br ,bstore
astore rtor irr,m1
  addn scalar,r19,s1
  bsave bstore,[q50]
  rtor irr,m2
  addn scalar,r21,s2
  bsave bstore,[q50]
  rtor irr,m3
  addn scalar,r19,s3
  bsave bstore,[q50]
  rtor irr,m4
  addn scalar,r21,s4
  bsave bstore,[q50]
  rtor irr,m5
  addn scalar,r19,s5
  bsave bstore,[q50]
  rtor irr,m6
  addn scalar,r21,s6
  bsave bstore,[q50]
  rtor irr,m7
  addn scalar,r19,s7
  bsave bstore,[q50]
  rtor irr,m8
  addn scalar,r21,s8
  bsave bstore,[q50]
  rtor irr,m9
  addn scalar,r19,s9
  bsave bstore,[q50]
  rtor irr,m10
  addn scalar,r21,s10
  bsave bstore,[q50]
  rtor irr,m11
  addn scalar,r19,s11
  bsave bstore,[q50]
  rtor irr,m12
  addn scalar,r21,s12
  bsave bstore,[q50]
  rtor irr,m13
  addn scalar,r19,s13
  bsave bstore,[q50]
  addn scalar,r21,scalar
  es it,plus14
  rtor qstore,bstore
  sto [r,m1],s1
  sto [r,m2],s2
  sto [r,m3],s3
  sto [r,m4],s4
  sto [r,m5],s5
  sto [r,m6],s6
  sto [r,m7],s7
  sto [r,m8],s8
  sto [r,m9],s9
  sto [r,m10],s10
  sto [r,m11],s11
  sto [r,m12],s12
  sto [r,m13],s13
  sto [r,irr],scalar
p50 rtor ir,irr
  addx ir,mcolr,ir
  addx ia,mcola,ia
  subx it,r16,it
  dbnz i,[dotter]
  addx b,mrowb,b
  addx ir,inc,ir
  dbnz nrow,[q1]
  ex qstore,cstore
  ibxeq,rel ,dotter,rot2,qot2,
  addn ,r20,scalar
  addn scalar,r21,scalar
  ibnz ,[qstore,it]
rot2 addn ,r18,scalar
  addn scalar,r19,scalar
  ibnz ,[qstore,it]
cstore sto [r,m13],s13
  sto [r,m12],s12
  sto [r,m11],s11
  sto [r,m10],s10
  sto [r,m9],s9
  sto [r,m8],s8
  sto [r,m7],s7
  sto [r,m6],s6
  sto [r,m5],s5
  sto [r,m4],s4
  sto [r,m3],s3
  sto [r,m2],s2
  sto [r,m1],s1
  sto [r,irr],scalar
  swap save,r14,
  barb,br ,r1a
  end
mymo ident
*
*    matrix multiply by outer product algorithm --- r = a * b
*
*      subroutine mymo(a,mrowa,b,mcolb,mrowb,r,mrowr,
*     1 ncol,link,nrow)
*      dimension a(1),b(1),r(1)
*
*      a(ncol,link) has col.,row spacing of 1,mrowa
*      b(link,nrow) has col.,row spacing of mcolb,mrowb
*      r(ncol,nrow) has col.,row spacing of 1,mrowr
*
*    coded    v.r. saunders ---- aug 1985.
*
  msec 2
a equ 3*64
mrowa equ 4*64
b equ 5*64
mcolb equ 6*64
mrowb equ 7*64
temp equ mrowb
r equ 8*64
qload equ 9*64
mrowr equ #a*64
ncol equ #b*64
link equ #c*64
nrow equ #d*64
vrs32n equ #e*64
nsize equ #f*64
inc equ #10*64
mcolb64 equ #11*64
qmin equ #12*64
bgath equ #13*64
*
r14 equ #14*64
r16 equ #16*64
r17 equ #17*64
q60 equ r17
r18 equ #18*64
ia equ #19*64
r1a equ #1a*64
r1b equ #1b*64
*
k equ #20*64
q1 equ #21*64
zgath equ #22*64
q50 equ #23*64
skip equ #24*64
qskip equ #25*64
sskip equ #26*64
wgath equ #27*64
*
plus17 equ 17
nsave equ plus17+plus17+6
*
scalar equ #28*64
s2 equ scalar+64
s3 equ s2+64
s4 equ s3+64
s5 equ s4+64
s6 equ s5+64
s7 equ s6+64
s8 equ s7+64
s9 equ s8+64
s10 equ s9+64
s11 equ s10+64
s12 equ s11+64
s13 equ s12+64
s14 equ s13+64
s15 equ s14+64
s16 equ s15+64
s17 equ s16+64
*
m3 equ plus17*64+scalar
m4 equ m3+64
m5 equ m4+64
m6 equ m5+64
m7 equ m6+64
m8 equ m7+64
m9 equ m8+64
m10 equ m9+64
m11 equ m10+64
m12 equ m11+64
m13 equ m12+64
m14 equ m13+64
m15 equ m14+64
m16 equ m15+64
m17 equ m16+64
*
*           fetch arguments
*
  entry mymo
mymo es mcolb,3
  lod [r17,mcolb],mcolb
  es ncol,7
  lod [r17,ncol],ncol
  es nrow,9
  lod [r17,nrow],nrow
  es link,8
  lod [r17,link],link
  es mrowr,6
  lod [r17,mrowr],mrowr
  es mrowb,4
  lod [r17,mrowb],mrowb
  lod [r17,r16],mrowa
  es r,5
  lod [r17,r],r
  es b,2
  lod [r17,b],b
  lod [r17],a
  lod [mcolb],mcolb
  lod [ncol],ncol
  lod [nrow],nrow
  lod [link],link
  lod [mrowr],mrowr
  lod [mrowb],mrowb
  lod [mrowa],mrowa
  rtor r1b,r18
  elen r18,nsave
  swap ,r14,r18
*
*            calc. indices for b-matrix fetch
*
  addx mcolb,mcolb,m3
  addx m3,mcolb,m4
  addx m4,mcolb,m5
  addx m5,mcolb,m6
  addx m6,mcolb,m7
  addx m7,mcolb,m8
  addx m8,mcolb,m9
  addx m9,mcolb,m10
  addx m10,mcolb,m11
  addx m11,mcolb,m12
  addx m12,mcolb,m13
  addx m13,mcolb,m14
  addx m14,mcolb,m15
  addx m15,mcolb,m16
  addx m16,mcolb,m17
*
*      inc=mrowb-mcolb*link
*      nsize=link
*      if(inc.eq.0)nsize=link*nrow
*      ir=0
*      ib=0
*      ix=0
*      do 1 j=1,nrow
*      ia=0
*      do 2 k=1,link
*      if(ix.ne.0)goto 50
*      ix=min0(32,nsize)
*      nsize=nsize-ix
*      if(nsize.le.0)nsize=link
*      call gathr(ix,reg,b(ib+1),mcolb)
*      it=1
*50    scalar=reg(it)
*      it=it+1
*      if(scalar.eq.(0.0))goto 60
*      call triad(ncol,scalar,r(ir+1),a(ia+1))
*60    ib=ib+mcolb
*      ix=ix-1
*2     ia=ia+mrowa
*      ib=ib+inc
*1     ir=ir+mrowr
*
  shifti mcolb,6,mcolb64
  shifti mrowb,6,mrowb
  shifti mrowa,6,mrowa
  shifti mrowr,6,mrowr
  mpyx mcolb64,link,inc
  mpyx link,nrow,nsize
  es vrs32n,-plus17
  ex wgath,pgath
  ex qmin,pmin
  ex q1,p1
  ex q50,p50
  ex q60,p60
  subx mrowb,inc,inc
  ex qload,pload
  ex zgath,ygath
  ex qskip,pskip
  ex sskip,rskip
  pack ncol,r,r
  ibxeq,rel ,inc,p13,,
  subx b,inc,b
  ex q1,p11
p11 rtor link,nsize
  addx b,inc,b
p13 rtor wgath,bgath
  es temp,0
p1 rtor link,k
  es skip,1
  pack ncol,a,ia
  barb,br ,bgath
ygath rtor s2,scalar
  bsave bgath,[q50]
  rtor s3,scalar
  bsave bgath,[q50]
  rtor s4,scalar
  bsave bgath,[q50]
  rtor s5,scalar
  bsave bgath,[q50]
  rtor s6,scalar
  bsave bgath,[q50]
  rtor s7,scalar
  bsave bgath,[q50]
  rtor s8,scalar
  bsave bgath,[q50]
  rtor s9,scalar
  bsave bgath,[q50]
  rtor s10,scalar
  bsave bgath,[q50]
  rtor s11,scalar
  bsave bgath,[q50]
  rtor s12,scalar
  bsave bgath,[q50]
  rtor s13,scalar
  bsave bgath,[q50]
  rtor s14,scalar
  bsave bgath,[q50]
  rtor s15,scalar
  bsave bgath,[q50]
  rtor s16,scalar
  bsave bgath,[q50]
  rtor s17,scalar
  bsave bgath,[q50]
pgath rtor zgath,bgath
  ibxgt vrs32n,nsize,[qmin],,nsize
  subx ,nsize,temp
pmin lod [b],scalar
  ibnz ,[qload,temp]
pload lod [b,m17],s17
  lod [b,m16],s16
  lod [b,m15],s15
  lod [b,m14],s14
  lod [b,m13],s13
  lod [b,m12],s12
  lod [b,m11],s11
  lod [b,m10],s10
  lod [b,m9],s9
  lod [b,m8],s8
  lod [b,m7],s7
  lod [b,m6],s6
  lod [b,m5],s5
  lod [b,m4],s4
  lod [b,m3],s3
  lod [b,mcolb],s2
p50 ibxeq ,scalar,[q60],,
  dbnz skip,[qskip]
  mpysv,a scalar,[ia],r,
  barb,br ,q60
pskip linkv,rb
  mpysv,a scalar,[ia],ia,
  addnv [r],[ia],r,
p60 addx b,mcolb64,b
  addx ia,mrowa,ia
  dbnz k,[bgath]
  dbnz skip,[sskip]
  vtov,a ,r,
rskip addx r,mrowr,r
  dbnz nrow,[q1]
  swap r18,r14,
  barb,br ,r1a
  end
mxmo ident
*
*    matrix multiply by outer product algorithm --- r = r + a * b
*
*      subroutine mxmo(a,mrowa,b,mcolb,mrowb,r,mrowr,
*     1 ncol,link,nrow)
*      dimension a(1),b(1),r(1)
*
*      a(ncol,link) has col.,row spacing of 1,mrowa
*      b(link,nrow) has col.,row spacing of mcolb,mrowb
*      r(ncol,nrow) has col.,row spacing of 1,mrowr
*
*    coded    v.r. saunders ---- aug 1985.
*
  msec 2
a equ 3*64
mrowa equ 4*64
b equ 5*64
mcolb equ 6*64
mrowb equ 7*64
temp equ mrowb
r equ 8*64
qload equ 9*64
mrowr equ #a*64
ncol equ #b*64
link equ #c*64
nrow equ #d*64
vrs32n equ #e*64
nsize equ #f*64
inc equ #10*64
mcolb64 equ #11*64
qmin equ #12*64
bgath equ #13*64
*
r14 equ #14*64
r16 equ #16*64
r17 equ #17*64
q60 equ r17
r18 equ #18*64
ia equ #19*64
r1a equ #1a*64
r1b equ #1b*64
*
k equ #20*64
q1 equ #21*64
zgath equ #22*64
q50 equ #23*64
*
plus19 equ 19
nsave equ plus19+plus19+2
*
scalar equ #24*64
s2 equ scalar+64
s3 equ s2+64
s4 equ s3+64
s5 equ s4+64
s6 equ s5+64
s7 equ s6+64
s8 equ s7+64
s9 equ s8+64
s10 equ s9+64
s11 equ s10+64
s12 equ s11+64
s13 equ s12+64
s14 equ s13+64
s15 equ s14+64
s16 equ s15+64
s17 equ s16+64
s18 equ s17+64
s19 equ s18+64
*
m3 equ plus19*64+scalar
m4 equ m3+64
m5 equ m4+64
m6 equ m5+64
m7 equ m6+64
m8 equ m7+64
m9 equ m8+64
m10 equ m9+64
m11 equ m10+64
m12 equ m11+64
m13 equ m12+64
m14 equ m13+64
m15 equ m14+64
m16 equ m15+64
m17 equ m16+64
m18 equ m17+64
m19 equ m18+64
*
*           fetch arguments
*
  entry mxmo
mxmo es mcolb,3
  lod [r17,mcolb],mcolb
  es ncol,7
  lod [r17,ncol],ncol
  es nrow,9
  lod [r17,nrow],nrow
  es link,8
  lod [r17,link],link
  es mrowr,6
  lod [r17,mrowr],mrowr
  es mrowb,4
  lod [r17,mrowb],mrowb
  lod [r17,r16],mrowa
  es r,5
  lod [r17,r],r
  es b,2
  lod [r17,b],b
  lod [r17],a
  lod [mcolb],mcolb
  lod [ncol],ncol
  lod [nrow],nrow
  lod [link],link
  lod [mrowr],mrowr
  lod [mrowb],mrowb
  lod [mrowa],mrowa
  rtor r1b,r18
  elen r18,nsave
  swap ,r14,r18
*
*            calc. indices for b-matrix fetch
*
  addx mcolb,mcolb,m3
  addx m3,mcolb,m4
  addx m4,mcolb,m5
  addx m5,mcolb,m6
  addx m6,mcolb,m7
  addx m7,mcolb,m8
  addx m8,mcolb,m9
  addx m9,mcolb,m10
  addx m10,mcolb,m11
  addx m11,mcolb,m12
  addx m12,mcolb,m13
  addx m13,mcolb,m14
  addx m14,mcolb,m15
  addx m15,mcolb,m16
  addx m16,mcolb,m17
  addx m17,mcolb,m18
  addx m18,mcolb,m19
*
*      inc=mrowb-mcolb*link
*      nsize=link
*      if(inc.eq.0)nsize=link*nrow
*      ir=0
*      ib=0
*      ix=0
*      do 1 j=1,nrow
*      ia=0
*      do 2 k=1,link
*      if(ix.ne.0)goto 50
*      ix=min0(32,nsize)
*      nsize=nsize-ix
*      if(nsize.le.0)nsize=link
*      call gathr(ix,reg,b(ib+1),mcolb)
*      it=1
*50    scalar=reg(it)
*      it=it+1
*      if(scalar.eq.(0.0))goto 60
*      call triad(ncol,scalar,r(ir+1),a(ia+1))
*60    ib=ib+mcolb
*      ix=ix-1
*2     ia=ia+mrowa
*      ib=ib+inc
*1     ir=ir+mrowr
*
  shifti mcolb,6,mcolb64
  shifti mrowb,6,mrowb
  shifti mrowa,6,mrowa
  shifti mrowr,6,mrowr
  mpyx mcolb64,link,inc
  mpyx link,nrow,nsize
  es vrs32n,-plus19
  ex qmin,pmin
  ex q1,p1
  ex q50,p50
  ex q60,p60
  subx mrowb,inc,inc
  ex qload,pload
  ex zgath,ygath
  pack ncol,r,r
  ibxeq,rel ,inc,p13,,
  subx b,inc,b
  ex q1,p11
p11 rtor link,nsize
  addx b,inc,b
p13 ex bgath,pgath
  es temp,0
p1 rtor link,k
  pack ncol,a,ia
  barb,br ,bgath
ygath rtor s2,scalar
  bsave bgath,[q50]
  rtor s3,scalar
  bsave bgath,[q50]
  rtor s4,scalar
  bsave bgath,[q50]
  rtor s5,scalar
  bsave bgath,[q50]
  rtor s6,scalar
  bsave bgath,[q50]
  rtor s7,scalar
  bsave bgath,[q50]
  rtor s8,scalar
  bsave bgath,[q50]
  rtor s9,scalar
  bsave bgath,[q50]
  rtor s10,scalar
  bsave bgath,[q50]
  rtor s11,scalar
  bsave bgath,[q50]
  rtor s12,scalar
  bsave bgath,[q50]
  rtor s13,scalar
  bsave bgath,[q50]
  rtor s14,scalar
  bsave bgath,[q50]
  rtor s15,scalar
  bsave bgath,[q50]
  rtor s16,scalar
  bsave bgath,[q50]
  rtor s17,scalar
  bsave bgath,[q50]
  rtor s18,scalar
  bsave bgath,[q50]
  rtor s19,scalar
  bsave bgath,[q50]
pgath rtor zgath,bgath
  ibxgt vrs32n,nsize,[qmin],,nsize
  subx ,nsize,temp
pmin lod [b],scalar
  ibnz ,[qload,temp]
pload lod [b,m19],s19
  lod [b,m18],s18
  lod [b,m17],s17
  lod [b,m16],s16
  lod [b,m15],s15
  lod [b,m14],s14
  lod [b,m13],s13
  lod [b,m12],s12
  lod [b,m11],s11
  lod [b,m10],s10
  lod [b,m9],s9
  lod [b,m8],s8
  lod [b,m7],s7
  lod [b,m6],s6
  lod [b,m5],s5
  lod [b,m4],s4
  lod [b,m3],s3
  lod [b,mcolb],s2
p50 ibxeq ,scalar,[q60],,
  linkv,rb
  mpysv,a scalar,[ia],ia,
  addnv [r],[ia],r,
p60 addx b,mcolb64,b
  addx ia,mrowa,ia
  dbnz k,[bgath]
  addx r,mrowr,r
  dbnz nrow,[q1]
  swap r18,r14,
  barb,br ,r1a
  end
mxmbn ident
*  call mxmbn(a,mcola,mrowa,b,mcolb,mrowb,r,mcolr,mrowr,ncol,nlink,nrow)
*  matrix multiply routine
*  r(ncol,nrow) = - a(ncol,nlink) * b(nlink,nrow) + r(ncol,nrow)
*  r   ****must****   be pre-initialized
*  r,a,b  stored with row elements spaced at mrowr,mrowa,mrowb
*  r,a,b  stored with column elements spaced at mcolr,mcola,mcolb
*
*    coded    v.r.saunders ---- jun 84.
*
scr205 msec 3
temp res,128 64
  msec 2
r3 equ 3*64
r4 equ 4*64
r5 equ 5*64
r6 equ 6*64
r7 equ 7*64
r8 equ 8*64
r9 equ 9*64
ra equ #a*64
rb equ #b*64
rc equ #c*64
rd equ #d*64
re equ #e*64
rf equ #f*64
r10 equ #10*64
r11 equ #11*64
r12 equ #12*64
r13 equ #13*64
r15 equ #15*64
r16 equ #16*64
r17 equ #17*64
r1a equ #1a*64
r1b equ #1b*64
r1c equ #1c*64
r1d equ #1d*64
r1e equ #1e*64
r1f equ #1f*64
r20 equ #20*64
r21 equ #21*64
r22 equ #22*64
r23 equ #23*64
r24 equ #24*64
r25 equ #25*64
r26 equ #26*64
r27 equ #27*64
r28 equ #28*64
r29 equ #29*64
r2a equ #2a*64
r2b equ #2b*64
r2c equ #2c*64
r2d equ #2d*64
r2e equ #2e*64
r2f equ #2f*64
con32 equ 28
scal equ #2f
scal1 equ (scal+1)*64
scal2 equ (scal+2)*64
scal3 equ (scal+3)*64
scal4 equ (scal+4)*64
scal5 equ (scal+5)*64
scal6 equ (scal+6)*64
scal7 equ (scal+7)*64
scal8 equ (scal+8)*64
scal9 equ (scal+9)*64
scal10 equ (scal+10)*64
scal11 equ (scal+11)*64
scal12 equ (scal+12)*64
scal13 equ (scal+13)*64
scal14 equ (scal+14)*64
scal15 equ (scal+15)*64
scal16 equ (scal+16)*64
scal17 equ (scal+17)*64
scal18 equ (scal+18)*64
scal19 equ (scal+19)*64
scal20 equ (scal+20)*64
scal21 equ (scal+21)*64
scal22 equ (scal+22)*64
scal23 equ (scal+23)*64
scal24 equ (scal+24)*64
scal25 equ (scal+25)*64
scal26 equ (scal+26)*64
scal27 equ (scal+27)*64
scal28 equ (scal+28)*64
rcal equ scal+con32-2
rcal3 equ (rcal+3)*64
rcal4 equ (rcal+4)*64
rcal5 equ (rcal+5)*64
rcal6 equ (rcal+6)*64
rcal7 equ (rcal+7)*64
rcal8 equ (rcal+8)*64
rcal9 equ (rcal+9)*64
rcal10 equ (rcal+10)*64
rcal11 equ (rcal+11)*64
rcal12 equ (rcal+12)*64
rcal13 equ (rcal+13)*64
rcal14 equ (rcal+14)*64
rcal15 equ (rcal+15)*64
rcal16 equ (rcal+16)*64
rcal17 equ (rcal+17)*64
rcal18 equ (rcal+18)*64
rcal19 equ (rcal+19)*64
rcal20 equ (rcal+20)*64
rcal21 equ (rcal+21)*64
rcal22 equ (rcal+22)*64
rcal23 equ (rcal+23)*64
rcal24 equ (rcal+24)*64
rcal25 equ (rcal+25)*64
rcal26 equ (rcal+26)*64
rcal27 equ (rcal+27)*64
rcal28 equ (rcal+28)*64
con12 equ 12
rreg equ rcal+con32-2
rreg3 equ (rreg+3)*64
rreg4 equ (rreg+4)*64
rreg5 equ (rreg+5)*64
rreg6 equ (rreg+6)*64
rreg7 equ (rreg+7)*64
rreg8 equ (rreg+8)*64
rreg9 equ (rreg+9)*64
rreg10 equ (rreg+10)*64
rreg11 equ (rreg+11)*64
rreg12 equ (rreg+12)*64
areg equ rreg+con12-2
areg3 equ (areg+3)*64
areg4 equ (areg+4)*64
areg5 equ (areg+5)*64
areg6 equ (areg+6)*64
areg7 equ (areg+7)*64
areg8 equ (areg+8)*64
areg9 equ (areg+9)*64
areg10 equ (areg+10)*64
areg11 equ (areg+11)*64
areg12 equ (areg+12)*64
rre equ areg+con12
rre1 equ (rre+1)*64
rre2 equ (rre+2)*64
rre3 equ (rre+3)*64
rre4 equ (rre+4)*64
rre5 equ (rre+5)*64
rre6 equ (rre+6)*64
rre7 equ (rre+7)*64
rre8 equ (rre+8)*64
rre9 equ (rre+9)*64
rre10 equ (rre+10)*64
rre11 equ (rre+11)*64
rre12 equ (rre+12)*64
are equ rre+con12
are1 equ (are+1)*64
are2 equ (are+2)*64
are3 equ (are+3)*64
are4 equ (are+4)*64
are5 equ (are+5)*64
are6 equ (are+6)*64
are7 equ (are+7)*64
are8 equ (are+8)*64
are9 equ (are+9)*64
are10 equ (are+10)*64
are11 equ (are+11)*64
are12 equ (are+12)*64
*
*  register assignments
*
*   r3   vl/a for  do 2
*   r4   vrsc
*   r5   vrsd   or  vrsdd    or  rstore+(con12-ncol)*32
*   r6   nrow for  do 1
*   r7   vl/r for  do 1 and do 2
*   r8   triav  or  triad+(con12-ncol)*32
*   r9   nlink for do 2
*   ra   nlink or nlink*nrow for do 1 and do 2
*   rb   mcolb
*   rc   #ffff   or   ncol remainder
*   rd   nlink   or   nlink*nrow   permanently
*   re   vrsb
*   rf   link register for  b-fetching  routine
*   r10  vrsa  or  triad+(con12-ncol)*64
*   r11  con12  or  con12-ncol
*   r12  0/b  for do 1 and do 2
*   r13  0/b  permanently
*   r1a  vl/a  for do 65535
*   r1b  mrowa*64  permanently
*   r1d  (mrowb-nlink*mcolb)*64    permanently
*   r1e  0/r  for do 65535
*   r1f  mrowr*64    permanently
*   r20  ncol  for do 65535 loop
*   r21  nlink    permanently
*   r22  nrow     permanently
*   r23  con32
*   r24  vrs5   or   vrs4     // joop for zero argument jump
*   r25  mcolb*64   permanently
*   r26  vl/mcolr
*   r27  mcolr*64
*   r28  temp
*   r29  vl/r  or  vl/temp
*   r2a  vrs44  or  vrs55   rfetch+(con12-ncol)*32
*   r2b  ncol/mcola
*   r2c  mcola*64
*   r2d  vl/a  or  vl/temp  for  do 1
*   r2e  temp or temp+ncol*nlink*64
*   r2f  mrowa*64  or  ncol*64
*   r30  onwards for con32 registers  --  b-fetch area
*   referenced as scal1 , scal2 , scal3     etc
*   r30+con32  onwards for con32-2 mcolb*2,mcolb*3, etc
*   referenced as rcal3 , rcal4 , rcal5     etc
*
  entry mxmbn
mxmbn swap ,r15,r1c
  es areg11,10
  lod [r17,areg11],r21
  es areg6,5
  lod [r17,areg6],r1d
  es areg5,4
  lod [r17,areg5],rb
  ex r2a,vrs44
  es areg12,11
  lod [r17,areg12],r22
  ex r28,temp
  es areg3,2
  lod [r17,areg3],r1b
  ex r8,triav
  es areg9,8
  lod [r17,areg9],r1f
  es areg10,9
  lod [r17,areg10],r20
  lod [r17],r1a
  es areg4,3
  lod [r17,areg4],r13
  es areg7,6
  lod [r17,areg7],r1e
  es areg8,7
  lod [r17,areg8],r26
  lod [r17,r16],r2b
  lod [r21],r21
  lod [r1d],ra
  lod [rb],rb
  lod [r22],r22
  lod [r1b],r1b
  lod [r1f],r1f
  lod [r20],r20
*  **********************************************************
*  return if nlink (r21) or nrow (r22) or ncol(r20)  .le. 0
*   if so jump to label joop , where registers are restored
*     jvl  feb 88
  ex r24,joop
  ibxle ,r21,[r24],,
  ibxle ,r22,[r24],,
  ibxle ,r20,[r24],,
*  *********************************************************
  lod [r26],r26
  lod [r2b],r2b
  es r11,con12
  es r23,con32
  ex rc,#ffff
  ex r4,vrsc
  ex r5,vrsd
  ex re,vrsb
  ex r10,vrsa
  ex r24,vrs4
  ibxne,rel ,r21,vrs88,r16,
  rtor ra,rb
vrs88 mpyx r21,rb,r1d
  mpyx r21,r22,rd
  shifti r1b,6,r1b
  shifti r1f,6,r1f
  ibxeq,rel ,r1d,vrs1,ra,
  shifti ra,6,ra
  shifti r1d,6,r1d
  ex r24,vrs5
  rtor r21,rd
  subx ra,r1d,r1d
vrs1 shifti r2b,6,r2c
  shifti rb,6,r25
  rtor rd,ra
  shifti r26,6,r27
  addx rb,rb,rcal3
  addx rb,rcal3,rcal4
  addx rb,rcal4,rcal5
  addx rb,rcal5,rcal6
  addx rb,rcal6,rcal7
  addx rb,rcal7,rcal8
  addx rb,rcal8,rcal9
  addx rb,rcal9,rcal10
  addx rb,rcal10,rcal11
  addx rb,rcal11,rcal12
  addx rb,rcal12,rcal13
  addx rb,rcal13,rcal14
  addx rb,rcal14,rcal15
  addx rb,rcal15,rcal16
  addx rb,rcal16,rcal17
  addx rb,rcal17,rcal18
  addx rb,rcal18,rcal19
  addx rb,rcal19,rcal20
  addx rb,rcal20,rcal21
  addx rb,rcal21,rcal22
  addx rb,rcal22,rcal23
  addx rb,rcal23,rcal24
  addx rb,rcal24,rcal25
  addx rb,rcal25,rcal26
  addx rb,rcal26,rcal27
  addx rb,rcal27,rcal28
  ibxeq,rel ,r26,top,r16,
  ex r2a,vrs55
  ex r5,vrsdd
*
*    end    of    entry code
*
*
*   begin   do 65535 loop=1,nseg
*
top ibxge,rel ,r20,vrs2,rc,
  rtor r20,rc
  ibxgt,rel ,r20,vrs2,r11,
  subx r11,r20,r11
  shifti r11,5,r5
  ex r2a,rfetch
  ex r8,triad
  addx r26,r26,rreg3
  addx r26,rreg3,rreg4
  addx r26,rreg4,rreg5
  addx r26,rreg5,rreg6
  addx r26,rreg6,rreg7
  addx r26,rreg7,rreg8
  addx r26,rreg8,rreg9
  addx r26,rreg9,rreg10
  addx r26,rreg10,rreg11
  addx r26,rreg11,rreg12
  addx r2a,r5,r2a
  addx r8,r5,r10
  addx r8,r5,r8
  ix r5,rstore
vrs2 rtor r28,r2e
  rtor r1b,r2f
  pack rc,r1a,r2d
  ibxeq,rel ,r2b,vrs22,r16,
*  transpose  a  matrix
  shifti rc,6,r2f
  pack rc,r2b,r2b
  rtor r21,r9
  rtor r28,r3
  ex r12,vrszz
vrszz vxtov,fia [r2b],r2d,r3
  addx r2d,r1b,r2d
  addx r3,r2f,r3
  dbnz r9,[r12]
vrsm mpyx r2f,r21,r2e
  pack rc,r28,r2d
  addx r2e,r28,r2e
vrs22 rtor r13,r12
  rtor re,rf
  rtor r22,r6
  pack rc,r1e,r7
  pack rc,r26,r26
*
*   begin   do 1 i=1,nrow
*
  rtor r7,r29
  barb,br ,r2a
*  gather a column of r to temp
vrs55 vxtov,fia [r26],r7,r2e
  pack rc,r2e,r29
vrs44 rtor r2d,r3
  rtor r21,r9
*
*   begin   do 2 j=1,nlink
*
  barb,br ,rf
*
*  fetch   con32  elements  of  b
*
vrsb subx ra,r23,ra
  rtor r4,rf
  es scal2,0
  ibxgt,rel ,ra,vrs7,,
  subx scal2,ra,scal2
  rtor rd,ra
vrs7 lod [r12],scal1
  subx r23,scal2,are1
  mpyx are1,r25,are1
  bim scal2,vrs77
vrs77 lod [r12,rcal28],scal28
  lod [r12,rcal27],scal27
  lod [r12,rcal26],scal26
  lod [r12,rcal25],scal25
  lod [r12,rcal24],scal24
  lod [r12,rcal23],scal23
  lod [r12,rcal22],scal22
  lod [r12,rcal21],scal21
  lod [r12,rcal20],scal20
  lod [r12,rcal19],scal19
  lod [r12,rcal18],scal18
  lod [r12,rcal17],scal17
  lod [r12,rcal16],scal16
  lod [r12,rcal15],scal15
  lod [r12,rcal14],scal14
  lod [r12,rcal13],scal13
  lod [r12,rcal12],scal12
  lod [r12,rcal11],scal11
  lod [r12,rcal10],scal10
  lod [r12,rcal9],scal9
  lod [r12,rcal8],scal8
  lod [r12,rcal7],scal7
  lod [r12,rcal6],scal6
  lod [r12,rcal5],scal5
  lod [r12,rcal4],scal4
  lod [r12,rcal3],scal3
  lod [r12,rb],scal2
  addx r12,are1,r12
vrsa ibxne ,scal1,[r8],,
  addx r3,r2f,r3
  dbnz r9,[rf]
  barb,br ,r5
*
*   end   of do 2 j=1,nlink
*
triav linkv,rb
  mpysv,a scal1,[r3],r3,
  subnv [r29],[r3],r29,
  addx r3,r2f,r3
  dbnz r9,[rf]
  barb,br ,r5
vrsc rtor scal2,scal1
  bsave rf,[r10]
  rtor scal3,scal1
  bsave rf,[r10]
  rtor scal4,scal1
  bsave rf,[r10]
  rtor scal5,scal1
  bsave rf,[r10]
  rtor scal6,scal1
  bsave rf,[r10]
  rtor scal7,scal1
  bsave rf,[r10]
  rtor scal8,scal1
  bsave rf,[r10]
  rtor scal9,scal1
  bsave rf,[r10]
  rtor scal10,scal1
  bsave rf,[r10]
  rtor scal11,scal1
  bsave rf,[r10]
  rtor scal12,scal1
  bsave rf,[r10]
  rtor scal13,scal1
  bsave rf,[r10]
  rtor scal14,scal1
  bsave rf,[r10]
  rtor scal15,scal1
  bsave rf,[r10]
  rtor scal16,scal1
  bsave rf,[r10]
  rtor scal17,scal1
  bsave rf,[r10]
  rtor scal18,scal1
  bsave rf,[r10]
  rtor scal19,scal1
  bsave rf,[r10]
  rtor scal20,scal1
  bsave rf,[r10]
  rtor scal21,scal1
  bsave rf,[r10]
  rtor scal22,scal1
  bsave rf,[r10]
  rtor scal23,scal1
  bsave rf,[r10]
  rtor scal24,scal1
  bsave rf,[r10]
  rtor scal25,scal1
  bsave rf,[r10]
  rtor scal26,scal1
  bsave rf,[r10]
  rtor scal27,scal1
  bsave rf,[r10]
  rtor scal28,scal1
  rtor re,rf
  barb,br ,r10
*
*   end  of   b-fetch   code
*
*  scatter a column of r from temp
vrsdd vtovx,fia [r26],r2e,r7
vrsd addx r7,r1f,r29
  addx r7,r1f,r7
  barb,br ,r24
vrs5 addx r12,r1d,r12
  rtor re,rf
vrs4 dbnz r6,[r2a]
*
*   end  of  do 1 i=1,nrow
*
  mpyx rc,r2c,r3
  mpyx rc,r27,r9
  subx r20,rc,r20
  addx r1a,r3,r1a
  addx r1e,r9,r1e
  ibxne,rel ,r20,top,,
*
*   end   of   do 65535 loop=1,nseg
*
joop  swap r1c,r15,
  barb,br ,r1a
*
*  short vector  r-fetch
*
rfetch lod [r7,rreg12],rre12
  lod [r7,rreg11],rre11
  lod [r7,rreg10],rre10
  lod [r7,rreg9],rre9
  lod [r7,rreg8],rre8
  lod [r7,rreg7],rre7
  lod [r7,rreg6],rre6
  lod [r7,rreg5],rre5
  lod [r7,rreg4],rre4
  lod [r7,rreg3],rre3
  lod [r7,r26],rre2
  lod [r7],rre1
  rtor r2d,r3
  rtor r21,r9
  barb,br ,rf
*
*  short vector triad code
*
triad lod [r3,areg12],are12
  lod [r3,areg11],are11
  lod [r3,areg10],are10
  lod [r3,areg9],are9
  lod [r3,areg8],are8
  lod [r3,areg7],are7
  lod [r3,areg6],are6
  lod [r3,areg5],are5
  lod [r3,areg4],are4
  lod [r3,areg3],are3
  lod [r3,r16],are2
  lod [r3],are1
  bim r11,triae
triae mpys are12,scal1,are12
  mpys are11,scal1,are11
  mpys are10,scal1,are10
  mpys are9,scal1,are9
  mpys are8,scal1,are8
  mpys are7,scal1,are7
  mpys are6,scal1,are6
  mpys are5,scal1,are5
  mpys are4,scal1,are4
  mpys are3,scal1,are3
  mpys are2,scal1,are2
  mpys are1,scal1,are1
  bim r11,triaf
triaf subn rre12,are12,rre12
  subn rre11,are11,rre11
  subn rre10,are10,rre10
  subn rre9,are9,rre9
  subn rre8,are8,rre8
  subn rre7,are7,rre7
  subn rre6,are6,rre6
  subn rre5,are5,rre5
  subn rre4,are4,rre4
  subn rre3,are3,rre3
  subn rre2,are2,rre2
  subn rre1,are1,rre1
  addx r3,r2f,r3
  dbnz r9,[rf]
  barb,br ,r5
*
*  short vector r-store
*
rstore sto [r7,rreg12],rre12
  sto [r7,rreg11],rre11
  sto [r7,rreg10],rre10
  sto [r7,rreg9],rre9
  sto [r7,rreg8],rre8
  sto [r7,rreg7],rre7
  sto [r7,rreg6],rre6
  sto [r7,rreg5],rre5
  sto [r7,rreg4],rre4
  sto [r7,rreg3],rre3
  sto [r7,r26],rre2
  sto [r7],rre1
  addx r7,r1f,r7
  barb,br ,r24
  end
mxmb ident
*  call mxmb(a,mcola,mrowa,b,mcolb,mrowb,r,mcolr,mrowr,ncol,nlink,nrow)
*  matrix multiply routine
*  r(ncol,nrow) = a(ncol,nlink) * b(nlink,nrow) + r(ncol,nrow)
*  r   ****must****   be pre-initialized
*  r,a,b  stored with row elements spaced at mrowr,mrowa,mrowb
*  r,a,b  stored with column elements spaced at mcolr,mcola,mcolb
*
*    coded     v.r.saunders ---- jun 84.
*
scr205 msec 3
temp res,128 64
  msec 2
r3 equ 3*64
r4 equ 4*64
r5 equ 5*64
r6 equ 6*64
r7 equ 7*64
r8 equ 8*64
r9 equ 9*64
ra equ #a*64
rb equ #b*64
rc equ #c*64
rd equ #d*64
re equ #e*64
rf equ #f*64
r10 equ #10*64
r11 equ #11*64
r12 equ #12*64
r13 equ #13*64
r15 equ #15*64
r16 equ #16*64
r17 equ #17*64
r1a equ #1a*64
r1b equ #1b*64
r1c equ #1c*64
r1d equ #1d*64
r1e equ #1e*64
r1f equ #1f*64
r20 equ #20*64
r21 equ #21*64
r22 equ #22*64
r23 equ #23*64
r24 equ #24*64
r25 equ #25*64
r26 equ #26*64
r27 equ #27*64
r28 equ #28*64
r29 equ #29*64
r2a equ #2a*64
r2b equ #2b*64
r2c equ #2c*64
r2d equ #2d*64
r2e equ #2e*64
r2f equ #2f*64
con32 equ 28
scal equ #2f
scal1 equ (scal+1)*64
scal2 equ (scal+2)*64
scal3 equ (scal+3)*64
scal4 equ (scal+4)*64
scal5 equ (scal+5)*64
scal6 equ (scal+6)*64
scal7 equ (scal+7)*64
scal8 equ (scal+8)*64
scal9 equ (scal+9)*64
scal10 equ (scal+10)*64
scal11 equ (scal+11)*64
scal12 equ (scal+12)*64
scal13 equ (scal+13)*64
scal14 equ (scal+14)*64
scal15 equ (scal+15)*64
scal16 equ (scal+16)*64
scal17 equ (scal+17)*64
scal18 equ (scal+18)*64
scal19 equ (scal+19)*64
scal20 equ (scal+20)*64
scal21 equ (scal+21)*64
scal22 equ (scal+22)*64
scal23 equ (scal+23)*64
scal24 equ (scal+24)*64
scal25 equ (scal+25)*64
scal26 equ (scal+26)*64
scal27 equ (scal+27)*64
scal28 equ (scal+28)*64
rcal equ scal+con32-2
rcal3 equ (rcal+3)*64
rcal4 equ (rcal+4)*64
rcal5 equ (rcal+5)*64
rcal6 equ (rcal+6)*64
rcal7 equ (rcal+7)*64
rcal8 equ (rcal+8)*64
rcal9 equ (rcal+9)*64
rcal10 equ (rcal+10)*64
rcal11 equ (rcal+11)*64
rcal12 equ (rcal+12)*64
rcal13 equ (rcal+13)*64
rcal14 equ (rcal+14)*64
rcal15 equ (rcal+15)*64
rcal16 equ (rcal+16)*64
rcal17 equ (rcal+17)*64
rcal18 equ (rcal+18)*64
rcal19 equ (rcal+19)*64
rcal20 equ (rcal+20)*64
rcal21 equ (rcal+21)*64
rcal22 equ (rcal+22)*64
rcal23 equ (rcal+23)*64
rcal24 equ (rcal+24)*64
rcal25 equ (rcal+25)*64
rcal26 equ (rcal+26)*64
rcal27 equ (rcal+27)*64
rcal28 equ (rcal+28)*64
con12 equ 12
rreg equ rcal+con32-2
rreg3 equ (rreg+3)*64
rreg4 equ (rreg+4)*64
rreg5 equ (rreg+5)*64
rreg6 equ (rreg+6)*64
rreg7 equ (rreg+7)*64
rreg8 equ (rreg+8)*64
rreg9 equ (rreg+9)*64
rreg10 equ (rreg+10)*64
rreg11 equ (rreg+11)*64
rreg12 equ (rreg+12)*64
areg equ rreg+con12-2
areg3 equ (areg+3)*64
areg4 equ (areg+4)*64
areg5 equ (areg+5)*64
areg6 equ (areg+6)*64
areg7 equ (areg+7)*64
areg8 equ (areg+8)*64
areg9 equ (areg+9)*64
areg10 equ (areg+10)*64
areg11 equ (areg+11)*64
areg12 equ (areg+12)*64
rre equ areg+con12
rre1 equ (rre+1)*64
rre2 equ (rre+2)*64
rre3 equ (rre+3)*64
rre4 equ (rre+4)*64
rre5 equ (rre+5)*64
rre6 equ (rre+6)*64
rre7 equ (rre+7)*64
rre8 equ (rre+8)*64
rre9 equ (rre+9)*64
rre10 equ (rre+10)*64
rre11 equ (rre+11)*64
rre12 equ (rre+12)*64
are equ rre+con12
are1 equ (are+1)*64
are2 equ (are+2)*64
are3 equ (are+3)*64
are4 equ (are+4)*64
are5 equ (are+5)*64
are6 equ (are+6)*64
are7 equ (are+7)*64
are8 equ (are+8)*64
are9 equ (are+9)*64
are10 equ (are+10)*64
are11 equ (are+11)*64
are12 equ (are+12)*64
*
*  register assignments
*
*   r3   vl/a for  do 2
*   r4   vrsc
*   r5   vrsd   or  vrsdd    or  rstore+(con12-ncol)*32
*   r6   nrow for  do 1
*   r7   vl/r for  do 1 and do 2
*   r8   triav  or  triad+(con12-ncol)*32
*   r9   nlink for do 2
*   ra   nlink or nlink*nrow for do 1 and do 2
*   rb   mcolb
*   rc   #ffff   or   ncol remainder
*   rd   nlink   or   nlink*nrow   permanently
*   re   vrsb
*   rf   link register for  b-fetching  routine
*   r10  vrsa  or  triad+(con12-ncol)*64
*   r11  con12  or  con12-ncol
*   r12  0/b  for do 1 and do 2
*   r13  0/b  permanently
*   r1a  vl/a  for do 65535
*   r1b  mrowa*64  permanently
*   r1d  (mrowb-nlink*mcolb)*64    permanently
*   r1e  0/r  for do 65535
*   r1f  mrowr*64    permanently
*   r20  ncol  for do 65535 loop
*   r21  nlink    permanently
*   r22  nrow     permanently
*   r23  con32
*   r24  vrs5   or   vrs4    // temp.  label joop for argument test
*   r25  mcolb*64   permanently
*   r26  vl/mcolr
*   r27  mcolr*64
*   r28  temp
*   r29  vl/r  or  vl/temp
*   r2a  vrs44  or  vrs55   rfetch+(con12-ncol)*32
*   r2b  ncol/mcola
*   r2c  mcola*64
*   r2d  vl/a  or  vl/temp  for  do 1
*   r2e  temp or temp+ncol*nlink*64
*   r2f  mrowa*64  or  ncol*64
*   r30  onwards for con32 registers  --  b-fetch area
*   referenced as scal1 , scal2 , scal3     etc
*   r30+con32  onwards for con32-2 mcolb*2,mcolb*3, etc
*   referenced as rcal3 , rcal4 , rcal5     etc
*
  entry mxmb
mxmb swap ,r15,r1c
  es areg11,10
  lod [r17,areg11],r21
  es areg6,5
  lod [r17,areg6],r1d
  es areg5,4
  lod [r17,areg5],rb
  ex r2a,vrs44
  es areg12,11
  lod [r17,areg12],r22
  ex r28,temp
  es areg3,2
  lod [r17,areg3],r1b
  ex r8,triav
  es areg9,8
  lod [r17,areg9],r1f
  es areg10,9
  lod [r17,areg10],r20
  lod [r17],r1a
  es areg4,3
  lod [r17,areg4],r13
  es areg7,6
  lod [r17,areg7],r1e
  es areg8,7
  lod [r17,areg8],r26
  lod [r17,r16],r2b
  lod [r21],r21
  lod [r1d],ra
  lod [rb],rb
  lod [r22],r22
  lod [r1b],r1b
  lod [r1f],r1f
  lod [r20],r20
*  **********************************************************
*  return if nlink (r21) or nrow (r22) or ncol(r20)  .le. 0
*   if so jump to label joop , where registers are restored
*     jvl  feb 88
  ex r24,joop
  ibxle ,r21,[r24],,
  ibxle ,r22,[r24],,
  ibxle ,r20,[r24],,
*  *********************************************************
  lod [r26],r26
  lod [r2b],r2b
  es r11,con12
  es r23,con32
  ex rc,#ffff
  ex r4,vrsc
  ex r5,vrsd
  ex re,vrsb
  ex r10,vrsa
  ex r24,vrs4
  ibxne,rel ,r21,vrs88,r16,
  rtor ra,rb
vrs88 mpyx r21,rb,r1d
  mpyx r21,r22,rd
  shifti r1b,6,r1b
  shifti r1f,6,r1f
  ibxeq,rel ,r1d,vrs1,ra,
  shifti ra,6,ra
  shifti r1d,6,r1d
  ex r24,vrs5
  rtor r21,rd
  subx ra,r1d,r1d
vrs1 shifti r2b,6,r2c
  shifti rb,6,r25
  rtor rd,ra
  shifti r26,6,r27
  addx rb,rb,rcal3
  addx rb,rcal3,rcal4
  addx rb,rcal4,rcal5
  addx rb,rcal5,rcal6
  addx rb,rcal6,rcal7
  addx rb,rcal7,rcal8
  addx rb,rcal8,rcal9
  addx rb,rcal9,rcal10
  addx rb,rcal10,rcal11
  addx rb,rcal11,rcal12
  addx rb,rcal12,rcal13
  addx rb,rcal13,rcal14
  addx rb,rcal14,rcal15
  addx rb,rcal15,rcal16
  addx rb,rcal16,rcal17
  addx rb,rcal17,rcal18
  addx rb,rcal18,rcal19
  addx rb,rcal19,rcal20
  addx rb,rcal20,rcal21
  addx rb,rcal21,rcal22
  addx rb,rcal22,rcal23
  addx rb,rcal23,rcal24
  addx rb,rcal24,rcal25
  addx rb,rcal25,rcal26
  addx rb,rcal26,rcal27
  addx rb,rcal27,rcal28
  ibxeq,rel ,r26,top,r16,
  ex r2a,vrs55
  ex r5,vrsdd
*
*    end    of    entry code
*
*
*   begin   do 65535 loop=1,nseg
*
top ibxge,rel ,r20,vrs2,rc,
  rtor r20,rc
  ibxgt,rel ,r20,vrs2,r11,
  subx r11,r20,r11
  shifti r11,5,r5
  ex r2a,rfetch
  ex r8,triad
  addx r26,r26,rreg3
  addx r26,rreg3,rreg4
  addx r26,rreg4,rreg5
  addx r26,rreg5,rreg6
  addx r26,rreg6,rreg7
  addx r26,rreg7,rreg8
  addx r26,rreg8,rreg9
  addx r26,rreg9,rreg10
  addx r26,rreg10,rreg11
  addx r26,rreg11,rreg12
  addx r2a,r5,r2a
  addx r8,r5,r10
  addx r8,r5,r8
  ix r5,rstore
vrs2 rtor r28,r2e
  rtor r1b,r2f
  pack rc,r1a,r2d
  ibxeq,rel ,r2b,vrs22,r16,
*  transpose  a  matrix
  shifti rc,6,r2f
  pack rc,r2b,r2b
  rtor r21,r9
  rtor r28,r3
  ex r12,vrszz
vrszz vxtov,fia [r2b],r2d,r3
  addx r2d,r1b,r2d
  addx r3,r2f,r3
  dbnz r9,[r12]
vrsm mpyx r2f,r21,r2e
  pack rc,r28,r2d
  addx r2e,r28,r2e
vrs22 rtor r13,r12
  rtor re,rf
  rtor r22,r6
  pack rc,r1e,r7
  pack rc,r26,r26
*
*   begin   do 1 i=1,nrow
*
  rtor r7,r29
  barb,br ,r2a
*  gather a column of r to temp
vrs55 vxtov,fia [r26],r7,r2e
  pack rc,r2e,r29
vrs44 rtor r2d,r3
  rtor r21,r9
*
*   begin   do 2 j=1,nlink
*
  barb,br ,rf
*
*  fetch   con32  elements  of  b
*
vrsb subx ra,r23,ra
  rtor r4,rf
  es scal2,0
  ibxgt,rel ,ra,vrs7,,
  subx scal2,ra,scal2
  rtor rd,ra
vrs7 lod [r12],scal1
  subx r23,scal2,are1
  mpyx are1,r25,are1
  bim scal2,vrs77
vrs77 lod [r12,rcal28],scal28
  lod [r12,rcal27],scal27
  lod [r12,rcal26],scal26
  lod [r12,rcal25],scal25
  lod [r12,rcal24],scal24
  lod [r12,rcal23],scal23
  lod [r12,rcal22],scal22
  lod [r12,rcal21],scal21
  lod [r12,rcal20],scal20
  lod [r12,rcal19],scal19
  lod [r12,rcal18],scal18
  lod [r12,rcal17],scal17
  lod [r12,rcal16],scal16
  lod [r12,rcal15],scal15
  lod [r12,rcal14],scal14
  lod [r12,rcal13],scal13
  lod [r12,rcal12],scal12
  lod [r12,rcal11],scal11
  lod [r12,rcal10],scal10
  lod [r12,rcal9],scal9
  lod [r12,rcal8],scal8
  lod [r12,rcal7],scal7
  lod [r12,rcal6],scal6
  lod [r12,rcal5],scal5
  lod [r12,rcal4],scal4
  lod [r12,rcal3],scal3
  lod [r12,rb],scal2
  addx r12,are1,r12
vrsa ibxne ,scal1,[r8],,
  addx r3,r2f,r3
  dbnz r9,[rf]
  barb,br ,r5
*
*   end   of do 2 j=1,nlink
*
triav linkv,rb
  mpysv,a scal1,[r3],r3,
  addnv [r29],[r3],r29,
  addx r3,r2f,r3
  dbnz r9,[rf]
  barb,br ,r5
vrsc rtor scal2,scal1
  bsave rf,[r10]
  rtor scal3,scal1
  bsave rf,[r10]
  rtor scal4,scal1
  bsave rf,[r10]
  rtor scal5,scal1
  bsave rf,[r10]
  rtor scal6,scal1
  bsave rf,[r10]
  rtor scal7,scal1
  bsave rf,[r10]
  rtor scal8,scal1
  bsave rf,[r10]
  rtor scal9,scal1
  bsave rf,[r10]
  rtor scal10,scal1
  bsave rf,[r10]
  rtor scal11,scal1
  bsave rf,[r10]
  rtor scal12,scal1
  bsave rf,[r10]
  rtor scal13,scal1
  bsave rf,[r10]
  rtor scal14,scal1
  bsave rf,[r10]
  rtor scal15,scal1
  bsave rf,[r10]
  rtor scal16,scal1
  bsave rf,[r10]
  rtor scal17,scal1
  bsave rf,[r10]
  rtor scal18,scal1
  bsave rf,[r10]
  rtor scal19,scal1
  bsave rf,[r10]
  rtor scal20,scal1
  bsave rf,[r10]
  rtor scal21,scal1
  bsave rf,[r10]
  rtor scal22,scal1
  bsave rf,[r10]
  rtor scal23,scal1
  bsave rf,[r10]
  rtor scal24,scal1
  bsave rf,[r10]
  rtor scal25,scal1
  bsave rf,[r10]
  rtor scal26,scal1
  bsave rf,[r10]
  rtor scal27,scal1
  bsave rf,[r10]
  rtor scal28,scal1
  rtor re,rf
  barb,br ,r10
*
*   end  of   b-fetch   code
*
*  scatter a column of r from temp
vrsdd vtovx,fia [r26],r2e,r7
vrsd addx r7,r1f,r29
  addx r7,r1f,r7
  barb,br ,r24
vrs5 addx r12,r1d,r12
  rtor re,rf
vrs4 dbnz r6,[r2a]
*
*   end  of  do 1 i=1,nrow
*
  mpyx rc,r2c,r3
  mpyx rc,r27,r9
  subx r20,rc,r20
  addx r1a,r3,r1a
  addx r1e,r9,r1e
  ibxne,rel ,r20,top,,
*
*   end   of   do 65535 loop=1,nseg
*
joop  swap r1c,r15,
  barb,br ,r1a
*
*  short vector  r-fetch
*
rfetch lod [r7,rreg12],rre12
  lod [r7,rreg11],rre11
  lod [r7,rreg10],rre10
  lod [r7,rreg9],rre9
  lod [r7,rreg8],rre8
  lod [r7,rreg7],rre7
  lod [r7,rreg6],rre6
  lod [r7,rreg5],rre5
  lod [r7,rreg4],rre4
  lod [r7,rreg3],rre3
  lod [r7,r26],rre2
  lod [r7],rre1
  rtor r2d,r3
  rtor r21,r9
  barb,br ,rf
*
*  short vector triad code
*
triad lod [r3,areg12],are12
  lod [r3,areg11],are11
  lod [r3,areg10],are10
  lod [r3,areg9],are9
  lod [r3,areg8],are8
  lod [r3,areg7],are7
  lod [r3,areg6],are6
  lod [r3,areg5],are5
  lod [r3,areg4],are4
  lod [r3,areg3],are3
  lod [r3,r16],are2
  lod [r3],are1
  bim r11,triae
triae mpys are12,scal1,are12
  mpys are11,scal1,are11
  mpys are10,scal1,are10
  mpys are9,scal1,are9
  mpys are8,scal1,are8
  mpys are7,scal1,are7
  mpys are6,scal1,are6
  mpys are5,scal1,are5
  mpys are4,scal1,are4
  mpys are3,scal1,are3
  mpys are2,scal1,are2
  mpys are1,scal1,are1
  bim r11,triaf
triaf addn rre12,are12,rre12
  addn rre11,are11,rre11
  addn rre10,are10,rre10
  addn rre9,are9,rre9
  addn rre8,are8,rre8
  addn rre7,are7,rre7
  addn rre6,are6,rre6
  addn rre5,are5,rre5
  addn rre4,are4,rre4
  addn rre3,are3,rre3
  addn rre2,are2,rre2
  addn rre1,are1,rre1
  addx r3,r2f,r3
  dbnz r9,[rf]
  barb,br ,r5
*
*  short vector r-store
*
rstore sto [r7,rreg12],rre12
  sto [r7,rreg11],rre11
  sto [r7,rreg10],rre10
  sto [r7,rreg9],rre9
  sto [r7,rreg8],rre8
  sto [r7,rreg7],rre7
  sto [r7,rreg6],rre6
  sto [r7,rreg5],rre5
  sto [r7,rreg4],rre4
  sto [r7,rreg3],rre3
  sto [r7,r26],rre2
  sto [r7],rre1
  addx r7,r1f,r7
  barb,br ,r24
  end
mxmtr ident
*
*    r(rectangle) = a(triangle) * b(rectangle) + r(rectangle)
*
*    coded    v.r. saunders ---- oct 1985.
*
*      subroutine mxmtr(a,b,r,ncol,nrow)
*      dimension a(1),b(ncol,nrow),r(ncol,nrow)
*
  msec 2
*
plus21 equ 21
*
a equ 3*64
b equ 4*64
r equ 5*64
ncol equ 6*64
nrow equ 7*64
ncol64 equ 8*64
vic equ 9*64
bget equ #a*64
q50 equ #c*64
q60 equ #d*64
temp equ #e*64
ia equ #f*64
moop equ #10*64
*
r15 equ #15*64
r16 equ #16*64
r17 equ #17*64
r1a equ #1a*64
r1c equ #1c*64
*
scalar equ #20*64
s2 equ scalar+64
s3 equ s2+64
s4 equ s3+64
s5 equ s4+64
s6 equ s5+64
s7 equ s6+64
s8 equ s7+64
s9 equ s8+64
s10 equ s9+64
s11 equ s10+64
s12 equ s11+64
s13 equ s12+64
s14 equ s13+64
s15 equ s14+64
s16 equ s15+64
s17 equ s16+64
s18 equ s17+64
s19 equ s18+64
s20 equ s19+64
s21 equ s20+64
m3 equ plus21*64+scalar
m4 equ m3+64
m5 equ m4+64
m6 equ m5+64
m7 equ m6+64
m8 equ m7+64
m9 equ m8+64
m10 equ m9+64
m11 equ m10+64
m12 equ m11+64
m13 equ m12+64
m14 equ m13+64
m15 equ m14+64
m16 equ m15+64
m17 equ m16+64
m18 equ m17+64
m19 equ m18+64
m20 equ m19+64
m21 equ m20+64
  entry mxmtr
*
*    fetch arguments
*
mxmtr swap ,r15,r1c
  es m4,3
  lod [r17,m4],ncol
  es m5,4
  lod [r17,m5],nrow
  lod [r17],a
  lod [r17,r16],b
  es m3,2
  lod [r17,m3],r
  es m6,5
  es m7,6
  es m8,7
  es m9,8
  es m10,9
  es m11,10
  es m12,11
  es m13,12
  es m14,13
  es m15,14
  es m16,15
  es m17,16
  es m18,17
  es m19,18
  es m20,19
  es m21,20
  lod [ncol],ncol
  lod [nrow],nrow
  ex vic,vrs
  ex q50,p50
  ex q60,p60
  shifti ncol,6,ncol64
*
*      ib=1
*      ir=1
*      do 1 loop=1,nrow
*      ia=1
*      do 2 moop=1,ncol
*      call triad(moop,b(ib),r(ir),a(ia))
*      ib=ib+1
*2     ia=ia+moop
*1     ir=ir+ncol
*
  rtor a,ia
  es moop,1
  bim ,alan
vrs rtor s2,scalar
  bsave bget,[q50]
  rtor s3,scalar
  bsave bget,[q50]
  rtor s4,scalar
  bsave bget,[q50]
  rtor s5,scalar
  bsave bget,[q50]
  rtor s6,scalar
  bsave bget,[q50]
  rtor s7,scalar
  bsave bget,[q50]
  rtor s8,scalar
  bsave bget,[q50]
  rtor s9,scalar
  bsave bget,[q50]
  rtor s10,scalar
  bsave bget,[q50]
  rtor s11,scalar
  bsave bget,[q50]
  rtor s12,scalar
  bsave bget,[q50]
  rtor s13,scalar
  bsave bget,[q50]
  rtor s14,scalar
  bsave bget,[q50]
  rtor s15,scalar
  bsave bget,[q50]
  rtor s16,scalar
  bsave bget,[q50]
  rtor s17,scalar
  bsave bget,[q50]
  rtor s18,scalar
  bsave bget,[q50]
  rtor s19,scalar
  bsave bget,[q50]
  rtor s20,scalar
  bsave bget,[q50]
  rtor s21,scalar
  bsave bget,[q50]
alan rtor vic,bget
  lod [b],scalar
  lod [b,r16],s2
  lod [b,m3],s3
  lod [b,m4],s4
  lod [b,m5],s5
  lod [b,m6],s6
  lod [b,m7],s7
  lod [b,m8],s8
  lod [b,m9],s9
  lod [b,m10],s10
  lod [b,m11],s11
  lod [b,m12],s12
  lod [b,m13],s13
  lod [b,m14],s14
  lod [b,m15],s15
  lod [b,m16],s16
  lod [b,m17],s17
  lod [b,m18],s18
  lod [b,m19],s19
  lod [b,m20],s20
  lod [b,m21],s21
p50 ibxeq ,scalar,[q60],,
  pack moop,r,r
  pack moop,ia,ia
  linkv,rb
  mpysv,a scalar,[ia],ia,
  addnv [r],[ia],r,
p60 shifti moop,6,temp
  is b,64
  addx ia,temp,ia
  ibxle r16,moop,[bget],ncol,moop
  addx r,ncol64,r
  rtor a,ia
  es moop,1
  dbnz nrow,[bget]
  swap r1c,r15,
  barb,br ,r1a
  end
mxmm ident
*
*    coded     v.r.saunders ---- jul 84.
*
*  call mxmm(a,mrowa,b,r,mrowr,ncol,nrow)
*  matrix multiply routine
*  r(ncol,nrow) = a(ncol,nrow) * b(nrow*(nrow+1)/2 u triangle)
*  + r(ncol,nrow)
*  -----warning----  ncol  ---must--- be .le. 65535
*  r   ****must****   be pre-initialized
*
*      subroutine mxmm(a,mrowa,b,r,mrowr,ncol,nrow)
*      dimension a(mrowa,1),b(1),r(mrowr,1)
*      m=1
*      do 1 i=1,nrow
*      do 2 j=1,i
*      top=b(m)
*      if(top)4,2,4
*4     do 5 loop=1,ncol
*5     r(loop,i)=r(loop,i)+a(loop,j)*top
*2     m=m+1
*1     continue
*      return
*      end
  msec 2
r3 equ 3*64
r4 equ 4*64
r5 equ 5*64
r6 equ 6*64
r7 equ 7*64
r9 equ 9*64
ra equ #a*64
rb equ #b*64
rc equ #c*64
rd equ #d*64
re equ #e*64
rf equ #f*64
r10 equ #10*64
r12 equ #12*64
r15 equ #15*64
r16 equ #16*64
r17 equ #17*64
r1a equ #1a*64
r1c equ #1c*64
plus22 equ 22
scal equ #1c
scal1 equ (scal+1)*64
scal2 equ (scal+2)*64
scal3 equ (scal+3)*64
scal4 equ (scal+4)*64
scal5 equ (scal+5)*64
scal6 equ (scal+6)*64
scal7 equ (scal+7)*64
scal8 equ (scal+8)*64
scal9 equ (scal+9)*64
scal10 equ (scal+10)*64
scal11 equ (scal+11)*64
scal12 equ (scal+12)*64
scal13 equ (scal+13)*64
scal14 equ (scal+14)*64
scal15 equ (scal+15)*64
scal16 equ (scal+16)*64
scal17 equ (scal+17)*64
scal18 equ (scal+18)*64
scal19 equ (scal+19)*64
scal20 equ (scal+20)*64
scal21 equ (scal+21)*64
scal22 equ (scal+22)*64
rcal equ scal+plus22-2
rcal3 equ (rcal+3)*64
rcal4 equ (rcal+4)*64
rcal5 equ (rcal+5)*64
rcal6 equ (rcal+6)*64
rcal7 equ (rcal+7)*64
rcal8 equ (rcal+8)*64
rcal9 equ (rcal+9)*64
rcal10 equ (rcal+10)*64
rcal11 equ (rcal+11)*64
rcal12 equ (rcal+12)*64
rcal13 equ (rcal+13)*64
rcal14 equ (rcal+14)*64
rcal15 equ (rcal+15)*64
rcal16 equ (rcal+16)*64
rcal17 equ (rcal+17)*64
rcal18 equ (rcal+18)*64
rcal19 equ (rcal+19)*64
rcal20 equ (rcal+20)*64
rcal21 equ (rcal+21)*64
rcal22 equ (rcal+22)*64
  entry mxmm
mxmm swap ,r15,r1c
  lod [r17,r16],rb
  es rcal5,4
  lod [r17,rcal5],ra
  es rcal6,5
  lod [r17,rcal6],r5
  es rcal7,6
  lod [r17,rcal7],r6
  lod [r17],rc
  es rcal3,2
  lod [r17,rcal3],r12
  es rcal4,3
  lod [r17,rcal4],r7
  ex r10,vrsa
  ex r4,vrsc
  ex re,vrs3
  lod [rb],rb
  lod [ra],ra
  lod [r5],r5
  lod [r6],r6
  es rd,1
  es rcal8,7
  es rcal9,8
  es rcal10,9
  es rcal11,10
  es rcal12,11
  es rcal13,12
  es rcal14,13
  es rcal15,14
  es rcal16,15
  es rcal17,16
  es rcal18,17
  es rcal19,18
  es rcal20,19
  es rcal21,20
  es rcal22,21
  shifti rb,6,rb
  shifti ra,6,ra
  pack r5,r7,r7
  pack r5,rc,r3
  rtor rd,r9
  bim ,vrsb
*  fetch plus22 elements of b
vrsc rtor scal2,scal1
  bsave rf,[r10]
  rtor scal3,scal1
  bsave rf,[r10]
  rtor scal4,scal1
  bsave rf,[r10]
  rtor scal5,scal1
  bsave rf,[r10]
  rtor scal6,scal1
  bsave rf,[r10]
  rtor scal7,scal1
  bsave rf,[r10]
  rtor scal8,scal1
  bsave rf,[r10]
  rtor scal9,scal1
  bsave rf,[r10]
  rtor scal10,scal1
  bsave rf,[r10]
  rtor scal11,scal1
  bsave rf,[r10]
  rtor scal12,scal1
  bsave rf,[r10]
  rtor scal13,scal1
  bsave rf,[r10]
  rtor scal14,scal1
  bsave rf,[r10]
  rtor scal15,scal1
  bsave rf,[r10]
  rtor scal16,scal1
  bsave rf,[r10]
  rtor scal17,scal1
  bsave rf,[r10]
  rtor scal18,scal1
  bsave rf,[r10]
  rtor scal19,scal1
  bsave rf,[r10]
  rtor scal20,scal1
  bsave rf,[r10]
  rtor scal21,scal1
  bsave rf,[r10]
  rtor scal22,scal1
  bsave rf,[r10]
vrsb rtor r4,rf
  lod [r12],scal1
  lod [r12,r16],scal2
  lod [r12,rcal3],scal3
  lod [r12,rcal4],scal4
  lod [r12,rcal5],scal5
  lod [r12,rcal6],scal6
  lod [r12,rcal7],scal7
  lod [r12,rcal8],scal8
  lod [r12,rcal9],scal9
  lod [r12,rcal10],scal10
  lod [r12,rcal11],scal11
  lod [r12,rcal12],scal12
  lod [r12,rcal13],scal13
  lod [r12,rcal14],scal14
  lod [r12,rcal15],scal15
  lod [r12,rcal16],scal16
  lod [r12,rcal17],scal17
  lod [r12,rcal18],scal18
  lod [r12,rcal19],scal19
  lod [r12,rcal20],scal20
  lod [r12,rcal21],scal21
  lod [r12,rcal22],scal22
vrsa ibxeq ,scal1,[re],,
  linkv,rb
  mpysv,a scal1,[r3],r3,
  addnv [r7],[r3],r7,
vrs3 addx r3,rb,r3
  is r12,64
  dbnz r9,[rf]
*  end of do 2
  is rd,1
  addx r7,ra,r7
  rtor rd,r9
  pack r5,rc,r3
  dbnz r6,[rf]
*  end of do 1
  swap r1c,r15,
  barb,br ,r1a
  end
mxmd ident
  msec 2
r3 equ 3*64
r4 equ 4*64
r5 equ 5*64
r6 equ 6*64
r7 equ 7*64
r8 equ 8*64
r17 equ #17*64
  extc mxmb
  entry mxmd
*   call mxmd(a,ia,ja,b,ib,jb,r,ir,jr,i,j,k)
*
*            coded   v.r.saunders     oct 1984
*
*      dimension a(1),b(1),c(1)
*      call szero(r,i*k)
*      call mxmb(a,ia,ja,b,ib,jb,r,ir,jr,i,j,k)
*
mxmd es r4,9
  lod [r17,r4],r4
  es r5,11
  lod [r17,r5],r5
  es r3,6
  lod [r17,r3],r3
  ex r6,#ffff
  ex r7,mxmb
  lod [r4],r4
  lod [r5],r5
  mpyx r4,r5,r4
mxmd1 pack r4,r3,r3
mxmd2 vtov,a ,r3,
  ltor r3,r8
  ibxeq ,r8,[r7],r4,
  shifti r8,6,r5
  subx r4,r8,r4
  addx r3,r5,r3
  ibxle,rel ,r4,mxmd1,r6,
  pack r6,r3,r3
  bim ,mxmd2
  end
dagger ident
*
*    coded     v.r.saunders ---- may 84.
*
  msec 2
r3 equ 3*64
r4 equ 4*64
r5 equ 5*64
r6 equ 6*64
r7 equ 7*64
r8 equ 8*64
r9 equ 9*64
rb equ #b*64
rc equ #c*64
rd equ #d*64
re equ #e*64
rf equ #f*64
r10 equ #10*64
r11 equ #11*64
r12 equ #12*64
r16 equ #16*64
r17 equ #17*64
r1a equ #1a*64
  entry dagger
*   call dagger(ncol,nrow,a,mrowa,r,mrowr)
*   r(nrow,col) = a(ncol,nrow)   transposed
*
*   fortran equivalent
*
*      subroutine dagger(ncol,nrow,a,mrowa,r,mrowr)
*      dimension a(mrowa,1),r(mrowr,1)
*      if(ncol.gt.nrow)goto 99
*      do 1 j=1,ncol
*      do 1 i=1,nrow
*1     r(i,j)=a(j,i)
*      return
*99    do 2 j=1,nrow
*      do 2 i=1,ncol
*2     r(j,i)=a(i,j)
*      return
*      end
dagger lod [r17,r16],r4
  lod [r17],r3
  es r6,3
  lod [r17,r6],r6
  es r8,5
  lod [r17,r8],r8
  es r5,2
  lod [r17,r5],r5
  es r7,4
  lod [r17,r7],r7
  lod [r4],r4
  lod [r3],r3
  lod [r6],r6
  lod [r8],r8
  ex rb,#ffff
  ex r9,koop1
  ex r12,koop2
  ibxgt,rel ,r3,p99,r4,
  shifti r8,6,rc
koop1 rtor r4,rf
  rtor r5,rd
  rtor r7,re
xrs1 pack rf,r6,r6
wrs1 vxtov,fia [r6],rd,re
  ltor r6,r10
  ibxeq,rel ,r10,vic1,rf,
  shifti r10,6,r11
  subx rf,r10,rf
  mpyx r11,r6,r12
  addx re,r11,re
  addx rd,r12,rd
  ibxle,rel ,rf,xrs1,rb,
  pack rb,r6,r6
  bim ,wrs1
vic1 addx r7,rc,r7
  is r5,64
  dbnz r3,[r9]
  barb,br ,r1a
*
p99 shifti r6,6,rc
koop2 rtor r3,rf
  rtor r5,rd
  rtor r7,re
xrs2 pack rf,r8,r8
wrs2 vtovx,fia [r8],rd,re
  ltor r8,r10
  ibxeq,rel ,r10,vic2,rf,
  shifti r10,6,r11
  subx rf,r10,rf
  mpyx r11,r8,r9
  addx rd,r11,rd
  addx re,r9,re
  ibxle,rel ,rf,xrs2,rb,
  pack rb,r8,r8
  bim ,wrs2
vic2 addx r5,rc,r5
  is r7,64
  dbnz r4,[r12]
  barb,br ,r1a
  end
leadz ident
*   bramo assai, poco spero, nulla chieggio.
*                tasso
*
*    coded     v.r.saunders nov 83 to nov 84.
*
  msec 1
ordv res,128 64
disc msec 3
irep res,128 128
npb res,64 64*51
discne msec 3
ipbpos res,128 16*64
mlen res,64 16*64
length res,64 16*64
iostat res,64 16*64
iobuff res,64 16*64
itab res,64 16*64
istabf res,64 9*64
ipribf res,64 9*64
junibf res,64 9*64
jpbbf res,64 9*64
junjpb res,64 9*64
ibfbas res,64 9*64
igetm res,64 9*64
iputm res,64 9*64
mbuff res,64 8*64
  msec 2
  entry leadz
r3 equ 3*64
r4 equ 4*64
r5 equ 5*64
r6 equ 6*64
r6h equ r6+32
r7 equ 7*64
r7h equ r7+32
r8 equ 8*64
r9 equ 9*64
ra equ #a*64
rb equ #b*64
rc equ #c*64
r16 equ #16*64
r17 equ #17*64
r18 equ #18*64
r18h equ r18+32
r1a equ #1a*64
leadz lod [r17],r3
  es r4,0
  elen r4,47
  es r5,47
  es r6,0
  lod [r3],r3
  elen r5,17
  extb r3,r4,r4
  ibxeq,rel ,r4,part2,,
  addn r4,,r4
  exp r4,r4
  subx r6,r4,r18
  barb,br ,r1a
part2 extb r3,r5,r5
  es r6,17
  ibxeq,rel ,r5,part3,,r5
  addn r5,,r5
  exp r5,r5
  subx r6,r5,r18
  barb,br ,r1a
part3 es r18,64
  barb,br ,r1a
  entry lead1
lead1 lod [r17],r3
  es r4,0
  elen r4,47
  es r5,47
  es r6,0
  lod [r3],r3
  shifti ,-63,rc
  elen r5,17
  rxor r3,rc,r3
  extb r3,r4,r4
  ibxeq,rel ,r4,qart2,,
  addn r4,,r4
  exp r4,r4
  subx r6,r4,r18
  barb,br ,r1a
qart2 extb r3,r5,r5
  es r6,17
  ibxeq,rel ,r5,qart3,,r5
  addn r5,,r5
  exp r5,r5
  subx r6,r5,r18
  barb,br ,r1a
qart3 es r18,64
  barb,br ,r1a
  entry popcnt
popcnt lod [r17],r3
  elen r3,64
  cnto [r3],r18
  barb,br ,r1a
  entry pack2
*  b=pack2(left32,right32)
pack2 lod [r17],r3
  lod [r17,r16],r4
  lodh [r3,r16],r18
  lodh [r4,r16],r18h
  barb,br ,r1a
  entry upack2
* call upack2(b,left32,right32)
upack2 lod [r17],r3
  lod [r17,r16],r4
  es r6,2
  lod [r17,r6],r5
  es r7,0
  lodh [r3],r6h
  lodh [r3,r16],r7h
  sto [r4],r6
  sto [r5],r7
  barb,br ,r1a
  entry shiftr
shiftr lod [r17,r16],r4
  lod [r17],r3
  es r5,64
  lod [r4],r4
  lod [r3],r18
  ibxeq ,r4,[r1a],,
  subx r5,r4,r7
  rtor r18,r3
  es r18,0
  pack r7,,r7
  ibxge ,r4,[r1a],r5,
  extb r3,r7,r18
  barb,br ,r1a
  entry shiftl
shiftl lod [r17,r16],r4
  lod [r17],r3
  es r5,64
  lod [r4],r4
  lod [r3],r18
  ibxeq ,r4,[r1a],,
  subx r5,r4,r7
  rtor r18,r3
  es r18,0
  pack r7,,r7
  ibxge ,r4,[r1a],r5,
  insb r3,r7,r18
  barb,br ,r1a
  entry fmove
*   call fmove(a,r,n)
*   r = a
fmove es r5,2
  lod [r17,r5],r5
  lod [r17,r16],r4
  lod [r17],r3
  lod [r5],r5
  ibxeq ,r3,[r1a],r4,
  elen r3,#ffff
  ltor r3,r7
fm1 pack r5,r4,r4
fm2 vtov [r3],r4,
  ltor r4,ra
  ibxeq ,ra,[r1a],r5,
  shifti ra,6,r9
  subx r5,ra,r5
  addx r4,r9,r4
  addx r3,r9,r3
  ibxle,rel ,r5,fm1,r7,
  pack r7,r4,r4
  bim ,fm2
  entry szero
*   call szero(r,n)
*   r = 0.0
szero lod [r17,r16],r4
  lod [r17],r3
  ex r7,#ffff
  lod [r4],r4
ze1 pack r4,r3,r3
ze2 vtov,a ,r3,
  ltor r3,ra
  ibxeq ,ra,[r1a],r4,
  shifti ra,6,r9
  addx r3,r9,r3
  subx r4,ra,r4
  ibxle,rel ,r4,ze1,r7,
  pack r7,r3,r3
  bim ,ze2
  entry vecsum
*   scalar=vecsum(a,b,n)   =   dot product of a,b
vecsum es r5,2
  lod [r17,r5],r5
  lod [r17],r3
  lod [r17,r16],r4
  rtor ,r18
  lod [r5],r5
  pack r5,r3,r3
  pack r5,r4,r4
  dotv [r3],[r4],ra,
  ibxeq ,r5,[r1a],,
  ltor r3,r9
  ibxne,rel ,r9,ve3,r5,
  addn ,ra,ra
  addn ra,rb,r18
  barb,br ,r1a
ve3 elen r3,#ffff
  elen r4,#ffff
  ltor r3,r7
ve2 shifti r9,6,rc
  subx r5,r9,r5
  addx r3,rc,r3
  addx r4,rc,r4
  ibxge,rel ,r5,ve1,r7,
  pack r5,r3,r3
  pack r5,r4,r4
ve1 addn ,ra,ra
  addn ra,rb,ra
  addn ra,r18,r18
  ibxeq ,r5,[r1a],,
  dotv [r3],[r4],ra,
  ltor r3,r9
  bim ,ve2
  entry sumup
*   scalar=sumup(a,n)   =   sum of elements of  a
sumup lod [r17,r16],r5
  lod [r17],r3
  rtor ,r18
  lod [r5],r5
  pack r5,r3,r3
  sum [r3],ra,
  ibxeq ,r5,[r1a],,
  ltor r3,r9
  ibxne,rel ,r9,sumup3,r5,
  addn ,ra,ra
  addn ra,rb,r18
  barb,br ,r1a
sumup3 elen r3,#ffff
  ltor r3,r7
sumup2 shifti r9,6,rc
  subx r5,r9,r5
  addx r3,rc,r3
  ibxge,rel ,r5,sumup1,r7,
  pack r5,r3,r3
sumup1 addn ,ra,ra
  addn ra,rb,ra
  addn ra,r18,r18
  ibxeq ,r5,[r1a],,
  sum [r3],ra,
  ltor r3,r9
  bim ,sumup2
  entry subvec
*    call subvec(r,a,b,n)
*    r = a - b
subvec es r5,3
  lod [r17,r5],r5
  lod [r17],rb
  lod [r17,r16],r3
  es r4,2
  lod [r17,r4],r4
  lod [r5],r5
  elen r3,#ffff
  elen r4,#ffff
  ltor r3,r7
su1 pack r5,rb,rb
su2 subnv [r3],[r4],rb,
  ltor rb,ra
  ibxeq ,ra,[r1a],r5,
  shifti ra,6,r9
  subx r5,ra,r5
  addx rb,r9,rb
  addx r3,r9,r3
  addx r4,r9,r4
  ibxle,rel ,r5,su1,r7,
  pack r7,rb,rb
  bim ,su2
_IFN(charmm)
  entry addvec
*    call addvec(r,a,b,n)
*    r = a + b
addvec es r5,3
  lod [r17,r5],r5
  lod [r17],rb
  lod [r17,r16],r3
  es r4,2
  lod [r17,r4],r4
  lod [r5],r5
  elen r3,#ffff
  elen r4,#ffff
  ltor r3,r7
ad1 pack r5,rb,rb
ad2 addnv [r3],[r4],rb,
  ltor rb,ra
  ibxeq ,ra,[r1a],r5,
  shifti ra,6,r9
  subx r5,ra,r5
  addx rb,r9,rb
  addx r3,r9,r3
  addx r4,r9,r4
  ibxle,rel ,r5,ad1,r7,
  pack r7,rb,rb
  bim ,ad2
_ENDIF
  entry uuau
*    call uuau(n,ir,ia,ib)
*    ir = ia + ib
uuau lod [r17],r5
  lod [r17,r16],rb
  es r3,2
  lod [r17,r3],r3
  es r4,3
  lod [r17,r4],r4
  lod [r5],r5
  elen r3,#ffff
  elen r4,#ffff
  ltor r3,r7
uuau1 pack r5,rb,rb
uuau2 addxv [r3],[r4],rb,
  ltor rb,ra
  ibxeq ,ra,[r1a],r5,
  shifti ra,6,r9
  subx r5,ra,r5
  addx rb,r9,rb
  addx r3,r9,r3
  addx r4,r9,r4
  ibxle,rel ,r5,uuau1,r7,
  pack r7,rb,rb
  bim ,uuau2
  entry uumu
*    call uumu(n,ir,ia,ib)
*    ir = ia - ib
uumu lod [r17],r5
  lod [r17,r16],rb
  es r3,2
  lod [r17,r3],r3
  es r4,3
  lod [r17,r4],r4
  lod [r5],r5
  elen r3,#ffff
  elen r4,#ffff
  ltor r3,r7
uumu1 pack r5,rb,rb
uumu2 subxv [r3],[r4],rb,
  ltor rb,ra
  ibxeq ,ra,[r1a],r5,
  shifti ra,6,r9
  subx r5,ra,r5
  addx rb,r9,rb
  addx r3,r9,r3
  addx r4,r9,r4
  ibxle,rel ,r5,uumu1,r7,
  pack r7,rb,rb
  bim ,uumu2
  entry triad
*            call triad(n,scalar,a,b)
*            a=a+scalar*b
triad lod [r17,r16],rc
  lod [r17],r5
  es r4,3
  lod [r17,r4],r4
  es r3,2
  lod [r17,r3],r3
  lod [rc],rc
  lod [r5],r5
  elen r4,#ffff
  ltor r4,r7
  ibxeq ,rc,[r1a],,
tr1 pack r5,r3,r3
tr2 linkv,rb
  mpysv,a rc,[r4],r3,
  addnv [r3],[r3],r3,
  ltor r3,ra
  ibxeq ,ra,[r1a],r5,
  shifti ra,6,r9
  subx r5,ra,r5
  addx r4,r9,r4
  addx r3,r9,r3
  ibxle,rel ,r5,tr1,r7,
  pack r7,r3,r3
  bim ,tr2
  entry gtriad
*            call gtriad(n,scalar,r,a,b)
*            r=a+scalar*b
gtriad lod [r17],r5
  lod [r17,r16],rc
  es r4,4
  lod [r17,r4],r4
  es r3,3
  lod [r17,r3],r3
  es r6,2
  lod [r17,r6],r6
  lod [r5],r5
  lod [rc],rc
  elen r3,#ffff
  elen r4,#ffff
  ltor r3,r7
gtr1 pack r5,r6,r6
gtr2 linkv,rb
  mpysv,a rc,[r4],r6,
  addnv [r3],[r6],r6,
  ltor r6,ra
  ibxeq ,ra,[r1a],r5,
  shifti ra,6,r9
  subx r5,ra,r5
  addx r4,r9,r4
  addx r3,r9,r3
  addx r6,r9,r6
  ibxle,rel ,r5,gtr1,r7,
  pack r7,r6,r6
  bim ,gtr2
  entry gtrian
*            call gtrian(n,scalar,r,a,b)
*            r=scalar*b - a
gtrian lod [r17],r5
  lod [r17,r16],rc
  es r4,4
  lod [r17,r4],r4
  es r3,3
  lod [r17,r3],r3
  es r6,2
  lod [r17,r6],r6
  lod [r5],r5
  lod [rc],rc
  elen r3,#ffff
  elen r4,#ffff
  ltor r3,r7
gtrn1 pack r5,r6,r6
gtrn2 linkv,ra
  mpysv,b [r4],rc,r6,
  subnv [r6],[r3],r6,
  ltor r6,ra
  ibxeq ,ra,[r1a],r5,
  shifti ra,6,r9
  subx r5,ra,r5
  addx r4,r9,r4
  addx r3,r9,r3
  addx r6,r9,r6
  ibxle,rel ,r5,gtrn1,r7,
  pack r7,r6,r6
  bim ,gtrn2
  entry gtriaa
*            call gtriaa(n,scalar,r,a,b)
*            r=(b+scalar)*a
gtriaa lod [r17],r5
  lod [r17,r16],rc
  es r4,4
  lod [r17,r4],r4
  es r3,3
  lod [r17,r3],r3
  es r6,2
  lod [r17,r6],r6
  lod [r5],r5
  lod [rc],rc
  elen r3,#ffff
  elen r4,#ffff
  ltor r3,r7
atr1 pack r5,r6,r6
atr2 linkv,rb
  addnv,a rc,[r4],r6,
  mpysv [r3],[r6],r6,
  ltor r6,ra
  ibxeq ,ra,[r1a],r5,
  shifti ra,6,r9
  subx r5,ra,r5
  addx r4,r9,r4
  addx r3,r9,r3
  addx r6,r9,r6
  ibxle,rel ,r5,atr1,r7,
  pack r7,r6,r6
  bim ,atr2
  entry gtriab
*            call gtriab(n,scalar,r,a,b)
*            r=(a+b)*scalar
gtriab lod [r17],r5
  lod [r17,r16],rc
  es r4,4
  lod [r17,r4],r4
  es r3,3
  lod [r17,r3],r3
  es r6,2
  lod [r17,r6],r6
  lod [r5],r5
  lod [rc],rc
  elen r3,#ffff
  elen r4,#ffff
  ltor r3,r7
btr1 pack r5,r6,r6
btr2 linkv,rb
  addnv [r3],[r4],r6,
  mpysv,a rc,[r6],r6,
  ltor r6,ra
  ibxeq ,ra,[r1a],r5,
  shifti ra,6,r9
  subx r5,ra,r5
  addx r4,r9,r4
  addx r3,r9,r3
  addx r6,r9,r6
  ibxle,rel ,r5,btr1,r7,
  pack r7,r6,r6
  bim ,btr2
  entry gtriac
*            call gtriac(n,scalar,r,a,b)
*            r=(a-b)*scalar
gtriac lod [r17],r5
  lod [r17,r16],rc
  es r4,4
  lod [r17,r4],r4
  es r3,3
  lod [r17,r3],r3
  es r6,2
  lod [r17,r6],r6
  lod [r5],r5
  lod [rc],rc
  elen r3,#ffff
  elen r4,#ffff
  ltor r3,r7
ctr1 pack r5,r6,r6
ctr2 linkv,rb
  subnv [r3],[r4],r6,
  mpysv,a rc,[r6],r6,
  ltor r6,ra
  ibxeq ,ra,[r1a],r5,
  shifti ra,6,r9
  subx r5,ra,r5
  addx r4,r9,r4
  addx r3,r9,r3
  addx r6,r9,r6
  ibxle,rel ,r5,ctr1,r7,
  pack r7,r6,r6
  bim ,ctr2
  entry vvtsas
*    call vvtsas(n,scalar1,scalar2,r,a)
*    r = a*scalar1 + scalar2
vvtsas lod [r17],r3
  lod [r17,r16],r4
  es r5,2
  lod [r17,r5],r5
  es r7,4
  lod [r17,r7],r7
  es r6,3
  lod [r17,r6],r6
  lod [r3],r3
  lod [r4],r4
  lod [r5],r5
  elen r7,#ffff
  ltor r7,r8
vvtsas1 pack r3,r6,r6
vvtsas2 linkv,rb
  mpysv,b [r7],r4,r6,
  addnv,a r5,[r6],r6,
  ltor r6,ra
  ibxeq ,ra,[r1a],r3,
  shifti ra,6,r9
  subx r3,ra,r3
  addx r6,r9,r6
  addx r7,r9,r7
  ibxle,rel ,r3,vvtsas1,r8,
  pack r8,r6,r6
  bim ,vvtsas2
  entry scaler
*            call scaler(n,scalar,r,b)
*            r=scalar*b
scaler lod [r17],r5
  lod [r17,r16],rc
  es r4,3
  lod [r17,r4],r4
  es r6,2
  lod [r17,r6],r6
  lod [r5],r5
  lod [rc],rc
  elen r4,#ffff
  ltor r4,r7
sca1 pack r5,r6,r6
sca2 mpysv,b [r4],rc,r6,
  ltor r6,ra
  ibxeq ,ra,[r1a],r5,
  shifti ra,6,r9
  subx r5,ra,r5
  addx r4,r9,r4
  addx r6,r9,r6
  ibxle,rel ,r5,sca1,r7,
  pack r7,r6,r6
  bim ,sca2
  entry uvts
*            call uvts(n,scalar,ir,b)
*            ir=scalar*b
uvts lod [r17],r5
  lod [r17,r16],rc
  es r4,3
  lod [r17,r4],r4
  es r6,2
  lod [r17,r6],r6
  lod [r5],r5
  lod [rc],rc
  elen r4,#ffff
  ltor r4,r7
uvts1 pack r5,r6,r6
uvts2 linkv,ra
  mpysv,b [r4],rc,r6,
  flrv [r6],r6,
  ltor r6,ra
  ibxeq ,ra,[r1a],r5,
  shifti ra,6,r9
  subx r5,ra,r5
  addx r4,r9,r4
  addx r6,r9,r6
  ibxle,rel ,r5,uvts1,r7,
  pack r7,r6,r6
  bim ,uvts2
  entry vuts
*            call vuts(n,scalar,r,ib)
*            r=scalar*ib
vuts lod [r17],r5
  lod [r17,r16],rc
  es r4,3
  lod [r17,r4],r4
  es r6,2
  lod [r17,r6],r6
  lod [r5],r5
  lod [rc],rc
  elen r4,#ffff
  ltor r4,r7
vuts1 pack r5,r6,r6
vuts2 linkv,ra
  addnv,a ,[r4],r6,
  mpysv,b [r6],rc,r6,
  ltor r6,ra
  ibxeq ,ra,[r1a],r5,
  shifti ra,6,r9
  subx r5,ra,r5
  addx r4,r9,r4
  addx r6,r9,r6
  ibxle,rel ,r5,vuts1,r7,
  pack r7,r6,r6
  bim ,vuts2
  entry vsav
*            call vsav(n,scalar,r,b)
*            r=scalar+b
vsav lod [r17],r5
  lod [r17,r16],rc
  es r4,3
  lod [r17,r4],r4
  es r6,2
  lod [r17,r6],r6
  lod [r5],r5
  lod [rc],rc
  elen r4,#ffff
  ltor r4,r7
vsav1 pack r5,r6,r6
vsav2 addnv,a rc,[r4],r6,
  ltor r6,ra
  ibxeq ,ra,[r1a],r5,
  shifti ra,6,r9
  subx r5,ra,r5
  addx r4,r9,r4
  addx r6,r9,r6
  ibxle,rel ,r5,vsav1,r7,
  pack r7,r6,r6
  bim ,vsav2
  entry vsavsq
*            call vsavsq(n,scalar,r,b)
*            r=(scalar+b) ** 2
vsavsq lod [r17],r5
  lod [r17,r16],rc
  es r4,3
  lod [r17,r4],r4
  es r6,2
  lod [r17,r6],r6
  lod [r5],r5
  lod [rc],rc
  elen r4,#ffff
  ltor r4,r7
vsavsq1 pack r5,r6,r6
vsavsq2 linkv,ra,rb
  addnv,a rc,[r4],r6,
  mpysv [r6],[r6],r6,
  ltor r6,ra
  ibxeq ,ra,[r1a],r5,
  shifti ra,6,r9
  subx r5,ra,r5
  addx r4,r9,r4
  addx r6,r9,r6
  ibxle,rel ,r5,vsavsq1,r7,
  pack r7,r6,r6
  bim ,vsavsq2
  entry vvmvsq
*      call vvmvsq(n,r,a,b)
*      r = (a - b) ** 2
vvmvsq lod [r17],r5
  es r3,2
  lod [r17,r3],r3
  es r4,3
  lod [r17,r4],r4
  lod [r17,r16],r6
  lod [r5],r5
  elen r3,#ffff
  elen r4,#ffff
  ltor r3,r7
vvmvsq1 pack r5,r6,r6
vvmvsq2 linkv,ra,rb
  subnv [r3],[r4],r6,
  mpysv [r6],[r6],r6,
  ltor r6,ra
  ibxeq ,ra,[r1a],r5,
  shifti ra,6,r9
  subx r5,ra,r5
  addx r6,r9,r6
  addx r3,r9,r3
  addx r4,r9,r4
  ibxle,rel ,r5,vvmvsq1,r7,
  pack r7,r6,r6
  bim ,vvmvsq2
  entry vsdv
*            call vsdv(n,scalar,r,b)
*            r=scalar/b
vsdv lod [r17],r5
  lod [r17,r16],rc
  es r4,3
  lod [r17,r4],r4
  es r6,2
  lod [r17,r6],r6
  lod [r5],r5
  lod [rc],rc
  elen r4,#ffff
  ltor r4,r7
vsdv1 pack r5,r6,r6
vsdv2 divsv,a rc,[r4],r6,
  ltor r6,ra
  ibxeq ,ra,[r1a],r5,
  shifti ra,6,r9
  subx r5,ra,r5
  addx r4,r9,r4
  addx r6,r9,r6
  ibxle,rel ,r5,vsdv1,r7,
  pack r7,r6,r6
  bim ,vsdv2
  entry ujau
*            call ujau(n,jscalar,ir,ib)
*            ir=jscalar+ib
ujau lod [r17],r5
  lod [r17,r16],rc
  es r4,3
  lod [r17,r4],r4
  es r6,2
  lod [r17,r6],r6
  lod [r5],r5
  lod [rc],rc
  elen r4,#ffff
  ltor r4,r7
ujau1 pack r5,r6,r6
ujau2 addxv,a rc,[r4],r6,
  ltor r6,ra
  ibxeq ,ra,[r1a],r5,
  shifti ra,6,r9
  subx r5,ra,r5
  addx r4,r9,r4
  addx r6,r9,r6
  ibxle,rel ,r5,ujau1,r7,
  pack r7,r6,r6
  bim ,ujau2
  entry absmax
*    dmax=absmax(n,tester,c)
absmax lod [r17],r3
  lod [r17,r16],r18
  es r4,2
  lod [r17,r4],r4
  lod [r3],r3
  lod [r18],r18
  ex ra,pipe
  elen r4,#ffff
  ltor r4,r5
  shifti r5,6,r8
pipe ibxle ,r3,[r1a],,
  ibxge,rel ,r3,drum,r5,
  pack r3,r4,r4
drum max,ma [r4],r6,r9,
  subx r3,r5,r3
  addx r4,r8,r8
  abs r9,r9
  blt r9,r18,ra
  rtor r9,r18
  barb,br ,ra
  entry absmin
*    dmin=absmin(n,tester,c)
absmin lod [r17],r3
  lod [r17,r16],r18
  es r4,2
  lod [r17,r4],r4
  lod [r3],r3
  lod [r18],r18
  ex ra,piano
  elen r4,#ffff
  ltor r4,r5
  shifti r5,6,r8
piano ibxle ,r3,[r1a],,
  ibxge,rel ,r3,flute,r5,
  pack r3,r4,r4
flute min,ma [r4],r6,r9,
  subx r3,r5,r3
  addx r4,r8,r8
  abs r9,r9
  bge r9,r18,ra
  rtor r9,r18
  barb,br ,ra
  entry minimum
*    index=minimum(n,c)
*    n  ***must*** be  .ge. 1  and  .le. 65535
minimum lod [r17],r3
  lod [r17,r16],r4
  lod [r3],r3
  pack r3,r4,r4
  min [r4],r18,r5,
  barb,br ,r1a
  entry maximum
*    index=maximum(n,c)
*    n  ***must*** be  .ge. 1  and  .le. 65535
maximum lod [r17],r3
  lod [r17,r16],r4
  lod [r3],r3
  pack r3,r4,r4
  max [r4],r18,r5,
  barb,br ,r1a
  entry minabs
*    index=minabs(n,c)
*    n  ***must*** be  .ge. 1  and  .le. 65535
minabs lod [r17],r3
  lod [r17,r16],r4
  lod [r3],r3
  pack r3,r4,r4
  min,ma [r4],r18,r5,
  barb,br ,r1a
  entry maxabs
*    index=maxabs(n,c)
*    n  ***must*** be  .ge. 1  and  .le. 65535
maxabs lod [r17],r3
  lod [r17,r16],r4
  lod [r3],r3
  pack r3,r4,r4
  max,ma [r4],r18,r5,
  barb,br ,r1a
  entry setsto
*            call setsto(n,scalar,r)
*            r=scalar
setsto lod [r17],r5
  lod [r17,r16],rc
  es r6,2
  lod [r17,r6],r6
  ex r7,#ffff
  lod [r5],r5
  lod [rc],rc
sets1 pack r5,r6,r6
sets2 vtov,a rc,r6,
  ltor r6,ra
  ibxeq ,ra,[r1a],r5,
  shifti ra,6,r9
  addx r6,r9,r6
  subx r5,ra,r5
  ibxle,rel ,r5,sets1,r7,
  pack r7,r6,r6
  bim ,sets2
  entry gather
*   call gather(n,r,a,map)
*   r(loop)=a(map(loop))
gather lod [r17],r5
  es r3,2
  lod [r17,r3],r3
  es r4,3
  lod [r17,r4],r4
  lod [r17,r16],r6
  ex r7,#ffff
  lod [r5],r5
  is r3,-64
gat1 pack r5,r4,r4
gat2 vxtov [r4],r3,r6
  ltor r4,ra
  ibxeq ,ra,[r1a],r5,
  shifti ra,6,r9
  subx r5,ra,r5
  addx r4,r9,r4
  addx r6,r9,r6
  ibxle,rel ,r5,gat1,r7,
  pack r7,r4,r4
  bim ,gat2
  entry gathr
*  call gathr(n,r,a,map)
*  r(loop)=a((loop-1)*map+1)
gathr lod [r17],r5
  es r4,3
  lod [r17,r4],r4
  lod [r17,r16],r3
  es r6,2
  lod [r17,r6],r6
  ex r7,#ffff
  lod [r5],r5
  lod [r4],r4
gatr1 pack r5,r4,r4
gatr2 vxtov,fia [r4],r6,r3
  ltor r4,ra
  ibxeq ,ra,[r1a],r5,
  shifti ra,6,r9
  mpyx r4,r9,rb
  subx r5,ra,r5
  addx r3,r9,r3
  addx r6,rb,r6
  ibxle,rel ,r5,gatr1,r7,
  pack r7,r4,r4
  bim ,gatr2
  entry scatter
*   call scatter(n,r,map,a)
*   r(map(loop)) = a(loop)
scatter lod [r17],r5
  lod [r17,r16],r3
  es r4,2
  lod [r17,r4],r4
  es r6,3
  lod [r17,r6],r6
  ex r7,#ffff
  lod [r5],r5
  es r9,64
  subx r3,r9,r3
scat1 pack r5,r4,r4
scat2 vtovx [r4],r6,r3
  ltor r4,ra
  ibxeq ,ra,[r1a],r5,
  shifti ra,6,r9
  subx r5,ra,r5
  addx r4,r9,r4
  addx r6,r9,r6
  ibxle,rel ,r5,scat1,r7,
  pack r7,r4,r4
  bim ,scat2
  entry scatt
*   call scatt(n,scalar,r,map)
*   r(map(loop)) = scalar
scatt lod [r17],r5
  lod [r17,r16],r6
  es r3,2
  lod [r17,r3],r3
  es r4,3
  lod [r17,r4],r4
  ex r7,#ffff
  lod [r5],r5
  lod [r6],r6
  es r9,64
  subx r3,r9,r3
scata1 pack r5,r4,r4
scata2 vtovx,b [r4],r6,r3
  ltor r4,ra
  ibxeq ,ra,[r1a],r5,
  shifti ra,6,r9
  subx r5,ra,r5
  addx r4,r9,r4
  ibxle,rel ,r5,scata1,r7,
  pack r7,r4,r4
  bim ,scata2
  entry scattr
*   call scattr(n,r,map,a)
*   r((loop-1)*map+1) = a(loop)
scattr lod [r17],r5
  es r4,2
  lod [r17,r4],r4
  lod [r17,r16],r3
  es r6,3
  lod [r17,r6],r6
  ex r7,#ffff
  lod [r5],r5
  lod [r4],r4
scatr1 pack r5,r4,r4
scatr2 vtovx,fia [r4],r6,r3
  ltor r4,ra
  ibxeq ,ra,[r1a],r5,
  shifti ra,6,r9
  mpyx r4,r9,rb
  subx r5,ra,r5
  addx r6,r9,r6
  addx r3,rb,r3
  ibxle,rel ,r5,scatr1,r7,
  pack r7,r4,r4
  bim ,scatr2
  entry ibasgn
*  call ibasgn(n,ia,ib,ir)
*  ir(loop)=(loop-1)*ib+ia
ibasgn lod [r17],r3
  lod [r17,r16],r4
  es r5,2
  lod [r17,r5],r5
  es r6,3
  lod [r17,r6],r6
  ex r7,#ffff
  lod [r3],r3
  lod [r4],r4
  lod [r5],r5
ibasg1 pack r3,r6,r6
ibasg2 interval r4,r5,r6,
  ltor r6,r8
  ibxeq ,r8,[r1a],r3,
  mpyx r5,r8,ra
  shifti r8,6,r9
  subx r3,r8,r3
  addx r4,ra,ra
  addx r6,r9,r6
  ibxle,rel ,r3,ibasg1,r7,
  pack r7,r6,r6
  bim ,ibasg2
  entry vvtv
*      call vvtv(n,r,a,b)
*      r = a * b
vvtv lod [r17],r5
  es r3,2
  lod [r17,r3],r3
  es r4,3
  lod [r17,r4],r4
  lod [r17,r16],r6
  lod [r5],r5
  elen r3,#ffff
  elen r4,#ffff
  ltor r3,r7
vvtv1 pack r5,r6,r6
vvtv2 mpysv [r3],[r4],r6,
  ltor r6,ra
  ibxeq ,ra,[r1a],r5,
  shifti ra,6,r9
  subx r5,ra,r5
  addx r6,r9,r6
  addx r3,r9,r3
  addx r4,r9,r4
  ibxle,rel ,r5,vvtv1,r7,
  pack r7,r6,r6
  bim ,vvtv2
  entry uvtv
*      call uvtv(n,ir,a,b)
*      ir = a * b
uvtv lod [r17],r5
  es r3,2
  lod [r17,r3],r3
  es r4,3
  lod [r17,r4],r4
  lod [r17,r16],r6
  lod [r5],r5
  elen r3,#ffff
  elen r4,#ffff
  ltor r3,r7
uvtv1 pack r5,r6,r6
uvtv2 linkv,ra,rb
  mpysv [r3],[r4],r6,
  flrv [r6],r6,
  ltor r6,ra
  ibxeq ,ra,[r1a],r5,
  shifti ra,6,r9
  subx r5,ra,r5
  addx r6,r9,r6
  addx r3,r9,r3
  addx r4,r9,r4
  ibxle,rel ,r5,uvtv1,r7,
  pack r7,r6,r6
  bim ,uvtv2
  entry vvdv
*      call vvdv(n,r,a,b)
*      r = a / b
vvdv lod [r17],r5
  es r3,2
  lod [r17,r3],r3
  es r4,3
  lod [r17,r4],r4
  lod [r17,r16],r6
  lod [r5],r5
  elen r3,#ffff
  elen r4,#ffff
  ltor r3,r7
vvdv1 pack r5,r6,r6
vvdv2 divsv [r3],[r4],r6,
  ltor r6,ra
  ibxeq ,ra,[r1a],r5,
  shifti ra,6,r9
  subx r5,ra,r5
  addx r6,r9,r6
  addx r3,r9,r3
  addx r4,r9,r4
  ibxle,rel ,r5,vvdv1,r7,
  pack r7,r6,r6
  bim ,vvdv2
  entry vsqrtv
*     call vsqrtv(n,r,a)
*     r = dsqrt(a)
vsqrtv lod [r17],r3
  lod [r17,r16],r4
  es r5,2
  lod [r17,r5],r5
  lod [r3],r3
  elen r5,#ffff
  ltor r5,r6
vsqr1 pack r3,r4,r4
vsqr2 sqrtv [r5],r4,
  ltor r4,r7
  ibxeq ,r7,[r1a],r3,
  shifti r7,6,r8
  subx r3,r7,r3
  addx r4,r8,r4
  addx r5,r8,r5
  ibxlt,rel ,r3,vsqr1,r6,
  pack r6,r4,r4
  bim ,vsqr2
  entry vfloatv
*     call vfloatv(n,r,ia)
*     r = dfloat(ia)
vfloatv lod [r17],r3
  lod [r17,r16],r4
  es r5,2
  lod [r17,r5],r5
  lod [r3],r3
  elen r5,#ffff
  ltor r5,r6
vfloat1 pack r3,r4,r4
vfloat2 addnv,a ,[r5],r4,
  ltor r4,r7
  ibxeq ,r7,[r1a],r3,
  shifti r7,6,r8
  subx r3,r7,r3
  addx r4,r8,r4
  addx r5,r8,r5
  ibxlt,rel ,r3,vfloat1,r6,
  pack r6,r4,r4
  bim ,vfloat2
  entry vfixv
*     call vfixv(n,ir,a)
*     ir = fix(a)
vfixv lod [r17],r3
  lod [r17,r16],r4
  es r5,2
  lod [r17,r5],r5
  lod [r3],r3
  elen r5,#ffff
  ltor r5,r6
vfix1 pack r3,r4,r4
vfix2 flrv [r5],r4,
  ltor r4,r7
  ibxeq ,r7,[r1a],r3,
  shifti r7,6,r8
  subx r3,r7,r3
  addx r4,r8,r4
  addx r5,r8,r5
  ibxlt,rel ,r3,vfix1,r6,
  pack r6,r4,r4
  bim ,vfix2
  entry locat1
*  n  =  locat1(list,nlist,text)
*  to locate element in list equal to text
*  if not located    locat1=0
*   ==== warning ====    o.k. only if nlist .le. 65535
locat1 lod [r17,r16],r4
  es r5,2
  lod [r17,r5],r5
  lod [r17],r6
  es r7,0
  shifti ,-63,rb
  lod [r4],r4
  lod [r5],r5
  pack r4,r6,r6
  mcmpw [r6,r7],r5,rb
  addx r16,r7,r18
  ibxne ,r7,[r1a],r4,
  es r18,0
  barb,br ,r1a
  entry nocat1
*  n  =  nocat1(list,nlist,text)
*  to locate element in list not equal to text
*  if not located    nocat1=0
*   ==== warning ====    o.k. only if nlist .le. 65535
nocat1 lod [r17,r16],r4
  es r5,2
  lod [r17,r5],r5
  lod [r17],r6
  es r7,0
  shifti ,-63,rb
  lod [r4],r4
  lod [r5],r5
  pack r4,r6,r6
  mcmpw,neq [r6,r7],r5,rb
  addx r16,r7,r18
  ibxne ,r7,[r1a],r4,
  es r18,0
  barb,br ,r1a
  entry revpri
revpri ex r4,mbuff
  lod [r4],r4
  ex r3,npb
  lod [r3],r3
  ex r6,ordv
  ex r7,ipribf
  ex r5,ipribf-64
  es rb,0
  lod [r5,r4],ra
  pack r3,r6,r6
  pack r3,r7,r7
  cmplt,b [r7],ra,r6
  sto [r5,r4],rb
  addxv,b [r7],r16,r7,r6
  barb,br ,r1a
  end
flip ident
*
*    call flip(n,r,s,x,y,a,b)
*
*    if(x.ge.y) r=a , s=b
*    if(x.lt.y) s=a , r=b
*
*    r,s must    ****not****   overwrite a,b
*
*    coded     v.r.saunders ---- nov 84.
*
maskvrs msec 3
temp res,128 2048*64
  msec 2
n equ 3*64
r equ 4*64
s equ 5*64
x equ 6*64
y equ 7*64
a equ 8*64
b equ 9*64
ffff equ #a*64
masker equ #b*64
cvl equ #c*64
r16 equ #16*64
r17 equ #17*64
r1a equ #1a*64
  entry flip
flip lod [r17],n
  lod [r17,r16],r
  es s,2
  lod [r17,s],s
  es x,3
  lod [r17,x],x
  es y,4
  lod [r17,y],y
  es a,5
  lod [r17,a],a
  es b,6
  lod [r17,b],b
  lod [n],n
  ex ffff,#ffff
  ex masker,temp
  pack ffff,x,x
  pack ffff,y,y
alan pack n,masker,masker
gary cmpge [x],[y],masker
  ltor masker,cvl
  ibxeq,rel ,cvl,neil,,
  subx n,cvl,n
  maskv a,b,r,masker
  maskv b,a,s,masker
neil ibxeq ,n,[r1a],,
  shifti cvl,6,cvl
  addx x,cvl,x
  addx y,cvl,y
  addx a,cvl,a
  addx b,cvl,b
  addx r,cvl,r
  addx s,cvl,s
  ibxle,rel ,n,alan,ffff,
  pack ffff,masker,masker
  bim ,gary
  end
symm1 ident
*
*    coded     v.r.saunders ---- apr 84.
*
*  call symm1(r,a,n)
*  r(triangle)= a + a(dagger)
*
*  fortran equivalent
*
*      subroutine symm1(r,a,n)
*      dimension r(1),a(n,1)
*      m=1
*      do 1 i=1,n
*      do 1 j=1,i
*      r(m)=a(i,j)+a(j,i)
*1     m=m+1
*      return
*      end
  msec 2
r3 equ 3*64
r4 equ 4*64
r5 equ 5*64
r6 equ 6*64
r7 equ 7*64
r8 equ 8*64
r9 equ 9*64
r16 equ #16*64
r17 equ #17*64
r1a equ #1a*64
  entry symm1
symm1 es r3,2
  lod [r17,r3],r3
  lod [r17,r16],r4
  lod [r17],r5
  es r8,1
  lod [r3],r3
  rtor r4,r6
  shifti r3,6,r7
back pack r8,r3,r3
  vxtov,fia [r3],r6,r5
  is r6,64
  shifti r8,6,r9
  pack r8,r5,r5
  pack r8,r4,r4
  addnv [r5],[r4],r5,
  addx r4,r7,r4
  addx r5,r9,r5
  ibxle,rel r16,r8,back,r3,r8
  barb,br ,r1a
  end
anti1 ident
*
*    coded     v.r.saunders ---- apr 84.
*
  msec 2
r3 equ 3*64
r4 equ 4*64
r5 equ 5*64
r6 equ 6*64
r7 equ 7*64
r8 equ 8*64
r9 equ 9*64
r16 equ #16*64
r17 equ #17*64
r1a equ #1a*64
  entry anti1
*  call anti1(r,a,n)
*  r(triangle with no diagonal) = a - a(dagger)
*
*  fortran equivalent
*
*      subroutine anti1(r,a,n)
*      dimension r(1),a(n,1)
*      m=1
*      do 1 i=2,n
*      im1=i-1
*      do 1 j=1,im1
*      r(m)=a(j,i)-a(i,j)
*1     m=m+1
*      return
*      end
anti1 es r3,2
  lod [r17,r3],r3
  lod [r17,r16],r4
  lod [r17],r5
  es r8,1
  lod [r3],r3
  rtor r4,r6
  shifti r3,6,r7
bak is r6,64
  pack r8,r3,r3
  vxtov,fia [r3],r6,r5
  pack r8,r4,r4
  pack r8,r5,r5
  shifti r8,6,r9
  addx r4,r7,r4
  subnv [r4],[r5],r5,
  addx r5,r9,r5
  ibxlt,rel r16,r8,bak,r3,r8
  barb,br ,r1a
  end
square ident
*
*    coded     v.r.saunders ---- mar 84.
*
  msec 2
r3 equ 3*64
r4 equ 4*64
r5 equ 5*64
r6 equ 6*64
r7 equ 7*64
r8 equ 8*64
r9 equ 9*64
ra equ #a*64
r16 equ #16*64
r17 equ #17*64
r1a equ #1a*64
  entry square
*  call square(r,a,mr,n)
*  r(square) =  symmetrized  a(triangle)
*
*  fortran equivalent
*
*      subroutine square(r,a,mr,n)
*      dimension r(mr,1),a(1)
*      m=0
*      do 1 loop=1,n
*      do 2 moop=1,loop
*      r(moop,loop)=a(m+moop)
*2     r(loop,moop)=a(m+moop)
*1     m=m+loop
*      return
*      end
square es r5,2
  lod [r17,r5],r5
  es r6,3
  lod [r17,r6],r6
  lod [r17],r3
  lod [r17,r16],r4
  es r8,1
  lod [r5],r5
  lod [r6],r6
  rtor r3,r7
  shifti r5,6,r9
squa1 pack r8,r4,r4
  pack r8,r3,r3
  vtov [r4],r3,
  ibxeq,rel ,r8,squa2,r16,
  subx r8,r16,ra
  pack ra,r5,r5
  vtovx,fia [r5],r4,r7
squa2 shifti r8,6,ra
  is r7,64
  addx r3,r9,r3
  addx r4,ra,r4
  ibxle,rel r16,r8,squa1,r6,r8
  barb,br ,r1a
  end
sqtrip ident
*
*    coded    v.r.saunders ---- mar 84.
*
  msec 2
r3 equ 3*64
r4 equ 4*64
r5 equ 5*64
r6 equ 6*64
r7 equ 7*64
r8 equ 8*64
r9 equ 9*64
ra equ #a*64
r16 equ #16*64
r17 equ #17*64
r1a equ #1a*64
  entry sqtrip
*  call sqtrip(r,a,n)
*  r(square) = anti-symmetrized a(triangle without diagonal)
*  upper triangle +ve    lower triangle -ve   of   r
*
*  fortran equivalent
*
*      subroutine sqtrip(r,a,n)
*      dimension r(n,1),a(1)
*      m=0
*      do 1 j=2,n
*      jm1=j-1
*      do 2 i=1,jm1
*      r(i,j)=a(m+i)
*2     r(j,i)=-a(m+i)
*1     m=m+jm1
*      do 3 i=1,n
*3     r(i,i)=0.0
*      return
*      end
sqtrip es r6,2
  lod [r17,r6],r6
  lod [r17],r3
  lod [r17,r16],r4
  es r8,1
  lod [r6],r6
  rtor r3,r7
  rtor r3,r5
  shifti r6,6,r9
squt1 pack r8,r3,r3
  addx r3,r9,r3
  pack r8,r4,r4
  subnv,a ,[r4],r3,
  is r7,64
  pack r8,r6,r6
  vtovx,fia [r6],r3,r7
  vtov [r4],r3,
  shifti r8,6,ra
  addx r4,ra,r4
  ibxlt,rel r16,r8,squt1,r6,r8
  pack r6,r6,r6
  is r6,1
  vtovx,b,fia [r6],,r5
  barb,br ,r1a
  end
locate ident
*
*    coded     v.r.saunders ---- feb 84.
*
scr205 msec 3
temp res,128 64
  msec 2
r3 equ 3*64
r4 equ 4*64
r5 equ 5*64
r6 equ 6*64
r7 equ 7*64
r8 equ 8*64
r9 equ 9*64
ra equ #a*64
rb equ #b*64
rc equ #c*64
rd equ #d*64
re equ #e*64
r16 equ #16*64
r17 equ #17*64
r18 equ #18*64
r1a equ #1a*64
  entry locate
*  n  =  locate(list,nlist,spacer,text)
*  to locate element in list equal to text
*  if not located    result=0
locate lod [r17,r16],r4
  es r3,2
  lod [r17,r3],r3
  es r5,3
  lod [r17,r5],r5
  lod [r17],r6
  es r18,1
  shifti ,-63,rb
  lod [r4],r4
  lod [r3],r3
  lod [r5],r5
  ex rd,#ffff
  ex rc,temp
  pack r4,r3,r9
vrs3 vxtov,fia [r9],r6,rc
  ltor r9,r8
  mpyx r8,r3,r9
  ibxeq,rel ,r4,vrs2,,
  shifti r9,6,r9
  pack r8,rc,rc
  es r7,0
  addx r6,r9,r6
  mcmpw [rc,r7],r5,rb
  subx r4,r8,r4
  pack rd,r3,r9
  ibxge,rel ,r4,vrs1,rd,
  pack r4,r3,r9
vrs1 addx r18,r7,r18
  ibxne ,r7,[r1a],r8,
  ibxne,rel ,r4,vrs3,,
vrs2 es r18,0
  barb,br ,r1a
  entry vivits
*  call vivits(n,scalar,r,mr,a,ma)
*  r((mr-1)*loop+1)=a((ma-1)*loop+1)*scalar
vivits lod [r17],r3
  es r8,5
  lod [r17,r8],r8
  es r6,3
  lod [r17,r6],r6
  lod [r17,r16],r4
  es r7,4
  lod [r17,r7],r7
  es r5,2
  lod [r17,r5],r5
  ex r9,#ffff
  ex ra,temp
  lod [r3],r3
  lod [r8],r8
  lod [r6],r6
  lod [r4],r4
  pack r3,r8,rb
  shifti r8,6,r8
  shifti r6,6,rc
civ3 vxtov,fia [rb],r7,ra
  ltor rb,rd
  mpyx rd,r8,re
  ibxeq,rel ,rd,civ1,,
  pack rd,ra,ra
  addx r7,re,r7
  mpysv,b [ra],r4,ra,
  subx r3,rd,r3
  pack rd,r6,r6
  mpyx rc,rd,re
  vtovx,fia [r6],ra,r5
  addx r5,re,r5
civ1 ibxeq ,r3,[r1a],,
  ibxge,rel ,r3,civ2,r9,
  rtor r3,r9
civ2 pack r9,rb,rb
  bim ,civ3
  entry vxvxts
*  call vxvxts(n,scalar,r,mr,a,ma)
*  r(mr(loop))=a(ma(loop))*scalar
vxvxts lod [r17],r3
  lod [r17,r16],r4
  es r7,4
  lod [r17,r7],r7
  es r5,2
  lod [r17,r5],r5
  es r8,5
  lod [r17,r8],r8
  es r6,3
  lod [r17,r6],r6
  ex r9,#ffff
  ex ra,temp
  lod [r3],r3
  lod [r4],r4
  is r7,-64
  is r5,-64
vic2 pack r3,r8,r8
vic3 vxtov [r8],r7,ra
  ltor r8,rd
  ibxeq,rel ,rd,vic1,,
  pack rd,ra,ra
  shifti rd,6,re
  addx r8,re,r8
  mpysv,b [ra],r4,ra,
  subx r3,rd,r3
  pack rd,r6,r6
  vtovx [r6],ra,r5
  addx r6,re,r6
vic1 ibxeq ,r3,[r1a],,
  ibxle,rel ,r3,vic2,r9,
  pack r9,r8,r8
  bim ,vic3
  entry rvtvar
*  call rvtvar(n,r,a,b)
*  r = a * b + r
rvtvar lod [r17],r3
  es r5,2
  lod [r17,r5],r5
  es r6,3
  lod [r17,r6],r6
  lod [r17,r16],r4
  ex ra,temp
  lod [r3],r3
  elen r5,#ffff
  elen r6,#ffff
  ltor r5,r9
srv2 pack r3,ra,ra
srv3 mpysv [r5],[r6],ra,
  ltor ra,rd
  ibxeq,rel ,rd,srv1,,
  shifti rd,6,re
  pack rd,r4,r4
  subx r3,rd,r3
  addx r5,re,r5
  addx r6,re,r6
  addnv [ra],[r4],r4,
  addx r4,re,r4
srv1 ibxeq ,r3,[r1a],,
  ibxle,rel ,r3,srv2,r9,
  pack r9,ra,ra
  bim ,srv3
  end
upak8d ident
*
*   gamess version   (changed jvl 1988)
*   disks3 => bufb  / ranio => blksiz (different contents)
*
scr205 msec 3
i205 res,128 8190*4*64
bufb   msec 3
nwb res,128 12288*64
blksiz msec 3
idt res,128 64*3
nsz170 res,64 64
  msec 2
r3 equ 3*64
r4 equ 4*64
r5 equ 5*64
r6 equ 6*64
r7 equ 7*64
r8 equ 8*64
r9 equ 9*64
ra equ #a*64
rb equ #b*64
rc equ #c*64
r14 equ #14*64
r16 equ #16*64
r17 equ #17*64
r1a equ #1a*64
  entry upak8d
*    call upak8d to unpack data from psort file in direct
*    so routines need common/scr205/ i205(8190),j205(8190),k205..
*
*    coded     v.r.saunders ---- aug 84.
*
upak8d ex r8,nsz170
  lod [r8],r8
  ex r3,nwb
  lod [r3],r7
  ex r4,i205
  ex r9,8190*64
  is r3,64*2
  es r5,40
  es r6,#ff
  rtor r8,rc
  shifti r8,7,rb
  shifti r8,6,ra
  addx r3,rb,r3
  ibxgt,rel ,r7,vrs1,r8,
  rtor r7,rc
vrs1 pack rc,r4,r4
  pack rc,r3,r3
* unpack first nsz170   i
  linkv,ra
  shiftv,b [r3],r5,r4,
  andv,b [r4],r6,r4,
  addx r4,ra,ra
  es r5,48
  addx r4,r9,r4
* unpack first nsz170   j
  linkv,ra
  shiftv,b [r3],r5,r4,
  andv,b [r4],r6,r4,
  es r5,56
  addx r4,r9,r4
* unpack first nsz170   k
  linkv,ra
  shiftv,b [r3],r5,r4,
  andv,b [r4],r6,r4,
  addx r4,r9,r4
* unpack first nsz170   l
  andv,b [r3],r6,r4,
  ibxle ,r7,[r1a],r8,
  subx r7,r8,r7
  es r5,8
  pack r7,ra,r4
* unpack second nsz170   i
  linkv,ra
  shiftv,b [r3],r5,r4,
  andv,b [r4],r6,r4,
  addx r4,r9,r4
  es r5,16
* unpack second nsz170   j
  linkv,ra
  shiftv,b [r3],r5,r4,
  andv,b [r4],r6,r4,
  addx r4,r9,r4
  es r5,24
* unpack second nsz170   k
  linkv,ra
  shiftv,b [r3],r5,r4,
  andv,b [r4],r6,r4,
  addx r4,r9,r4
* unpack second nsz170   l
  linkv,ra
  shiftv,b [r3],r14,r4,
  andv,b [r4],r6,r4,
  barb,br ,r1a
  end
upak4d ident
*
*  call upak4d(nwb,msz408,ia,ir)
*
*  unpacks  nwb  16-bit integers from ia, result to ir
*  ia is of length msz408 (divisible by 4)
*
*    coded    v.r. saunders ---- aug 85.
*
  msec 2
msz102 equ 3*64
nwb equ 4*64
ia equ 5*64
ir equ 6*64
msz816 equ 7*64
msz408 equ 8*64
maskk equ 9*64
shif equ #a*64
r16 equ #16*64
r17 equ #17*64
r1a equ #1a*64
  entry upak4d
upak4d lod [r17,r16],msz408
  lod [r17],nwb
  es ia,2
  lod [r17,ia],ia
  es ir,3
  lod [r17,ir],ir
  lod [msz408],msz408
  lod [nwb],nwb
  elen ia,#ffff
  es shif,48
  ltor ia,maskk
  shifti msz408,62,msz102
  shifti msz408,4,msz816
  pack msz102,ir,ir
* unpack first msz102 integers
  andv,b [ia],maskk,ir,
* unpack second,third,fourth msz102 integers
top subx nwb,msz102,nwb
  ibxle ,nwb,[r1a],,
  addx ir,msz816,ir
  ibxge,rel ,nwb,middle,msz102,
  pack nwb,ir,ir
middle linkv,ra
  shiftv,b [ia],shif,ir,
  andv,b [ir],maskk,ir,
  is shif,-16
  bim ,top
  end
upack ident
*
*      subroutine upack(iconf,nint,iorb)
*c     unpack occupation scheme in 'conf' into 'iorb('nint')'
*c     conf has 2 bits/orbital ; 1 bit/electron
*c     so room for 32 orbitals / word
*      dimension i0102(4),iconf(2),iorb(2)
*      data i0102/0,1,0,2/
*      mm=min0(32,nint)
*      joker=iconf(1)
*      ii=1
*11    do 10 i=ii,mm
*      joker=shift(joker,2)
*10    iorb(i)=i0102(and(joker,3)+1)
*      if(mm.eq.nint)return
*      mm=nint
*      joker=iconf(2)
*      ii=33
*      goto 11
*      end
*
*    coded    v.r. saunders ---- aug 85.
*
  msec 2
iconf equ 3*64
nint equ 4*64
iorb equ 5*64
joker equ 6*64
temp equ 7*64
tem equ 8*64
ii equ 9*64
mm equ #a*64
r16 equ #16*64
r17 equ #17*64
r1a equ #1a*64
  entry upack
upack lod [r17,r16],nint
  lod [r17],iconf
  es iorb,2
  lod [r17,iorb],iorb
  lod [nint],nint
  lod [iconf],joker
  es mm,32
  es ii,1
  is iorb,-64
  ibxgt,rel ,nint,p10,mm,
p11 rtor nint,mm
p10 shifti joker,1,temp
  shifti joker,2,joker
  rand temp,r16,temp
  rand joker,r16,tem
  addx temp,tem,tem
  sto [iorb,ii],tem
  ibxle,rel r16,ii,p10,mm,ii
  ibxeq ,mm,[r1a],nint,
  lod [iconf,r16],joker
  bim ,p11
  end
mixup ident
maskvrs msec 3
temp res,128 1024*64
tem res,64 1024*64
  msec 2
a equ 3*64
b equ 4*64
c equ 5*64
n equ 6*64
masker equ 7*64
i equ 8*64
j equ 9*64
k equ #a*64
l equ #b*64
maske equ #c*64
mx equ #d*64
my equ #e*64
r16 equ #16*64
r17 equ #17*64
r18 equ #18*64
r1a equ #1a*64
  entry mixup
*
*    call mixup(n,a,b)
*    if j.lt.l a=b  or if j.ge.l a=a
*
*    called by gabcd in sort in direct
*
*     coded    v.r.saunders ---- nov 84.
*
mixup lod [r17],n
  lod [r17,r16],a
  es b,2
  lod [r17,b],b
  lod [n],n
  ex masker,temp
  pack n,masker,masker
  maskv a,b,a,masker
  barb,br ,r1a
  entry logicq
*
*    nn=logicq(n,a,b,c,i,k,j,l)
*   if i.eq.k or j.eq.l lose corresponding elements of a,b,c
*   resultant length of a,b,c lists = function value
*
*    called by gabcd in sort in direct
*
logicq lod [r17],n
  es i,4
  lod [r17,i],i
  es k,5
  lod [r17,k],k
  es j,6
  lod [r17,j],j
  es l,7
  lod [r17,l],l
  lod [r17,r16],a
  es b,2
  lod [r17,b],b
  es c,3
  lod [r17,c],c
  lod [n],n
  ex masker,temp
  ex maske,tem
  pack n,i,i
  pack n,k,k
  pack n,masker,masker
  cmpeq [i],[k],masker
  pack n,j,j
  pack n,l,l
  pack n,maske,maske
  cmpeq [j],[l],maske
  is n,63
  shifti n,-6,n
  pack n,masker,mx
  pack n,maske,my
  iorv [mx],[my],mx,
  cpsv,z a,a,masker
  cpsv,z b,b,masker
  ltor a,r18
  ibxle ,r18,[r1a],,
  cpsv,z c,c,masker
  barb,br ,r1a
  end
stackx ident
scra  msec 3
wv res,128 3680*64
wib res,64 3680*64
wi res,64 3680*64
wj res,64 3680*64
wk res,64 3680*64
wl res,64 3680*64
  msec 2
nstak equ 3*64
i equ 5*64
j equ 6*64
k equ 7*64
l equ 8*64
r12 equ #12*64
r13 equ #13*64
*
r16 equ #16*64
r17 equ #17*64
r1a equ #1a*64
*
*    coded     v.r.saunders ---- oct 84.
*
*     changed local to scra to fit with gamess-direct  routine stack
*
*      subroutine stackx(nstak)
*      common/scra/v(3680),ib(3680),
*     *i(3680),j(3680),k(3680),l(3680)
*      do 1 n=1,nstak
*1     i(n)=or(shift(i(n),24),or(shift(j(n),16),or(shift(k(n),8),
*     *l(n))))
*      return
*      end
*
  entry stackx
stackx lod [r17],nstak
  ex k,wk
  ex r13,3680*64
  es r12,8
  lod [nstak],nstak
  pack nstak,k,k
  addx k,r13,l
*     pack k and l
  linkv,ra
  shiftv,b [k],r12,k,
  iorv [k],[l],k,
*     pack j and k,l
  subx k,r13,j
  es r12,16
  linkv,ra
  shiftv,b [j],r12,j,
  iorv [j],[k],j,
*     pack i and j,k,l
  subx j,r13,i
  es r12,24
  linkv,ra
  shiftv,b [i],r12,i,
  iorv [i],[j],i,
  barb,br ,r1a
  end
ld340 ident
  msec 2
wmstak equ 3*64
nstak equ 4*64
nsz170 equ 5*64
nwbuf equ 6*64
lnumb equ 7*64
v equ 8*64
ib equ 9*64
ai equ #a*64
buf equ #b*64
mstak equ #c*64
nsz340 equ #d*64
nw equ #e*64
vv equ #f*64
aa equ #10*64
aah equ aa+32
ln equ #11*64
r13 equ #13*64
*
r16 equ #16*64
r17 equ #17*64
r18 equ #18*64
r1a equ #1a*64
*
*    coded     v.r.saunders ---- oct 84.
*
*      function ld340(mstack,nstak,nsz170,nwbuf,lnumb,
*     *v,ib,ai,buf)
*      dimension nwbuf(1),lnumb(1),v(1),ib(1),ai(1),buf(1)
*      mstak=mstack
*      nsz340=nsz170+nsz170
*888   ld340=ib(mstak)
*      vv=v(mstak)
*      aa=ai(mstak)
*      nw=nwbuf(ld340)+1
*      ln=lnumb(ld340)+nw
*      buf(ln)=vv
*      if(nw.gt.nsz170)goto 777
*      buf(ln+nsz340)=aa
*555   nwbuf(ld340)=nw
*      mstak=mstak+1
*      if(mstak.le.nstak)goto 888
*      ld340=0
*      return
*777   buf(ln+nsz170)=pack2(aa,buf(ln+nsz170))
*      if(nw.lt.nsz340)goto 555
*      nwbuf(ld340)=0
*      mstack=mstak
*      return
*      end
*
*    used in sort1,sort2 of 4-index transformer
*    and in presort of direct-ci
*
  entry ld340
ld340 lod [r17],wmstak
  lod [r17,r16],nstak
  es nsz170,2
  lod [r17,nsz170],nsz170
  es nwbuf,3
  lod [r17,nwbuf],nwbuf
  es lnumb,4
  lod [r17,lnumb],lnumb
  es v,5
  lod [r17,v],v
  es ib,6
  lod [r17,ib],ib
  es ai,7
  lod [r17,ai],ai
  es buf,8
  lod [r17,buf],buf
  lod [wmstak],mstak
  lod [nstak],nstak
  lod [nsz170],nsz170
  is nwbuf,-64
  is lnumb,-64
  is v,-64
  is ai,-64
  subx mstak,r16,r13
  lod [ib,r13],r18
  addx nsz170,nsz170,nsz340
p888 lod [nwbuf,r18],nw
  lod [lnumb,r18],ln
  lod [v,mstak],vv
  lod [ai,mstak],aa
  addx ln,nw,ln
  addx ln,nsz340,r13
  sto [buf,ln],vv
  ibxgt,rel r16,nw,p777,nsz170,nw
  sto [buf,r13],aa
p555 sto [nwbuf,r18],nw
  lod [ib,mstak],r18
  ibxle,rel r16,mstak,p888,nstak,mstak
  es r18,0
  barb,br ,r1a
p777 addx r13,ln,ln
  stoh [buf,ln],aah
  ibxlt,rel ,nw,p555,nsz340,
  es r13,0
  sto [wmstak],mstak
  sto [nwbuf,r18],r13
  barb,br ,r1a
  end
load408 ident
*      function load408(mstack,nstak,msz102,mwbuf,mnumb,v,ib,ai,buf)
*      dimension mwbuf(1),mnumb(1),v(1),ib(1),ai(1),buf(1)
*      mstak=mstack
*      msz204=msz102+msz102
*      msz306=msz204+msz102
*      msz408=msz204+msz204
*888   load408=ib(mstak)
*      vv=v(mstak)
*      aa=ai(mstak)
*      nw=mwbuf(load408)+1
*      ln=mnumb(load408)+nw
*      buf(ln)=vv
*      if(nw.gt.msz102)goto 2
*      buf(ln+msz408)=aa
*      goto 101
*2     if(nw.gt.msz204)goto 3
*      buf(ln+msz306)=or(shift(aa,16),buf(ln+msz306))
*      goto 101
*3     if(nw.gt.msz306)goto 4
*      buf(ln+msz204)=or(shift(aa,32),buf(ln+msz204))
*      goto 101
*4     buf(ln+msz102)=or(shift(aa,48),buf(ln+msz102))
*      if(nw.eq.msz408)goto 777
*101   mwbuf(load408)=nw
*      mstak=mstak+1
*      if(mstak.le.nstak)goto 888
*      load408=0
*      return
*777   mwbuf(load408)=0
*      mstack=mstak
*      return
*      end
*
*    coded   v.r. saunders ---- aug 85.
*
  msec 2
mstack equ 3*64
nstak equ 4*64
msz102 equ 5*64
msz204 equ 6*64
msz306 equ 7*64
msz408 equ 8*64
mwbuf equ 9*64
mnumb equ #a*64
v equ #b*64
ib equ #c*64
ai equ #d*64
buf equ #e*64
mstak equ #f*64
vv equ #10*64
aa equ #11*64
nw equ #12*64
ln equ #13*64
r16 equ #16*64
r17 equ #17*64
iborig equ r17
r18 equ #18*64
r1a equ #1a*64
  entry load408
load408 lod [r17],mstack
  es msz102,2
  lod [r17,msz102],msz102
  lod [r17,r16],nstak
  es ib,6
  lod [r17,ib],ib
  es v,5
  lod [r17,v],v
  es ai,7
  lod [r17,ai],ai
  es mwbuf,3
  lod [r17,mwbuf],mwbuf
  es mnumb,4
  lod [r17,mnumb],mnumb
  es buf,8
  lod [r17,buf],buf
  lod [mstack],mstak
  lod [msz102],msz102
  lod [nstak],nstak
  rtor ib,iborig
  is ib,-64
  is v,-64
  is ai,-64
  is mwbuf,-64
  is mnumb,-64
  is buf,-64
  lod [ib,mstak],r18
  addx msz102,msz102,msz204
  addx msz204,msz102,msz306
  shifti msz102,2,msz408
p888 lod [mwbuf,r18],nw
  lod [mnumb,r18],ln
  lod [v,mstak],vv
  lod [ai,mstak],aa
  ibxgt,rel r16,nw,p2,msz102,nw
  addx ln,nw,ln
  sto [buf,ln],vv
  addx ln,msz408,ln
  bim ,p102
p2 addx ln,nw,ln
  sto [buf,ln],vv
  ibxgt,rel ,nw,p3,msz204,
  addx ln,msz306,ln
  lod [buf,ln],vv
  shifti aa,16,aa
  bim ,p101
p3 ibxgt,rel ,nw,p4,msz306,
  addx ln,msz204,ln
  lod [buf,ln],vv
  shifti aa,32,aa
  bim ,p101
p4 addx ln,msz102,ln
  lod [buf,ln],vv
  shifti aa,48,aa
  ibxeq,rel ,nw,p777,msz408,
p101 rior vv,aa,aa
p102 sto [mwbuf,r18],nw
  lod [iborig,mstak],r18
  sto [buf,ln],aa
  ibxle,rel r16,mstak,p888,nstak,mstak
  es r18,0
  barb,br ,r1a
p777 rior vv,aa,vv
  es aa,0
  sto [mstack],mstak
  sto [buf,ln],vv
  sto [mwbuf,r18],aa
  barb,br ,r1a
  end
upak8w ident
*
*    coded     v.r.saunders ---- jun 84.
*
scr205 msec 3
ijklt res,128 340*64*4
mi res,64 340*4*64
vrs res,64 340*2
  msec 2
r3 equ 3*64
r4 equ 4*64
r5 equ 5*64
r6 equ 6*64
r8 equ 8*64
r9 equ 9*64
ra equ #a*64
rb equ #b*64
rc equ #c*64
rd equ #d*64
re equ #e*64
rf equ #f*64
r10 equ #10*64
r11 equ #11*64
r12 equ #12*64
r13 equ #13*64
r14 equ #14*64
r16 equ #16*64
r17 equ #17*64
r1a equ #1a*64
  entry upak8w
*  call upak8w(ijkl(170),nwb;i205(340),...l205(340);mapper)
*
*   fortran equivalent
*
*      subroutine upak8w(ijkl,i205,mapper)
*      dimension ijkl(170),i205(340),j205(340),k205(340),l205(340)
*      dimension mapper(255)
*      nwb=ijkl(171)
*      do 1 iw=1,nwb
*      call unp8(ijkl,iw,i,j,k,l)
*      i=mapper(i)
*      j=mapper(j)
*      k=mapper(k)
*      l=mapper(l)
*      if(i.ge.j)goto 2
*      m=i
*      i=j
*      j=m
*2     if(k.ge.l)goto 3
*      m=k
*      k=l
*      l=m
*3      if((i*256+j).ge.(k*256+l))goto 4
*      m=i
*      i=k
*      k=m
*      m=j
*      j=l
*      l=m
*4     i205(iw)=i
*      i205(iw+340)=j
*      i205(iw+680)=k
*1     i205(iw+1020)=l
*      return
*      end
upak8w lod [r17],r3
  lod [r17,r16],r5
  es r13,2
  lod [r17,r13],r13
  es r4,170
  ex r12,ijklt
  es r6,#ff
  pack r4,r12,r12
  lod [r3,r4],re
*  unpack first 170   i
  pack r4,r3,r3
  es ra,40
  linkv,ra
  shiftv,b [r3],ra,r3,
  andv,b [r3],r6,r12,
  ibxge,rel ,re,vrs1,r4,
  pack re,r12,r12
*  unpack first 170   k
vrs1 shifti re,6,rd
  addx r12,rd,r11
  es ra,56
  linkv,ra
  shiftv,b [r3],ra,r11,
  andv,b [r11],r6,r11,
*  unpack first 170   j
  addx r11,rd,r8
  es ra,48
  linkv,ra
  shiftv,b [r3],ra,r8,
  andv,b [r8],r6,r8,
*  unpack first 170   l
  addx r8,rd,rc
  andv,b [r3],r6,rc,
  ibxle,rel ,re,vrs2,r4,
*  unpack second 170   i
  subx re,r4,rf
  es ra,8
  pack rf,r12,r10
  is r10,170*64
  linkv,ra
  shiftv,b [r3],ra,r10,
  andv,b [r10],r6,r10,
*  unpack second 170   k
  addx r10,rd,r10
  es ra,24
  linkv,ra
  shiftv,b [r3],ra,r10,
  andv,b [r10],r6,r10,
*  unpack second 170   j
  addx r10,rd,r10
  es ra,16
  linkv,ra
  shiftv,b [r3],ra,r10,
  andv,b [r10],r6,r10,
*  unpack second 170   l
  addx r10,rd,r10
  linkv,ra
  shiftv,b [r3],r14,r10,
  andv,b [r10],r6,r10,
*  gather  mapper(i,k,j,l)
vrs2 shifti re,2,rf
  es r10,64
  subx r13,r10,r13
  ex r10,mi
  pack rf,r12,r12
  vxtov [r12],r13,r10
*  compare  mapi,mapk  with  mapj,mapl
  addx re,re,r13
  pack r13,r10,r10
  ex rc,vrs
  pack r13,rc,rc
  shifti r13,6,r9
  addx r10,r9,rb
  cmpge [r10],[rb],rc
*  select larger/smaller as i,k/j,l
  maskv r10,rb,r12,rc
  maskv rb,r10,r8,rc
*  generate shift32(i,k) .or. (j,l) = ij,kl
  linkv,ra
  shiftv,b [r12],r14,r12,
  iorv [r12],[r8],r12,
*  compare ij with kl
  pack re,r12,r12
  pack re,r11,r11
  pack re,rc,rc
  cmpge [r12],[r11],rc
*  select larger/smaller as ij/kl
  maskv r12,r11,r10,rc
  maskv r11,r12,rb,rc
*  unpack   i205
  pack re,r5,r5
  linkv,ra
  shiftv,b [r10],r14,r10,
  andv,b [r10],r6,r5,
*  unpack   j205
  is r5,340*64
  andv,b [r10],r6,r5
*  unpack   k205
  is r5,340*64
  linkv,ra
  shiftv,b [rb],r14,rb,
  andv,b [rb],r6,r5,
*  unpack   l205
  is r5,340*64
  andv,b [rb],r6,r5,
  barb,br ,r1a
  end
isort1 ident
*
*    coded     v.r.saunders ---- jan 84.
*
junk msec 3
nwbuck res,128 1500*64
itx res,64 3400*64
ktx res,64 3400*64
gtx res,64 3400*64
blkin msec 3
gin res,128 340*64
ijkl res,64 170*64
mword res,64 64
junke msec 3
maxt res,128 64
ires res,64 64
ipass res,64 64
nteff res,64 64
npass1 res,64 64
npass2 res,64 64
lentri res,64 64
nbuck res,64 64
mloww res,64 64
mhi res,64 64
ntri res,64 64
mapper msec 3
iky res,128 256*64
scr205 msec 3
ikyi res,128 1700*64
i res,64 1700*64
j res,64 3400*64
stak msec 3
btri res,128 64
mlow res,64 64
nstack res,64 64
iblock res,64 64
mstack res,64 64
  msec 2
r3 equ 3*64
r4 equ 4*64
r5 equ 5*64
r6 equ 6*64
r7 equ 7*64
r8 equ 8*64
r9 equ 9*64
ra equ #a*64
rb equ #b*64
rc equ #c*64
rd equ #d*64
re equ #e*64
rf equ #f*64
r10 equ #10*64
r11 equ #11*64
r12 equ #12*64
r13 equ #13*64
r14 equ #14*64
r16 equ #16*64
r17 equ #17*64
r18 equ #18*64
r1a equ #1a*64
  entry isort1
*   n = isort1(itri(680);ktri(680);gtri(680))
*
*   fortran  equivalent
*
*      function isort1(itri,ktri,gtri)
*      dimension itri(680),ktri(680),gtri(680)
*      common/blkin/gin(340),ijkl(170),mword
*      common/mapper/iky(256)
*      common/junke/maxt,ires,ipass,nteff,npass1,npass2,
*     *lentri,nbuck,mloww,mhi,ntri
*      n=0
*      do 1 loop=1,mword
*      call unp8(ijkl,loop,i,j,k,l)
*      itx=iky(i)+j
*      ktx=iky(k)+l
*      gtx=gin(loop)
*      if(mloww.ge.itx.or.mhi.lt.itx)goto 2
*      n=n+1
*      itri(n)=itx
*      ktri(n)=ktx
*      gtri(n)=gtx
*2     if(mloww.ge.ktx.or.mhi.lt.ktx.or.itx.eq.ktx)goto 1
*      n=n+1
*      itri(n)=ktx
*      ktri(n)=itx
*      gtri(n)=gtx
*1     continue
*      isort1=n
*      return
*      end
isort1 ex r7,mword
  lod [r7],r7
  ex r10,mloww
  lod [r10],r10
  ex r11,mhi
  lod [r11],r11
  lod [r17],r9
  lod [r17,r16],re
  es r13,2
  lod [r17,r13],r13
  ex r3,ijkl
  ex r4,i
  es r8,170
  es r5,40
  es r6,#ff
  elen r3,170
  elen r4,170
* unpack first 170   i
  linkv,ra
  shiftv,b [r3],r5,r3,
  andv,b [r3],r6,r4,
  shifti r7,6,ra
  addx r4,ra,rb
  es r5,56
  ex rc,j
  elen rc,170
  ibxle,rel ,r8,up81,r7,
  pack r7,rb,rb
  pack r7,rc,rc
* unpack first 170   k
up81 linkv,ra
  shiftv,b [r3],r5,rb,
  andv,b [rb],r6,rb,
  es r5,48
* unpack first 170   j
  linkv,ra
  shiftv,b [r3],r5,rc,
  andv,b [rc],r6,rc,
  addx rc,ra,rd
* unpack first 170   l
  andv,b [r3],r6,rd,
  ibxle,rel ,r7,up82,r8,
  subx r7,r8,r8
  es r5,8
  pack r8,r4,rf
  is rf,170*64
* unpack second 170   i
  linkv,ra
  shiftv,b [r3],r5,rf,
  andv,b [rf],r6,rf,
  addx rf,ra,rf
  es r5,24
* unpack second 170   k
  linkv,ra
  shiftv,b [r3],r5,rf,
  andv,b [rf],r6,rf,
  pack r8,rc,rf
  is rf,170*64
  es r5,16
* unpack second 170   j
  linkv,ra
  shiftv,b [r3],r5,rf,
  andv,b [rf],r6,rf,
  addx rf,ra,rf
* unpack second 170   l
  linkv,ra
  shiftv,b [r3],r14,rf,
  andv,b [rf],r6,rf,
*   gather the iky(i) and iky(k)
up82 addx r7,r7,r5
  ex rf,iky-64
  ex r6,ikyi
  pack r5,r4,r4
  vxtov [r4],rf,r6
*   form  iky(i)+j   and   iky(k)+l
  pack r5,r6,r6
  pack r5,rc,rc
  addxv [r6],[rc],rc,
*   if  mloww  .ge.  itri/ktri set 1 mask
  cmpge,a r10,[rc],r4
*   if mhi .lt. itri/ktri set 1 mask
  cmplt,a r11,[rc],r6
  es r12,11
  pack r12,r4,rf
  pack r12,r6,r12
*   merge the  mloww  and  mhi  masks
  iorv [rf],[r12],rf,
  pack r7,rd,rd
  pack r7,rc,rc
  pack r7,r6,r6
*   if itri .eq.  ktri set 1 mask
  cmpeq [rc],[rd],r6
  pack r7,r4,r4
  addx r4,r7,r3
*  merge mloww/mhi and itri.eq.ktri masks
  ior [r3],[r6],[r6]
*   select   itri/ktri and gtri   on zeros in the mask
  cpsv,z rc,r9,r4
  cpsv,z rd,re,r4
  ltor r9,r3
  ex r12,gin
  ibxeq,rel ,r3,vrsj,,
  cpsv,z r12,r13,r4
  shifti r3,6,r8
  addx r9,r8,r9
  addx re,r8,re
  addx r13,r8,r13
*   select   ktri/itri and gtri on zeros in the mask
vrsj cpsv,z rd,r9,r6
  cpsv,z rc,re,r6
  ltor r9,r4
  addx r3,r4,r18
  ibxeq ,r4,[r1a],,
  cpsv,z r12,r13,r6
  barb,br ,r1a
  entry jsort1
*
*   fortran equivalent
*
*      subroutine jsort1
*      common/stak/btri,mlow,nstack,iblock,mstack
*      common/scr205/ibuk(3400),itxktx(3400)
*      common/scrp/nwbuck(1500),itx(3400),ktx(3400),gtx(3400)
*      mstack=1
*      do 1 i=1,nstack
*      ibuk(i)=(itx(i)-mlow)*btri+1
*1     itxktx(i)=shift(itx(i),16).or.ktx(i)
*      return
*      end
*
jsort1 ex r3,nstack
  lod [r3],r3
  ex r4,mlow
  lod [r4],r4
  ex r9,ntri
  lod [r9],r9
  ex r5,btri
  lod [r5],r5
  ex ra,j
  ex rb,ktx
  es r12,16
  ex rf,mstack
  sto [rf],r16
  ex r6,itx
  pack r3,r6,r6
  ex r7,ikyi
  pack r3,r7,r7
  subx r9,r4,r4
*   generate  (itx-mlow)*btri + 1  =  ibuck
  linkv,ra
  addnv,b [r6],r4,r6,
  mpysv,b [r6],r5,r7,
  truv [r7],r7,
*   generate itx (left shift 16) . or . ktx
  pack r3,ra,ra
  pack r3,rb,rb
  linkv,ra
  shiftv,b [r6],r12,r6,
  iorv [r6],[rb],ra,
  barb,br ,r1a
  entry isort2
*   n = isort2(itri(340);ktri(340);gtri(340))
*
*  fortran equivalent
*
*      function isort2(itri,ktri,gtri)
*      dimension itri(340),ktri(340),gtri(340)
*      common/blkin/gin(340),ijkl(170),mword
*      common/junke/maxt,ires,ipass,nteff,npass1,npass2,
*     *lentri,nbuck,mloww,mhi,ntri
*      n=0
*      do 1 loop=1,mword
*      call unp16(ijkl,loop,kl,ij)
*      if(mloww.ge.ij.or.mhi.lt.ij)goto 1
*      n=n+1
*      gtri(n)=gin(loop)
*      itri(n)=ij
*      ktri(n)=kl
*1     continue
*      isort2=n
*      return
*      end
isort2 ex r7,mword
  lod [r7],r7
  ex r10,mloww
  lod [r10],r10
  ex r11,mhi
  lod [r11],r11
  lod [r17],r9
  lod [r17,r16],re
  es r13,2
  lod [r17,r13],r13
  ex rc,j
  es r8,170
  es r5,48
  ex r6,#ffff
  ex r3,ijkl
  elen r3,170
  elen rc,170
*   unpack first 170   kl
  linkv,ra
  shiftv,b [r3],r5,r3,
  andv,b [r3],r6,rc
*   unpack first 170   ij
  shifti r7,6,ra
  addx rc,ra,rd
  ibxle,rel ,r8,vp81,r7,
  pack r7,rd,rd
vp81 andv,b [r3],r6,rd
  ibxle,rel ,r7,vp82,r8,
*   unpack second 170   kl
  subx r7,r8,r8
  pack r8,rc,rf
  is rf,170*64
  expv [r3],rf,
*   unpack second 170   ij
  addx rf,ra,rf
  linkv,ra
  shiftv,b [r3],r14,rf,
  andv,b [rf],r6,rf,
*   if mloww.ge.ij or mhi.lt.ij set 1 mask
vp82 ex r4,i
  pack r7,r4,r4
  pack r7,rd,rd
  cmpge,a r10,[rd],r4
  rtor r4,r3
  is r3,8192
  cmplt,a r11,[rd],r3
  es r12,6
  pack r12,r3,rf
  pack r12,r4,r12
  iorv [rf],[r12],r12,
*   compress to itri ktri gtri
  cpsv,z rd,r9,r4
  cpsv,z rc,re,r4
  ltor r9,r18
  ibxeq ,r18,[r1a],,
  ex r11,gin
  cpsv,z r11,r13,r4
  barb,br ,r1a
  end
pak8x ident
*
*    coded     v.r.saunders ---- jun 84.
*
*   call pak8x(i(340),j(340),k(340),l(340);ijkl(170),nwb)
*   packs i,j,k,l permuting to ensure i.ge.j k.ge.l ij.ge.kl
scr205 msec 3
ij res,128 340*64
ordv res,64 64*6
  msec 2
r3 equ 3*64
r4 equ 4*64
r5 equ 5*64
r6 equ 6*64
r7 equ 7*64
r8 equ 8*64
r9 equ 9*64
ra equ #a*64
rb equ #b*64
rc equ #c*64
rd equ #d*64
r14 equ #14*64
r16 equ #16*64
r17 equ #17*64
r1a equ #1a*64
  entry pak8x
pak8x lod [r17,r16],r4
  lod [r17],r3
  es r8,170
  es r5,8
  ex ra,ij
  shifti r8,7,r7
  shifti r8,6,r6
  lod [r4,r8],r9
*   generate   ij
  pack r9,ra,ra
  addx ra,r7,rd
  pack r9,r3,r3
  addx r3,r7,rc
  cmpge [r3],[rc],rd
  linkv,ra
  shiftv,b [r3],r5,ra,rd
  iorv [ra],[rc],ra,rd
  linkv,ra
  shiftv,b,z [rc],r5,ra,rd
  iorv,z [ra],[r3],ra,rd
*   generate   kl
  addx rc,r7,r3
  addx r3,r7,rb
  cmpge [r3],[rb],rd
  linkv,ra
  shiftv,b [r3],r5,rc,rd
  iorv [rc],[rb],rc,rd
  linkv,ra
  shiftv,b,z [rb],r5,rc,rd
  iorv,z [rc],[r3],rc,rd
*   generate   ijkl
  cmpge [ra],[rc],rd
  es r5,16
  ibxgt,rel ,r9,vrs1,r8,
  pack r9,r4,r3
vrs1 linkv,ra
  shiftv,b [ra],r5,r3,rd
  iorv [r3],[rc],r3,rd
  linkv,ra
  shiftv,b,z [rc],r5,r3,rd
  iorv,z [r3],[ra],r3,rd
  ibxle ,r9,[r1a],r8,
*   generate ijklijkl
  subx r9,r8,r9
  ibxeq,rel ,r9,vrs2,r8,
  shifti r9,6,rd
  subx r8,r9,rb
  addx r4,rd,ra
  addx r3,rd,r5
  pack rb,ra,ra
  vtov [r5],ra,
vrs2 addx r3,r6,r6
  pack r9,r4,r4
  linkv,ra
  shiftv,b [r6],r14,r4,
  iorv [r4],[r3],r4,
  barb,br ,r1a
  end
pak4v ident
*
*    coded     v.r.saunders ---- nov 83.
*
*   call pak4v(ij205(340),kl205(340);ijkl(170),nwb)
  msec 2
  entry pak4v
r3 equ 3*64
r4 equ 4*64
r5 equ 5*64
r6 equ 6*64
r7 equ 7*64
r8 equ 8*64
r14 equ #14*64
r16 equ #16*64
r17 equ #17*64
r1a equ #1a*64
pak4v lod [r17,r16],r3
  lod [r17],r4
  es r5,170
  es r6,16
  lod [r3,r5],r7
  pack r5,r3,r3
  pack r5,r4,r4
  shifti r5,7,r8
  addx r4,r8,r8
*   pack first 170  ij   and   kl
  ibxge,rel ,r7,pak1,r5,
  pack r7,r3,r3
pak1 linkv,ra
  shiftv,b [r4],r6,r3,
  iorv [r3],[r8],r3,
  ibxle ,r7,[r1a],r5,
  subx r7,r5,r7
  is r4,170*64
  is r8,170*64
*   pack second 170  ij
  pack r7,r3,r3
  packv [r4],[r3],r3,
*   pack second 170  kl
  linkv,ra
  shiftv,b [r8],r14,r3,
  iorv [r3],[r3],r3,
  barb,br ,r1a
  entry pak6v
*   call pak6v(ij205(204),kl205(204);ijkl(102),nwb)
pak6v lod [r17,r16],r3
  lod [r17],r4
  es r5,102
  es r6,16
  lod [r3,r5],r7
  pack r5,r3,r3
  pack r5,r4,r4
*   pack first 102  ij
  shiftv,b [r4],r6,r3,
  ibxge,rel ,r7,pak11,r5,
  pack r7,r3,r3
  pack r7,r4,r4
*   pack first 102  kl
pak11 es r8,204*64
  addx r4,r8,r8
  iorv [r3],[r8],r3,
  ibxle ,r7,[r1a],r5,
  subx r7,r5,r7
  is r4,102*64
  is r8,102*64
*   pack second 102  ij
  pack r7,r3,r3
  pack r7,r4,r4
  packv [r4],[r3],r3,
*   pack second 102  kl
  pack r7,r8,r8
  linkv,ra
  shiftv,b [r8],r14,r8,
  iorv [r8],[r3],r3,
  barb,br ,r1a
  entry pak2v
*   call pak2v(ijkl205(340);ijkl(170),nwb)
pak2v lod [r17,r16],r4
  lod [r17],r3
  es r5,170
  shifti r5,6,r6
  shifti r5,1,r8
  lod [r4,r5],r7
  pack r5,r4,r4
  pack r5,r3,r3
  addx r3,r6,r6
  ibxeq,rel ,r7,vrs,r8,
  vtov [r3],r4,
  ibxle ,r7,[r1a],r5,
  subx r7,r5,r8
  pack r8,r4,r4
vrs linkv,ra
  shiftv,b [r6],r14,r4,
  iorv [r4],[r3],r4,
  barb,br ,r1a
  end
transc ident
*
*    coded     v.r.saunders ---- jan 84.
*
blksiz msec 3
nsz res,128 3*64
nsz170 res,64 4*64
bufb msec 3
nkk res,128 128
gin res,64 64
scr205 msec 3
ij res,128 8000*64
junke msec 3
maxt res,128 8*64
mloww res,64 3*64
junkf msec 3
lword res,128 7*64
lenbas res,64 8*64
  msec 2
  entry transc
r3 equ 3*64
r4 equ 4*64
r5 equ 5*64
r6 equ 6*64
r7 equ 7*64
r8 equ 8*64
r9 equ 9*64
ra equ #a*64
rb equ #b*64
rc equ #c*64
rd equ #d*64
re equ #e*64
rf equ #f*64
r10 equ #10*64
r11 equ #11*64
r14 equ #14*64
r17 equ #17*64
r1a equ #1a*64
*        call transc(qq)  in calc1,calc2 of 4-index code
transc ex r3,nsz170
  lod [r3],r3
  ex r4,nkk
  lod [r4],r5
  lod [r17],r6
  ex rc,mloww
  lod [rc],rc
  ex rd,lenbas
  lod [rd],rd
  es r7,48
  ex r8,#ffff
  is r4,128
  addx r3,r3,r9
  shifti r3,7,re
  addx re,r4,ra
  pack r3,ra,ra
  ex rb,ij
  pack r3,rb,rb
*   unpack first  nsz170  ij
  linkv,ra
  shiftv,b [ra],r7,ra,
  andv,b [ra],r8,rb,
  addx rb,re,rf
  ibxle,rel ,r3,vrsa,r5,
  pack r5,rf,rf
*   unpack first  nsz170  kl
vrsa andv,b [ra],r8,rf,
  ibxle,rel ,r5,vrsb,r3,
  subx r5,r3,r10
  shifti r3,6,r11
  addx rb,r11,r11
  pack r10,r11,r11
*   unpack second  nsz170  ij
  expv [ra],r11,
  addx r11,re,r11
*   unpack second  nsz170  kl
  linkv,ra
  shiftv,b [ra],r14,r11,
  andv,b [r11],r8,r11,
vrsb pack r5,rb,rb
  pack r5,rf,rf
  addn rd,,rd
*   generate  (ij - mloww) * lenbas
  linkv,ra
  subnv,b [rb],rc,rb,
  mpysv,b [rb],rd,rb,
*   add on  kl
  truv [rb],rb,
  addxv [rb],[rf],rb,
  is r6,-64
*   scatter  gin  to  qq
  vtovx [rb],r4,r6
  barb,br ,r1a
  end
indcmp ident
*
*      function indcmp(n,t,r,s,ind,jnd)
*      dimension r(1),s(1),ind(1),jnd(1)
*      nn=0
*      do 1 i=1,n
*      if(dabs(r(i)).lt.t)goto 1
*      nn=nn+1
*      s(nn)=r(i)
*      jnd(nn)=ind(i)
*1     continue
*      indcmp=nn
*      return
*      end
*
*      coded    v.r. saunders ----- nov 1985
*
maskvrs msec 3
vmask res,128 2048*64
  msec 2
n equ 3*64
t equ 4*64
r equ 5*64
s equ 6*64
ind equ 7*64
jnd equ 8*64
vm equ 9*64
*
r16 equ #16*64
r17 equ #17*64
r18 equ #18*64
r1a equ #1a*64
  entry indcmp
indcmp lod [r17],n
  lod [r17,r16],t
  es r,2
  lod [r17,r],r
  es s,3
  lod [r17,s],s
  es ind,4
  lod [r17,ind],ind
  es jnd,5
  lod [r17,jnd],jnd
  lod [n],n
  lod [t],t
  ex vm,vmask
  pack n,r,r
  arithcps,ma,b [r],t,s,vm
  cpsv ind,jnd,vm
  ltor s,r18
  barb,br ,r1a
  end
upak8v ident
*
*    coded     v.r.saunders ---- nov 83.
*
  msec 2
r3 equ 3*64
r4 equ 4*64
r5 equ 5*64
r6 equ 6*64
r7 equ 7*64
r8 equ 8*64
ra equ #a*64
r14 equ #14*64
r16 equ #16*64
r17 equ #17*64
r1a equ #1a*64
  entry upak8v
* call upak8v(ijkl(170),nword;i(340),j(340),k(340),l(340))
upak8v lod [r17],r3
  lod [r17,r16],r4
  es r8,170
  es r5,40
  es r6,#ff
  elen r3,170
  elen r4,170
  lod [r3,r8],r7
* unpack first 170   i
  linkv,ra
  shiftv,b [r3],r5,r3,
  andv,b [r3],r6,r4,
  es ra,170*64
  addx r4,ra,ra
  es r5,48
  is r4,340*64
  ibxle,rel ,r8,up81,r7,
  pack r7,r4,r4
* unpack first 170   j
up81 linkv,ra
  shiftv,b [r3],r5,r4,
  andv,b [r4],r6,r4,
  es r5,56
  is r4,340*64
* unpack first 170   k
  linkv,ra
  shiftv,b [r3],r5,r4,
  andv,b [r4],r6,r4,
  is r4,340*64
* unpack first 170   l
  andv,b [r3],r6,r4,
  ibxle ,r7,[r1a],r8,
  subx r7,r8,r7
  es r5,8
  pack r7,ra,r4
* unpack second 170   i
  linkv,ra
  shiftv,b [r3],r5,r4,
  andv,b [r4],r6,r4,
  is r4,340*64
  es r5,16
* unpack second 170   j
  linkv,ra
  shiftv,b [r3],r5,r4,
  andv,b [r4],r6,r4,
  is r4,340*64
  es r5,24
* unpack second 170   k
  linkv,ra
  shiftv,b [r3],r5,r4,
  andv,b [r4],r6,r4,
  is r4,340*64
* unpack second 170   l
  linkv,ra
  shiftv,b [r3],r14,r4,
  andv,b [r4],r6,r4,
  barb,br ,r1a
  entry upak4v
* call upak4v(ijkl(170),nword;ij(340),kl(340))
upak4v lod [r17],r3
  lod [r17,r16],r4
  es r8,170
  es r5,48
  ex r6,#ffff
  elen r3,170
  elen r4,170
  lod [r3,r8],r7
* unpack first 170   ij
  linkv,ra
  shiftv,b [r3],r5,r3,
  andv,b [r3],r6,r4,
  es ra,170*64
  addx r4,ra,ra
  is r4,340*64
  ibxle,rel ,r8,up41,r7,
  pack r7,r4,r4
* unpack first 170   kl
up41 andv,b [r3],r6,r4,
  ibxle ,r7,[r1a],r8,
  subx r7,r8,r7
  es r5,16
  pack r7,ra,r4
* unpack second 170   ij
  linkv,ra
  shiftv,b [r3],r5,r4,
  andv,b [r4],r6,r4,
  is r4,340*64
* unpack second 170   kl
  linkv,ra
  shiftv,b [r3],r14,r4,
  andv,b [r4],r6,r4,
  barb,br ,r1a
  entry upak6v
* call upak6v(ijkl(102),nword;ij(204),kl(204))
upak6v lod [r17],r3
  lod [r17,r16],r4
  es r8,102
  es r5,48
  ex r6,#ffff
  elen r3,102
  elen r4,102
  lod [r3,r8],r7
* unpack first 102   ij
  linkv,ra
  shiftv,b [r3],r5,r3,
  andv,b [r3],r6,r4,
  es ra,102*64
  addx r4,ra,ra
  is r4,204*64
  ibxle,rel ,r8,up42,r7,
  pack r7,r3,r3
  pack r7,r4,r4
* unpack first 102   kl
up42 andv,b [r3],r6,r4,
  ibxle ,r7,[r1a],r8,
  subx r7,r8,r7
  es r5,16
  pack r7,r3,r3
  pack r7,ra,r4
* unpack second 102   ij
  linkv,ra
  shiftv,b [r3],r5,r3,
  andv,b [r3],r6,r4,
  is r4,204*64
* unpack second 102   kl
  linkv,ra
  shiftv,b [r3],r14,r3,
  andv,b [r3],r6,r4,
  barb,br ,r1a
  entry upak2v
*   call upak2v(ijkl,nwb;ijkl205)
upak2v lod [r17],r3
  lod [r17,r16],r4
  es r5,170
  ex r6,#ffffffff
  pack r5,r3,r3
  pack r5,r4,r4
  lod [r3,r5],r7
  andv,b [r3],r6,r4,
  ibxle ,r7,[r1a],r5,
  is r4,170*64
  subx r7,r5,r5
  pack r5,r4,r4
  linkv,ra
  shiftv,b [r3],r14,r4,
  andv,b [r4],r6,r4,
  barb,br ,r1a
  end
srot ident
*
*    coded     v.r.saunders ---- mar 84.
*
scr205 msec 3
temp res,128 64
  msec 2
r3 equ 3*64
r4 equ 4*64
r5 equ 5*64
r6 equ 6*64
r7 equ 7*64
r8 equ 8*64
r9 equ 9*64
r16 equ #16*64
r17 equ #17*64
r1a equ #1a*64
  entry srot
*  call drot(n,a,b,cos,sin)
*  a=a*cos+b*sin
*  b=b*cos-a*sin
srot lod [r17],r3
  es r7,4
  lod [r17,r7],r7
  es r6,3
  lod [r17,r6],r6
  lod [r17,r16],r4
  es r5,2
  lod [r17,r5],r5
  ex r8,temp
  lod [r3],r3
  lod [r7],r7
  lod [r6],r6
  shifti r3,6,r9
  pack r3,r4,r4
  pack r3,r8,r8
*    temp=sin*a
  mpysv,b [r4],r7,r8,
  pack r3,r5,r5
  addx r8,r9,r9
*     tem=cos*a
  mpysv,b [r4],r6,r9,
*     a=sin*b+tem
  linkv,ra
  mpysv,b [r5],r7,r8,
  addnv [r8],[r9],r4,
*     b=cos*b-temp
  linkv,ra
  mpysv,b [r5],r6,r9,
  subnv [r9],[r8],r5,
  barb,br ,r1a
  end
sgmat ident
*
*    coded     v.r.saunders ---- dec 83.
*
savl equ 42
  msec 1
len res,128 6*64
* lenbas,nshell,iky,jbas,kbas,pbas
blkin msec 3
g res,128 340*64
ijkl res,64 170*64
nwb res,64 64
scr205 msec 3
i res,128 340*4*64
ij res,64 340*6*64
pij res,64 340*6*64
*  ij,jk,kl, ik,jl,  il
g2 res,64 340*3*64
sav res,128 savl*64
ordv res,64 340
  msec 2
r3 equ 3*64
r4 equ 4*64
r5 equ 5*64
r6 equ 6*64
r7 equ 7*64
r8 equ 8*64
r9 equ 9*64
ra equ #a*64
rb equ #b*64
rc equ #c*64
rd equ #d*64
re equ #e*64
rf equ #f*64
r10 equ #10*64
r11 equ #11*64
r12 equ #12*64
r13 equ #13*64
*
r14 equ #14*64
r16 equ #16*64
r17 equ #17*64
r1a equ #1a*64
*
r20 equ #20*64
r21 equ #21*64
r22 equ #22*64
r23 equ #23*64
r24 equ #24*64
r25 equ #25*64
r26 equ #26*64
r27 equ #27*64
r28 equ #28*64
r29 equ #29*64
r2a equ #2a*64
r2b equ #2b*64
r2c equ #2c*64
r2d equ #2d*64
r2e equ #2e*64
r2f equ #2f*64
*
r30 equ #30*64
r31 equ #31*64
r32 equ #32*64
r33 equ #33*64
r34 equ #34*64
r35 equ #35*64
r36 equ #36*64
r37 equ #37*64
r38 equ #38*64
r39 equ #39*64
r3a equ #3a*64
r3b equ #3b*64
r3c equ #3c*64
r3d equ #3d*64
r3e equ #3e*64
r3f equ #3f*64
*
r40 equ #40*64
r41 equ #41*64
r42 equ #42*64
r43 equ #43*64
r44 equ #44*64
r45 equ #45*64
r46 equ #46*64
r47 equ #47*64
r48 equ #48*64
  entry sgmat
*
*   generates  nshell  coulomb and exchange matrices
*   from integral list + density matrices
*
sgmat ex rf,sav
  elen rf,savl
  ex r3,len
  elen r3,6
  swap r3,r14,rf
  ex r7,nwb
  lod [r7],r7
  ex r3,ijkl
  elen r3,170
  ex ra,i
  elen ra,170
  es r5,40
  es r6,#ff
*  unpack first 170  i
  linkv,ra
  shiftv,b [r3],r5,r3,
  andv,b [r3],r6,ra,
  shifti r7,6,r13
  es r8,170
  es r5,48
  addx ra,r13,rb
  ibxle,rel ,r8,sgm1,r7,
  pack r7,rb,rb
*  unpack first 170  j
sgm1 linkv,ra
  shiftv,b [r3],r5,rb,
  andv,b [rb],r6,rb,
  addx rb,r13,rc
  es r5,56
*  unpack first 170  k
  linkv,ra
  shiftv,b [r3],r5,rc,
  andv,b [rc],r6,rc,
  addx rc,r13,rd
*  unpack first 170  l
  andv,b [r3],r6,rd,
  ibxle,rel ,r7,sgm2,r8,
  subx r7,r8,r8
  ex r4,i+170*64
  es r5,8
  pack r8,r4,r4
*  unpack second 170  i
  linkv,ra
  shiftv,b [r3],r5,r4,
  andv,b [r4],r6,r4,
  addx r4,r13,r4
  es r5,16
*  unpack second 170  j
  linkv,ra
  shiftv,b [r3],r5,r4,
  andv,b [r4],r6,r4,
  addx r4,r13,r4
  es r5,24
*  unpack second 170  k
  linkv,ra
  shiftv,b [r3],r5,r4,
  andv,b [r4],r6,r4,
  addx r4,r13,r4
*  unpack second 170  l
  linkv,ra
  shiftv,b [r3],r14,r4,
  andv,b [r4],r6,r4,
sgm2 shifti r7,2,r3
  ex r4,pij
  pack r3,ra,ra
*  gather the  iky(i  j  k  l)
  vxtov [ra],r22,r4
*  generate iky(i)+j   iky(j)+k   iky(k)+l
  subx r3,r7,r6
  pack r6,r4,r4
  pack r6,rb,rb
  ex r8,ij
  pack r6,r8,r8
  addxv [r4],[rb],r8,
*  generate iky(i)+k   iky(j)+l
  shifti r6,6,re
  subx r6,r7,r12
  addx r8,re,r9
  pack r12,rc,rc
  pack r12,r9,r9
  addxv [r4],[rc],r9,
*  generate iky(i)+l
  shifti r12,6,r11
  addx r11,r9,r10
  pack r7,r10,r10
  pack r7,r4,r4
  pack r7,rd,rd
  addxv [r4],[rd],r10,
*  generate iky(k)+j  in  g2
  ex r26,g2
  pack r7,r26,r26
  addx r4,r11,r5
  addxv [r5],[rb],r26,
*  merge kj with jk
  ex r27,ordv
  pack r7,r27,r27
  cmpge [rb],[rc],r27
  addx r8,r13,r28
  maskv r28,r26,r28,r27
*  generate iky(l)+j  in  g2
  addx r5,r13,r29
  addxv [r29],[rb],r26,
*  merge lj with jl
  cmpge [rb],[rd],r27
  subx r10,r13,r2a
  maskv r2a,r26,r2a,r27
*  generate 2*g in g2
  ex r22,g
  pack r7,r22,r22
  addnv [r22],[r22],r26,
*  generate  gik
  pack r7,ra,ra
  cmpeq [ra],[rc],r27
  addx r26,r13,r2b
  addx r2b,r13,r2c
  cmpeq [rb],[rd],r2c
  ior [r27],[r2c],[r27]
  maskv r26,r22,r2b,r27
*  generate  gil
  cmpne [rb],[rc],r27
  maskv r22,r26,r2c,r27
* setup loop over shells
  addx r6,r6,r6
  pack r6,r8,r8
  addx r4,r13,r33
  pack r12,r5,r5
  pack r12,r26,r27
  addx r5,r11,r6
  pack r12,r2b,ra
  ex r3,sgm3
  addx r28,r13,r11
  addx r6,r13,r34
  ex r35,sgm4
*  gather the p elements
sgm3 vxtov [r8],r25,r4
  addx r25,r20,r25
*  generate  g2*pij
  mpysv [r26],[r4],r4,
*  generate  gil*pjk
  mpysv [r2c],[r33],r33,
*  generate  g2*pkl  and  gik*pik
  mpysv [r27],[r5],r5,
*  generate  gik*pjl  and  gil*pil
  mpysv [ra],[r6],r6,
*
*  jij=jij+g2*pkl
*  jkl=jkl+g2*pij
*  kjl=kjl+gik*pik
*  kjk=kjk+gil*pil
*  kil=kil+gil*pjk
*  kik=kik+gik*pjl
*
  es rb,0
  ibxeq,rel r16,rb,sgm5,r7,rc
  shifti r7,-1,r48
*  fetch  ij,kl,jl,jk,il,ik
sgm4 lod [r8,rb],rd
  lod [r11,rb],re
  lod [r2a,rb],r12
  lod [r28,rb],r13
  lod [r10,rb],r22
  lod [r9,rb],r2b
*  fetch  pkl,pij,pik,pil,pjk,pjl
  lod [r5,rb],r2d
  lod [r4,rb],r2e
  lod [r29,rb],r2f
  lod [r34,rb],r30
  lod [r33,rb],r31
  lod [r6,rb],r32
*  fetch  ij,kl,jl,jk,il,ik  (2)
  lod [r8,rc],r3c
  lod [r11,rc],r3d
  lod [r2a,rc],r3e
  lod [r28,rc],r3f
  lod [r10,rc],r40
  lod [r9,rc],r41
*  fetch  jij,jkl,kjl,kjk,kil,kik
  lod [r23,rd],r36
  lod [r23,re],r37
  lod [r24,r12],r38
  lod [r24,r13],r39
  lod [r24,r22],r3a
  lod [r24,r2b],r3b
*  fetch  pkl,pij,pik,pil,pjk,pjl  (2)
  lod [r5,rc],r42
  lod [r4,rc],r43
  lod [r29,rc],r44
  lod [r34,rc],r45
  lod [r33,rc],r46
  lod [r6,rc],r47
*  add updates to  j  and  k
  addn r36,r2d,r2d
  addn r37,r2e,r2e
  addn r38,r2f,r2f
  addn r39,r30,r30
*  store  jij
  sto [r23,rd],r2d
  addn r3a,r31,r31
  addn r3b,r32,r32
*  store  jkl
  sto [r23,re],r2e
*  fetch jij,jkl  (2)
  lod [r23,r3c],r36
  lod [r23,r3d],r37
*  store kjl,kjk,kil,kik
  sto [r24,r12],r2f
  sto [r24,r13],r30
  sto [r24,r22],r31
  sto [r24,r2b],r32
*  fetch  kjl,kjk,kil,kik  (2)
  lod [r24,r3e],r38
  lod [r24,r3f],r39
  lod [r24,r40],r3a
  lod [r24,r41],r3b
*  update  j  (2)
  addn r36,r42,r42
  is rb,2
  addn r37,r43,r43
  is rc,2
  sto [r23,r3c],r42
  sto [r23,r3d],r43
*  update  k  (2)
  addn r38,r44,r44
  addn r39,r45,r45
  addn r3a,r46,r46
  addn r3b,r47,r47
  sto [r24,r3e],r44
  sto [r24,r3f],r45
  sto [r24,r40],r46
  sto [r24,r41],r47
  dbnz r48,[r35]
  ibxeq,rel ,rb,sgm6,r7,
*  fetch  ij,kl,jl,jk,il,ik
sgm5 lod [r8,rb],rd
  lod [r11,rb],re
  lod [r2a,rb],r12
  lod [r28,rb],r13
  lod [r10,rb],r22
  lod [r9,rb],r2b
*  fetch  pkl,pij,pik,pil,pjk,pjl
  lod [r5,rb],r2d
  lod [r4,rb],r2e
  lod [r29,rb],r2f
  lod [r34,rb],r30
  lod [r33,rb],r31
  lod [r6,rb],r32
*  fetch  jij,jkl,kjl,kjk,kil,kik
  lod [r23,rd],r36
  lod [r23,re],r37
  lod [r24,r12],r38
  lod [r24,r13],r39
  lod [r24,r22],r3a
  lod [r24,r2b],r3b
*  add updates to  j  and  k
  addn r36,r2d,r2d
  addn r37,r2e,r2e
  addn r38,r2f,r2f
  addn r39,r30,r30
  sto [r23,rd],r2d
  addn r3a,r31,r31
  addn r3b,r32,r32
  sto [r23,re],r2e
  sto [r24,r12],r2f
  sto [r24,r13],r30
  sto [r24,r22],r31
  sto [r24,r2b],r32
sgm6 addx r23,r20,r23
  addx r24,r20,r24
  dbnz r21,[r3]
  swap rf,r14,
  barb,br ,r1a
*
*
  entry sgmat1
*  call sgmat1(lenbas,nshell,iky,j,k,p) to initialize sgmat
sgmat1 lod [r17],r3
  lod [r17,r16],r4
  es rb,2
  lod [r17,rb],r5
  es rc,3
  lod [r17,rc],r6
  es rd,4
  lod [r17,rd],r7
  es re,5
  lod [r17,re],r8
  lod [r3],r3
  lod [r4],r4
  is r5,-64
  is r6,-64
  is r7,-64
  is r8,-64
  ex ra,len
  sto [ra,rb],r5
  sto [ra,rc],r6
  sto [ra,rd],r7
  sto [ra,re],r8
  shifti r3,6,r3
  sto [ra,r16],r4
  sto [ra],r3
  barb,br ,r1a
  end
proc2 ident
*
*    coded     v.r.saunders ---- dec 83.
*
savl equ 34
  msec 1
len res,128 6*64
* bbas,qbas,iky,jbas,kbas,pbas
blkin msec 3
g res,128 340*64
ijkl res,64 170*64
nwb res,64 64
scr205 msec 3
i res,128 340*4*64
ij res,64 340*6*64
pij res,64 340*6*64
*  ij,jk,kl, ik,jl,  il
g2 res,64 340*3*64
sav res,128 savl*64
ordv res,64 340
  msec 2
r3 equ 3*64
r4 equ 4*64
r5 equ 5*64
r6 equ 6*64
r7 equ 7*64
r8 equ 8*64
r9 equ 9*64
ra equ #a*64
rb equ #b*64
rc equ #c*64
rd equ #d*64
re equ #e*64
rf equ #f*64
r10 equ #10*64
r11 equ #11*64
r12 equ #12*64
r13 equ #13*64
*
r14 equ #14*64
r16 equ #16*64
r17 equ #17*64
r1a equ #1a*64
*
r20 equ #20*64
r21 equ #21*64
r22 equ #22*64
r23 equ #23*64
r24 equ #24*64
r25 equ #25*64
r26 equ #26*64
r27 equ #27*64
r28 equ #28*64
r29 equ #29*64
r2a equ #2a*64
r2b equ #2b*64
r2c equ #2c*64
r2d equ #2d*64
r2e equ #2e*64
r2f equ #2f*64
*
r30 equ #30*64
r31 equ #31*64
r32 equ #32*64
r33 equ #33*64
r34 equ #34*64
r35 equ #35*64
r36 equ #36*64
r37 equ #37*64
r38 equ #38*64
r39 equ #39*64
r3a equ #3a*64
r3b equ #3b*64
r3c equ #3c*64
r3d equ #3d*64
r3e equ #3e*64
r3f equ #3f*64
*
r40 equ #40*64
r41 equ #41*64
  entry proc2
proc2 ex rf,sav
  elen rf,savl
  ex r3,len
  elen r3,6
  swap r3,r14,rf
  ex r7,nwb
  lod [r7],r7
  ex r3,ijkl
  elen r3,170
  ex ra,i
  elen ra,170
  es r5,40
  es r6,#ff
*  unpack first 170  i
  linkv,ra
  shiftv,b [r3],r5,r3,
  andv,b [r3],r6,ra,
  shifti r7,6,r13
  es r8,170
  es r5,48
  addx ra,r13,rb
  ibxle,rel ,r8,sgm1,r7,
  pack r7,rb,rb
*  unpack first 170  j
sgm1 linkv,ra
  shiftv,b [r3],r5,rb,
  andv,b [rb],r6,rb,
  addx rb,r13,rc
  es r5,56
*  unpack first 170  k
  linkv,ra
  shiftv,b [r3],r5,rc,
  andv,b [rc],r6,rc,
  addx rc,r13,rd
*  unpack first 170  l
  andv,b [r3],r6,rd,
  ibxle,rel ,r7,sgm2,r8,
  subx r7,r8,r8
  ex r4,i+170*64
  es r5,8
  pack r8,r4,r4
*  unpack second 170  i
  linkv,ra
  shiftv,b [r3],r5,r4,
  andv,b [r4],r6,r4,
  addx r4,r13,r4
  es r5,16
*  unpack second 170  j
  linkv,ra
  shiftv,b [r3],r5,r4,
  andv,b [r4],r6,r4,
  addx r4,r13,r4
  es r5,24
*  unpack second 170  k
  linkv,ra
  shiftv,b [r3],r5,r4,
  andv,b [r4],r6,r4,
  addx r4,r13,r4
*  unpack second 170  l
  linkv,ra
  shiftv,b [r3],r14,r4,
  andv,b [r4],r6,r4,
sgm2 shifti r7,2,r3
  ex r4,pij
  pack r3,ra,ra
*  gather the  iky(i  j  k  l)
  vxtov [ra],r22,r4
*  generate iky(i)+j   iky(j)+k   iky(k)+l
  subx r3,r7,r6
  pack r6,r4,r4
  pack r6,rb,rb
  ex r8,ij
  pack r6,r8,r8
  addxv [r4],[rb],r8,
*  generate iky(i)+k   iky(j)+l
  shifti r6,6,re
  subx r6,r7,r12
  addx r8,re,r9
  pack r12,rc,rc
  pack r12,r9,r9
  addxv [r4],[rc],r9,
*  generate iky(i)+l
  shifti r12,6,r11
  addx r11,r9,r10
  pack r7,r10,r10
  pack r7,r4,r4
  pack r7,rd,rd
  addxv [r4],[rd],r10,
*  generate iky(k)+j  in  g2
  ex r26,g2
  pack r7,r26,r26
  addx r4,r11,r5
  addxv [r5],[rb],r26,
*  merge kj with jk
  ex r27,ordv
  pack r7,r27,r27
  pack r7,rc,r2f
  cmpge [rb],[r2f],r27
  addx r8,r13,r28
  maskv r28,r26,r28,r27
*  generate iky(l)+j  in  g2
  addx r5,r13,r29
  addxv [r29],[rb],r26,
*  merge lj with jl
  cmpge [rb],[rd],r27
  subx r10,r13,r2a
  maskv r2a,r26,r2a,r27
*  generate 2*g in g2
  ex r22,g
  pack r7,r22,r22
  addnv [r22],[r22],r26,
*  generate  gik
  pack r7,ra,ra
  cmpeq [ra],[r2f],r27
  addx r26,r13,r2b
  addx r2b,r13,r2c
  cmpeq [rb],[rd],r2c
  ior [r27],[r2c],[r27]
  maskv r26,r22,r2b,r27
*  generate  gil
  cmpne [rb],[r2f],r27
  maskv r22,r26,r2c,r27
*  gather the p elements
  addx r6,r6,r2f
  pack r2f,r8,r8
  vxtov [r8],r25,r4
*  generate  g2*pij
  mpysv [r26],[r4],r4,
*  generate  gil*pjk
  addx r4,r13,r33
  mpysv [r2c],[r33],r33,
*  generate  g2*pkl  and  gik*pik
  pack r12,r5,r5
  pack r12,r26,r27
  mpysv [r27],[r5],r5,
*  gather  the  qjk
  vxtov [r28],r21,ra
*  gather  the  qik,qjl,qil
  pack r6,r9,r2f
  vxtov [r2f],r21,rb
*  generate  gik*pjl  and  gil*pil
  addx r5,r11,r6
  pack r12,r2b,r2f
  mpysv [r2f],[r6],r6,
*  generate  gil*qjk
  mpysv [r2c],[ra],ra,
*  generate  gik*qik
  mpysv [r2b],[rb],rb,
*  generate  gik*qjl  and  gil*qil
  mpysv [r2f],[rc],rc,
*
*  jij=jij+g2*pkl
*  jkl=jkl+g2*pij
*
*  kjl=kjl+gik*pik
*  kjk=kjk+gil*pil
*  kil=kil+gil*pjk
*  kik=kik+gik*pjl
*
  es r40,64
  ex r35,sgm4
  subx r8,r40,r8
  subx r2a,r40,r2a
  subx r28,r40,r28
  subx r9,r40,r9
  subx r5,r40,r5
  subx r10,r40,r10
  subx r4,r40,r4
  subx r29,r40,r29
  subx r33,r40,r33
  subx r6,r40,r6
  subx rb,r40,rb
  subx rd,r40,rd
  subx ra,r40,ra
  subx rc,r40,rc
  addx r28,r13,r11
  addx r6,r13,r34
*  fetch  ij,kl,jl,jk,il,ik
sgm4 lod [r8,r7],r3
  lod [r11,r7],re
  lod [r2a,r7],r12
  lod [r28,r7],r13
  lod [r10,r7],r22
  lod [r9,r7],r2b
*  fetch  pkl,pij,qik,qil,qjk,qjl
  lod [r5,r7],r2d
  lod [r4,r7],r2e
  lod [rb,r7],r2f
  lod [rd,r7],r30
  lod [ra,r7],r31
  lod [rc,r7],r32
*  fetch  pik,pil,pjk,pjl
  lod [r29,r7],r3c
  lod [r34,r7],r3d
  lod [r33,r7],r3e
  lod [r6,r7],r3f
*  fetch  jij,jkl,bjl,bjk,bil,bik
  lod [r23,r3],r36
  lod [r23,re],r37
  lod [r20,r12],r38
  lod [r20,r13],r39
  lod [r20,r22],r3a
  lod [r20,r2b],r3b
*  fetch kjl,kjk,kil,kik
  lod [r24,r12],r40
  lod [r24,r13],r41
  lod [r24,r22],r26
  lod [r24,r2b],r27
*  add updates to  j  ,  b  and  k
  addn r36,r2d,r2d
  addn r37,r2e,r2e
  addn r38,r2f,r2f
  addn r39,r30,r30
  addn r3a,r31,r31
  addn r3b,r32,r32
  subn r40,r3c,r3c
  subn r41,r3d,r3d
  subn r26,r3e,r3e
  subn r27,r3f,r3f
  sto [r23,r3],r2d
  sto [r23,re],r2e
  sto [r20,r12],r2f
  sto [r20,r13],r30
  sto [r20,r22],r31
  sto [r20,r2b],r32
  sto [r24,r12],r3c
  sto [r24,r13],r3d
  sto [r24,r22],r3e
  sto [r24,r2b],r3f
  dbnz r7,[r35]
  swap rf,r14,
  barb,br ,r1a
*
*
  entry proc21
* call proc21(b,q,iky,j,k,p) to initialize proc2
proc21 lod [r17],r3
  lod [r17,r16],r4
  es rb,2
  lod [r17,rb],r5
  es rc,3
  lod [r17,rc],r6
  es rd,4
  lod [r17,rd],r7
  es re,5
  lod [r17,re],r8
  is r3,-64
  is r4,-64
  is r5,-64
  is r6,-64
  is r7,-64
  is r8,-64
  ex ra,len
  sto [ra],r3
  sto [ra,r16],r4
  sto [ra,rb],r5
  sto [ra,rc],r6
  sto [ra,rd],r7
  sto [ra,re],r8
  barb,br ,r1a
  end
 finis
_ENDIF
_IF(fps)
          $title dbuild
          $entry dbuild
          $insert '=tm.key'
"
" apal fock matrix builder.
"
" equvialent fortran is
"      subroutine dbuildf(fock,dmat)
"      implicit REAL  (a-h,p-w),integer * 4(i-n),logical *4 (o)
"      implicit character *8 (z),character *1 (x)
"      implicit character *4 (y)
"$insert vp.bigscf.apftn64(sizes)
"      REAL  fock(1),dmat(1)
"      common/craypk/ijkl(4,340)
"      common/mapper/iky(maxorb),ikyp(4,maxorb),i4096(maxorb)
"      common/shlt/tol,cutoff,icount,ic4,out
"      common/blkin/goutx(510),nword
"c ***
"       do 10 i=1,icount-1
"        i1=ijkl(1,i)
"        i2=ijkl(2,i)
"        i3=ijkl(3,i)
"        i4=ijkl(4,i)
"        if(i1.lt.i2) then
"          m=i1
"          i1=i2
"          i2=m
"        endif
"        if(i3.lt.i4) then
"          m=i3
"          i3=i4
"          i4=m
"        endif
"        if(i1-i3) 1,2,3
"1         m=i1
"          i1=i3
"          i3=m
"          m=i2
"          i2=i4
"          i4=m
"        goto 3
"2       if(i2.eq.i4) goto 1
"3       itr12=iky(i1)+i2
"        itr13=iky(i1)+i3
"        itr14=iky(i1)+i4
"        itr34=iky(i3)+i4
"        itr23=iky(max0(i2,i3))+min0(i2,i3)
"        itr24=iky(max0(i2,i4))+min0(i2,i4)
"        val=goutx(i)
"        val2=val+val
"        val4=val2+val2
"        val13=val
"        val14=val
"        if(i1.eq.i3 .or. i2.eq.i4) val13=val2
"        if(i2.eq.i3) val14=val2
"        f12 = val4*dmat(itr34) + fock(itr12)
"        fock(itr34) = val4*dmat(itr12) + fock(itr34)
"        fock(itr12) = f12
"        f23 = fock(itr23) - val14*dmat(itr14)
"        f14 = fock(itr14) - val14*dmat(itr23)
"        f13 = fock(itr13) - val13*dmat(itr24)
"        fock(itr24) = fock(itr24) - val13*dmat(itr13)
"        fock(itr23) = f23
"        fock(itr14) = f14
"        fock(itr13) = f13
"10      continue
"        icount=1
"        return
"        end
"
" spad definitions
"
          cap       $equ sp(0)
          fock      $equ sp(1)
          dmat      $equ sp(2)
          iky       $equ sp(3)
          goutx     $equ sp(4)
          i         $equ sp(5)
          j         $equ sp(6)
          k         $equ sp(7)
          l         $equ sp(8)
          ij        $equ sp(9)
          ik        $equ sp(10)
          il        $equ sp(11)
          jk        $equ sp(12)
          jl        $equ sp(13)
          kl        $equ sp(14)
          ijkl      $equ sp(15)
          icount    $equ sp(16)
"
          tmp       $equ sp(10)
          ki        $equ sp(10)
          li        $equ sp(11)
          kj        $equ sp(12)
          lj        $equ sp(13)
"
" data pad definitions
"
          xval      $equ dpx(0)
          tmpi      $equ dpx(1)
          tmpj      $equ dpx(2)
          two       $equ dpx(-1)
          four      $equ dpx(-2)
"
          yval      $equ dpy(0)
          tmpk      $equ dpy(1)
          tmpl      $equ dpy(2)
          val       $equ dpy(3)
          val2      $equ dpy(-1)
          val4      $equ dpy(-2)
"
          $pentf dbuild
"
" load the two arguments and set up spads for common blocks
"
          inc cap; setma
          inc cap; setma
          ldspi iky; db=iky_-1
          ldspi fock; db=md
          ldspi dmat; db=md
          dec fock
          dec dmat
          ldma; db=icount_
          ldspi ijkl; db=ijkl_-4
          ldtma; db=tm$two
          ldspi icount; db=md
          two<db; db=tm
          ldtma; db=tm$four
          ldspi goutx; db=goutx_-1
          four<db; db=tm
"
          dec icount
          bgt loop
home:     jmp finished
"
" now start the loop. 49 to 58 cycles per integral.
"
loop:     addi 4,ijkl; setma
"
" get the unordered i,j,k,l and put into partial order
"
          incma; dec icount
          incma; blt home
          incma; ldspi i; db=md; tmpi<db
                 ldspi j; db=md; tmpj<db; ficmp tmpi,md
                 ldspi k; db=md; tmpk<db; fapush
                 ldspi l; db=md; tmpl<db; ficmp tmpk,md
          inc goutx; setma; fapush;       fapush; bfge i_ge_j
swap_ij:  ldspi i; db=tmpj
          ldspi j; db=tmpi; bfge k_ge_l
swap_kl:  ldspi k; db=tmpl
          ldspi l; db=tmpk; br k_ge_l
i_ge_j:   nop
          bflt swap_kl
k_ge_l:   val<db; db=md; fmul two,md
"
" here have i>=j and k>=l. get triangle indices and then
" check for ij>=kl
"
          add# i,iky; setma; fmul four,md
          add# k,iky; setma; fmpush
          add# j,iky; setma; fmpush; val2<fm
          ldspi ij; db=md; tmpi<db;  val4<fm
          ldspi kl; db=md; tmpk<db
                    db=md; tmpj<db; add# l,iky; setma
          add j,ij
          add l,kl
          cmp ij,kl;db=md; tmpl<db
          mov i,tmp; ble ij_ge_kl
"
" swap ij,kl.
"
          mov k,i
          mov tmp,k;    fiadd zero,tmpi; tmpi<tmpk
          mov j,tmp;    fiadd zero,tmpj; tmpj<tmpl
          mov l,j;      fapush;          tmpk<fa
          mov tmp,l;                     tmpl<fa
"
" for rest of this section have (i>=j) >= (k>=l).
" do the coulomb parts in parallel with other address computations.
" triangle indices are used also for tests on i,j,k,l.
"
ij_ge_kl: add# kl,dmat; setma
          add# ij,dmat; setma               ; ficmp tmpj,tmpk
          nop            ; ldspi jl; db=tmpj; ficmp tmpj,tmpl
          fmul val4,md; add# ij,fock; setma ; nop
          fmul val4,md; add  kl,fock; setma ; inoload; bfge j_ge_k
j_lt_k:   fmpush                    ; fapush; ldspi jk; db=tmpk
          fmpush; fadd fm,md        ; nop   ; add j,jk
          fadd fm,md           ; bfge j_ge_l; add l,jl
j_lt_l:   fapush; add# ij,fock; setma; mi<fa
          add# kl,fock; setma;  mi<fa
          ldspi jl; db=tmpl
          add j,jl; br cont
"
j_ge_k:   fmpush                            ; ldspi jk; db=tmpj
          fmpush; fadd fm,md                ; add k,jk
          fadd fm,md                        ; add l,jl
j_ge_l:   fapush; add# ij,fock; setma; mi<fa
          add# kl,fock; setma;  mi<fa
"
" have jk,jl set up. set up ik,il, do the exchange parts and
" also test for (i=k .or. j=l) and also j=k
"
cont:     ldspi ik; db=tmpi
          add k,ik
          ldspi il; db=tmpi                ; ficmp tmpi,tmpk
          add l,il                         ; ficmp tmpj,tmpl
"
          add# ik,dmat; setma              ; ficmp tmpj,tmpk
"
                    add jl,dmat; setma     ; bfne i_ne_k;
                                inoload
"
"
i_eq_k:                      add il,dmat; setma; bfeq j_eq_l;
                                          inoload
"
j_eq_l:   xval<db;  db=md;
                                        add# jk,dmat; setma;
                                             fapush
"
          fmul xval,val2; add# jl,fock; setma;
                    xval<db; db=md
"
          fmpush;
                    fmul xval,val2; add ik,fock; setma; inoload;
                              xval<db; db=md; bfne j_ne_k
"
j_eq_k:   fmpush;
                    fmpush;
                              fmul xval,val2; add# jk,fock; setma;
                                        xval<db; db=md
"
          fsubr fm,md;
                    fmpush;
                              fmpush;
                                        fmul xval,val2; add# il,fock;
                                                            setma
"
          fapush;
                    fsubr fm,md;
                              fmpush;
                                        fmpush; br cont2
"
"
i_ne_k:                      add il,dmat; setma; bfeq j_eq_l;
                                          inoload
"
j_ne_l:   xval<db;  db=md;
                                        add# jk,dmat; setma;
                                             fapush
"
          fmul xval,val; add# jl,fock; setma;
                    xval<db; db=md
"
          fmpush;
                    fmul xval,val; add ik,fock; setma; inoload;
                              xval<db; db=md; bfeq j_eq_k
"
"
j_ne_k:   fmpush;
                    fmpush;
                              fmul xval,val; add# jk,fock; setma;
                                        xval<db; db=md
"
          fsubr fm,md;
                    fmpush;
                              fmpush;
                                        fmul xval,val; add# il,fock;
                                                            setma
"
          fapush;
                    fsubr fm,md;
                              fmpush;
                                        fmpush
"
"
cont2:    add# jl,fock; setma; mi<fa;
                    fapush;
                              fsubr fm,md;
                                        fmpush
"
          nop;
                    add# ik,fock; setma; mi<fa;
                              fapush;
                                        fsubr fm,md
"
          nop;
                    nop;
                              add# jk,fock; setma; mi<fa;
                                        fapush
"
          nop;
                    nop;
                              nop;
                                        add# il,fock; setma;
                                                      mi<fa
"
"
          jmp loop
"
"
finished: ldtma; db=tm$i1
          ldma; db=icount_-1
          incma; mi<db; db=tm
          dpx(0)<zero
          $pexit dbuild
"
          $psect craypk,md,ovl
"      common/craypk/ijkl(4,340)
ijkl_:    $rs 1360
"
          $psect mapper,md,ovl
"      common/mapper/iky(maxorb),ikyp(4,maxorb),i4096(maxorb)
iky_:     $rs 1
"
          $psect shlt,md,ovl
"      common/shlt/tol,cutoff,icount,ic4,out
tol_:     $rs 2
icount_:  $rs 1
          $psect blkin,md,ovl
"      common/blkin/goutx(510),nword
goutx_:   $rs 510
          $end
          $title gsup    " first card in routine gsup
          $entry gsup
          $ext upak4v
"
" apal equivalent to the following fortran. inner loop
" is about twice as fast as fortran and setup is substantially
" reduced.
"
"      subroutine gsup(h,p)
"      dimension h(1),p(1)
"      common/blkin/g(340),gij(170),mword
"      common/craypk/ij205(340),kl205(340)
"      dimension iij(1)
"      equivalence (iij(1),gij(1))
"      call upak4v(iij,ij205)
"       do 1 iw=1,mword
"       ij=ij205(iw)
"       kl=kl205(iw)
"      h(ij)=p(kl)*g(iw)+h(ij)
" 1    h(kl)=p(ij)*g(iw)+h(kl)
"      return
"      end
"
" register definitions
"
          cap    $equ sp(0)
          h      $equ sp(1)
          p      $equ sp(2)
          ij     $equ sp(3)
          kl     $equ sp(4)
          ij205  $equ sp(5)
          kl205  $equ sp(7)
          g      $equ sp(8)
          mword  $equ sp(9)
"
          ijd    $equ dpx(1)
          kld    $equ dpy(1)
          gg     $equ dpy(2)
"
         $pentf gsup,save_area,rap
          ldma; db=save_area
          incma; mov cap,cap; mi<spfn
         $callf upak4v,upak4v_parm      " call upak4v
"
          ldma; db=save_area+1
          nop
          nop
          ldspi cap; db=md
          inc cap; setma         " now get the arguments
          inc cap; setma
          ldspi ij205; db=craypk-1
          ldspi h; db=md
          ldspi p; db=md
          ldma; db=blkin+510        " get mword off end of /blkin/
          ldspi kl205; db=craypk+340-1
          ldspi g; db=blkin-1
          ldspi mword; db=md
          dec h
          dec p
"
" now start the loop. not pipelined 'cos of all the spad operations
"
          dec# mword
          jmplt finished           " mword le.0 so return
loop:     inc kl205; setma         " get kl
          inc g; setma             " get g
          inc ij205; setma         " get ij
          ldspi kl; kld<db; db=md  " save kl
          add# kl,p; setma;        " get p(kl)
              gg<md                "   save g
          ldspi ij; ijd<db; db=md  " save ij
          add# ij,p; setma;        " get p(ij)
              ficmp ijd,kld        "   see if ij=kl
          fmul gg,md; fapush;      " g*p(kl)
              add# ij,h; setma     "   get h(ij)
          fmpush
          fmul gg,md;              " g*p(ij)
              add kl,h; setma; inoload;    "   get h(kl)
                  bfeq ij_eq_kl    "     jump if ij=kl
          fadd fm,md; fmpush       " g*p(kl)+hij
          fapush; fmpush
          add# ij,h; setma; mi<fa; " write hij
              fadd fm,md           "   g*p(ij)+h(kl)
          fapush; dec mword
          add kl,h; setma; mi<fa;  " write h(kl)
              inoload; bgt loop
           jmp finished
"
ij_eq_kl: fadd fm,md               " ij=kl so just add on twice
          fapush
          fadd fm,fa
          fapush
          add# kl,h; setma; mi<fa  " write h(kl=ij)
          dec mword
          jmpgt loop
"
finished:  dpx(0)<zero
           $pexit gsup
"
"
         $psect .data.,md,cat
upak4v_parm: $word 2           " paramters for call to upak4v
             $word blkin+340
             $word craypk
"
save_area:   $rs 2
"
          $psect blkin,md,ovl    " common/blkin/
blkin:    $rs 511
"
          $psect craypk,md,ovl   " common/craypk/
craypk:   $rs 680
"
"
          $end     " last card in routine gsup
          $title gsupa    " first card in routine gsupa
          $entry gsupa
          $ext upak6v
"
" apal equivalent to the following fortran. inner loop
" is about twice as fast as fortran and setup is substantially
" reduced.
"
"      subroutine gsupa(h,p)
"      dimension h(1),p(1)
"      common/blkin/g(408),gij(102),mword
"      common/craypk/ij205(204),kl205(204)
"      dimension iij(1)
"      equivalence (iij(1),gij(1))
"      call upak6v(iij,ij205)
"       do 1 iw=1,mword
"       ij=ij205(iw)
"       kl=kl205(iw)
"      h(ij)=p(kl)*g(iw)+h(ij)
" 1    h(kl)=p(ij)*g(iw)+h(kl)
"      return
"      end
"
" register definitions
"
          cap    $equ sp(0)
          h      $equ sp(1)
          p      $equ sp(2)
          ij     $equ sp(3)
          kl     $equ sp(4)
          ij205  $equ sp(5)
          kl205  $equ sp(7)
          g      $equ sp(8)
          mword  $equ sp(9)
"
          ijd    $equ dpx(1)
          kld    $equ dpy(1)
          gg     $equ dpy(2)
"
         $pentf gsupa,save_area,rap
          ldma; db=save_area
          incma; mov cap,cap; mi<spfn
         $callf upak6v,upak6v_parm      " call upak6v
"
          ldma; db=save_area+1
          nop
          nop
          ldspi cap; db=md
          inc cap; setma         " now get the arguments
          inc cap; setma
          ldspi ij205; db=craypk-1
          ldspi h; db=md
          ldspi p; db=md
          ldma; db=blkin+510        " get mword off end of /blkin/
          ldspi kl205; db=craypk+204-1
          ldspi g; db=blkin-1
          ldspi mword; db=md
          dec h
          dec p
"
" now start the loop. not pipelined 'cos of all the spad operations
"
          dec# mword
          jmplt finished           " mword le.0 so return
loop:     inc kl205; setma         " get kl
          inc g; setma             " get g
          inc ij205; setma         " get ij
          ldspi kl; kld<db; db=md  " save kl
          add# kl,p; setma;        " get p(kl)
              gg<md                "   save g
          ldspi ij; ijd<db; db=md  " save ij
          add# ij,p; setma;        " get p(ij)
              ficmp ijd,kld        "   see if ij=kl
          fmul gg,md; fapush;      " g*p(kl)
              add# ij,h; setma     "   get h(ij)
          fmpush
          fmul gg,md;              " g*p(ij)
              add kl,h; setma; inoload;    "   get h(kl)
                  bfeq ij_eq_kl    "     jump if ij=kl
          fadd fm,md; fmpush       " g*p(kl)+hij
          fapush; fmpush
          add# ij,h; setma; mi<fa; " write hij
              fadd fm,md           "   g*p(ij)+h(kl)
          fapush; dec mword
          add kl,h; setma; mi<fa;  " write h(kl)
              inoload; bgt loop
           jmp finished
"
ij_eq_kl: fadd fm,md               " ij=kl so just add on twice
          fapush
          fadd fm,fa
          fapush
          add# kl,h; setma; mi<fa  " write h(kl=ij)
          dec mword
          jmpgt loop
"
finished:  dpx(0)<zero
           $pexit gsupa
"
"
         $psect .data.,md,cat
upak6v_parm: $word 2           " paramters for call to upak6v
             $word blkin+408
             $word craypk
"
save_area:   $rs 2
"
          $psect blkin,md,ovl    " common/blkin/
blkin:    $rs 511
"
          $psect craypk,md,ovl   " common/craypk/
craypk:   $rs 408
"
"
          $end     " last card in routine gsupa
          $title gsuppk    " first card in routine gsup
          $entry gsuppk
          $insert '=tm.key'
"
"      subroutine gsupp(h,p,g,ij,kl205,nij)
"      dimension h(1),p(1),g(1),kl205(1)
"      temp=0.0
"      pij = p(ij)
"cdir$ ivdep
"      do 20 iw = 1,nij
"          kl = kl205(iw)
"          gg = g(iw)
"          h(kl) = h(kl) + pij * gg
"20        temp = temp + p(kl) * gg
"      h(ij) = h(ij) + temp
"      return
"      end
"
"
" register definitions
"
          cap    $equ sp(0)
          h      $equ sp(1)
          p      $equ sp(2)
          ij     $equ sp(3)
          kl     $equ sp(4)
          kl205  $equ sp(5)
          g      $equ sp(6)
          oldkl  $equ sp(7)
"
          gg     $equ dpx(0)
          pij    $equ dpy(0)
          nij    $equ dpx(1)
          hkl    $equ dpy(1)
          temp   $equ dpy(2)
"
gsuppk: inc# cap; setma
        incma
        incma
        incma; ldspi h; db=md
        incma; ldspi p; db=md
        incma; ldspi g; db=md
        sldma; db=md
        ldspi kl205; db=md
        sldma; db=md
        ldspi ij; db=md
        sldtma; db=tm$i1
        nij < db; db=md; dec p
        dec h
        add# ij,p; setma
        temp<db; db=zero
"
        mov kl205,kl205; setma
        mov g,g; setma; pij<md
        nop
        ldspi kl; db=md
        add# kl,h; setma; gg<md
        add# kl,p; setma; fmul gg,pij
        fmpush
"
loop:   inc kl205; setma;                hkl<md; fmpush
        inc g; setma;                    fmul gg,md; fadd fm,hkl
        nop;                             fmpush; fapush; mov kl,oldkl
        ldspi kl; db=md;                 fmpush; hkl<fa; fisubr tm,nij
        add# kl,h; setma; gg<md;         fadd fm,temp
        add# kl,p; setma; fmul gg,pij;   fapush; nij < fa
        fmpush;                          temp<fa; add h,oldkl; setma;
                                           mi<db; db=hkl; bfgt loop
"
        nop " have to wait in case last kl=ij
        nop
        nop
        add# ij,h; setma
        nop
        nop
        fadd temp,md
        fapush
        add# ij,h; setma; mi<fa
finished:  dpx(0)<zero
"
"
        return
        $end     " last card in routine gsuppk
          $title k2
          $entry k2
          $insert '=tm.key'
"
" following fortran has the same functionality and
" operates similarly. note that on fps it is actually more
" efficient to pass j and k as separate matrices and isolate
" the nshell=1 case. the former would destory call compatibility
" with the vectorised xmp and c2 versions as is no advantage
" to a handwritten version.
"
" the apal source below may be improved by at least 30%
" by using the adder unit to do some of the logic and
" address generation. the inner loop of nshell is about 7 cycles
" slower than it need be! feel free to play ... this was only
" written in one day!
"
"      subroutine k1(l2,nshell,dm1,h2,p,iq,jq,kq,lq,gq,dm2,dm3,ic)
"      implicit real*8(a-h,o-z)
"      common/mapper/iky(256)
"      dimension h2(*),p(*),gq(*)
"      integer iq(*),jq(*),kq(*),lq(*)
"      icoff = l2*nshell
"cdir$ align
"      do 10 iw = 1,ic
"          i = iq(iw)
"          j = jq(iw)
"          k = kq(iw)
"          l = lq(iw)
"          ikyi = iky(i)
"          ikyk = iky(k)
"          ikyj = iky(j)
"          ij = ikyi + j
"          ik = ikyi + k
"          il = ikyi + l
"          jk = ikyj + k
"          jl = ikyj + l
"          kl = ikyk + l
"          if(j.ge.k) goto 1
"          jk = ikyk + j
"          if(j.ge.l) goto 1
"          jl = iky(l) + j
"1         continue
"          g = gq(iw)
"          xik = g
"          if(i.eq.k .or. j.eq.l) xik = xik+xik
"          xil = g
"          if(j.eq.k) xil = xil+xil
"c$dir scalar
"          do 15 ishell = 1,nshell
"              cij = h2(ij+icoff) + g*p(kl)
"              h2(kl+icoff) = h2(kl+icoff) + g*p(ij)
"              h2(ij+icoff) = cij
"              hik = h2(ik) + xik * p(jl)
"              hil = h2(il) + xil * p(jk)
"              hjk = h2(jk) + xil * p(il)
"              h2(jl) = h2(jl) + xik * p(ik)
"              h2(jk) = hjk
"              h2(il) = hil
"              h2(ik) = hik
"              ij = ij + l2
"              ik = ik + l2
"              il = il + l2
"              jk = jk + l2
"              jl = jl + l2
"              kl = kl + l2
"15    continue
"10    continue
"      return
"      end
"
"
" spad definitions
"
          cap       $equ sp(0)
          l2        $equ sp(1)
          nshell    $equ sp(2)
          coul      $equ sp(3)
          exch      $equ sp(4)
          p         $equ sp(5)
          iq        $equ sp(6)
          jq        $equ sp(7)
          kq        $equ sp(8)
          lq        $equ sp(9)
          gq        $equ sp(10)
          nint      $equ sp(11)
          iky       $equ sp(12)
          i         $equ sp(13)
          j         $equ sp(14)
          k         $equ sp(15)
          l         $equ sp(16)
          ij        $equ sp(17)
          ik        $equ sp(18)
          il        $equ sp(19)
          jk        $equ sp(20)
          jl        $equ sp(21)
          kl        $equ sp(22)
          ishell    $equ sp(23)
"
" data pad definitions
"
          two       $equ dpx(0)
          g         $equ dpx(1)
          hjk       $equ dpx(2)
          pjl       $equ dpx(3)
          pjk       $equ dpx(-1)
          pil       $equ dpx(-2)
          pik       $equ dpx(-3)
          xik       $equ dpy(0)
          xil       $equ dpy(1)
          hik       $equ dpy(2)
          hil       $equ dpy(3)
          $pentf k2
"
" load arguments and address of iky
"
          inc# cap; setma
          incma
          incma
          incma; ldspi l2; db<md
          incma; ldspi nshell; db<md
          incma; nop
          incma; ldspi coul; db=md
          incma; ldspi p; db=md
          incma; ldspi iq; db=md
          incma; ldspi jq; db=md
          incma; ldspi kq; db=md
          incma; ldspi lq; db=md
          incma; ldspi gq; db=md
          mov l2,l2; setma
          mov nshell,nshell; setma
          ldma; db=md
          ldspi l2; db=md
          ldspi nshell; db=md
          ldspi nint; db=md
"
" decrement array address for use as base + index
" and generate correct offset for coulomb matrix
"
          dec coul
          mov coul,exch
          mov nshell,ishell
init:     add l2,coul
          dec ishell
          bgt init
          ldtma; db=tm$two
          dec p
          dec iq; two<db; db=tm
          dec jq
          dec kq
          dec lq
          dec gq
          ldspi iky; db=iky_-1
"
" now start loop. check iteration count and then
" start forming ij,ik,il,jk,jl,kl. this could be done
" faster by using the data pads and the adder unit.
" get it working first!
"
loop:     dec nint
          jmplt finished
          inc iq; setma
          nop              ; inc jq; setma
          nop              ; nop                 ; inc kq; setma
          ldspi i; db=md   ; nop                 ; nop
          nop              ; ldspi j; db=md      ; nop
          nop              ; nop                 ; ldspi k; db=md
          add# i,iky; setma; nop                 ; nop
          nop              ; add# j,iky; setma   ; nop
          nop              ; nop                 ; add# k,iky; setma
          ldspi ij; db=md  ; nop                 ; nop
          nop              ; ldspi jk; db=md     ; nop
          nop              ; nop                 ; ldspi kl; db=md
"
" got i,j,k and iky of these ... get l, iky(l) and form triangles
"
          inc lq; setma
          nop                ; mov ij,ik
          nop                ; mov ij,il
          ldspi l; db=md
          add# l,iky; setma
          nop                ; add j,ij
          nop                ; add k,ik
          ldspi jl; db=md
          nop                ; add l,il
" compare j and k
          cmp j,k
          bgt j_lt_k
          mov jk,jl
          add k,jk
          add l,jl
          br j_done
"
j_lt_k:   cmp j,l
          bgt j_lt_l
          mov jk,jl
          add l,jl
          mov kl,jk
          add j,jk
          br j_done
"
j_lt_l:   add j,jl
          mov kl,jk
          add j,jk
"
j_done:   add l,kl
"
" now have all the triangle indices. get the integral
" and form xik,xil.
"
          inc gq; setma
          nop                        ; cmp i,k
          nop                        ; beq i_eq_k
          g<db; db=md; fmul two,md   ; cmp j,l
          fmpush; xik<g              ; beq j_eq_l
          fmpush; xil<g              ; cmp j,k
          nop                        ; beq j_eq_k
          br got_g
"
i_eq_k:   g<db; db=md; fmul two,md
          fmpush
j_eq_l:   fmpush; xil<g              ; cmp j,k
          xik<fm                     ; bne got_g
"
j_eq_k:   xil<fm
got_g:    mov nshell,ishell
"
" now loop over nshell ... do the coulomb part and pre-fetch the
" exchange density matrix elements
"
loopn:    add# kl,p; setma
"
          nop                          ; add# ij,p; setma
"
          nop                          ; nop     ; add# jl,p;setma
"
          fmul g,md; add# ij,coul;setma; nop
"
          fmpush                       ; fmul g,md; add# kl,coul; setma
"
          fmpush                       ; fmpush  ; pjl<db; db=md;
                                                   add# jk,p; setma
"
          fadd fm,md                   ; fmpush  ; add# il,p; setma
"
          fapush                       ; fadd fm,md; add# ik,p; setma
"
          add# ij,coul; setma; mi<fa   ; fapush  ; pjk<db; db=md
"
          nop                          ; add# kl,coul; setma; mi<fa;
                                                   pil<db; db=md
"
" now do the exchange bits
"
  add# ik,exch; setma; fmul xik,pjl ; pik<db; db=md
"
  fmpush ;
          add# il,exch; setma; fmul xil,pjk
"
  fmpush ;
          fmpush ;
                  add# jk,exch; setma; fmul xil,pil
"
  fadd fm,md ;
          fmpush ;
                  fmpush;
                         add# jl,exch; setma; fmul xik,pik
"
  fapush ;
          fadd fm,md ;
                  fmpush;
                         fmpush
"
  add# ik,exch; setma; mi<fa;
          fapush;
                  fadd fm,md;
                         fmpush
"
          add# il,exch; setma; mi<fa;
                  fapush;
                         fadd fm,md
"
                  add# jk,exch; setma; mi<fa;
                         fapush
"
                         add# jl,exch; setma; mi<fa
"
  dec ishell
  jmple loop
  add l2,ij
  add l2,ik
  add l2,il
  add l2,jk
  add l2,jl
  add l2,kl
  jmp loopn
"
finished:  dpx(0)<zero
          $pexit k2
"
          $psect mapper,md,ovl
"      common/mapper/iky(maxorb),ikyp(4,maxorb),i4096(maxorb)
iky_:     $rs 1
"
          $end
                  $title popcnt
                  $entry popcnt
                  $hfunc popcnt,popcnt/u,idum/u/inout
                  $psect .code.,ps ,cat
"
                  dpmask $equ dpy(1)
                  num    $equ dpx(0)
                  value  $equ dpx(1)
                  cnt    $equ sp(4)
                  adrs   $equ sp(5)
                  ptr    $equ sp(6)
"
"
"
                  $pentf popcnt
                  inc 0; setma
                  ldspi ptr; db=-64
                  sldma; db=mask
                  fadd zero,zero; ldspi adrs; db=md
                  fapush; mov adrs,adrs; setma
                   dpmask<md; num<fa
                   ldspi cnt; db=0
                    value<md; inc ptr; flsh md
                    fapush
                  fland dpmask,fa
                  flsh value;inc ptr
                  fapush
                  fland dpmask,fa; bfeq loop2
                  inc cnt
"
loop1:            flsh value;inc ptr
                  fapush;beq end
                  fland dpmask,fa; bfeq loop2
                  br loop1;inc cnt
"
loop2:            flsh value; inc ptr
                  fapush;beq end
                  fland dpmask,fa; bfeq loop2
                  br loop1;inc cnt
"
end:              fland dpmask,fa;bfeq end1
                  inc cnt
"
end1:             fapush
                  nop
                  bfeq end2;mov cnt,cnt
                  inc cnt
                  num<db;db=spfn;spdbr;return
"
end2:             num<db;db=spfn;spdbr;return
"
"
                  $pexit popcnt
"
"
pop:              $psect .data., md, cat
"
lcl:              $rs 1
                  $loc lcl + 0
"
mask:             $word 1
"
                  $end
                  $title leadz
                  $entry leadz
                  $hfunc leadz,leadz/u,idum/u/inout
                  $psect .code.,ps ,cat
"
                  dpmask $equ dpy(1)
                  num    $equ dpx(0)
                  value  $equ dpx(1)
                  cnt    $equ sp(4)
                  adrs   $equ sp(5)
                  ptr    $equ sp(6)
"
"
"
                  $pentf leadz
                  inc 0; setma
                  ldspi ptr; db=-64
                  sldma; db=mask
                  fadd zero,zero; ldspi adrs; db=md
                  fapush; mov adrs,adrs; setma
                  dpmask<md; num<fa
                  ldspi cnt; db=0
                  value<md; inc ptr; flsh md
                  fapush
                  fland dpmask,fa
                  flsh value;inc ptr
                  fapush
                  fland dpmask,fa; bfgt end
                  inc cnt
"
loop:             flsh value; inc ptr
                  fapush; beq end1
                  fland dpmask,fa; bfgt end
                  br loop; inc cnt
"
end:              mov cnt,cnt
                  num<db; db=spfn; spdbr; return
"
end1:             fland dpmask,fa; bfgt end2; mov cnt,cnt
                  inc cnt; fapush
                  nop
                  bfgt end2; mov cnt,cnt
                  inc cnt
end11:            num<db;db=spfn;spdbr;return
"
end2:             num<db; db=spfn; spdbr; return
"
                  $pexit leadz
"
"
pop:              $psect  .data., md, cat
"
lcl:              $rs 1
                  $loc lcl + 0
"
mask:             $word 1
"
                  $end
                  $title save_arg,'apal_ici.01'
"                  $insert'=ps.key'
"                  $insert'=sys.key'
"                  $insert'=tm.key'
                  $entry save_arg
                  $entry get_arg
                  $hsubr save_arg
                  $psect .code.,ps ,cat
"
save_arg:
                  addi# 1,0; setma
                  nop
                  nop
                  ldma; db=md
                  incma
                  mov 60,60; decma; mi<spfn; return
get_arg:
                  addi# 1,0; setma
                  nop
                  nop
                  ldma; db=md
                  nop
                  nop
                  ldspi 60; db=md; return
                  $psect .data.,md ,cat
save:             $rs       1
                  $end
          $title upak4v
          $entry upak4v
          $insert '=tm.key'
"
" apal equivalent to following fortran. apal runs twice as
" fast as opt(3) fortran.
"
"      subroutine upak4v(ix,intij)
"      word ix
"      dimension ix(1),intij(1)
"      int=1
"      do 1 i=1,170
"      intij(int  )=extract(ix(i),0,16)
"      intij(int+1)=extract(ix(i),16,16)
"      intij(int+2)=extract(ix(i),32,16)
"      intij(int+3)=extract(ix(i),48,16)
" 1    int=int+4
"      return
"      end
"
"
" definitions
"
           mask $equ dpx(0)
           save $equ dpy(0)
"
           cap  $equ sp(0)
           ix   $equ sp(1)
           i    $equ sp(2)
           int4 $equ sp(3)
           count $equ sp(4)
"
upak4v:    inc cap; setma
           inc cap; setma
           movi 4,int4
           ldspi ix; db=md
           ldspi i; db=md
           ldtma; db=tm$qw3msk
           sub int4,i
           mask<db; db=tm
           movi 170,count
"
           mov ix,ix; setma
           nop
           nop
           flshi -48,md
           flshi -32,md
           flshi -16,md; save<fa
           fland mask,fa; dec count
"
loop:      inc ix; setma;        fland mask,fa; ble done
           nop;                  save<fa; add int4,i; fland mask,md;
                                          setma; mi<db; db=save
           nop;                  save<fa; incma; mi<db; db=save; fapush
           flshi -48,md;         save<fa; incma; mi<db; db=save
           flshi -32,md;                  incma; mi<db; db=save
           flshi -16,md; save<fa
           fland mask,fa; dec count; br loop
"
done:      save<fa; add int4,i; fland mask,md;
                    setma; mi<db; db=save
           save<fa; incma; mi<db; db=save; fapush
           save<fa; incma; mi<db; db=save
                    incma; mi<db; db=save
           dpx(0)<zero
           return
          $end
          $title upak6v
          $entry upak6v
          $insert '=tm.key'
"
" apal equivalent to following fortran. apal runs twice as
" fast as opt(3) fortran.
"
"      subroutine upak6v(ix,intij)
"      word ix
"      dimension ix(1),intij(1)
"      int=1
"      do 1 i=1,102
"      intij(int  )=extract(ix(i),0,16)
"      intij(int+1)=extract(ix(i),16,16)
"      intij(int+2)=extract(ix(i),32,16)
"      intij(int+3)=extract(ix(i),48,16)
" 1    int=int+4
"      return
"      end
"
"
" definitions
"
           mask $equ dpx(0)
           save $equ dpy(0)
"
           cap  $equ sp(0)
           ix   $equ sp(1)
           i    $equ sp(2)
           int4 $equ sp(3)
           count $equ sp(4)
"
upak6v:    inc cap; setma
           inc cap; setma
           movi 4,int4
           ldspi ix; db=md
           ldspi i; db=md
           ldtma; db=tm$qw3msk
           sub int4,i
           mask<db; db=tm
           movi 102,count
"
           mov ix,ix; setma
           nop
           nop
           flshi -48,md
           flshi -32,md
           flshi -16,md; save<fa
           fland mask,fa; dec count
"
loop:      inc ix; setma;        fland mask,fa; ble done
           nop;                  save<fa; add int4,i; fland mask,md;
                                          setma; mi<db; db=save
           nop;                  save<fa; incma; mi<db; db=save; fapush
           flshi -48,md;         save<fa; incma; mi<db; db=save
           flshi -32,md;                  incma; mi<db; db=save
           flshi -16,md; save<fa
           fland mask,fa; dec count; br loop
"
done:      save<fa; add int4,i; fland mask,md;
                    setma; mi<db; db=save
           save<fa; incma; mi<db; db=save; fapush
           save<fa; incma; mi<db; db=save
                    incma; mi<db; db=save
           dpx(0)<zero
           return
          $end

"****** viindx          vector integer index     rel 1.1         jan 82
"
" v e c t o r   i n t e g e r   i n d e x
"
"  history:
"       ap120 original  feb 77          s. berkowitz (vindex)
"       x64 conversion  nov 80          g. lee
"       revision        jan 82          t. skinner
"                                       tar #3336
"                                       make indexing like fortran
"       integer version oct 82          r. bair
"
"  purpose:
"       to form a vector by using the elements of one vector as
"       the subscripts by which to select the elements from a
"       second vector.
"
"  calling sequence:
"       call viindx(a,b,j,c,k,n)
"
"  input parameters:
"       a               REAL  array
"                       source vector (values)
"       b               integer         array
"                       source vector (indices)
"       j               integer         scalar
"                       md element step for vector b
"       k               integer         scalar
"                       md element step for vector c
"       n               integer         scalar
"                       element count for vectors a, b and c
"
"  output parameters:
"       c               REAL  array
"                       result vector
"
"  subprograms used:
"       mth$pu_reslve   (psrom)
"
"  error conditions:
"       none
"
"  description:
"       c(m)=a( b(m) )
"
"  global registers modified:
"       none
"
"  scratch registers modified:
"       sp(1,3,5), dpx(0,1), dpy(0)
"       fa, md, spfn
"
"  speed (cycles/element):
"       best:    5
"       typical: 5
"       worst:   6
"       (procedure setup: 71 cycles)
"       (mlsp setup:      13 cycles)
"
"  ps size:
"       34
"
"  md size:
"       1
"
"
        $title viindx
        $radix d'8'
        $hsubr viindx /udc, a,b,j,c,k,n
        $entry mth$sp_viindx
        $insert '=tm.key'
"
        $pentf viindx /proc, glosav, aps3
        ldspi 16; db=hval(z'f0ffffff')  "exception disable mask
        inc 0; setma                    "begin s-pad parameter setup
        incma; ldspi 17;raps3           "read apstat3
        incma; and 16,17;ldaps3;db=spfn "disable exception interrupts
        incma; ldspi  0; db=md
        incma; ldspi  1; db=md
        incma; ldspi  2; db=md
               ldspi  3; db=md
               ldspi  4; db=md
               ldspi  5; db=md
        mov  2, 2; setma
        mov  4, 4; setma
        mov  5, 5; setma
        ldspi  2; db=md
        ldspi  4; db=md
        ldspi  5; db=md
        jsr mth$sp_viindx               "jsr to mlsp entry
        ldspi 16; db=hval(z'ff781fff')  "exception indicator mask
        ldspi 17;raps                   "read apstatus
        and 16,17;ldaps;db=spfn         "clear exception indicators
        $pexit viindx
"
        $psect glosav,md,ovl
glosav: $rs 1
        $psect
"
"       spad mnemonics
        a $equ 0
        b $equ 1
        j $equ 2
        c $equ 3
        k $equ 4
        n $equ 5
"
        $pentf mth$sp_viindx/nonstd     "mlsp entry
        dec a;                                  "initially decrement a,
                                                "so that indexing will
                                                "go from 1 to n
                 dpx(1)<spfn                    "save base address of a
                      mov b,b; setma            "fetch first b
        nop
        nop
        add j,b;setma;                              "1,fetch b
                fiadd dpx(1),md                     "1,compute address
        sub k,c;                                    "2,back up pointer
                fapush                              "2,push
        mov n,n;                                 "is element count = 0?
                fiadd zero,fa                       "3,fix address
        fiadd dpx(1),md;                            "4,compute address
          beq done                                  "exit if count = 0
        add j,b; setma;                             "1,get index
                fapush;                             "1,push
                        dpx<fa                      "1,save address
                        ldma; db=dpx                "2,fetch value
                fiadd zero,fa                       "3,fix address
        fiadd dpx(1),md                             "4,compute address
"here is loop
loop:   add j,b; setma;                             "1, get index
                fapush;                             "1, push
                        dpx<fa;                     "1, save address
                                dpy<md              "1, save value
                        ldma; db=dpx                "2, fetch value
                fiadd zero,fa;                      "3, fix address
                                dec n               "3, decrement count
        fiadd dpx(1),md;                            "4, compute address
                                add k,c; setma;mi<dpy; "4, store value
                                bne loop         "4, back if more values
done:   return                                        "exit
        $end

"****** vipk32           vect 32-bit integr pack rel 1.0         aug 81
"
"  v e c t o r   32 - b i t   i n t e g e r   p a c k
"
"  history:
"       164 original    dec 80          d. davis (vpk32)
"       integer version may 82          r. bair
"
"  purpose:
"       to pack sets of 2 integer words from an array into
"       two integers of single 64-bit words.
"
"  calling sequence:
"       call vipk32(a,i,c,k,n)
"
"  input parameters:
"       a               word (integer)  array
"                       source vector
"       i               integer         scalar
"                       md element step for vector a
"       k               integer         scalar
"                       md element step for vector c
"       n               integer         scalar
"                       element count for vectors a and c
"
"       output parameters:
"       c               word            array
"                       result vector
"
"
"  subprograms used:
"       myh$pu_reslve
"
"  error conditions:
"       none
"
"  description:
"       from an array, sets of  2  integer  words  are
"       packed into single 64 bit words.
"       the integers are packed from left to right.
"
"  global registers modified:
"       none
"
"  scratch registers modified:
"       sp(0,2,4),  dpx(0),  dpy(0),  fa, fm, md
"
"  speed (cycles/integer):
"       best:           2
"       typical:        2.25
"       worst:          3
"       (procedure setup:       80 cycles)
"       (mlsp setup:            10 cycles)
"
"  ps size:
"       46
"
"  md size:
"       1
"
"
"
        $title vipk32
        $radix d'8'
        $hsubr vipk32/udc, a,i,c,k,n
        $entry mth$sp_vipk32
        $insert '=tm.key'
"
        $pentf vipk32/proc, glosav, aps3        "procedure entry
        ldspi 16; db=hval(z'f0ffffff')  "exception disable mask
        inc 0; setma                    "begin s-pad parameter setup
        incma; ldspi 17;raps3           "read apstat3
        incma; and 16,17;ldaps3;db=spfn "disable exception interrupts
        incma; ldspi  0; db=md
        incma; ldspi  1; db=md
        incma; ldspi  2; db=md
               ldspi  3; db=md
               ldspi  4; db=md
        mov  1, 1; setma
        mov  3, 3; setma
        mov  4, 4; setma
        ldspi  1; db=md
        ldspi  3; db=md
        ldspi  4; db=md
        jsr mth$sp_vipk32                "jsr to mlsp entry
        fapush
        fapush
        ldspi 16; db=hval(z'ff781fff')  "exception indicator mask
        ldspi 17;raps                   "read apstatus
        and 16,17;ldaps;db=spfn         "clear exception indicators
        $pexit vipk32
"
        $psect glosav,md,ovl
glosav: $rs 1
        $psect
"
"s-pads
        a $equ 0                "source vector address
        i $equ 1                "source vector increment
        c $equ 2                "destination vector address
        k $equ 3                "destination vector increment
        n $equ 4                "destination vector length
"-----------------------------------------------------------------------
        $pentf mth$sp_vipk32/nonstd
        mov a,a;setma           "fetch first source vector
        add i,a;setma            "fetch 2nd source vector
        sub k,c                 "initialize destination
        dpy<md                  "save 1st vector
        dpx<md                  "save 2nd source vector
        add i,a;setma           "fetch 3rd source vector
        add i,a;setma           "fetch 4th source vector
        nop
        fiadd zero,md            "save 3rd source vector
"-----------------------------------------------------------------------
loop1:  add i,a;setma;          "fetch 2n+3 source vector
                fiadd zero,md   "save 2n+2 source vector
        add i,a;setma;          "fetch 2n+4 source vector
                dpy<fa;          "save 2n+1 saved source vector
                        iwrtdlr;dpx<dpy
                                "write left half from right half
                                "of 2n-1 saved word
                fadd;
                        dec n   "dec destination words
        fiadd zero,md;          "save 2n+3 source word
                dpx<fa;         "save 2n+2 saved word
                        add k,c;setma;mi<dpx;bne loop1
                                "store n concatenated word
"-----------------------------------------------------------------------
        return
        $end

"****** vipk16           vect 16-bit integr pack rel 1.0         aug 81
"
"  v e c t o r   16 - b i t   i n t e g e r   p a c k
"
"  history:
"       164 original    dec 80          d. davis (vpk16)
"       integer version may 82          r. bair
"
"  purpose:
"       to pack sets of 4 floating point words from an array into
"       four 16-bit integers of single 64-bit words.
"
"  calling sequence:
"       call vipk16(a,i,c,k,n)
"
"  input parameters:
"       a               word (integer)  array
"                       source vector
"       i               integer         scalar
"                       md element step for vector a
"       k               integer         scalar
"                       md element step for vector c
"       n               integer         scalar
"                       element count for vectors a and c
"
"       output parameters:
"       c               word            array
"                       result vector
"
"
"  subprograms used:
"       myh$pu_reslve
"
"  error conditions:
"       none
"
"  description:
"       from an array, sets of  4  floating  point  words  are
"       packed into single 64 bit words.
"       the 16-bit integers are packed from left to right.
"
"  global registers modified:
"       none
"
"  scratch registers modified:
"       sp(0,2,4,13,14),  dpx(0),  dpy(0,1),  fa, fm, md, tm
"
"  speed (cycles/integer):
"       best:           4
"       typical:        4.1
"       worst:          4.25
"       (procedure setup:       86 cycles)
"       (mlsp setup:            16 cycles)
"
"  ps size:
"       82
"
"  md size:
"       1
"
"
"
        $title vipk16
        $radix d'8'
        $hsubr vipk16/udc, a,i,c,k,n
        $entry mth$sp_vipk16
        $insert '=tm.key'
"
        $pentf vipk16/proc, glosav, aps3        "procedure entry
        ldspi 16; db=hval(z'f0ffffff')  "exception disable mask
        inc 0; setma                    "begin s-pad parameter setup
        incma; ldspi 17;raps3           "read apstat3
        incma; and 16,17;ldaps3;db=spfn "disable exception interrupts
        incma; ldspi  0; db=md
        incma; ldspi  1; db=md
        incma; ldspi  2; db=md
               ldspi  3; db=md
               ldspi  4; db=md
        mov  1, 1; setma
        mov  3, 3; setma
        mov  4, 4; setma
        ldspi  1; db=md
        ldspi  3; db=md
        ldspi  4; db=md
        jsr mth$sp_vipk16                "jsr to mlsp entry
        fapush
        fapush
        ldspi 16; db=hval(z'ff781fff')  "exception indicator mask
        ldspi 17;raps                   "read apstatus
        and 16,17;ldaps;db=spfn         "clear exception indicators
        $pexit vipk16
"
        $psect glosav,md,ovl
glosav: $rs 1
        $psect
"
"
"       s-pads
        a $equ 0                        "source vector address
        i $equ 1                        "source vector increment
        c $equ 2                        "destination vector address
        k $equ 3                        "destination vector increment
        n $equ 4                        "destination vector length
        frt8 $equ 15
        sxt $equ 16
"       data pads
        m $equ 1
        mask $equ tm$lpwmsk+60
"---------------------------------------------------------
        $pentf mth$sp_vipk16/nonstd
        mov a,a;setma                   "fetch i'st word
        ldtma;db=mask            "fetch mask  (0:47) = 0;  (48:63) = 1
        ldspi sxt;db=16.                 "load sp(sxt) = 16.
"----------------------------------------------------------------
f1:     dpy(m)<tm;                      "save mask
        add i,a;setma;                  "fetch 2'nd word
          fiadd zero,md                  "save 1'st word
        fadd;
          sub k,c               "initialize destination
        fland dpy(m),fa    "mask out only low 16 bit integer of 1'st wor
        ldspi frt8;db=48.;              "load sp(frt8) = 48.
          fadd
        add i,a;setma;          "fetch 3'rd word
          fiadd zero,md;                 "save 2'nd word
            dpx<fa         "save masked low integer(16 bits) of 1'st wor
        flsh dpx;mov frt8,frt8          "shift 1'st word to far left
        fland dpy(m),fa          "mask out only low integer of 2'nd word
        flor zero,fa                    "begin concatenation
        add i,a;setma;                  "fetch 4'th word
         fiadd zero,md;                  "save 3'rd word
            dpx<fa                    "save only low integer of 2'nd wor
        dpy<fa;                 "continue concatenation
          flsh dpx;sub sxt,frt8     "shift 2'nd word to next left intege
        fland dpy(m),fa               "mask out only low integer
        flor dpy,fa
"-----------------------------------------------------------
loop1:   add i,a;setma;          "fetch (4n+1)'th word
                        fiadd zero,md;dpx<fa   "save 4n'th word
                        flsh dpx;sub sxt,frt8;dpy<fa   "shift 4n-1 word
                        fland dpy(m),fa        "mask out low 4n word
                        flor dpy,fa     "concatenate 3 integers for n'th
                                                "result
        add i,a;setma;                         "fetch (4n+2)'th word
          fiadd zero,md;                 "save 4n+1 word
                        dpx<fa
                        dpy<fa;
                        flsh dpx;sub sxt,frt8   "shift 4n word
        fland dpy(m),fa                     "mask out 4n+1 low integer
        ldspi frt8;db=48.;                      "load sp(frt8) = 48.
                        flor dpy,fa     "concatenate 4 integers
        add i,a;setma;          "fetch (4n+3)'th word
          fiadd zero,md;                 "save 4n+3 word
            dpx<fa
        flsh dpx;mov frt8,frt8;         "shift 4n+1 word
                        dpy<fa  "save concatenated 4 integers
        fland dpy(m),fa         "mask low integer of 4n+2 word
        flor zero,fa;           "begin concatenation
                        add k,c;setma;mi<dpy    "store resultin n'th
                                                        "word
        add i,a;setma;                    "fetch 4n+4 word
         fiadd zero,md;                  "save 4n+3 word
            dpx<fa
        dpy<fa;
          flsh dpx;sub sxt,frt8         "shift 4n+2 word
        fland dpy(m),fa;                "mask 4n+3 low integer
                        dec n           "decrement count
        flor dpy,fa;                    "concatenate 2 words
                        bne loop1        "branch to loop if not done
"-----------------------------------------------------------
        return
        $end

"****** vipk8            vector 8-bit byte pack  rel 1.0         aug 81
"
"  v e c t o r   8 - b i t   b y t e   p a c k
"
"  history:
"       164 original    dec 80          d. davis (vpk8)
"       integer version may 82          r. bair
"
"  purpose:
"       to pack sets of 8 floating point words from an array into
"       eight 8-bit integers of single 64-bit words.
"
"  calling sequence:
"       call vipk8(a,i,c,k,n)
"
"  input parameters:
"       a               word (integer)  array
"                       source vector
"       i               integer         scalar
"                       md element step for vector a
"       k               integer         scalar
"                       md element step for vector c
"       n               integer         scalar
"                       element count for vectors a and c
"
"       output parameters:
"       c               word            array
"                       result vector
"
"
"  subprograms used:
"       myh$pu_reslve
"
"  error conditions:
"       none
"
"  description:
"       from an array, sets of  8  integer words  are
"       packed into single 64 bit words.
"       the 8-bit integers are packed from left to right.
"
"  global registers modified:
"       none
"
"  scratch registers modified:
"       sp(0,2,4,13,14),  dpx(0),  dpy(0,1),  fa, fm, md, tm
"
"  speed (cycles/integer):
"       best:           4
"       typical:        4
"       worst:          4.1
"       (procedure setup:       72 cycles)
"       (mlsp setup:            12 cycles)
"
"  ps size:
"       146
"
"  md size:
"       1
"
"
"
        $title vipk8
        $radix d'8'
        $hsubr vipk8/udc, a,i,c,k,n
        $entry mth$sp_vipk8
        $insert '=tm.key'
"
        $pentf vipk8/proc, glosav, aps3        "procedure entry
        ldspi 16; db=hval(z'f0ffffff')  "exception disable mask
        inc 0; setma                    "begin s-pad parameter setup
        incma; ldspi 17;raps3           "read apstat3
        incma; and 16,17;ldaps3;db=spfn "disable exception interrupts
        incma; ldspi  0; db=md
        incma; ldspi  1; db=md
        incma; ldspi  2; db=md
               ldspi  3; db=md
               ldspi  4; db=md
        mov  1, 1; setma
        mov  3, 3; setma
        mov  4, 4; setma
        ldspi  1; db=md
        ldspi  3; db=md
        ldspi  4; db=md
        jsr mth$sp_vipk8                 "jsr to mlsp entry
        fapush
        fapush
        ldspi 16; db=hval(z'ff781fff')  "exception indicator mask
        ldspi 17;raps                   "read apstatus
        and 16,17;ldaps;db=spfn         "clear exception indicators
        $pexit vipk8
"
        $psect glosav,md,ovl
glosav: $rs 1
        $psect
"
"       s-pads
        a $equ 0                        "source vector address
        i $equ 1                        "source vector increment
        c $equ 2                        "destination vector address
        k $equ 3                        "destination vector increment
        n $equ 4                        "destination vector length
        fvt6 $equ 15
        egt $equ 16
"       data pads
        m $equ 1
        mask $equ tm$lpwmsk+70
"---------------------------------------------------------
        $pentf mth$sp_vipk8/nonstd               "mlsp entry
         mov a,a;setma            "fetch 1'st vector
        ldtma;db=mask  "fetch mask       (0:55) = 0; (56:63) = 1
        ldspi egt;db=8.                  "load sp(egt) = 8.
"----------------------------------------------------------------
f1:     dpy(m)<tm;                      "save mask
        add i,a;setma;                  "fetch 2'nd word
          fiadd zero,md                  "save 1'st word
        fadd;
          sub k,c                          "initialize result address
        fland dpy(m),fa                "mask out low byte of 1'st word
        ldspi fvt6;db=56.;                      "load sp(fvt6) = 56.
          fadd
        add i,a;setma;                       "fetch 3'rd word
          fiadd zero,md;                 "save 2'nd word
            dpx<fa
        flsh dpx;mov fvt6,fvt6          "shift 1'st byte to far left
        fland dpy(m),fa         "mask out low byte of 2'nd word
        flor zero,fa            "start concatenation
        add i,a;setma;                  "fetch 4'th word
         fiadd zero,md;                  "save 3'rd word
            dpx<fa
        dpy<fa;
          flsh dpx;sub egt,fvt6         "shift 2'nd byte to next left
        fland dpy(m),fa         "mask out low byte of3'rd word
        flor dpy,fa             "concatenate 2 bytes
       add i,a;setma;                   "fetch 5'th word
          fiadd zero,md;                 "save 4'th word
            dpx<fa
        dpy<fa;
          flsh dpx;sub egt,fvt6         "shift 3'rd byte to next left
        fland dpy(m),fa         "mask out low byte of 4'th word
        flor dpy,fa                     "concatenate 3 bytes
        add i,a;setma;                  "fetch 6'th word
          fiadd zero,md;                 "save 5'th word
            dpx<fa
        flsh dpx;sub egt,fvt6;          "shift 4'th byte to next left
          dpy<fa
        fland dpy(m),fa         "mask out low byte of 5'th word
        flor dpy,fa             "concatenate 4 bytes
        add i,a;setma;          "fetch 7'th word
          fiadd zero,md;         "save 6'th word
            dpx<fa
        dpy<fa;
        flsh dpx;sub egt,fvt6           "shift 5'th byte to next left
        fland dpy(m),fa         "mask out low byte of 6'th word
        flor dpy,fa                     "concatenate 5 bytes
        add i,a;setma;          "fetch 8'th word
         fiadd zero,md;                  "save 7'th word
            dpx<fa
        dpy<fa;
          flsh dpx;sub egt,fvt6         "shift 6'th byte to next left
        fland dpy(m),fa         "mask out low byte of 7'th word
        flor dpy,fa;br loop1                     "concatenate 6 bytes
"----------------------------------------------------------------
out1:    return
"-----------------------------------------------------------
loop1:   add i,a;setma;          "fetch 8n+1 word
                        fiadd zero,md;dpx<fa "save 8n word
                        flsh dpx;sub egt,fvt6;dpy<fa
                                "shift 8n-1 byte to next left
                        fland dpy(m),fa "mask out low byte of 8n word
                        flor dpy,fa     "conctenate 7 bytes
        add i,a;setma;                  "fetch 8n+2 word
          fiadd zero,md;        "save 8n+1 word
                        dpx<fa
                        dpy<fa;
                        flsh dpx;sub egt,fvt6   "leave 8n byte as is
        fland dpy(m),fa                 "mask low byte of 8n+1 word
                        flor dpy,fa;    "concatenate 8 bytes
                        ldspi fvt6;db=56.       "load sp(fvt6) = 56.
        add i,a;setma;                  "fetch 8n+3 word
          fiadd zero,md;                 "save 8n+2 word
            dpx<fa
        flsh dpx;mov fvt6,fvt6;           "shift 8n+1 byte to far left
                        dpy<fa
        fland dpy(m),fa;               "mask out low byte of 8n+2 word
                        dec n           "decrement count
        flor zero,fa;           "begin concatenation
                        add k,c;setma;mi<dpy;beq out1
                                "store n'th result, branch out if done
        add i,a;setma;                  "fetch 8n+4 word
         fiadd zero,md;                  "save 8n+3 word
            dpx<fa
        dpy<fa;
          flsh dpx;sub egt,fvt6        "shift 8n+2 byte to next left
        fland dpy(m),fa         "mask out low byte of 8n+3 word
        flor dpy,fa             "concatenate 2 words
        add i,a;setma;          "fetch 8n+5 word
        fiadd zero,md;dpx<fa     "save 8n+4 word
        flsh dpx;sub egt,fvt6;dpy<fa    "shift 8n+3 byte to next left
        fland dpy(m),fa         "mask out low byte of 8n+4 word
        flor dpy,fa             "concatenate 3 bytes
        add i,a;setma;          "fetch 8n+6 word
          fiadd zero,md;         "save 8n+5 word
        dpx<fa
        dpy<fa;
        flsh dpx;sub egt,fvt6   "shift 8n+4 byte to next left
        fland dpy(m),fa         "mask out low byte of 8n+5 word
        flor dpy,fa             "concatenate 4 bytes
        add i,a;setma;          "fetch 8n+7 word
          fiadd zero,md;                 "save 8n+6 word
            dpx<fa
        flsh dpx;sub egt,fvt6;         "shift 8n+5 byte to next left
        dpy<fa
        fland dpy(m),fa               "mask out low byte of 8n+6 word
        flor dpy,fa                     "concatenate 5 bytes
        add i,a;setma;                  "fetch 8n+8 word
         fiadd zero,md;                  "save 8n+7 word
            dpx<fa
        dpy<fa;
          flsh dpx;sub egt,fvt6       "shift 8n+6 byte to next left
        fland dpy(m),fa               "mask out low byte of 8n+7 word
        flor dpy,fa;                    "concatenate 7 bytes
                        jmp loop1
"------------------------------------------------------------------
        $end

"****** viup32           vec 32-bt integr unpack rel 1.0         aug 81
"
"  v e c t o r   3 2 - b i t   i n t e g e r   u n p a c k
"
"  history:
"       164 original    dec 80          d. davis (vup32)
"       integer version may 82          r. bair
"
"  purpose:
"       to unpack two 32-bit integers from each single 64 bit word
"       of an array.
"
"  calling sequence:
"       call viup32(a,i,c,k,n)
"
"  input parameters:
"       a               word            array
"                       source vector
"       i               integer         scalar
"                       md element step for vector a
"       k               integer         scalar
"                       md element step for vector c
"       n               integer         scalar
"                       element count for vectors a and c
"
"       output parameters:
"       c               word (integer)  array
"                       result vector
"
"  subprograms used:
"       myh$pu_reslve
"
"  error conditions:
"       none
"
"  description:
"       two 32-bit integers from each 64-bit word of an array
"       are unpacked. the integers are unpacked from left to
"       right.
"
"  global registers modified:
"       none
"
"  scratch registers modified:
"       sp(0,2,4),  dpx(0),  dpy(0),  fa, fm, md
"
"  speed (cycles/integer):
"       best:           2
"       typical:        2.25
"       worst:          3
"       (procedure setup:       68 cycles)
"       (mlsp setup:            8 cycles)
"
"  ps size:
"       32
"
"  md size:
"       1
"
"
"
        $title viup32
        $radix d'8'
        $hsubr viup32/udc, a,i,c,k,n
        $entry mth$sp_viup32
        $insert '=tm.key'
"
        $pentf viup32/proc, glosav, aps3        "procedure entry
        ldspi 16; db=hval(z'f0ffffff')  "exception disable mask
        inc 0; setma                    "begin s-pad parameter setup
        incma; ldspi 17;raps3           "read apstat3
        incma; and 16,17;ldaps3;db=spfn "disable exception interrupts
        incma; ldspi  0; db=md
        incma; ldspi  1; db=md
               ldspi  2; db=md
               ldspi  3; db=md
               ldspi  4; db=md
        mov  1, 1; setma
        mov  3, 3; setma
        mov  4, 4; setma
        ldspi  1; db=md
        ldspi  3; db=md
        ldspi  4; db=md
        jsr mth$sp_viup32                "jsr to mlsp entry
        fapush
        fapush
        ldspi 16; db=hval(z'ff781fff')  "exception indicator mask
        ldspi 17;raps                   "read apstatus
        and 16,17;ldaps;db=spfn         "clear exception indicators
        $pexit viup32
"
        $psect glosav,md,ovl
glosav: $rs 1
        $psect
"
"
"s-pads
        a $equ 0                "source vector address
        i $equ 1                "source vector increment
        c $equ 2                "destination vector address
        k $equ 3                "destination vector increment
        n $equ 4                "source word count
"-----------------------------------------------------------------------
        $pentf mth$sp_viup32/nonstd
        mov a,a;setma                   "fetch first source vector
        sub k,c;dpx<zero                 "initialize c
        dpy<zero
"-----------------------------------------------------------------------
                dpy<md;iwrtdrl          "save left half of first into
                                        "right half
        add i,a;setma;                  "fetch 2nd source vector
                dpx<md;iwrtr          "save right half of first into
                                        "right half
                fiadd zero,dpy         "save left half of first
                fiadd zero,dpx         "save right half of first
"-----------------------------------------------------------------------
loop:           dpy<md;iwrtdrl;          "save left half to right half
                                                "of n+1
                        add k,c;setma;mi<fa
                                        "store left half of n
        add i,a;setma;                  "fetch  source vector n+2
                dpx<md;iwrtr           "save right half of n+1
                fiadd zero,dpy;        "save left half of n+1
                        dec n           "decrement count
                fiadd zero,dpx;        "save right half of n+1
                        add k,c;setma;mi<fa;bne loop
                                        "store packed n'th word
"-----------------------------------------------------------------------
        return
        $end

"****** viup16           vec 16-bt integr unpack rel 1.0         aug 81
"
"  v e c t o r   1 6 - b i t   i n t e g e r   u n p a c k
"
"  history:
"       164 original    dec 80          d. davis (vup16)
"       integer version may 82          r. bair
"
"  purpose:
"       to unpack four 16-bit integers from each single 64 bit word
"         of an array.
"
"  calling sequence:
"       call viup16(a,i,c,k,n)
"
"  input parameters:
"       a               word            array
"                       source vector
"       i               integer         scalar
"                       md element step for vector a
"       k               integer         scalar
"                       md element step for vector c
"       n               integer         scalar
"                       element count for vectors a and c
"
"       output parameters:
"       c               word (integer)  array
"                       result vector
"
"  subprograms used:
"       myh$pu_reslve
"
"  error conditions:
"       none
"
"  description:
"       four 16-bit integers from each 64-bit word of an array
"       are unpacked. the integers are unpacked from left to
"       right.
"
"  global registers modified:
"       none
"
"  scratch registers modified:
"       sp(0,2,4,13,14),  dpx(0,1,2),  dpy(0,1),  fa, fm, md, tm
"
"  speed (cycles/integer):
"       best:           3
"       typical:        3.1
"       worst:          3.25
"       (procedure setup:       79 cycles)
"       (mlsp setup:            19 cycles)
"
"  ps size:
"       51
"
"  md size:
"       1
"
"
"
        $title viup16
        $radix d'8'
        $hsubr viup16/udc, a,i,c,k,n
        $entry mth$sp_viup16
        $insert '=tm.key'
"
        $pentf viup16/proc, glosav, aps3        "procedure entry
        ldspi 16; db=hval(z'f0ffffff')  "exception disable mask
        inc 0; setma                    "begin s-pad parameter setup
        incma; ldspi 17;raps3           "read apstat3
        incma; and 16,17;ldaps3;db=spfn "disable exception interrupts
        incma; ldspi  0; db=md
        incma; ldspi  1; db=md
               ldspi  2; db=md
               ldspi  3; db=md
               ldspi  4; db=md
        mov  1, 1; setma
        mov  3, 3; setma
        mov  4, 4; setma
        ldspi  1; db=md
        ldspi  3; db=md
        ldspi  4; db=md
        jsr mth$sp_viup16                "jsr to mlsp entry
        fapush
        fapush
        ldspi 16; db=hval(z'ff781fff')  "exception indicator mask
        ldspi 17;raps                   "read apstatus
        and 16,17;ldaps;db=spfn         "clear exception indicators
        $pexit viup16
"
        $psect glosav,md,ovl
glosav: $rs 1
        $psect
"
"
"       s-pads
        a $equ 0                "source vector address
        i $equ 1                "source vector increment
        c $equ 2                "destination vector address
        k $equ 3                "destination vector increment
        n $equ 4                "source vector length
        frt8 $equ 15
        sxt $equ 16
"   data pads
        m $equ 0
        mask $equ tm$lpwmsk+60
"-----------------------------------------------------------------------
        $pentf mth$sp_viup16/nonstd
        ldtma;db=mask   "fetch mask     (0:47) = 0; (48:63) = 1
        sub k,c         "initialize output address
        mov a,a;setma;dpy(m)<tm         "fetch 1'st word, save mask
        ldspi sxt;db=16.                "load sp(sxt) = 16.
        ldspi frt8;db=-48.              "load sp(frt8) = 48.
         flsh md;mov frt8,frt8;dpy(1)<md       "shift left most 16 bit
                                                "byte to far right
        flsh dpy(1);add sxt,frt8        "shift 2'nd left byte to
                                        "far right
        fland dpy(m),fa                 "mask out only 16 bit byte
        fland dpy(m),fa                 "mask out only 16 bit byte
        fiadd zero,fa                   "save left most byte
        fiadd zero,fa;                  "save next left byte
           add i,a;setma                "fetch 2'nd word
        flsh dpy(1);add sxt,frt8;dpx(1)<fa    "shift 3'rd left byte to
                                                "far right
        flsh dpy(1);add sxt,frt8;dpx(2)<fa      "right most byte,
                                                "leave as is
        fland dpy(m),fa;add k,c;setma;mi<dpx(1) "mask out only byte
                                                "store left most byte
        fland dpy(m),fa;ldspi frt8;db=-48.      "mask out only byte
                                        "load sp(frt8) = 48.
        fiadd zero,fa;dec n     "save 3'rd byte
        fiadd zero,fa;add k,c;setma;mi<dpx(2);  "save 4'th byte
                                        "store 2'nd left byte
                                beq done        "done if n=1
"--------------------------------------------------
loop:   flsh md;mov frt8,frt8;dpy(1)<md;              "shift left most
                                                "byte to far right
                        dpx(1)<fa
        flsh dpy(1);add sxt,frt8;       "shift next left byte to
                                                "far right
                        dpx(2)<fa
        fland dpy(m),fa;                "mask out only byte
                        add k,c;setma;mi<dpx(1)      "store 3'rd byte
                                                "0f n'th word
        fland dpy(m),fa                     "mask out only byte
        fiadd zero,fa;                          "save left most byte
                                                "of (n+1)'th word
                        add k,c;setma;mi<dpx(2)    "store 4'th byte
                                                "of n'th word
        fiadd zero,fa;                  "save 2'nd left byte n+1 word
           add i,a;setma                "fetch (n+2)'th word
        flsh dpy(1);add sxt,frt8;dpx(1)<fa   "shift 3'rd byte far right
        flsh dpy(1);add sxt,frt8;dpx(2)<fa     "4'th byte, leave as is
        fland dpy(m),fa;add k,c;setma;mi<dpx(1)    "mask out only byte
                                        "store 4'th byte of n'th word
        fland dpy(m),fa;ldspi frt8;db=-48.      "mask out only byte
                                                "load sp(frt8) = 48.
        fiadd zero,fa;              "save 3'rd byte of n+1 word
                        dec n           "decrement count
        fiadd zero,fa;add k,c;setma;mi<dpx(2); "save 4'th byte of n+1
                                 "word, store 2'nd byte of n+1 word
                        bne loop        "branch to loop if not done
"--------------------------------------------------
done:   fadd;add k,c;setma;mi<fa        "store 3'rd byte of n'th word
        add k,c;setma;mi<fa;return      "store 4'th byte of n'th word
        $end

"****** viup8            vect 8-bit byte unpack  rel 1.0         aug 81
"
"  v e c t o r   8 - b i t   b y t e   u n p a c k
"
"  history:
"       164 original    dec 80          d. davis (vup8)
"       integer version may 82          r. bair
"
"  purpose:
"       to unpack eight 8-bit bytes from each single 64 bit word
"         of an array.
"
"  calling sequence:
"       call viup8(a,i,c,k,n)
"
"  input parameters:
"       a               word            array
"                       source vector
"       i               integer         scalar
"                       md element step for vector a
"       k               integer         scalar
"                       md element step for vector c
"       n               integer         scalar
"                       element count for vectors a and c
"
"       output parameters:
"       c               word (integer)  array
"                       result vector
"
"  subprograms used:
"       myh$pu_reslve
"
"  error conditions:
"       none
"
"  description:
"       eight 8-bit bytes from each 64-bit word of an array
"       are unpacked. the bytes are unpacked from left to
"       right.
"
"  global registers modified:
"       none
"
"  scratch registers modified:
"       sp(0,2,4,13,14),  dpx(0,1,2),  dpy(0,1),  fa, fm, md, tm
"
"  speed (cycles/integer):
"       best:           3
"       typical:        3
"       worst:          3.1
"       (procedure setup:       70 cycles)
"       (mlsp setup:            10 cycles)
"
"  ps size:
"       73
"
"  md size:
"       1
"
"
"
        $title viup8
        $radix d'8'
        $hsubr viup8/udc, a,i,c,k,n
        $entry mth$sp_viup8
        $insert '=tm.key'
"
        $pentf viup8/proc, glosav, aps3        "procedure entry
        ldspi 16; db=hval(z'f0ffffff')  "exception disable mask
        inc 0; setma                    "begin s-pad parameter setup
        incma; ldspi 17;raps3           "read apstat3
        incma; and 16,17;ldaps3;db=spfn "disable exception interrupts
        incma; ldspi  0; db=md
        incma; ldspi  1; db=md
               ldspi  2; db=md
               ldspi  3; db=md
               ldspi  4; db=md
        mov  1, 1; setma
        mov  3, 3; setma
        mov  4, 4; setma
        ldspi  1; db=md
        ldspi  3; db=md
        ldspi  4; db=md
        jsr mth$sp_viup8                 "jsr to mlsp entry
        fapush
        fapush
        ldspi 16; db=hval(z'ff781fff')  "exception indicator mask
        ldspi 17;raps                   "read apstatus
        and 16,17;ldaps;db=spfn         "clear exception indicators
        $pexit viup8
"
        $psect glosav,md,ovl
glosav: $rs 1
        $psect
"
"
"  s-pad values
        a $equ 0                "source vector address
        i $equ 1                "source vector increment
        c $equ 2                "destination vector address
        k $equ 3                "destination vector increment
        n $equ 4                "source vector length
        fvt6 $equ 15
        egt $equ 16
"   data pads
        m $equ 0
        mask $equ tm$lpwmsk+70
"----------------------------------------------------------
        $pentf mth$sp_viup8/nonstd
        ldtma;db=mask;sub k,c         "fetch mask (0:55)=0; (56:63)=1
        mov a,a;setma           "fetch 1'st word
        ldspi fvt6;db=-56.              "load sp(fvt6) = 56.
        ldspi egt;db=8.         "load sp(egt) = 8.
"--------------------------------------------------
        flsh md;mov fvt6,fvt6;dpy(1)<md     "shift left most byte to
                                                "far right
        flsh dpy(1);add egt,fvt6;dpy(m)<tm       "shift 2'nd left byte
                                                "to far right
        fland dpy(m),fa                "maskout only 1'st byte
        fland dpy(m),fa                "maskout only 2'nd byte
        fiadd zero,fa                           "save 1'st byte
        fiadd zero,fa                           "save 2'nd byte
        flsh dpy(1);add egt,fvt6;dpx(1)<fa      "shift 3'rd left byte
                                                "to far right
        flsh dpy(1);add egt,fvt6;dpx(2)<fa      "shift 4'th left byte
                                                 "to far right
        fland dpy(m),fa;add k,c;setma;mi<dpx(1)    "mask out only 3'rd
                                                        "byte
        fland dpy(m),fa                     "mask out only 4'th byte
        fiadd zero,fa                           "save 3'rd byte
        fiadd zero,fa;add k,c;setma;mi<dpx(2)   "save 4'th byte
        flsh dpy(1);add egt,fvt6;       "shift 5'th byte to far right
               dpx(1)<fa
        flsh dpy(1);add egt,fvt6;       "shift 6'th byte to far right
               dpx(2)<fa
        fland dpy(m),fa;                "mask out only 5'th byte
                add k,c;setma;mi<dpx(1)         "store 3'th byte
        fland dpy(m),fa                 "mask out only 6'th byte
        fiadd zero,fa;                  "save 5'th byte
                add k,c;setma;mi<dpx(2)         "store 4'th byte
        fiadd zero,fa;                  "save 6'th byte
                add i,a;setma           "fetch 2'nd word
        flsh dpy(1);add egt,fvt6;dpx(1)<fa      "shift 7'th byte to far
                                                        "right
        flsh dpy(1);add egt,fvt6;dpx(2)<fa      "shift 8'th byte to
                                                "far right
        fland dpy(m),fa;add k,c;setma;mi<dpx(1)      "mask out only
                                 "7'th byte, store 5'th byte
        fland dpy(m),fa;ldspi fvt6;db=-56.      "mask out only
                                                        "8'th byte
        fiadd zero,fa                   "save 7'th byte
        fiadd zero,fa;add k,c;setma;mi<dpx(2);br loop "save 8'th byte
                                                "store 6'th byte
"------------------------------------------------------------
out:    return
"--------------------------------------------------------------
loop:   flsh md;mov fvt6,fvt6;dpy(1)<md;        "shift left most byte
                                       "of n+1 word to far right
                                dpx(1)<fa
        flsh dpy(1);add egt,fvt6;                 "shift 2'nd byte
                                        "of n+1 word to far right
                                dpx(2)<fa
        fland dpy(m),fa;         "mask out only 1'st byte of n+1 word
                                add k,c;setma;mi<dpx(1)
                                        "store 7'th byte of n'th word
        fland dpy(m),fa;         "mask out only 2'nd byte of n+1 word
                                dec n   "decrement count
        fiadd zero,fa;                    "save 1'st byte
                                add k,c;setma;mi<dpx(2);beq out
                                        "store 8'th byte of n'th word
        fiadd zero,fa                   "save 2'nd byte
        flsh dpy(1);add egt,fvt6;dpx(1)<fa      "shift 3'rd byte
                                        "of n+1 word to far right
        flsh dpy(1);add egt,fvt6;dpx(2)<fa      "shift 4'th byte
                                        "of n+1 word to far right
        fland dpy(m),fa;add k,c;setma;mi<dpx(1)        "mask out only
                  "3'rd byte of n+1 word, store 1'st byte of n+1 word
        fland dpy(m),fa         "mask out only 4'th byte of n+1 word
        fiadd zero,fa                 "save 3'rd byte of n+1 word
        fiadd zero,fa;add k,c;setma;mi<dpx(2)   "save 4'th word
                                        "store 2'nd byte of n+1 word
        flsh dpy(1);add egt,fvt6;       "shift 5'th byte of n+1 word
                                                "to far right
                dpx(1)<fa
        flsh dpy(1);add egt,fvt6;       "shift 6'th byte of n+1 word
                                                        "to far right
                dpx(2)<fa
        fland dpy(m),fa;        "mask out only 5'th byte of n+1 word
                add k,c;setma;mi<dpx(1)   "store 3'rd byte of n+1 byte
        fland dpy(m),fa         "mask out 6'th byte of n+1 word
        fiadd zero,fa;                  "save 5'th byte of n+1 word
                add k,c;setma;mi<dpx(2)   "store 4'th byte of n+1 word
        fiadd zero,fa;                  "save 6'th byte of n+1 word
                add i,a;setma           "fetch n+1 word
        flsh dpy(1);add egt,fvt6;dpx(1)<fa      "shift 7'th byte of n+1
        flsh dpy(1);add egt,fvt6;dpx(2)<fa      "shift 8'th byte of n+1
        fland dpy(m),fa;add k,c;setma;mi<dpx(1) "mask out 7'th byte
                                "of n+1 word, store 5'th byte of n+1
        fland dpy(m),fa;ldspi fvt6;db=-56.      "mask out 8'th byte
                                "of n+1 word, load sp(fvt6) = 56.
        fiadd zero,fa;add k,c;setma;mi<dpx(2)   "save 7'th byte of n+1,
                                                "store 6'th byte of n+1
        fiadd zero,fa;jmp loop  "save 8'th byte of n+1, jump to loop
"----------------------------------------------------------------
        $end

"****** vius32          vec 32-bt int unpak sgn rel 1.0         aug 81
"
"  v e c t o r   32 - b i t   i n t e g e r   u n p a c k   s i g n e d
"
"  history:
"       164 original    dec 80          d. davis (vups32)
"       integer version may 82          r. bair
"
"  purpose:
"       to unpack two 32-bit integers signed from each single
"       64 bit word of an array.
"
"  calling sequence:
"       call vius32(a,i,c,k,n)
"
"  input parameters:
"       a               word            array
"                       source vector
"       i               integer         scalar
"                       md element step for vector a
"       k               integer         scalar
"                       md element step for vector c
"       n               integer         scalar
"                       element count for vectors a and c
"
"       output parameters:
"       c               word (integer)  array
"                       result vector
"
"  subprograms used:
"       myh$pu_reslve
"
"  error conditions:
"       none
"
"  description:
"       two 32-bit integers from each 64-bit word of an array
"       are unpacked with sign extend. the integers are unpacked
"       from left to right.
"
"  global registers modified:
"       none
"
"  scratch registers modified:
"      sp(0,2,4,5,15),  dpx(0,1,2,3),  dpy(0,1),  fa, fm, md, tm
"
"  speed (cycles/integer):
"       best:           3
"       typical:        3
"       worst:          3
"       (procedure setup:       75 cycles)
"       (mlsp setup:            15 cycles)
"
"  ps size:
"       78
"
"  md size:
"       1
"
"
"
        $title vius32
        $radix d'8'
        $hsubr vius32/udc, a,i,c,k,n
        $entry mth$sp_vius32
        $insert '=tm.key'
"
        $pentf vius32/proc, glosav, aps3        "procedure entry
        ldspi 16; db=hval(z'f0ffffff')  "exception disable mask
        inc 0; setma                    "begin s-pad parameter setup
        incma; ldspi 17;raps3           "read apstat3
        incma; and 16,17;ldaps3;db=spfn "disable exception interrupts
        incma; ldspi  0; db=md
        incma; ldspi  1; db=md
               ldspi  2; db=md
               ldspi  3; db=md
               ldspi  4; db=md
        mov  1, 1; setma
        mov  3, 3; setma
        mov  4, 4; setma
        ldspi  1; db=md
        ldspi  3; db=md
        ldspi  4; db=md
        jsr mth$sp_vius32               "jsr to mlsp entry
        fapush
        fapush
        ldspi 16; db=hval(z'ff781fff')  "exception indicator mask
        ldspi 17;raps                   "read apstatus
        and 16,17;ldaps;db=spfn         "clear exception indicators
        $pexit vius32
"
        $psect glosav,md,ovl
glosav: $rs 1
        $psect
"
"
        a = 0
        i = 1
        c = 2
        k = 3
        n = 4
"
        lo = 5
        m = 3
"
        $pentf mth$sp_vius32/nonstd
        ldtma;db=tm$lpwmsk+53           "fetch mask (0-42)=0: (43-63)=1
        ldspi 17;db=40          "load 32 into sp(17)
        dpx(m)<tm;mov a,a;setma         "save mask, fetch first word
        flsh dpx(m);mov 17,17           "shift mask left 32 places
        fadd;sldtma;db=tm$slide1+40   "mask (0-31)=0: (32)=1: (33-63)=0
        add i,a;setma;dpy<zero           "fetch 2'nd word
        dpx(m)<fa;ldspi lo;dpy<md;iwrtdrl      "save right half of 1'st
                                      "word in lo and left half in dpy
                                        "save sign extend
        fland tm,dpy;sub k,c            "test for sign of left half
                                                "initialize c
        fadd
                        mov lo,lo;dpx<spfn;     "put right half of
                                                "1'st word in dpx
                        flor dpy,dpx(m)       "mask sign extend in
                                                "left half 1'st in case
        add i,a;setma;                          "fetch 3'rd word
                        bfeq pos1;      "branch to pos1 if 1'st left>0
                        dpx(2)<dpy;fadd      "move 1'st left to dpx(2)
                 ldspi lo;dpy<md;        "save right half of 2'nd word
                iwrtdrl;          "in lo and left half in right of dpy
                        dpx(2)<fa     "save neg left of 1'st in dpx(2)
cont1:          fland tm,dpy;  "test for sign of left half of 2'nd word
                        dpx(1)<dpx        "transfer right half of 1'st
                                                "into dpx(1)
                        fiadd zero,dpx(2);   "save left half of 1'st
                        br loop                 "branch to loop
"-------------------------------------------------------------------
pos1:            ldspi lo;dpy<md;        "save right half of 2'nd word
                iwrtdrl;br cont1   "in lo and left half in right of dpy
"--------------------------------------------------------------
loop:                   mov lo,lo;dpx<spfn;    "put right half of (n+1)
                                                "word in dpx
                        flor dpy,dpx(m)      "mask sign extend into
                                               "left of (n+1) in case
        add i,a;setma;                  "fetch (n+3)'th word
                        bfeq pos;     "branch to pos if (n+1) left > 0
                        dpx(2)<dpy;         "move (n+1) left to dpx(2)
                                fadd;dpy(1)<fa  "save (n) right
                ldspi lo;dpy<md;        "save right half of (n+2) word
                iwrtdrl;          "in lo and left half in right of dpy
                        dpx(2)<fa               "save neg left (n+1)
cont:                           fiadd zero,dpx(1); "save right half (n)
                                add k,c;setma;mi<dpy(1)
                        "store left half (n)
                fland tm,dpy;           "test sign of left half (n+2)
                        dpx(1)<dpx;             "move right half (n+1)
                                                "in dpx(1)
                                fadd;dec n      "decrement count
                        fiadd zero,dpx(2);           "save left (n+1)
                                add k,c;setma;mi<fa;
                                        "store right half (n)
                                bne loop    "branch to loop if not done
"-----------------------------------------------------------
        return
"-------------------------------------------------------------------
pos:            ldspi lo;dpy<md;
                iwrtdrl;br cont
"--------------------------------------------------------------
        $end

"****** vius16          vec 16-bt int unpak sgn rel 1.0         aug 81
"
"  v e c t o r   16 - b i t   i n t e g e r   u n p a c k   s i g n e d
"
"  history:
"       164 original    dec 80          d. davis (vups16)
"       integer version may 82          r. bair
"
"  purpose:
"       to unpack four 16-bit integers signed from each single
"       64 bit word of an array.
"
"  calling sequence:
"       call vius16(a,i,c,k,n)
"
"  input parameters:
"       a               word            array
"                       source vector
"       i               integer         scalar
"                       md element step for vector a
"       k               integer         scalar
"                       md element step for vector c
"       n               integer         scalar
"                       element count for vectors a and c
"
"       output parameters:
"       c               word (integer)  array
"                       result vector
"
"  subprograms used:
"       myh$pu_reslve
"
"  error conditions:
"       none
"
"  description:
"       four 16-bit integers from each 64-bit word of an array
"       are unpacked with sign extend. the integers are unpacked
"       from left to right.
"
"  global registers modified:
"       none
"
"  scratch registers modified:
"      sp(0,2,4,13,14),  dpx(-2,-1,0,1),  dpy(0,1),  fa, fm, md, tm
"
"  speed (cycles/integer):
"       best:           5
"       typical:        5
"       worst:          5
"       (procedure setup:       80 cycles)
"       (mlsp setup:            20 cycles)
"
"  ps size:
"       78
"
"  md size:
"       1
"
"
"
        $title vius16
        $radix d'8'
        $hsubr vius16/udc, a,i,c,k,n
        $entry mth$sp_vius16
        $insert '=tm.key'
"
        $pentf vius16/proc, glosav, aps3        "procedure entry
        ldspi 16; db=hval(z'f0ffffff')  "exception disable mask
        inc 0; setma                    "begin s-pad parameter setup
        incma; ldspi 17;raps3           "read apstat3
        incma; and 16,17;ldaps3;db=spfn "disable exception interrupts
        incma; ldspi  0; db=md
        incma; ldspi  1; db=md
               ldspi  2; db=md
               ldspi  3; db=md
               ldspi  4; db=md
        mov  1, 1; setma
        mov  3, 3; setma
        mov  4, 4; setma
        ldspi  1; db=md
        ldspi  3; db=md
        ldspi  4; db=md
        jsr mth$sp_vius16               "jsr to mlsp entry
        fapush
        fapush
        ldspi 16; db=hval(z'ff781fff')  "exception indicator mask
        ldspi 17;raps                   "read apstatus
        and 16,17;ldaps;db=spfn         "clear exception indicators
        $pexit vius16
"
        $psect glosav,md,ovl
glosav: $rs 1
        $psect
"
"
"       s-pads
        a $equ 0
        i $equ 1
        c $equ 2
        k $equ 3
        n $equ 4
        frt8 $equ 15
        sxt $equ 16
"       data pads
        m $equ 1
        m1 $equ -1
        m2 $equ -2
        mask $equ tm$lpwmsk+60
        mask1 $equ tm$slide1+60
        mask2 $equ tm$lpwmsk+33
"--------------------------------------------------------
        $pentf mth$sp_vius16/nonstd
         ldtma;db=mask                 "mask:  (0:47)=0;   (48:63)=1
        ldtma;db=mask1             "mask1:(0:47)=0;(48:63)=100000
        dpx(m)<tm;sub k,c               "save mask ,initialize c
        ldtma;db=mask2                 "(0:26)=0; (27:63)=1   mask2
        mov a,a;setma;dpx(m1)<tm             "fetch first word
        dpx(m2)<tm                        "save mask2
        ldspi sxt;db=16.                "load sp(sxt) = 16.
        flsh dpx(m2);mov sxt,sxt        "shift mask2 16 bits to left
        fadd;ldspi frt8;db=-48.         "load sp(frt8) = -48.
        dpx(m2)<fa                "save (0:10)=0; (11:47)=1; (48:63)=0
"------------------------------------------------------------
        flsh md;mov frt8,frt8;dpy<md    "shift left most 16 bit
                            "byte to far right, save entire 1'st word
        fadd
        fland dpx(m),fa   "mask out only 16 far right bits of 1'st word
        fadd
        fland dpx(m1),fa;dpy(1)<fa      "test for sign of far right
                                             "16 bit byte of 1'st word
        flsh dpy;add sxt,frt8           "shift next left 16 bit byte
                                        "of 1'st word to far right
        add i,a;setma;                          "fetch 2'nd word
        fadd
        fland dpx(m),fa;            "mask out second byte of 17st word
        bfeq nosign10                "branch if 1'st byte is positive
        flor dpx(m2),dpy(1)             "mask in sign extend for neg
                                               "1'st byte of 1'st word
        fland dpx(m1),fa;dpy(1)<fa         "test for sign of 2'nd byte
                                                "of 1'st word
        flsh dpy;add sxt,frt8;dpx<fa    "shift 3'rd left byte to right
                                        "save negative 1'st byte
                                                "of 1'st word
cont10:  fiadd zero,dpx;                "save 1'st left byte
                                                "of 1'st word
                                jmp loop                "jump to loop
"-------------------------------------------------------------
"--------------------------------------------------------------------
nosign10:       fadd;dpx<dpy(1)         "save positive 1'st left byte
                                                "of 1'st word
        fland dpx(m1),fa;dpy(1)<fa        "test for sign of 2'rd byte
                                                "of 1'st word
        flsh dpy;add sxt,frt8;          "shift 3'rd left byte
                                                "of 1'st word to right
          br cont10
"------------------------------------------------------------------
nosign2:       fadd;dpx<dpy(1);
        add k,c;setma;mi<fa          "store left most byte of n'th word
        fland dpx(m1),fa;dpy(1)<fa      "test for sign of 3'rd word
                                                "of n'th byte
        flsh dpy;add sxt,frt8;       "4'th byte of n'th word at right
                                                "no shift
          br cont2
nosign3:       fadd;dpx<dpy(1);
        add k,c;setma;mi<fa             "store 2'nd byte of n'th word
        fland dpx(m1),fa;dpy(1)<fa            "test sign of 4'th byte
                                                "of n'th word
        flsh md;mov frt8,frt8;dpy<md;         "shift left most byte
                                             "of n'th word to far right
                                                "save entire word
          br cont3
"-----------------------------------------------------------
"-----------------------------------------------------------
loop:                   fland dpx(m),fa;bfeq nosign2
                                "mask out 3'rd byte of n'th word
                                "branch if 2'nd byte > 0
                        flor dpx(m2),dpy(1);
                "mask in sign extend bits for 2'nd byte of n'th word
                        add k,c;setma;mi<fa
                                "store left most byte of n'th word
                        fland dpx(m1),fa;dpy(1)<fa
                              "test for sign of 3'rd byte of n'th word
                        flsh dpy;add sxt,frt8;dpx<fa
                                "4'th byte at right, no shift
cont2:                  fiadd zero,dpx; "save 2'rd byte of n'th word
                        ldspi frt8;db=-48.      "load sp(frt8) = -48.
                        fland dpx(m),fa;bfeq nosign3
                    "test sign of 4'th byte, branch if 3'rd byte > 0
                        flor dpx(m2),dpy(1);
                        "mask in sign extend of 3'rd byte of n'th word
                        add k,c;setma;mi<fa     "store 2'nd byte
                                                        "of n'th word
                        fland dpx(m1),fa;dpy(1)<fa
                                "test sign of 4'th byte of n'th word
        flsh md;mov frt8,frt8;dpy<md;       "shift left most byte
                                        "of (n+1)'th word to far right
                        dpx<fa
cont3:                  fiadd zero,dpx   "save 3'rd byte of n'th word
        fland dpx(m),fa;                "mask out left most byte
                                                "of (n+1)'th word
                        bfeq nosign4    "branch if 4'th byte > 0
                        flor dpx(m2),dpy(1);
                        "mask in sign extend of 4'th byte of n'th word
                        add k,c;setma;mi<fa     "store 3'rd byte of
                                                "n'th word
        fland dpx(m1),fa;dpy(1)<fa      "test for sign of left most
                                        "byte of (n+1)'th word
        flsh dpy;add sxt,frt8;        "shift 2'nd byte of (n+1)'th
                                        "word to right
                        dpx<fa
cont4:  add i,a;setma;                  "fetch (n+1)'th word
                        fiadd zero,dpx  "save 4'th byte of n'th word
        fland dpx(m),fa;           "mask out 2'nd byte of (n+1)'th word
        bfeq nosign1;   "branch if left most byte of (n+1)'th > 0
                                dec n           "decrement count
        flor dpx(m2),dpy(1);    "mask in sign extend of 1'st byte of
                                        "(n+1)'th word
                        add k,c;setma;mi<fa;beq out
                                "store 4'th byte of n'th word
        fland dpx(m1),fa;dpy(1)<fa      "test sign of 2'nd byte of
                                        "(n+1)'th word
        flsh dpy;add sxt,frt8;dpx<fa    "shift 3'rd byte of (n+1)'th
                                                "word
cont1:  fiadd zero,dpx;         "save 1'st byte of (n+1)'th word
                                jmp loop
"-------------------------------------------------------
nosign4:       fadd;dpx<dpy(1);
        add k,c;setma;mi<fa     "store 3'rd byte of n'th word
        fland dpx(m1),fa;dpy(1)<fa      "test for sign of 1'st byte of
                                                "(n+1)'th word
        flsh dpy;add sxt,frt8;  "shift 2'nd byte of (n+1)'th
                                        "word to right
          br cont4
nosign1:       fadd;dpx<dpy(1);
        add k,c;setma;mi<fa;beq out     "store $'th byte of n'th word
        fland dpx(m1),fa;dpy(1)<fa      "test for sign of 2'nd byte
                                        "of (n+1)'th word
        flsh dpy;add sxt,frt8;  "shift 3'ed byte of (n+1)'th
                                "word to right
          br cont1
out:    return
        $end

"****** vius8           vec 8-bt byte unpak sgn rel 1.0         aug 81
"
"  v e c t o r   8 - b i t   b y t e   u n p a c k   s i g n e d
"
"  history:
"       164 original    dec 80          d. davis (vups8)
"       integer version may 82          r. bair
"
"  purpose:
"       to unpack eight 8-bit bytes signed from each single
"       64 bit word of an array.
"
"  calling sequence:
"       call vius8(a,i,c,k,n)
"
"  input parameters:
"       a               word            array
"                       source vector
"       i               integer         scalar
"                       md element step for vector a
"       k               integer         scalar
"                       md element step for vector c
"       n               integer         scalar
"                       element count for vectors a and c
"
"       output parameters:
"       c               word (integer)  array
"                       result vector
"
"  subprograms used:
"       myh$pu_reslve
"
"  error conditions:
"       none
"
"  description:
"       eight 8-bit bytes from each 64-bit word of an array
"       are unpacked with sign extend. the bytes are unpacked
"       from left to right.
"
"  global registers modified:
"       none
"
"  scratch registers modified:
"      sp(0,2,4,13,14,15),  dpx(-2,-1,0,1),  dpy(0,1),  fa, fm, md, tm
"
"  speed (cycles/integer):
"       best:           5
"       typical:        5
"       worst:          5.1
"       (procedure setup:       82 cycles)
"       (mlsp setup:            22 cycles)
"
"  ps size:
"       70
"
"  md size:
"       1
"
"
"
        $title vius8
        $radix d'8'
        $hsubr vius8/udc, a,i,c,k,n
        $entry mth$sp_vius8
        $insert '=tm.key'
"
        $pentf vius8/proc, glosav, aps3        "procedure entry
        ldspi 16; db=hval(z'f0ffffff')  "exception disable mask
        inc 0; setma                    "begin s-pad parameter setup
        incma; ldspi 17;raps3           "read apstat3
        incma; and 16,17;ldaps3;db=spfn "disable exception interrupts
        incma; ldspi  0; db=md
        incma; ldspi  1; db=md
               ldspi  2; db=md
               ldspi  3; db=md
               ldspi  4; db=md
        mov  1, 1; setma
        mov  3, 3; setma
        mov  4, 4; setma
        ldspi  1; db=md
        ldspi  3; db=md
        ldspi  4; db=md
        jsr mth$sp_vius8                "jsr to mlsp entry
        fapush
        fapush
        ldspi 16; db=hval(z'ff781fff')  "exception indicator mask
        ldspi 17;raps                   "read apstatus
        and 16,17;ldaps;db=spfn         "clear exception indicators
        $pexit vius8
"
        $psect glosav,md,ovl
glosav: $rs 1
        $psect
"
"
"
"       s-pads
        a $equ 0
        i $equ 1
        c $equ 2
        k $equ 3
        n $equ 4
        fvt6 $equ 15
        egt $equ 16
        sav $equ 17
"       data pads
        m $equ 1
        m1 $equ -1
        m2 $equ -2
        mask $equ tm$lpwmsk+70
        mask1 $equ tm$slide1+70
        mask2 $equ tm$lpwmsk+23
"-----------------------------------------------------------------------
        $pentf mth$sp_vius8/nonstd      "mlsp entry
        ldtma;db=mask                 "mask:  (0:55)=0;   (56:63)=1
        ldtma;db=mask1                "mask1:(0:55)=0;(56:63)=200
        dpx(m)<tm;sub k,c         "save mask, initiate result address
        ldtma;db=mask2           "mask2:(0:18)=0;(19:63)=1
        mov a,a;setma;dpx(m1)<tm          "fetch 1'st word, save mask1
        dpx(m2)<tm                      "save mask2
        ldspi egt;db=8.                 "load sp(egt) = 8.
        flsh dpx(m2);mov egt,egt                    "shift mask2
                                       "(0:10)=0; (11:55)=1; (56:63)=0
        fadd;ldspi sav;db=-70           "load sp(sac) = -56.
        dpx(m2)<fa                      "save new mask2
"-----------------------------------------------------------------------
        flsh md;mov sav,fvt6;dpy<md            "shift 1'st byte right
                                                "save entire word
        fadd
"-----------------------------------------------------------------------
                fland dpx(m),fa                    "mask out 1'st byte
                fadd
                fland dpx(m1),fa;dpy(1)<fa     "test for sign of
                                                        "1'st byte
        flsh dpy;add egt,fvt6                   "shift 2'nd byte right
        fadd
"-----------------------------------------------------------------------
"-----------------------------------------------------------------------
                fland dpx(m),fa;                 "mask out 2'nd byte
                        bfeq nos         "branch if 1'st byte >0
                        flor dpx(m2),dpy(1)       "add sign extend
                                                "bits to 1'st byte
                fland dpx(m1),fa;dpy(1)<fa               "test for sign
                                                        "of 2'nd byte
        flsh dpy;add egt,fvt6;             "shift 3'rd byte right
                        dpx<fa;br inloop
"-----------------------------------------------------------------------
nos:                    fadd;dpx<dpy(1)       "save positive 1'st byte
                fland dpx(m1),fa;dpy(1)<fa    "test sign of 2'nd byte
        flsh dpy;add egt,fvt6                   "shift 3'rd byte right
"-----------------------------------------------------------------------
"-----------------------------------------------------------------------
inloop:                         fiadd zero,dpx;          "save result
                                beq outloop     "branch if end of word
                                                "reached
                fland dpx(m),fa;                    "mask out byte
                        bfeq nosign                "branch if byte >0
                        flor dpx(m2),dpy(1);     "add sign extend bits
                                add k,c;setma;mi<fa     "store result
                fland dpx(m1),fa;dpy(1)<fa      "test for sign of byte
        flsh dpy;add egt,fvt6;                 "shift next byte right
                        dpx<fa;
                                br inloop       "branch back to loop
"-----------------------------------------------------------------------
outloop:  add i,a;setma;                        "fetch (n+3)'th word
                fland dpx(m),fa;                "mask out byte
                        bfeq outnos             "branch for byte > 0
                        flor dpx(m2),dpy(1);    "add sign extend bits
                                add k,c;setma;mi<fa     "store result
                fland dpx(m1),fa;dpy(1)<fa;     "test for sign of byte
                                dec n           "decrement count
        flsh md;mov sav,fvt6;dpy<md;         "shift 1'st byte right
                        dpx<fa;
                                bne inloop         "branch if not done
"-----------------------------------------------------------------------
        br out                          "branch to conclusion
"-------------------------------------------------------------------
nosign:                 fadd;dpx<dpy(1);        "save positive byte
                                add k,c;setma;mi<fa     "store result
                fland dpx(m1),fa;                   "test sign of byte
                dpy(1)<fa
        flsh dpy;add egt,fvt6;                  "shift next byte right
                                br inloop            "branch to inloop
"-----------------------------------------------------------------------
outnos:                 fadd;dpx<dpy(1);        "save positive byte
                                add k,c;setma;mi<fa     "store result
                fland dpx(m1),fa;                  "test sign of byte
                dpy(1)<fa;
                                dec n           "decrement count
        flsh md;mov sav,fvt6;dpy<md;            "shift 1'st byte right
                                bne inloop   "branch inloop if not done
"-----------------------------------------------------------------------
out:                            fiadd zero,dpx          "save result
                        bfeq nos1;          "branch for positive byte
                                fadd
                        flor dpx(m2),dpy(1);    "add sign extend bits
                                add k,c;setma;mi<fa     "store result
                        fadd
                        fiadd zero,fa           "save result
                        "------------------------------
                                fadd
                                add k,c;setma;mi<fa;    "store result
                                return
                        "-----------------------------------
nos1:                   fiadd zero,dpy(1);      "save result
                                add k,c;setma;mi<fa     "store result
                        fadd
                                add k,c;setma;mi<fa;    "store result
                                return
"----------------------------------------------------------------
        $end
_IF(max)
"
"***** pakmat4 *****  first card in routine
"***** pakmat4 *****
"                    routine packs matrix inserting skip factors
"                    into last 4 bits of the REAL  word.
"              note: previous versions used last byte of word.
"
" call pakmat4(b,mcolb,mrowb,ncol,nrow,lenpak,lenofl,pak,lenofp,npackd)
"
" b = input matrix
" mcolb = skip between column elements
" mrowb = skip between rows
" ncol = dimension of columns
" nrow = dimension of rows
" lenpak = output array for length of packed columns
" lenofl = length of input array lenpak
" pak = output array for packed matrix
" lenofp = length of input array pak
" npackd = number of vectors packed
"
           $title pakmat4
           $entry pakmat4
           $insert '=tm.key'
           $insert '=cpu.key'
"
" spad definitions
"
           cap       $equ sp(0)
           count     $equ sp(1)
           mcolb     $equ sp(2)
           pak       $equ sp(3)
           test      $equ sp(4)
           ib        $equ sp(5)
           ibb       $equ sp(6)
           mrowb     $equ sp(7)
           ncol      $equ sp(8)
           nrow      $equ sp(9)
           lenpak    $equ sp(10)
           paksav    $equ sp(11)
           lenofl    $equ sp(12)
           lenofp    $equ sp(13)
           npackd    $equ sp(14)
           nrow_sav  $equ sp(15)
"
" data pad definitions
"
           val       $equ dpx(0)
           skip      $equ dpx(1)
           mask      $equ dpx(2)
           rmcolb    $equ dpy(-2)
           two       $equ dpy(-1)
           one       $equ dpy(0)
           i253      $equ dpy(1)
           sk1       $equ dpy(2)
           sk2       $equ dpy(3)
           $pentf pakmat4
"
"
" get parameters
"
      inc cap; setma
      inc cap; setma
      inc cap; setma
      ldspi ib; db=md       " save addr(b)
      ldma; db=md
      ldma; db=md; mov ib,ibb
      inc cap; setma
      ldspi mcolb; db=md   " save mcolb
      ldspi mrowb; db=md   " save mrowb
      ldma; db=md
      inc cap; setma
      inc cap; setma
      ldspi ncol; db=md    " save ncol
      ldma; db=md
      ldspi lenpak; db=md  " save addr(lenpak)
      inc cap; setma
      ldspi nrow; db=md    " save nrow
      inc cap; setma
      ldma; db=md
      inc cap; setma
      ldspi pak; db=md     " save addr(pak)
      ldspi lenofl; db=md  " save lenofl
      ldma; db=md
      inc cap; setma
      ldtma; db=tm$i1
      ldspi lenofp; db=md
      ldspi npackd; db=md  " save addr(npackd)
      one<tm
      ldtma; db=tm$i2
      ldma; db=mask_4_bits  " note last version used tm$by7imsk
      two<tm
      val<db; db=13
      mask<md
      i253<val  " note that last version used 253 not 13 here.
"
" loop over number of rows
"
          mov nrow,nrow_sav
          dec lenpak
          dec pak
          sub mrowb,ib
next_row: add mrowb,ib
          mov ib,ibb
          sub mcolb,ibb
          mov pak,paksav
          skip<one
          mov ncol,count
          andi# 1,count
"
          beq shift_count              " test if count is even.
"
" do first odd element so final loop is even
"
      add mcolb,ibb; setma
      nop
      nop
      ldspl test; db=md; fland mask,md
      fapush; mov test,test
      flor skip,fa; beq is_zero
      fapush; skip<one
      inc pak; setma; mi<fa ; br dec_count
"
is_zero: skip<two; dec# count        " check if on last one.
         bgt dec_count
         inc pak; setma; mi<one      " is last one so write out.
"
dec_count: dec count
           jmple get_len
"
shift_count: movr count,count; br loop               "count/2
"
" loophead is loop
"
" have to write out last element of vector ... be lazy and write last tw
     +o out.
"
done:                          fland mask,md
      flor skip,fa
      fapush; skip<one        ;flor one,fa
      inc pak; setma;
                      mi<fa   ;fapush
                               inc pak; setma; mi<fa
      jmp get_len

"
zero12:                        skip<sk2; br loop
"
"
"      element b(1)             element b(2)                   counters
"     ____________________     _________________________      __________
     +____
"
loop: add mcolb,ibb; setma    ;nop                   ;fiadd skip,one
      nop                     ;add mcolb,ibb; setma  ;ficmp skip,i253
      nop                     ;nop                   ;fiadd skip,two;
                                                       dec count; sk1<fa
      ldspl test; db=md;
             fland mask,md    ;nop                   ;beq done; fapush
      mov test,test           ;fland mask,md         ;sk2<fa; bfge lgskp
      flor skip,fa;
                 beq zero1    ;ldspl test; db=md
      fapush; skip<one        ;flor one,fa; mov test,test
      inc pak; setma;
                      mi<fa   ;fapush; beq zero2
                               inc pak; setma; mi<fa   ;br loop
"
"
zero1:                         flor sk1,fa; mov test,test
                               fapush; skip<one; beq zero12
                               inc pak; setma; mi<fa   ;br loop
"
"
zero2:                         skip<two; br loop
"
"
lgskp: flor skip,fa           ;ldspl test; db=md
      fapush; skip<one        ;flor one,fa; mov test,test
      inc pak; setma;
                      mi<fa   ;fapush; beq zero2
                               inc pak; setma; mi<fa   ;br loop

"
" now work out packed length = (pak-paksav)
"
get_len: mov pak,count
         sub paksav,count; db=spfn; val<db
         inc lenpak; setma; mi<db; db=val   " write out packed length
"
" check on length of pak and lenpak to see if have enough room to
" store another vector
"
         sub count,lenofp
         sub# ncol,lenofp
         dec lenofl; blt finished
         dec nrow; beq finished
         jmpgt next_row
"
"
finished: sub nrow,nrow_sav; db=spfn; val<db
          mov npackd,npackd; setma; mi<val  " write npackd
          dpx(0)<zero
          $pexit pakmat4
"
          $psect .data.,md,cat
mask_4_bits: $word z'fffffffffffffff0'
      $end                    " last card in routine
"***** pkmxmm4 *****   first card in routine
"***** pkmxmm4 *****
"                   performs matrix multiply between a column packed
"                   matrix in md and one row of vectors in max
"                   r(i,j)=r(i,j) + amax(i,k)*bmd(k,j)
"                   b is assumed to have been packed into bmd
"                   with routine pakmat4 and a to have been loaded onto
"                   the max with prdtld.
"             note: skip factors are kept in last 4 bits of word not
"                   the last byte as in previous versions
"
           $title pkmxmm4
           $entry pkmxmm4
           $insert '=tm.key'
           $insert '=maxa.key'
"
" call pkmxmm4(b,lenb,len,n_bvec,nmax_vec,r,mcolr,mrowr,
"             istart,ifun,ierr)
"
" b = packed input array
" lenb = array containing packed lengths of b vectors
" len = unpacked length of b vectors
" nb_vec = number of packed b vectors
" nmax_vec = number of vectors stored on max (<=8*nmax_brd)
"
" r = output matrix
" mcolr = skip between elements of a column of r
" mrowr = skip between rows of r
" istart = pointer into max where a is stored
" ifun = 0 gives addition, all else subtraction
" ierr = error flag = 0 normal return
"                   = 1 wrong length when vector is unpacked
"                   = 2 packed vector length < 0
"                   = 3 unpacked length + istart > 2048
"                   = 4 unpacked length < 0
"                   = 5 packed length > unpacked length
"
" spad definitions
"
           cap              $equ sp(0)
           tmp_b            $equ sp(1)
           tmp_indx_reg     $equ sp(2)
           tmp_brd_cst      $equ sp(3)
           b                $equ sp(4)
           lenb             $equ sp(5)
           nb_vec           $equ sp(6)
           r                $equ sp(7)
           mcolr            $equ sp(8)
           mrowr            $equ sp(9)
           rr               $equ sp(10)
           rr_w             $equ sp(11)
           ierr             $equ sp(12)
           max_vec_inc      $equ sp(13)
           max_tbl          $equ sp(14)
           max_res          $equ sp(15)
           max_vec_len      $equ sp(16)
"
" dp definitions
"
           max_8            $equ dpx(-2)
           val              $equ dpx(0)
           last             $equ dpx(1)
           plen             $equ dpx(2)
           max_vec_cnt      $equ dpx(3)
"
           vec_cnt          $equ dpy(-4)
           lasty            $equ dpy(-3)
           len              $equ dpy(-2)
           nmax_vec         $equ dpy(-1)
           tmp_indx         $equ dpy(0)
           tmp_mask         $equ dpy(1)
           tmp_count        $equ dpy(2)
           indx             $equ dpy(3)
"
"
" get the parameters
"
pkmxmm4:    ldspi 1; db=hval(z'feffffff')
           ldspi 2; raps3
           and 2,1; ldaps3; db=spfn      " turn off floating interrupts.
           inc cap; setma
           inc cap; setma
           inc cap; setma
           ldspi b; db=md              " addr(packed input array b)
           ldspi lenb; db=md           " addr(array of packed lengths)
           ldma; db=md
           inc cap; setma
           inc cap; setma
           len<db; db=md; inc cap; setma " length of unpacked vectors
           ldma; db=md
           ldma; db=md
           ldspi r; db=md               " addr(output matrix)
           ldspi nb_vec; db=md          " no. packed b vectors
           nmax_vec<db; db=md; inc cap; setma  " no. vecs in max
           inc cap; setma
           inc cap; setma
           ldma; db=md
           ldma; db=md
           ldma; db=md
           ldspi mcolr; db=md      " mcolr
           ldspi mrowr; db=md      " mrowr
           indx<db; db=md; inc cap; setma " istart
           inc cap; setma
           last<db; db=len
           ldma; db=md
           ldspi ierr; db=md              " addr(ierr)
           nop
           val<db; db=md                 " ifun
"
" test if use addition or subtraction
"
           ficmp val,zero
           fapush; val<db; db=max$add+max$forw+max$vdot
           ldma; db=max$cntl+max$brd-1
           bfeq cntl
           val<db; db=max$sub+max$forw+max$vdot
cntl:      incma; mi<db; db=val
"
"
" vector is addressed from 0. decrement istart by two so just
" need to add on the skip factors. also should always end up
" with the index being last=istart+len-2
" at same time get first packed vector length, and set up tmp_b,tmp_indx
" and tmp_count
"
           ldtma; db=tm$i2
           fiadd indx,last; mov lenb,lenb; setma
           fisubr tm,indx; mov b,tmp_b
           fisubr tm,fa
           fapush; indx<fa ;plen<db; db=md      " packed length of b vec
           last<fa; tmp_indx<db; db=indx
           tmp_count<db; db=plen
" some paramter checks for the more sane amoungst us
           ficmp len,zero
           ficmp plen,zero
           ficmp len,plen
           fapush; jmpfle len_le_zero
           jmpfle plen_le_zero
           jmpflt len_lt_plen
           val<db; db=max$v_len-2
           lasty<db; db=val
           ficmp last,lasty
           fapush
           nop
           jmpfgt last_gt_2046
"
"
" now start loop over b vectors
"
loop_b_vec:  jsr packd_md_max_dot        " call kernel packed dot routin
"
" whilst waiting for the sum collapse do as much as possible
"
           ficmp tmp_indx,last      " check for mismatch
" set up values needed to get results
           fapush; ldspi max_vec_inc; db=max$b_v-max$a_v  " inc btwn max
           ldspi max_tbl; db=maxtbl-1  " addr(common/maxtbl/)-1
           jmpfne mismatch
           vec_cnt<db; db=nmax_vec      " temp counter vecs in max
           mov r,rr                     " rr = addr(to read r from)
           mov rr,rr_w                  " rr_w = addr(to write r to)
           sub mcolr,rr_w
           ldtma; db=tm$i1              " get one onto tm output
           ldspi max_vec_len; db=max$v_len-1
" now update pointers for the next b vector
           ldspi tmp_b; db=plen
           add tmp_b,b                  " b points to next packed vector
           mov b,tmp_b
           inc lenb; setma              " get next packed length
           tmp_indx<db; db=indx
           add mrowr,r                  " r points to next row
           plen<db; tmp_count<db; db=md
" check the packed length again
           ficmp plen,zero; dec# nb_vec
           ficmp len,plen; beq on_last_one
           fapush
           jmpfle plen_le_zero
           jmpflt len_lt_plen
           br nxtbrd
on_last_one: nop
           nop
           nop
           nop
nxtbrd:       inc max_tbl; setma
              val<db; db=8
              max_8<val
              ldspi max_res; db=md
"
              add max_vec_len,max_res; setma     " max_res is addr(max(2
              nop
              mov rr,rr; setma
"
              add max_vec_inc,max_res; setma;
                  val<db; db=md; fisubr tm,max_8
              nop;
                  nop; fisubr tm,vec_cnt
              add mcolr,rr; setma;
                  fadd val,md; max_8<fa
"
res_loop:     add max_vec_inc,max_res; setma;
                  val<db; db=md; fisubr tm,max_8;
                      bfeq last_el; vec_cnt<fa
"
              nop;
                  nop; fisubr tm,vec_cnt;
                      add mcolr,rr_w; setma; mi<fa; bfeq done
"
              add mcolr,rr; setma;
                  fadd val,md; max_8<fa; br res_loop
"
last_el:      add mcolr,rr_w; setma; mi<fa; bfeq done
              jmp nxtbrd
"
done:      dec nb_vec
           jmpgt loop_b_vec
"
" successful exit.
ierr_0:    val<db; db=zero
finish:    mov ierr,ierr; setma; mi<db; db=val   " set ierr to return co
           return
"
len_le_zero: ficmp len,zero
             fapush
             val<db; db=4
             jmpfeq ierr_0         " zero len is ok. return with zero co
             jmp finish        " -ve len not ok
"
plen_le_zero: ficmp plen,zero
           fapush
           val<db; db=2
           jmpfeq done        " if packed length is zero then just do ne
           jmp finish   " -ve plen not a healthy sign !!
"
len_lt_plen: val<db; db=5
             jmp finish
"
mismatch:  val<db; db=1
           jmp finish
"
last_gt_2046: val<db; db=3
              jmp finish
"
"
"
"
"
           $entry packd_md_max_dot
"
"
"
" kernel of the packed max dot product
"
" local definitions. on arrival tmp_b (address of packed vector),
" tmp_count (length of packed vector), tmp_indx (istart-2)
" must be set up.
"
" tmp_indx returns the final max index which should be equal to
" istart+unpaked_len-2 if unpacking has worked.
" any register named tmp_ is changed on return from previous value.
" routine also changes md and tm output.
"
"          tmp_b        $equ sp(1)
"          tmp_indx_reg $equ sp(2)
"          tmp_brd_cst  $equ sp(3)
"
           tmp_val      $equ dpx(0)
"
"          tmp_indx     $equ dpy(0)
"          tmp_mask     $equ dpy(1)
"          tmp_count    $equ dpy(2)
"
" zero the scalar registers
"
packd_md_max_dot: ldma; db=max$scl+max$brd-1
           ldspi tmp_brd_cst; db=16
           dec tmp_brd_cst
sloop:     incma; mi<db; db=zero; dec tmp_brd_cst; bgt sloop
"
" start setup for the loop
"
""""""           ldtma; db=tm$by7msk
          ldma; db=mask_4_bits
          nop
           ldtma; db=tm$i1
" 1)
           mov tmp_b,tmp_b; setma; tmp_mask<db; db=md  " was db=tm
" 2)
           ldspi tmp_indx_reg; db=max$indx+max$brd
" 3)
           ldspi tmp_brd_cst; db=max$dp+max$brd
" 4)
           fland tmp_mask,md
" 1)
           inc tmp_b; setma;
               fapush;
                   fisubr tm,tmp_count
" 2)
           nop;
               fiadd tmp_indx,fa
" 3)
           nop;
               fapush; tmp_val<md;
                   tmp_count<fa
" 4)
           fland tmp_mask,md;
              mov tmp_indx_reg,tmp_indx_reg; setma; mi<fa;
                 tmp_indx<fa; bfeq neq1
" 1)
loop:      inc tmp_b; setma;
               fapush;
                   fisubr tm,tmp_count
" 2)
           nop;
               fiadd tmp_indx,fa;
                   mov tmp_brd_cst,tmp_brd_cst; setma; mi<tmp_val;
                       fapush
" 3)
           nop;
               fapush; tmp_val<md;
                   tmp_count<fa
" 4)
           fland tmp_mask,md;
               mov tmp_indx_reg,tmp_indx_reg; setma; mi<fa;
                   tmp_indx<fa; bfgt loop
"
"
neq1:      nop
           mov tmp_brd_cst,tmp_brd_cst; setma; mi<tmp_val
"
           tmp_val<db; db=zero
           nop
           nop
           mov tmp_brd_cst,tmp_brd_cst; setma; mi<tmp_val
"
" initiate the collapse so that the final sums will be in max(2048)
"
           tmp_val<db; db=max$v_len-1
           mov tmp_indx_reg,tmp_indx_reg; setma; mi<db; db=tmp_val
           ldma; db=max$clps+max$brd
           statma; mi<db; db=zero
"
" results will be available 29 cycles from now.
"
           return
"
"
"
" definition of common /maxtbl/ set up in the first call of the job
" to routine prdtld
"
           $psect maxtbl,md,ovl
maxtbl:    $rs 15
n_max_brd: $rs 1
"
           $psect .data.,md,cat
mask_4_bits: $word z'000000000000000f'
           $end    " last card of pkmxmm4, packd_md_max_dot
"***** prdtld *****  first card in routine
"***** prdtld *****
"                    routine loads max for prdtmx
"
"  written by r. j. harrison, december, 1985.
"
           $title prdtld
           $entry prdtld
           $insert '=tm.key'
           $insert '=maxa.key'
           $ext sys$rdmaxtbl
"
" call prdtld(s,is,iss,len,ns,istart,ierr)
" s = input array
" is = skip distance between elements of a vector
" iss = skip distance between vectors
" len = length of vector (less than 2048)
" ns = number of vectors (no more than 8*nmax_brd)
" istart =  starting pointer into max.
" ierr = error flag : =0 normal return.
"                     =1 len < 0
"                     =2 ns <= 0
"                     =3 ns > 8*nmax_brd
"                     =4 is <= 0
"                     =5 iss < len (disabled in this version so can tran
     +spose)
"                     =6 len+istart>2048
"
"
" scratch pad definitions
"
           cap       $equ sp(0)
           max_tbl   $equ sp(1)
           vec_inc   $equ sp(2)
           max_indx  $equ sp(3)
           s         $equ sp(4)
           is        $equ sp(5)
           iss       $equ sp(6)
           ns        $equ sp(7)
           ierr      $equ sp(8)
           ipt       $equ sp(9)
           istart    $equ sp(10)
           vec_cnt   $equ sp(11)
           spt       $equ sp(12)
"
" data pad definitions
"
           val       $equ dpx(0)
           len       $equ dpx(1)
           num_brd   $equ dpx(2)
           maxlen    $equ dpy(0)
           temp      $equ dpy(1)
           last      $equ dpy(2)
           diss      $equ dpy(3)
"
           $pentf prdtld
"
" determine whether sys$rdmxtbl needs to be called
"
          ldma; db=n_max_brd
          nop
          nop
          ficmp zero,md; num_brd<md
          fapush
          nop
          jmpfne load_par       " jump if max table information is prese
     +nt
"
" call sys$rdmxtbl to get configuration table and number of active board
     +s (dpx)
"
          ldma; db=sav_area-1
          incma; mov cap,cap; mi<db; db=spfn
          ldspi cap; db=maxtbl_par
          jsr sys$rdmaxtbl
          num_brd<dpx(0)
"
" compress the table so inactive boards are not referenced.
"
          dpy(0)<db; db=max$no_addr
          ldspi 1; db=maxtbl-1              " addr to read from table
          ldspi 2; db=15
          ldspi 3; db=maxtbl-1              " addr to write to table
tblmak:   inc 1; setma                      " get element of max table
          nop
          nop
          fcmp dpy(0),md
          fapush
          nop
          bfeq dec_cnt
          inc 3; setma; mi<md               " write address of configure
     +d board
dec_cnt:  dec 2
          bne tblmak
          inc 1; mi<db; db=num_brd; setma  " write num. boards to mxtbl(
     +16)
"
" load parameters after restoring cap
"
          ldma; db=sav_area
          nop
          nop
          ldspi cap; db=md
load_par: inc cap; setma                               " get addr(s)
          inc cap; setma                               " get addr(is)
          inc cap; setma                               " get addr(iss)
          ldspi s; db=md                               " save addr(s)
          ldma; db=md                                  " get is
          ldma; db=md                                  " get iss
          inc cap; setma                               " get addr(len)
          ldspi is; db=md                              " save is
          ldspi iss; db=md; diss<db                   " save iss
          ldma; db=md                                  " get len
          inc cap; setma                               " get addr(ns)
          inc cap; setma                               " get addr(istart
     +)
          len<db; db=md                                " save len
          ldma; db=md                                  " get ns
          ldma; db=md                                  " get istart
          inc cap; setma                               " get addr(ierr)
          ldspi ns; db=md                              " save ns
          ldspi istart; db=md; temp<db                 " save istart
          ldspi ierr; db=md                            " save addr(ierr)
"
" check that parameters are reasonable
"
          maxlen<db; db=max$v_len
          fishi 3,num_brd
          fiadd temp,len; val<db; db=spfn; mov ns,ns
          fisub val,fa
          ficmp maxlen,fa
          ficmp len,zero
          jmpfgt ns_too_big; fapush                " jump if ns>24
          jmpflt len_gt_max                        " jump if len+istart>
     +2048
          jmpfle len_le_0                          " jump if len <= 0
          mov ns,ns
          jmple ns_le_0                            " jump if ns <= 0
          mov is,is
"
" not change to is < 0 not is <= 0
"
          jmplt is_le_0                            " jump if is <= 0
"          fisub diss,len
"          fapush
"          nop
"          jmpflt iss_lt_len
"
" set up pointers and addresses
"
          ldspi max_tbl; db=maxtbl-1               " addr(addr next brd)
          ldtma; db=tm$i1
          ldspi vec_inc; db=max$b_v-max$a_v        " max vector incremen
     +t
          subi 2,istart
"
" do the transfer
"
nxtbrd:   inc max_tbl; setma                           " get addr of nex
     +t board
          ldspi vec_cnt; db=8                          " reset num vec t
     +o 8
          nop
          ldspi ipt; db=md
          add istart,ipt                               " add on istart
          mov ipt,max_indx
"
nxtvec:   mov s,spt; setma; fisubr tm,len              " new vector
          fapush
veclop:   add is,spt; setma; fisubr tm,fa
          inc max_indx; setma; mi<md; fapush; bfgt veclop
"
          dec ns                                       " check if finish
     +ed
          jmpeq done
"
          add iss,s                                    " point to next m
     +d vector
          add vec_inc,ipt                              " point to next m
     +ax vect.
          mov ipt,max_indx
          dec vec_cnt
          bgt nxtvec
          jmp nxtbrd                                    " next board
"
len_le_0: ficmp len,zero
          fapush; val<db; db=1
          mov ierr,ierr; setma; mi<db; db=val
          jmpfeq done
          jmp finish
ns_le_0: val<db; db=2
         mov ierr,ierr; setma; mi<db; db=val
         jmp finish
ns_too_big: val<db; db=3
            mov ierr,ierr; setma; mi<db; db=val
            jmp finish
is_le_0: val<db; db=4
         mov ierr,ierr; setma; mi<db; db=val
         jmp finish
"iss_lt_len: val<db; db=5
"            mov ierr,ierr; setma; mi<db; db=val
"            jmp finish
len_gt_max: val<db; db=6
            mov ierr,ierr; setma; mi<db; db=val
            jmp finish
done:     mov ierr,ierr; setma; mi<db; db=zero
finish:   $pexit prdtld
          $psect .data.,md,cat
sav_area:    $rs 1
maxtbl_par:  $word 1
             $word maxtbl
"
" definition of common block /maxtbl/
"
          $psect maxtbl,md,ovl
maxtbl:      $rs 15
n_max_brd:   $word 0
          $end                    " last card in routine
_ENDIF
_ENDIF
_IF(ibm)
*member name = asslog
         title '64-bit logical operations - fortran callable functions'
***********************************************************************
*                                                                     *
***  zlogic64                    version 0   release 0.0      mar 85  *
*                                                                     *
*      fortran callable functions to perform various logical          *
*      operations on 64-bit quantities.
*                                                                     *
*                                                                     *
*                                                                     *
** copyright (c) 1985, science & engineering research council         *
*                                                                     *
***********************************************************************
zlogic64 csect
         dc    al1(dlh-*),c'zlogic64 (c) 1985 serc '
dlh      equ   *
r0       equ   0
r1       equ   1
r2       equ   2
r3       equ   3
r4       equ   4
r5       equ   5
r6       equ   6
r7       equ   7
r8       equ   8
r9       equ   9
r10      equ   10
r11      equ   11
r12      equ   12
r13      equ   13
r14      equ   14
r15      equ   15
         eject
***********************************************************************
*                                                                     *
**  popcnt:                                                           *
*     to return the number of bits set to 1 in the argument           *
*                                                                     *
**  leadz :                                                           *
*     to return the number of bits set to 0 prior to the first 1 bit  *
*                                                                     *
*                                                                     *
**  calling sequence:                                                 *
*        integer popcnt,leadz                                         *
*        ... = ... popcnt(z1) ...                                     *
*        ... = ... leadz (z2) ...                                     *
*   - z1 and z2:  treated as 8 byte bit-strings.                      *
*                                                                     *
***********************************************************************
         space 3
* leadz  - count number of 0 bits to the left of the leftmost 1 bit
*
         entry leadz
leadz    ds    0h                  entry point
         using *,r15               establish addressability
         stm   r2,r3,28(r13)       save regs
         l     r2,0(,r1)           address of argument
         lm    r2,r3,0(r2)         load argument
         sr    r0,r0               zero
         bctr  r0,0                minus one, all bits 1
         xr    r2,r0               reverse all bits
         xr    r3,r0                both regs
         sr    r0,r0               initialise count to zero
         space
lnxbit   ltr   r2,r2               test most sig bit
         bc    11,exit             zero, found first 1 bit of argument
         sldl  r2,1                shift up all bits, fill with zeros
         bct   r0,lnxbit           increment count (neg), next
         space 3
* popcnt - count number of 1 bits in argument
*
         entry popcnt
popcnt   ds    0h                  entry point for pop count function
         using *,r15               establish addressability
         stm   r2,r3,28(r13)       save regs
         l     r2,0(,r1)           address of argument
         lm    r2,r3,0(r2)         load argument
         la    r1,2                number of words in argument
         sr    r0,r0               initialise count to zero
         space
pnxbit   alr   r2,r2               shift word left, set cc shifted bit
         bc    12,pdone            shifted bit 0
         bctr  r0,0                shifted bit 1, incr count (neg)
pdone    bc    5,pnxbit            remainder not 0, repeat
         lr    r2,r3               next word
         bct   r1,pnxbit           go count its bits
         space 2
exit     lcr   r0,r0               make count positive
         sr    r15,r15             set return code (reqd by fortran)
         lm    r2,r3,28(r13)       restore regs
         br    r14                 return
         eject
***********************************************************************
*                                                                     *
**   and:                                                             *
*     to return the logical and of the two arguments                  *
*                                                                     *
**   or:                                                              *
*     to return the logical or of the two arguments                   *
*                                                                     *
**   xor:                                                             *
*     to return the exclusive or of the two arguments                 *
*                                                                     *
*                                                                     *
**  calling sequence:                                                 *
*        REAL  and, or, xor, z1,z2                                  *
*        ... = ...  and(z1,z2) ...                                    *
*        ... = ...  or(z1,z2) ...                                     *
*        ... = ...  xor(z1,z2) ...                                    *
*   - z1 and z2:  treated as 8 byte bit strings
*                                                                     *
***********************************************************************
         space 3
*  and - logical and of two 8-byte arguments
*
         entry and
and      ds    0h                  entry point
         using *,r15               establish addressability
         stm   r2,r2,28(r13)       save regs
         lm    r1,r2,0(r1)         load parameter addresses
         mvc   dwork,0(r1)         copy first parm
         nc    dwork,0(r2)         and in second
         ld    r0,dwork            put into fl.pt.reg 0 for return
         b     return
         space 3
* or - logical or of two 8-byte arguments
*
         entry or
or       ds    0h                  entry point
         using *,r15               establish addressability
         stm   r2,r2,28(r13)       save regs
         lm    r1,r2,0(r1)         load parameter addresses
         mvc   dwork,0(r1)         copy first parm
         oc    dwork,0(r2)         or in second
         ld    r0,dwork            put into fl.pt.reg 0 for return
         b     return
         space 3
* xor - exclusive or of two 8-byte arguments
*
         entry xor
xor      ds    0h                  entry point
         using *,r15               establish addressability
         stm   r2,r2,28(r13)       save regs
         lm    r1,r2,0(r1)         load parameter addresses
         mvc   dwork,0(r1)         copy first parm
         xc    dwork,0(r2)         xor in second
         ld    r0,dwork            put into fl.pt.reg 0 for return
         b     return
         eject
***********************************************************************
*                                                                     *
**   eq :                                                             *
*     to return the logical value of the equality of the arguments    *
*                                                                     *
**   pad:                                                             *
*     to return an 8-byte quantity from the 4-byte argument           *
*                                                                     *
**   ipad:                                                            *
*     to return a 4-byte quantity from the 8-byte argument            *
*                                                                     *
*                                                                     *
**  calling sequence:                                                 *
*        REAL  pad, z1,z2                                           *
*        integer*4  ipad, i1                                          *
*        logical*4 eq                                                 *
*        ... = ...  eq (z1,z2) ...                                    *
*        ... = ...  pad(i1) ...                                       *
*        ... = ...  ipad(z1) ...                                      *
*   - z1 and z2:  treated as 8 byte bit strings
*                                                                     *
***********************************************************************
         space 3
* eq  - return logical equality of the two arguments
*
         entry eq
eq       ds    0h                  entry point
         using *,r15               establish addressability
         stm   r2,r2,28(r13)       save regs
         lm    r1,r2,0(r1)         load parm addresses
         sr    r0,r0               preset result false
         mvc   dwork,0(r1)         copy first parm
         xc    dwork,0(r2)         xor in second - will be 0 if equal
         bnz   *+8                 non zero result - false (=x'00')
         la    r0,1                zero result - true (=x'01')
         b     return
         space 3
*  pad - return an 8-byte quantity from the 4-byte argument,
*          padded on the left with zeros
*
         entry pad
pad      ds    0h                  entry point
         using *,r15               establish addressability
         stm   r2,r2,28(r13)       entry point
         l     r1,0(,r1)           load parm address
         xc    dwork,dwork         set to zero
         mvc   dwork+4(4),0(r1)    copy in integer
         ld    r0,dwork            load into result reg
         b     return
         space 3
* ipad - return the lower 4 bytes of the 8-byte argument
*
         entry ipad
ipad     ds    0h                  entry point
         using *,r15               establish addressability
         stm   r2,r2,28(r13)       save regs
         l     r1,0(,r1)           load parm address
         l     r0,4(,r1)           load low half of parm
         b     return
         eject
***********************************************************************
*                                                                     *
**  shiftr:                                                           *
*     to return the argument right shifted                            *
*                                                                     *
**  shiftl:                                                           *
*     to return the argument left shifted                             *
*                                                                     *
**   shift:                                                           *
*     to return the argument circularly shifted                       *
*                                                                     *
*  in all cases the number of bits by which the argument is shifted   *
*     is calculated modulo 64.
*                                                                     *
**  calling sequence:                                                 *
*        REAL  shiftr,shiftl,shift,  z1                              *
*        integer*4 i                                                  *
*        ... = ... shiftl(z1,i) ...                                   *
*        ... = ... shiftr(z1,i) ...                                   *
*        ... = ...  shift(z1,i) ...                                   *
*   - z1 and z2:  treated as 8 byte bit strings
*                                                                     *
***********************************************************************
         space 3
* drshift - shift argument right, end off zero fill
*
         entry shiftr
shiftr   ds    0h                  entry point
         using *,r15               establish addressability
         stm   r2,r2,28(r13)       save regs
         lm    r1,r2,0(r1)         load parm addresses
         lm    r0,r1,0(r1)         load dword to be shifted
         l     r2,0(,r2)           load amount to shift by
         srdl  r0,0(r2)            perform shift
         stm   r0,r1,dwork         store result
         ld    r0,dwork            and load back into result reg
         b     return
         space 3
* shiftl - shift argument left, end off zero fill
*
         entry shiftl
shiftl   ds    0h                  entry point
         using *,r15               establish addressability
         stm   r2,r2,28(r13)       save regs
         lm    r1,r2,0(r1)         load parm addresses
         lm    r0,r1,0(r1)         load dword to be shifted
         l     r2,0(,r2)           load amount to shift by
         sldl  r0,0(r2)            perform shift
         stm   r0,r1,dwork         store result
         ld    r0,dwork            and load back into result reg
         b     return
         space 3
* shift - shift argument circularly (left)
*
         entry shift
shift    ds    0h                  entry point
         using *,r15               establish addressability
         stm   r2,r6,28(r13)       save regs
         lm    r1,r2,0(r1)         load parm addresses
         la    r6,63               mask for shift value
         n     r6,0(r2)            load shift amount
         lm    r2,r3,0(r1)         load dword to be shifted
         lm    r4,r5,0(r1)         and a second copy
         sldl  r2,0(r6)            shift right by requested amount
         la    r0,64               total number of bits
         sr    r0,r6               number of bits to shift other way
         lr    r6,r0               set complement shift
         srdl  r4,0(r6)            shift it
         or    r2,r4               or together
         or    r3,r5               the two halves
         stm   r2,r3,dwork         store
         ld    r0,dwork            and load back into result reg
         lm    r3,r6,32(r13)       restore extra regs
         eject
         drop  r15                 end using base reg
return   lm    r2,r2,28(r13)       restore regs used
         sr    r15,r15             set return code
         br    r14                 return
         space 2
dwork    ds    d                   double word work area
         print nogen
*        regdef
         end
*member name = nag77a
         macro
         aafmac &version,&hep
         lclc  &hop
&hop     setc  '&hep'
*
         aif   ('&version' eq 'watfiv').l1
x03aaf   csect
         ago   .l2
.l1      anop
x03aay   csect
.l2      anop
*
*        the following is the interface:
*            subroutine x03aaf(a,na,b,nb,n,ia,ib,c1,c2,d1,d2,sw,ifail)
*
         using *,r15
         b     entry
         dc    x'06',cl6'x03aaf'   routine name for tracebacks
savea    ds    18f
entry    stm   r14,r12,12(r13)     save registers
         lr    r14,r13
         la    r13,savea           point to new save area
         drop  r15
         using savea,r13           establish base register
         st    r13,8(r14)
         st    r14,4(r13)          chain save areas
*
*        copy the arguments as appropriate.
*
         lm    r14,r10,0(r1)       pick up all addresses
         st    r14,a
         mvc   na,0(r15)
         st    r0,b
         mvc   nb,0(r1)
         mvc   n,0(r2)
         mvc   ia,0(r3)
         mvc   ib,0(r4)
         mvc   c1,0(r5)
         mvc   c2,0(r6)
         st    r7,d1
         st    r8,d2
         mvc   sw,3(r9)            move last byte only
         st    r10,ifail
         aif   ('&version' ne 'watfiv').l3
*
*        the watfiv version must be called from an interface routine to
*        type check and to standardise the array arguments.
*
         l     r1,a
         mvc   a,4(r1)             pick up address of first elem. of a
         l     r1,b
         mvc   b,4(r1)             pick up address of first elem. of b
.l3      anop
*
*        now check and manipulate the arguments.
*
         l     r1,ifail            copy value of ifail
         mvc   ifail1,0(r1)
         l     r8,n                pick up n
         lr    r7,r8
         sh    r7,=h'1'            get (n-1)
         bnm   skip0
         slr   r7,r7               rounded up to zero
skip0    l     r9,ia
         ltr   r9,r9
         bnp   error1              ia <= 0
         l     r10,ib
         ltr   r10,r10
         bnp   error1              ib <=  0
         lr    r1,r7
         mr    r0,r9               (n-1)*ia
         c     r1,na
         bnl   error2              na <= (n-1)*ia
         lr    r1,r7
         mr    r0,r10              (n-1)*ib
         c     r1,nb
         bnl   error2              nb <= (n-1)*ib
         l     r6,a
         l     r7,b                pick up addresses
         sll   r9,3
         sll   r10,3               convert steps to bytes
         slr   r0,r0
         st    r0,ierror           set ierror to zero for exit
         tm    sw,x'ff'
         bz    basic               use basic precision if sw is false
*
*        evaluate the sum in extended precision.
*
         aif   ('&hep' ne 'both').m0
&hop     setc  'yes'
         l     r1,16               find cvt
         tm    182(r1),x'01'
         bz    nohep               hep is not present
.m0      anop
         ld    fp0,c1              pick up initial value
         sdr   fp2,fp2
         ld    fp4,c2              pick up second half
         sdr   fp6,fp6
         aif   ('&hop' ne 'yes').m1
         axr   fp0,fp4             c1+c2 in extended precision
         ago   .m2
.m1      anop
         bal   r14,axr04           simulate axr 0,4
.m2      anop
         ltr   r8,r8
         bnp   s1&hop              n is not positive
*
*        form the inner product.
*
l1&hop   ld    fp4,0(r6)           pick up a(i) and b(i)
         aif   ('&hop' ne 'yes').m3
         mxd   fp4,0(r7)
         axr   fp0,fp4             sum a(i)*b(i)
         ago   .m4
.m3      anop
         std   fp0,save            save the fp registers across call
         std   fp2,save1
         ld    fp0,0(r7)
         bal   r14,mxdr04          simulate mxdr 0,4
         ld    fp4,save
         ld    fp6,save1
         bal   r14,axr04           simulate axr 0,4
.m4      anop
         alr   r6,r9
         alr   r7,r10              update pointers and repeat
         bct   r8,l1&hop
*
*        round and return the results.
*
s1&hop   l     r6,d1               pick up argument addresses
         l     r7,d2
         aif   ('&hop' ne 'yes').m5
         lrdr  fp4,fp0             get d1 (rounded)
         ago   .m6
.m5      anop
         std   fp0,save            save fp registers across call
         std   fp2,save1
         bal   r14,lrdr00          simulate lrdr 0,0
         ldr   fp4,fp0
         ld    fp0,save
         ld    fp2,save1           restore fp registers
.m6      anop
         std   fp4,0(r6)
         sdr   fp6,fp6             calculate result - d1
         aif   ('&hop' ne 'yes').m7
         sxr   fp0,fp4             do subtraction
         ago   .m8
.m7      anop
         lcdr  fp4,fp4
         bal   r14,axr04           simulate axr 0,4
.m8      anop
         std   fp0,0(r7)           store in d2
         b     exit
         aif   ('&hep' ne 'both' or '&hop' ne 'yes').m99
&hop     setc  'no'
*
*        software extended precision code.
*
nohep    ds    0h                  come here if no hep
         ago   .m0
.m99     anop
*
*        evaluate the sum in double precision.
*
basic    ld    fp0,c1              pick up initial value
         ltr   r8,r8
         bnp   skip2               n is not positive
loop2    ld    fp4,0(r6)
         md    fp4,0(r7)           multiply a(i) by b(i) and sum
         adr   fp0,fp4
         alr   r6,r9               update pointers
         alr   r7,r10
         bct   r8,loop2            repeat
skip2    l     r6,d1
         l     r7,d2               pick up argument addresses
         std   fp0,0(r6)
         sdr   fp0,fp0             clear d2 to zero
         std   fp0,0(r7)
         b     exit                set ifail and return
*
*        error exits for invalid arguments.
*
error1   la    r0,1                ifail=1 for ia or ib negative
         b     error
error2   la    r0,2                ifail=2 for na <= (n-1)*ia etc.
error    st    r0,ierror
         aif   ('&version' eq 'watfiv').l4
         la    r1,arglist
         l     r15,=v(nagmsg)
         balr  r14,r15             call nagmsg for error handling
         st    r0,ierror
.l4      anop
*
*        set ifail and return.
*
exit     l     r0,ierror
         l     r1,ifail            set ifail as appropriate
         st    r0,0(r1)
         sdr   fp6,fp6             for watfiv only, in fact
         l     r13,savea+4
         lm    r14,r12,12(r13)     restore registers
         mvi   12(r13),x'ff'
         br    r14                 return
*
*        argument list.
*
         aif   ('&version' eq 'watfiv').l5
arglist  dc    a(ifail1,ierror)    ibm fortran argument list
         dc    x'80',al3(srname)
srname   ds    0d                  routine name kept in a real
         dc    cl8'x03aaf'
.l5      anop
*        sepcode &hep
         ltorg
*
*        copies of arguments, addresses thereof and workspace.
*
a        ds    a                   address of a
na       ds    f                   value of na
b        ds    a                   address of b
nb       ds    f                   value of nb
n        ds    f                   value of n
ia       ds    f                   value of ia
ib       ds    f                   value of ib
c1       ds    d                   value of c1
c2       ds    d                   value of c2
d1       ds    a                   address of d1
d2       ds    a                   address of d2
sw       ds    x                   value of sw
ifail    ds    a                   address of ifail
ifail1   ds    f                   arguments to nagmsg
ierror   ds    f
r0       equ   0
r1       equ   1
r2       equ   2
r3       equ   3
r4       equ   4
r5       equ   5
r6       equ   6
r7       equ   7
r8       equ   8
r9       equ   9
r10      equ   10
r11      equ   11
r12      equ   12
r13      equ   13
r14      equ   14
r15      equ   15
*
fp0      equ   0
fp2      equ   2
fp4      equ   4
fp6      equ   6
*
         mend
         aafmac ibm,yes
         end
*member name = tdaa
indvor csect
*               i*4 function indvor(i1,i2,i3,i4)
 using *,15
 stm 2,7,space
 lm 2,5,0(1)
 l 2,0(0,2)
 l 3,0(0,3)
 l 4,0(0,4)
 l 5,0(0,5)
 cr 2,3
 bc 11,n10
 lr 6,3
 lr 3,2
 lr 2,6
n10 cr 4,5
 bc 11,n20
 lr 6,5
 lr 5,4
 lr 4,6
n20 cr 2,4
 bc 2,n100
 bc 8,n30
 lr 6,4
 lr 4,2
 lr 2,6
 lr 6,5
 lr 5,3
 lr 3,6
 bc 15,n100
n30 cr 3,5
 bc 11,n100
 lr 6,5
 lr 5,3
 lr 3,6
 bc 15,n100
 entry indvx
*                i*4 function indvx(i,j,k,l)
 using *,15
indvx stm 2,7,space
 lm 2,5,0(1)
 l 2,0(0,2)         2=i
 l 3,0(0,3)         3=j
 l 4,0(0,4)         4=k
 l 5,0(0,5)         5=l
n100 balr 7,0
 using *,7
 sll 2,2            2=4i
 l 6,knd0
 a 3,0(2,6)         3=kndvec(i)+j = ib
 sll 3,2            3=4ib
 l 6,jnd0
 a 4,0(3,6)         4=k+jndvec(ib) = ia
 sll 4,2            4=4ia
 sll 5,2            5=4l
 l 6,ind0
 l 0,0(4,6)         0=indvec(ia)
 l 6,lnd0
 a 0,0(5,6)         0=indvec(ia)+lndvec(l)
 lm 2,7,space
 bcr 15,14
             entry linkor
*                            call linkor(indvec,jndvec,kndvec,lndvec)
*                            to initialize addresses
 using *,15
linkor stm 2,6,space
 lm 2,5,0(1)
 la 6,4
 sr 2,6
 sr 3,6
 sr 4,6
 sr 5,6
 stm 2,5,ind0
 lm 2,6,space
 bcr 15,14
space ds 6f
ind0  dc xl4'00000000'
jnd0  dc xl4'00000000'
knd0  dc xl4'00000000'
lnd0  dc xl4'00000000'
 end
*member name = gfa
indv csect
*                i*4 function indv(i,j,k,l)
 using *,15
 stm 2,6,space
 lm 2,5,0(1)
 l 2,0(0,2)         2=i
 l 3,0(0,3)         3=j
 l 4,0(0,4)         4=k
 l 5,0(0,5)         5=l
 sll 2,2            2=4i
 l 6,knd0
 a 3,0(2,6)         3=kndvec(i)+j = ib
 sll 3,2            3=4ib
 l 6,jnd0
 a 4,0(3,6)         4=k+jndvec(ib) = ia
 sll 4,2            4=4ia
 sll 5,2            5=4l
 l 6,ind0
 l 0,0(4,6)         0=indvec(ia)
 l 6,lnd0
 a 0,0(5,6)         0=indvec(ia)+lndvec(l)
 lm 2,6,space
 bcr 15,14
               entry indd
*              i*4 function indd(i,j)
*              orders i and j, returns triangle index
 using *,15
indd stm 2,4,space
 lm 2,3,0(1)
 l  0,0(0,3)                 indd=j
 l  4,knd0
 c  0,0(0,2)
 bc 12,n200                  branch if i.ge.j
 lr 3,0                      3=j
 l  0,0(0,2)                 indd=i
 sll 3,2                     3=4j
 a  0,0(3,4)                 indd=i+kndvec(j)
 lm 2,4,space
 bcr 15,14
n200 l 2,0(0,2)              2=i
 sll 2,2                     2=4i
 a  0,0(2,4)                 indd=j+kndvec(i)
 lm 2,4,space
 bcr 15,14
             entry linkin
*                            call linkin(indvec,jndvec,kndvec,lndvec)
*                            to initialize addresses
 using *,15
linkin stm 2,6,space
 lm 2,5,0(1)
 la 6,4
 sr 2,6
 sr 3,6
 sr 4,6
 sr 5,6
 stm 2,5,ind0
 lm 2,6,space
 bcr 15,14
space ds 5f
ind0  dc xl4'00000000'
jnd0  dc xl4'00000000'
knd0  dc xl4'00000000'
lnd0  dc xl4'00000000'
 end
_IFN(3090vf)
*member name = mxmbass
mxmb csect
 using *,15
 stm 0,14,sarea
 lm 1,12,0(1)
 l 10,0(0,10)
 l 11,0(0,11)
 l 13,0(0,12)
 ltr 10,10
 bz p999
 ltr 11,11
 bz p999
 ltr 0,13
 bz p999
 la 13,1(0,0)
 l 2,0(0,2)
 l 3,0(0,3)
 l 5,0(0,5)
 l 6,0(0,6)
 l 8,0(0,8)
 l 9,0(0,9)
 sldl 2,3
 sll 5,3
 sll 6,3
 sldl 8,3
 st 1,ia1
 st 11,nlink
 st 6,mrowb
 st 9,mrowr
p1 l 1,ia1
 st 4,ib
 l 11,nlink
p2 ld 2,0(0,4)
 ltdr 2,2
 be p22
 lcr 13,13
 bnl p44
 ldr 0,2
 lr 14,1
 b p22
p44 lr 12,1
 lr 6,7
 lr 9,10
p3 ld 6,0(0,12)
 ld 4,0(0,14)
 mdr 6,2
 mdr 4,0
 ar 12,2
 ar 14,2
 ad 6,0(0,6)
 adr 6,4
 std 6,0(0,6)
 ar 6,8
 bct 9,p3
p22 ar 4,5
 ar 1,3
 bct 11,p2
 ltr 13,13
 bnl p88
 lcr 13,13
 lr 6,7
 lr 9,10
p99 ld 6,0(0,14)
 mdr 6,0
 ar 14,2
 ad 6,0(0,6)
 std 6,0(0,6)
 ar 6,8
 bct 9,p99
p88 l 4,ib
 a 7,mrowr
 a 4,mrowb
 bct 0,p1
p999 lm 0,14,sarea
 br 14
sarea ds 15f
mrowb ds 1f
mrowr ds 1f
ia1 ds 1f
ib ds 1f
nlink ds 1f
 end
*member name = mxmbn
mxmbn csect
 using *,15
 stm 0,14,sarea
 lm 1,12,0(1)
 la 13,1(0,0)
 l 2,0(0,2)
 l 3,0(0,3)
 l 5,0(0,5)
 l 6,0(0,6)
 l 8,0(0,8)
 l 9,0(0,9)
 l 10,0(0,10)
 l 11,0(0,11)
 l 0,0(0,12)
 sldl 2,3
 sll 5,3
 sll 6,3
 sldl 8,3
 st 1,ia1
 st 11,nlink
 st 6,mrowb
 st 9,mrowr
p1 l 1,ia1
 st 4,ib
 l 11,nlink
p2 ld 2,0(0,4)
 lcdr 2,2
 be p22
 lcr 13,13
 bnl p44
 ldr 0,2
 lr 14,1
 b p22
p44 lr 12,1
 lr 6,7
 lr 9,10
p3 ld 6,0(0,12)
 ld 4,0(0,14)
 mdr 6,2
 mdr 4,0
 ar 12,2
 ar 14,2
 ad 6,0(0,6)
 adr 6,4
 std 6,0(0,6)
 ar 6,8
 bct 9,p3
p22 ar 4,5
 ar 1,3
 bct 11,p2
 ltr 13,13
 bnl p88
 lcr 13,13
 lr 6,7
 lr 9,10
p99 ld 6,0(0,14)
 mdr 6,0
 ar 14,2
 ad 6,0(0,6)
 std 6,0(0,6)
 ar 6,8
 bct 9,p99
p88 l 4,ib
 a 7,mrowr
 a 4,mrowb
 bct 0,p1
 lm 0,14,sarea
 br 14
sarea ds 15f
mrowb ds 1f
mrowr ds 1f
ia1 ds 1f
ib ds 1f
nlink ds 1f
 end
*member name = daxpy
***********************************************************************
*                                                                     *
*      subroutine daxpy(n,sa,sx,incx,sy,incy)
*c
*c     overwrite single precision sy with single precision sa*sx +sy.
*c
*      dimension sx(1),sy(1)
*      if(n.le.0.or.sa.eq.0.) return
*      if(incx.eq.incy) if(incx-1) 10,30,70
*10    continue
*c
*c        code for nonequal or nonpositive increments.
*c
*      ix = 1
*      iy = 1
*      if(incx.lt.0)ix = (-n+1)*incx + 1
*      if(incy.lt.0)iy = (-n+1)*incy + 1
*      do 20 i = 1,n
*        sy(iy) = sy(iy) + sa*sx(ix)
*        ix = ix + incx
*        iy = iy + incy
*20    continue
*      return
*c
*c        code for both increments equal to 1
*c
*c
*c        clean-up loop so remaining vector length is a multiple of 4.
*c
*30    m = n - (n/4)*4
*      if( m .eq. 0 ) go to 50
*      do 40 i = 1,m
*        sy(i) = sy(i) + sa*sx(i)
*40    continue
*      if( n .lt. 4 ) return
*50    mp1 = m + 1
*      do 60 i = mp1,n,4
*        sy(i) = sy(i) + sa*sx(i)
*        sy(i + 1) = sy(i + 1) + sa*sx(i + 1)
*        sy(i + 2) = sy(i + 2) + sa*sx(i + 2)
*        sy(i + 3) = sy(i + 3) + sa*sx(i + 3)
*60    continue
*      return
*c
*c        code for equal, positive, nonunit increments.
*c
*70    continue
*      ns = n*incx
*          do 80 i=1,ns,incx
*          sy(i) = sa*sx(i) + sy(i)
*80        continue
*      return
*      end
daxpy    csect
         using daxpy,15
         stm   2,6,28(13)
         lm    1,6,0(1)
         l     1,0(1)   n
         ltr   1,1
         bz    daxpy9
         ld    0,0(2)  sa
         ltdr  0,0
         bz    daxpy9
         lr    2,3  sx
         lr    3,5    sy
         l     4,0(4)     incx
         l     5,0(6)   incy
         lh    0,=h'1'
         cr    4,0
         bne   daxpy5
         cr    5,0
         bne   daxpy5
* increments are both 1 so we use fast code
         la    0,32          increment (relies on n=reg1)
         la    4,0           index
         sll   1,3
         sh    1,=h'32'
         bm    daxpy2
daxpy1   ld    2,0(4,2)
         ld    4,8(4,2)
         mdr   2,0
         mdr   4,0
         ad    2,0(4,3)
         ad    4,8(4,3)
         std   2,0(4,3)
         std   4,8(4,3)
         ld    2,16(4,2)
         ld    4,24(4,2)
         mdr   2,0
         mdr   4,0
         ad    2,16(4,3)
         ad    4,24(4,3)
         std   2,16(4,3)
         std   4,24(4,3)
         bxle  4,0,daxpy1
daxpy2   ah    1,=h'24'
         la    0,8
         sr    4,0           step index back one cos bxh at top
daxpy3   bxh   4,0,daxpy9
         ld    2,0(4,2)
         mdr   2,0
         ad    2,0(4,3)
         std   2,0(4,3)
         b     daxpy3
*
* code for at least one increment .ne. 1
* for -ve increments, not standard f77 - see cray manual
daxpy4   ds    f
daxpy5   ltr   4,4
         bh    daxpy6
         la    0,1
         sr    0,1
         st    4,daxpy4
         mh    0,daxpy4+2
         ar    2,0
daxpy6   ltr   5,5
         bh    daxpy7
         la    0,1
         sr    0,1
         st    5,daxpy4
         mh    0,daxpy4+2
         ar    3,0
daxpy7   sla   4,3
         sla   5,3
daxpy8   ld    2,0(2)
         mdr   2,0
         ad    2,0(3)
         std   2,0(3)
         ar    2,4
         ar    3,5
         bct   1,daxpy8
*
daxpy9   lm    2,6,28(13)
         br    14
         ltorg
         end
*member name = ddot
*      function ddot(n,sx,incx,sy,incy)
*c
*c     returns the dot product of REAL  sx and sy.
*c
*      dimension sx(1),sy(1)
*      ddot = 0.
*      if(n.le.0)return
*      if(incx.eq.incy) if(incx-1) 10,30,40
*10    continue
*c
*c         code for unequal or nonpositive increments.
*c
*      ix = 1
*      iy = 1
*      if(incx.lt.0)ix = (-n+1)*incx + 1
*      if(incy.lt.0)iy = (-n+1)*incy + 1
*      do 20 i = 1,n
*         ddot = ddot + sx(ix)*sy(iy)
*        ix = ix + incx
*        iy = iy + incy
*20    continue
*      return
*c
*c        code for both increments equal to 1.
*30    ddot = vecsum(sx,sy,n)
*      return
*c
*c         code for positive equal increments .ne.1.
*c
*40    continue
*      ns = n*incx
*          do 50 i=1,ns,incx
*          ddot = ddot + sx(i)*sy(i)
*50        continue
*      return
*      end
ddot     csect
         using ddot,15
         stm   2,6,28(13)
         sdr   0,0
         lm    1,5,0(1)
         l     1,0(1)   n
         ltr   1,1
         bz    ddot9
         lr    0,4
         l     4,0(3)   incx
         lr    3,0   sy
* sx equ 2
         l     5,0(5)   incy
         lh    0,=h'1'
         cr    4,0
         bne   ddot5
         cr    5,0
         bne   ddot5
* increments are both 1 so we use fast code
         la    0,32          increment (relies on n=reg1)
         la    6,0           index
         sll   1,3
         sh    1,=h'32'
         bm    ddot2
ddot1    ld    2,0(6,2)
         ld    4,8(6,2)
         md    2,0(6,3)
         md    4,8(6,3)
         adr   0,2
         adr   0,4
         ld    2,16(6,2)
         ld    4,24(6,2)
         md    2,16(6,3)
         md    4,24(6,3)
         adr   0,2
         adr   0,4
         bxle  6,0,ddot1
ddot2    ah    1,=h'24'
         la    0,8
         sr    6,0           step index back one cos bxh at top
ddot3    bxh   6,0,ddot9
         ld    2,0(6,2)
         md    2,0(6,3)
         adr   0,2
         b     ddot3
*
* code for at least one increment .ne. 1
* for -ve increments, not standard f77 - see cray manual
ddot4    ds    f
ddot5    ltr   4,4
         bh    ddot6
         la    0,1
         sr    0,1
         st    4,ddot4
         mh    0,ddot4+2
         ar    2,0
ddot6    ltr   5,5
         bh    ddot7
         la    0,1
         sr    0,1
         st    5,ddot4
         mh    0,ddot4+2
         ar    3,0
ddot7    sla   4,3
         sla   5,3
ddot8    ld    2,0(2)
         md    2,0(3)
         ar    2,4
         ar    3,5
         bct   1,ddot8
*
ddot9    lm    2,6,28(13)
         br    14
         ltorg
         end
_ENDIF
_IFN(fortio)
*member name = fmove
fmove csect
*
*              call fmove(a(1),b(1),nword)
*            a,b on double word boundaries
*            nword is integer*2
*            b=a for nword double words
*
 using *,15
 stm 2,6,28(13)
 l 3,8(0,1)
 l 4,0(0,3)
 srdl 4,2
 ltr 4,4
 l 6,0(0,1)
 l 1,4(0,1)
 bz e31
 bct 4,*+4
 sll 4,5
 lr 3,1
 ar 3,4
 l 2,i32
l32 ld 0,0(0,6)
 ld 2,8(0,6)
 ld 4,16(0,6)
 ld 6,24(0,6)
 ar 6,2
 std 0,0(0,1)
 std 2,8(0,1)
 std 4,16(0,1)
 std 6,24(0,1)
 bxle 1,2,l32
 sr 4,4
e31 sldl 4,2
 ltr 4,4
 bz e33
l34 ld 0,0(0,6)
 a 6,i8
 std 0,0(0,1)
 a 1,i8
 bct 4,l34
e33 lm 2,6,28(13)
 br 14
  entry srch
*
*               call srch(a(1),iblk)
*               for bsam only
*
  using *,15
srch stm 2,7,sarea+4
 drop 15
 balr 2,0
 using *,2
 l 3,vblock
 using ddisc,3
 l 6,isel
 l  4,4(0,1)
 sll 6,2
  l 4,0(0,4)
 l 5,ipos-4(6)
 st 4,ipos-4(6)
    la 15,iopen-36(6)
 sr 4,5
 bz qq4
 st 14,sarea
 l 6,dcbtab-4(6)
  l 7,0(0,1)
 bp qq1
 st 6,mask
 la 1,mask
 mvi mask,144
 svc 23   issue tclose
 ar 4,5
 b qq3
qq1 ts 0(15)
 bm qq2
 st 6,mask
 la 1,mask
 mvi mask,131
 svc 19 open (inout option)
 using ihadcb,6
 tm dcboflgs,x'10'
  bz qq5
qq2  st  7,ecbr2+12
   st 6,ecbr2+8
  la 7,ecbr2
qq7     lr 1,7
 l 15,48(0,6)
 balr 14,15   bsam read
 lr 1,7
 l 15,52(0,6)
   balr 14,15      call check routine
 bct 4,qq7
qq6 l 14,sarea
qq4 lm 2,7,sarea+4
   br 14
qq5   mvi reply,x'ff'
  b qq6
qq3  bct 4,qq2
   b qq6
  drop 2,3,6
 entry wtbsam
*               call wtbsam(ibuf) ,ibuf is i*4 buffer number
*                  writes next sequential block from buffer
*
 using *,15
wtbsam st 14,sarea
 stm 2,6,sarea+4
 drop 15
 balr 2,0
 using *,2
 l  5,0(0,1)
 l 3,vblock
 using ddisc,3
 l 5,0(0,5)
 l 4,isel
 s 5,hone
 st 4,iselw
 mh 5,i32760
 sll 4,2
 a 5,vbufa
 la 15,iopen-36(4)
 ts 0(15)
 l 6,dcbtab-4(4)
 bm q1
 st 6,mask
 la 1,mask
 mvi mask,135          outin option
 svc 19 issue open
 using ihadcb,6
 tm dcboflgs,x'10'
 bz q2
q1 la 1,ecbw2
 st 6,ecbw2+8
 st 5,ecbw2+12
 l 15,48(0,6)
 st 1,curwecb
 balr 14,15 bsam write
q3 l 6,ipos-4(4)
 l 14,sarea
 a 6,hone
 st 6,ipos-4(4)
 lm 2,6,sarea+4
 br 14
q2 mvi reply,x'ff'
 b q3
 drop 2,3,6
 entry rdbsam
*
*               call rdbsam(ibuf) ,ibuf is i*4 buffer number
*               reads next sequential block to correct buffer
*
 using *,15
rdbsam st 14,sarea
 stm 2,6,sarea+4
 drop 15
 balr 2,0
 using *,2
 l  5,0(0,1)
 l 3,vblock
 using ddisc,3
 l 5,0(0,5)
 l 4,isel
 s 5,hone
 st 4,iselr
 mh 5,i32760
 sll 4,2
 a 5,vbufa
 la 15,iopen-36(4)
 ts 0(15)
 l 6,dcbtab-4(4)
 bm p1
 st 6,mask
 la 1,mask
 mvi mask,131           inout option
 svc 19 issue open
 using ihadcb,6
 tm dcboflgs,x'10'
 bz p2
p1 la 1,ecbr2
 st 6,ecbr2+8
 st 5,ecbr2+12
 l 15,48(0,6)
 st 1,currecb
 balr 14,15 bsam read
p3 l 6,ipos-4(4)
 l 14,sarea
 a 6,hone
 st 6,ipos-4(4)
 lm 2,6,sarea+4
 br 14
p2 mvi reply,x'ff'
 b p3
 drop 2,3,6
 entry shut
*
*               call shut
*               closes all atmol data sets
*
 using *,15
shut stm 2,6,sarea+4
 xr 1,1
 balr 2,0
 using *,2
 drop 15
 la 3,dcbtab
 l 6,vblock
 using ddisc,6
 la 4,4(0,0)
 la 5,dcbtab+28
l61 tm ipos,x'80'
 bnz e62
 l 15,0(0,3)
 ar 1,4
 st 15,warea-4(1)
e62 ar 6,4
 bxle 3,4,l61
 la 5,dcbtab+60
 la 6,iopen
l63 tm 0(6),x'80'
 bz e64
 l 15,0(0,3)
 ar 1,4
 st 15,warea-4(1)
e64 ar 6,4
 bxle 3,4,l63
 ltr 1,1
 bz e65
 st 14,sarea
 la 6,warea-4(1)
 mvi 0(6),x'80'
 la 1,warea
 svc 20 issue close for all atmol data sets
 l 14,sarea
e65 lm 2,6,sarea+4
 br 14
 drop 2,6
 entry rdtpnm
*
*               call rdtpnm(title(1))
*               title is real*4 11 items long
*
 using *,15
rdtpnm st 14,sarea
 stm 2,7,sarea+4
 drop 15
 balr 2,0
 using *,2
 l 3,0(0,1)
 l 4,vblock
 using ddisc,4
 l 6,isel
 mvi warea,x'00'
 sll 6,2
 l 7,dcbtab-4(6)
 st 7,mask
 la 1,mask
 mvi mask,128
 svc 64 issue rdjfcb
 ar 4,6
 cli warea,x'00'
 be e71
 mvc 0(44,3),warea
 mvc ipos-4(4),hone
 c 6,i32
 bh e71
  st 7,mask
 la 1,mask
 mvi mask,143
 svc 19 issue open for bdam data set with output option
 using ihadcb,7
 tm dcboflgs,x'10'
 bnz e71
 mvi ipos-4,x'ff'
e71 mvi 12(13),x'ff'
 l 14,blksiz-4
 srl 14,3
 st 14,blksiz-4
 l 14,sarea
 lm 2,7,sarea+4
 br 14
dcbout mvc blksiz-2(2),dcbblksi
 br 14
 drop 2,4,7
 entry wtbdam
*
*               call wtbdam(ibuf,ipb)
*               ibuf=i*4,buffer number
*               buffer length 4095 REAL  words=32760 bytes
*               ipb=i*4,physical block number,start counting at zero
*
 using *,15
wtbdam st 14,sarea
 stm 2,5,sarea+4
 drop 15
 balr 2,0
 using *,2
 l 4,0(0,1)
 l 15,vblock
 using ddisc,15
 l 3,4(0,1)
 l 5,isel
 l 4,0(0,4)
 a 3,hone
 st 5,iselw
 bct 4,*+4
 sll 5,2
 mh 4,i32760
 l 5,dcbtab-4(5)
 a 4,vbufa
 la 1,ecbw1
 st 5,ecbw1+8
 st 4,ecbw1+12
 st 3,ecbw1+24
 l 15,48(0,5)
 st 1,curwecb
 balr 14,15 bdam write
 l 14,sarea
 lm 2,5,sarea+4
 br 14
 drop 2
 entry rdbdam
*
*               call rdbdam(ibuf,ipb)
*               ibuf=i*4,buffer number
*               buffer length 4095 REAL  words=32760 bytes
*               ipb=i*4,physical block number,start counting at zero
*
 using *,15
rdbdam st 14,sarea
 stm 2,5,sarea+4
 drop 15
 balr 2,0
 using *,2
 l 4,0(0,1)
 l 15,vblock
 using ddisc,15
 l 3,4(0,1)
 l 5,isel
 l 4,0(0,4)
 a 3,hone
 st 5,iselr
 bct 4,*+4
 sll 5,2
 mh 4,i32760
 l 5,dcbtab-4(5)
 a 4,vbufa
 la 1,ecbr1
 st 5,ecbr1+8
 st 4,ecbr1+12
 st 3,ecbr1+24
 l 15,48(0,5)
 st 1,currecb
 balr 14,15 bdam read
 l 14,sarea
 lm 2,5,sarea+4
 br 14
 drop 2
   entry setbf
 using *,15
setbf stm 6,7,sarea+4
 l 6,4(0,1)
 l 1,0(0,1)
 st 1,ecbwr+12
 st 14,sarea
 balr 2,0
 drop 15
 using *,2
 la 7,sort
 st 7,bask
 mvi bask,143
 la 1,bask
 svc 19
 using ihadcb,7
 tm dcboflgs,x'10'
 bnz f7171
 l 7,vblock
 using ddisc,7
 mvi reply,x'ff'
f7171 l 14,sarea
 lm 6,7,sarea+4
 br 14
 using ihadcb,7
sdcb lh 1,dcbblksi
 srl 1,3
 st 1,0(0,6)
 br 14
 drop 2,7
 entry srtrd1
 using *,15
srtrd1 l 0,recb
 st 2,sarea+4
 b f13
 entry srtpt1
 using *,15
srtpt1 l 0,wecb
 st 2,sarea+4
f13 balr 2,0
 using *,2
 drop 15
 st 14,sarea
 st 0,ecbwr+4
 l 1,0(0,1)
 a 1,hone
 st 1,ecbwr+24
 l 15,sort+48
 la 1,ecbwr
 balr 14,15
 mvi curr,x'00'
 l 14,sarea
 l 2,sarea+4
 br 14
 drop 2
 entry srtst1
 using *,15
srtst1 ts curr
 st 2,sarea
 bm  f12
 la 1,ecbwr
 l 0,hone
 svc 1
 balr 2,0
 using *,2
 drop 15
 clc 1(2,1),hzero
 be f122
 l 1,vblock
 using ddisc,1
 mvi reply,x'ff'
f122 l 2,sarea
 drop 1,2
f12 br 14
 entry closbf
 using *,15
closbf st 14,sarea
 st 2,sarea+4
 drop 15
 balr 2,0
 using *,2
 l 1,xarea
 svc 20  close sort file
 l 14,sarea
 l 2,sarea+4
 br 14
 drop 2
 entry checkw
*
*                  call checkw
*               checks current write order if no write going on
*               effect is null
*                  if i/o error sets reply in common/disc/ -ve
*
 using *,15
checkw stm 2,4,sarea+4
 l 1,curwecb
 ts curwecb
 l 3,vblock
 using ddisc,3
 l 4,iselw
 mvi iselw+3,x'00'
 b e11
 entry checkr
*
*               call checkr
*               checks current read order,if no read going on
*               effect is null
*                  if i/o error sets reply in common/disc/ -ve
*
 using *,15
checkr stm 2,4,sarea+4
 l 1,currecb
 ts currecb
 l 3,vblock
 l 4,iselr
 mvi iselr+3,x'00'
e11 mvi 12(13),x'ff'
 drop 15
 balr 2,0
 using *,2
 bm e12
 c 4,i8
 st 14,sarea
 bh e13
 l 0,hone
 svc 1
 clc 1(2,1),hzero
 be e14
 mvi reply,x'ff'
 b e14
e13 check (1)
e14 l 14,sarea
e12 lm 2,4,sarea+4
 br 14
syn mvi reply,x'ff'
 br 14
 drop 2,3
   entry setbf2
 using *,15
setbf2 stm 6,7,sarea+4
 l 6,4(0,1)
 l 1,0(0,1)
 st 1,ecbws+12
 st 14,sarea
 balr 2,0
 drop 15
 using *,2
 la 7,port
 st 7,cask
 mvi cask,143
 la 1,cask
 svc 19
 using ihadcb,7
 tm dcboflgs,x'10'
 bnz f8181
 l 7,vblock
 using ddisc,7
 mvi reply,x'ff'
f8181 l 14,sarea
 lm 6,7,sarea+4
 br 14
 using ihadcb,7
sdcc lh 1,dcbblksi
 srl 1,3
 st 1,0(0,6)
 br 14
 drop 2,7
 entry srtrd2
 using *,15
srtrd2 l 0,recc
 st 2,sarea+4
 b g13
 entry srtpt2
 using *,15
srtpt2 l 0,wecc
 st 2,sarea+4
g13 balr 2,0
 using *,2
 drop 15
 st 14,sarea
 st 0,ecbws+4
 l 1,0(0,1)
 a 1,hone
 st 1,ecbws+24
 l 15,port+48
 la 1,ecbws
 balr 14,15
 mvi curs,x'00'
 l 14,sarea
 l 2,sarea+4
 br 14
 drop 2
 entry srtst2
 using *,15
srtst2 ts curs
 st 2,sarea
 bm  g12
 la 1,ecbws
 l 0,hone
 svc 1
 balr 2,0
 using *,2
 drop 15
 clc 1(2,1),hzero
 be g122
 l 1,vblock
 using ddisc,1
 mvi reply,x'ff'
g122 l 2,sarea
 drop 1,2
g12 br 14
 entry closb2
 using *,15
closb2 st 14,sarea
 st 2,sarea+4
 drop 15
 balr 2,0
 using *,2
 l 1,yarea
 svc 20  close sort file
 l 14,sarea
 l 2,sarea+4
 br 14
 drop 2
   entry setbf3
 using *,15
setbf3 stm 6,7,sarea+4
 l 6,4(0,1)
 l 1,0(0,1)
 st 1,ecbwt+12
 st 14,sarea
 balr 2,0
 drop 15
 using *,2
 la 7,data
 st 7,dask
 mvi dask,143
 la 1,dask
 svc 19
 using ihadcb,7
 tm dcboflgs,x'10'
 bnz f9181
 l 7,vblock
 using ddisc,7
 mvi reply,x'ff'
f9181 l 14,sarea
 lm 6,7,sarea+4
 br 14
 using ihadcb,7
sdcd lh 1,dcbblksi
 srl 1,3
 st 1,0(0,6)
 br 14
 drop 2,7
 entry srtrd3
 using *,15
srtrd3 l 0,recd
 st 2,sarea+4
 b h13
 entry srtpt3
 using *,15
srtpt3 l 0,wecd
 st 2,sarea+4
h13 balr 2,0
 using *,2
 drop 15
 st 14,sarea
 st 0,ecbwt+4
 l 1,0(0,1)
 a 1,hone
 st 1,ecbwt+24
 l 15,data+48
 la 1,ecbwt
 balr 14,15
 mvi curt,x'00'
 l 14,sarea
 l 2,sarea+4
 br 14
 drop 2
 entry srtst3
 using *,15
srtst3 ts curt
 st 2,sarea
 bm  h12
 la 1,ecbwt
 l 0,hone
 svc 1
 balr 2,0
 using *,2
 drop 15
 clc 1(2,1),hzero
 be h122
 l 1,vblock
 using ddisc,1
 mvi reply,x'ff'
h122 l 2,sarea
 drop 1,2
h12 br 14
 entry closb3
 using *,15
closb3 st 14,sarea
 st 2,sarea+4
 drop 15
 balr 2,0
 using *,2
 l 1,zarea
 svc 20  close sort file
 l 14,sarea
 l 2,sarea+4
 br 14
 drop 2
 entry check
 using *,15
check stm 2,5,sarea
 l 3,4(0,1)
 l 3,0(0,3)
 ltr 3,3
 l 4,8(0,1)
 l 4,0(0,4)
 bnp x99
 s 3,hone
 l 1,0(0,1)
 sll 3,3
 l 2,i8
 xr 5,5
 ar 3,1
x98 x 4,0(0,1)
 x 5,4(0,1)
 bxle 1,2,x98
 xr 4,5
x99 l 5,vblock
 using ddisc,5
 st 4,ichek
 lm 2,5,sarea
 br 14
 drop 5,15
hzero dc f'0'
hone dc f'1'
iopen dc f'0,0,0,0,0,0,0,0'
i32760 dc h'32760'
i8 dc f'8'
i32 dc f'32'
 ds 0d
sarea ds 8f
warea ds 44f
ecbw2 dc f'0'
 dc x'00200000'
 dc 3f'0'
ecbr2 dc f'0'
 dc x'00800000'
 dc 3f'0'
ecbw1 dc f'0'
wecb dc x'00400000'
 dc 5f'0'
wecc dc x'00400000'
 dc 5f'0'
wecd dc x'00400000'
 dc 5f'0'
ecbr1 dc f'0'
recb dc x'00480000'
 dc 5f'0'
recc dc x'00480000'
 dc 5f'0'
recd dc x'00480000'
 dc 5f'0'
exlst dc x'87'
 dc al3(warea)
exlsd dc x'07'
 dc al3(warea)
 dc x'85'
 dc al3(dcbout)
exl1 dc x'07'
 dc al3(warea)
 dc x'85'
 dc al3(sdcb)
exl2 dc x'07'
 dc al3(warea)
 dc x'85'
 dc al3(sdcc)
exl3 dc x'07'
 dc al3(warea)
 dc x'85'
 dc al3(sdcd)
bask ds f
cask ds f
dask ds f
mask dc xl4'0c000000'
dcbtab dc a(ed0)
 dc a(ed1)
 dc a(ed2)
 dc a(ed3)
 dc a(ed4)
 dc a(ed5)
 dc a(ed6)
 dc a(ed7)
 dc a(mt0)
 dc a(mt1)
 dc a(mt2)
 dc a(mt3)
 dc a(mt4)
 dc a(mt5)
 dc a(mt6)
 dc a(mt7)
vblock dc v(disc)
vbufa dc v(bufa)
currecb dc f'-1'
curwecb dc f'-1'
curr dc f'-1'
curs dc f'-1'
curt dc f'-1'
ecbwr dc 2f'0'
 dc a(sort)
 dc 4f'0'
ecbws dc 2f'0'
 dc a(port)
 dc 4f'0'
ecbwt dc 2f'0'
 dc a(data)
 dc 4f'0'
xarea dc a(sort)
 dc xl1'80'
yarea dc a(port)
 dc xl1'80'
zarea dc a(data)
 dc xl1'80'
ed0 dcb ddname=ed0,dsorg=da,macrf=(ri,wi),optcd=r,exlst=exlsd
ed1 dcb ddname=ed1,dsorg=da,macrf=(ri,wi),optcd=r,exlst=exlsd
ed2 dcb ddname=ed2,dsorg=da,macrf=(ri,wi),optcd=r,exlst=exlsd
ed3 dcb ddname=ed3,dsorg=da,macrf=(ri,wi),optcd=r,exlst=exlsd
ed4 dcb ddname=ed4,dsorg=da,macrf=(ri,wi),optcd=r,exlst=exlsd
ed5 dcb ddname=ed5,dsorg=da,macrf=(ri,wi),optcd=r,exlst=exlsd
ed6 dcb ddname=ed6,dsorg=da,macrf=(ri,wi),optcd=r,exlst=exlsd
ed7 dcb ddname=ed7,dsorg=da,macrf=(ri,wi),optcd=r,exlst=exlsd
mt0      dcb   ddname=ft90f001,dsorg=ps,macrf=(r,w),exlst=exlst
mt1      dcb   ddname=ft91f001,dsorg=ps,macrf=(r,w),exlst=exlst
mt2      dcb   ddname=ft92f001,dsorg=ps,macrf=(r,w),exlst=exlst
mt3      dcb   ddname=ft93f001,dsorg=ps,macrf=(r,w),exlst=exlst
mt4      dcb   ddname=ft94f001,dsorg=ps,macrf=(r,w),exlst=exlst
mt5      dcb   ddname=ft95f001,dsorg=ps,macrf=(r,w),exlst=exlst
mt6      dcb   ddname=ft96f001,dsorg=ps,macrf=(r,w),exlst=exlst
mt7      dcb   ddname=ft97f001,dsorg=ps,macrf=(r,w),exlst=exlst
sort dcb ddname=sort,dsorg=da,macrf=(ri,wi),optcd=r,exlst=exl1
port dcb ddname=port,dsorg=da,macrf=(ri,wi),optcd=r,exlst=exl2
data dcb ddname=data,dsorg=da,macrf=(ri,wi),optcd=r,exlst=exl3
ddisc dsect
isel ds f
iselr ds f
iselw ds f
reply ds f
ichek ds f
ipos ds 16f
blksiz ds 16f
 dcbd dsorg=(ps,da),devd=(da,ta)
 end
_ENDIF
_IF(assem)
*member name = gmake
gmake csect
a equ 12
p equ 13
ij equ 1
ik equ 2
il equ 3
kl equ 4
jk equ 5
jl equ 10
i equ 6
j equ 7
k equ 8
l equ 9
ikyi equ il
ikyk equ jk
ikyj equ jl
gg1 equ 6
g2 equ 6
gg3 equ 4
gg4 equ 0
gg9 equ 2
g8 equ 4
 using *,15
 stm 2,14,sarea
 using dumg,14
 lm 11,14,zapg
 l 0,mword
p6000 lm i,l,zero
 lr 1,0
 ar 1,1
 ic i,iijj-2(1)
 ic j,iijj-1(1)
 ic k,kkll-2(1)
 ic l,kkll-1(1)
 sll 1,2
 ld gg4,g-8(1)
 sldl i,3
 sldl k,3
 ldr gg9,gg4
 lcdr gg4,gg4
 std gg4,gg2
 l ikyi,iky-8(i)
 l ikyk,iky-8(k)
 l ikyj,iky-8(j)
 std gg4,gg5
 lr ij,j
 ar ij,ikyi
 lr ik,k
 ar ik,ikyi
 std gg4,gg6
 lr kl,ikyk
 ar kl,l
 adr gg9,gg9
 cr i,j
 std gg4,gg7
 ldr gg3,gg4
 ldr gg1,gg4
 std gg9,gg8
 bne p1
 swr gg3,gg3 gg3=0
 std gg3,gg6 gg6=0
 swr gg4,gg4 gg4=0
 lcdr gg9,gg1 gg9=gg
 std gg3,gg7 gg7=0
p1 cr ij,kl
 bne p2
 swr gg9,gg9 gg9=0
 std gg9,gg5
 std gg9,gg6
 std gg9,gg7
p2 cr k,l
 bne p3
 lcdr gg4,gg1
 std gg4,gg8 gg8=gg
 swr gg4,gg4 gg4=0
 std gg4,gg2 gg2=0
 std gg4,gg7 gg7=0
p3 sr i,k
 bxh i,11,p33
 ad gg1,gg5 gg1=gg1+gg5
p33 cr j,k
 bl p4
 bh p6
 ad gg3,gg6 gg3=gg3+gg6
 b p7
p4 ld gg3,gg6 gg3=gg6
p7 ar jk,j
 cr j,l
 bl p8
 bh p10
 ad gg4,gg7 gg4=gg4+gg7
 b p10
p8 ld gg4,gg7  gg4=gg7
 l jl,iky-8(l)
 ar jl,j
 b p11
p6 lr jk,ikyj
 ar jk,k
p10 ar jl,l
p11 md gg1,0(p,jl)
 adr gg9,gg9
 ar il,l
 md gg3,0(p,il)
 ad gg1,0(a,ik)
 std gg1,0(a,ik)
 md gg9,0(p,ij)
 ad gg3,0(a,jk)
 ld g2,gg2
 std gg3,0(a,jk)
 md g2,0(p,jk)
 ad gg9,0(a,kl)
 ld g8,gg8
 std gg9,0(a,kl)
 md gg4,0(p,ik)
 ad g2,0(a,il)
 adr g8,g8
 std g2,0(a,il)
 md g8,0(p,kl)
 ad gg4,0(a,jl)
 std gg4,0(a,jl)
 ad g8,0(a,ij)
 std g8,0(a,ij)
 bct 0,p6000
 lm 2,14,sarea
 br 14
 entry gmake1
 using *,15
gmake1 stm 2,3,sarea
 lm 0,1,0(1)
 stm 0,1,abas
 lm 0,3,mmm
q1 st 0,iky(1)
 ar 1,2
 ar 0,1
 bct 3,q1
 lm 2,3,sarea
 br 14
 extrn blkin
gg2 ds d
gg5 ds d
gg6 ds d
gg7 ds d
gg8 ds d
sarea ds 7d
mmm dc f'-8,0,8,256'
zapg dc f'0'
abas ds f
pbas ds f
gbas dc a(blkin-8)
zero dc 4f'0'
iky ds 256d
dumg dsect
gspa ds d
g ds 340d
iijj ds 340h
kkll ds 340h
mword ds f
 end
*member name = gsup
gsup csect
h equ 2
p equ 3
ij equ 4
kl equ 5
 using *,15
 stm 2,6,sarea
 lm 1,3,base
 using dumg,1
 l 0,mword
p6000 lr 6,0
 ar 6,6
 lh 4,iijj-2(6)
 lh 5,kkll-2(6)
 sldl 4,3
 sll 6,2
 ld 0,g-8(6)
 ldr 2,0
 md 0,0(p,kl)
 md 2,0(p,ij)
 ad 0,0(h,ij)
 std 0,0(h,ij)
 ad 2,0(h,kl)
 std 2,0(h,kl)
 bct 0,p6000
 lm 2,6,sarea
 br 14
 entry gsup1
 using *,15
gsup1 lm 0,1,0(1)
 sh 0,m8
 sh 1,m8
 stm 0,1,hp
 br 14
 extrn blkin
base dc a(blkin-8)
hp ds 2f
m8 dc h'8'
sarea ds 3d
dumg dsect
gspa ds d
g ds 340d
iijj ds 340h
kkll ds 340h
mword ds f
 end
*member name = gsupa
gsupa csect
h equ 2
p equ 3
ij equ 4
kl equ 5
 using *,15
 stm 2,6,sarea
 lm 1,3,base
 using dumg,1
 l 0,mword
p6000 lr 6,0
 ar 6,6
 lh 4,iijj-2(6)
 lh 5,kkll-2(6)
 sldl 4,3
 sll 6,2
 ld 0,g-8(6)
 ldr 2,0
 md 0,0(p,kl)
 md 2,0(p,ij)
 ad 0,0(h,ij)
 std 0,0(h,ij)
 ad 2,0(h,kl)
 std 2,0(h,kl)
 bct 0,p6000
 lm 2,6,sarea
 br 14
 entry gsupa1
 using *,15
gsupa1 lm 0,1,0(1)
 sh 0,m8
 sh 1,m8
 stm 0,1,hp
 br 14
 extrn blkin
base dc a(blkin-8)
hp ds 2f
m8 dc h'8'
sarea ds 3d
dumg dsect
gspa ds d
g ds 408d
iijj ds 204h
kkll ds 204h
mword ds f
 end
*member name = jksupe
jksupe csect
iw equ 1
ij equ 2
kl equ 3
hc equ 4
h1 equ 5
h2 equ 6
pc equ 7
p equ 8
nshell equ 9
lenbas equ 10
basb equ 11
ishell equ 12
gj equ 4
gk equ 6
  using *,15
  extrn blkin
  stm 2,13,sarea
  lm hc,basb,hcbas
  using dumg,basb
  l iw,nint+8
p999 lr 13,iw
  ar 13,13
  lh ij,intij+6(13)
  lh kl,intkl+6(13)
  sll 13,2
  sldl ij,3
  ld gj,gjj(13)
  ld 0,0(pc,ij)
  ld 2,0(pc,kl)
  mdr 0,gj
  mdr 2,gj
  ad 0,0(hc,kl)
  ld gk,gkk(13)
  std 0,0(hc,kl)
  ad 2,0(hc,ij)
  lr ishell,nshell
  std 2,0(hc,ij)
p888 ld 2,0(p,ij)
  ldr 0,gj
  mdr 0,2
  mdr 2,gk
  ad 0,0(h1,kl)
  ad 2,0(h2,kl)
  std 0,0(h1,kl)
  ld 0,0(p,kl)
  std 2,0(h2,kl)
  ldr 2,0
  mdr 0,gj
  mdr 2,gk
  ad 0,0(h1,ij)
  ad 2,0(h2,ij)
  ar kl,lenbas
  std 0,0(h1,ij)
  std 2,0(h2,ij)
  ar ij,lenbas
  bct ishell,p888
  bct iw,p999
  lm 2,13,sarea
  br 14
  entry jksup1
  using *,15
jksup1 stm 2,8,sarea
  lm 2,8,0(1)
  la 1,8(0,0)
  l 7,0(0,7)
  l 8,0(0,8)
  sr 2,1
  sr 3,1
  sr 4,1
  sr 5,1
  sr 6,1
  sll 8,3
  stm 2,8,hcbas
  lm 2,8,sarea
  br 14
sarea ds 12f
hcbas ds f
h1bas ds f
h2bas ds f
pcbas ds f
pbas ds f
nshel ds f
lenb ds f
bas dc a(blkin-8)
dumg dsect
gjj ds 204d
gkk ds 204d
intij ds 204h
intkl ds 204h
nint ds f
  end
*member name = proc2
proc2 csect
a equ 12
p equ 13
b equ 14
q equ 6
ij equ 1
ik equ 2
il equ 3
kl equ 4
jk equ 5
jl equ 10
i equ 6
j equ 7
k equ 8
l equ 9
ikyi equ il
ikyk equ jk
ikyj equ jl
gg1 equ 6
g2 equ 6
gg3 equ 4
g8 equ 4
gg4 equ 0
ggg3 equ 0
ggg2 equ 0
gg9 equ 2
ggg4 equ  2
ggg1 equ 2
 using *,15
 stm 2,14,sarea
 using dumg,11
 lm 11,14,gapb
 l 0,mword
p6000 lm i,l,zero
 lr 1,0
 ar 1,1
 ic i,iijj-2(1)
 ic j,iijj-1(1)
 ic k,kkll-2(1)
 ic l,kkll-1(1)
 sll 1,2
 ld gg4,g-8(1)
 sldl i,3
 sldl k,3
 ldr gg9,gg4
 lcdr gg4,gg4
 std gg4,gg2
 l ikyi,iky-8(i)
 l ikyk,iky-8(k)
 l ikyj,iky-8(j)
 std gg4,gg5
 lr ij,j
 ar ij,ikyi
 lr ik,k
 ar ik,ikyi
 std gg4,gg6
 lr kl,ikyk
 ar kl,l
 adr gg9,gg9
 cr i,j
 std gg4,gg7
 ldr gg3,gg4
 ldr gg1,gg4
 std gg9,gg8
 bne p1
 swr gg3,gg3 gg3=0
 std gg3,gg6 gg6=0
 swr gg4,gg4 gg4=0
 lcdr gg9,gg1 gg9=gg
 std gg3,gg7 gg7=0
p1 cr ij,kl
 bne p2
 swr gg9,gg9 gg9=0
 std gg9,gg5
 std gg9,gg6
 std gg9,gg7
p2 cr k,l
 bne p3
 lcdr gg4,gg1
 std gg4,gg8 gg8=gg
 swr gg4,gg4 gg4=0
 std gg4,gg2 gg2=0
 std gg4,gg7 gg7=0
p3 cr i,k
 l q,qbas
 bne p33
 ad gg1,gg5 gg1=gg1+gg5
p33 cr j,k
 bl p4
 bh p6
 ad gg3,gg6 gg3=gg3+gg6
 b p7
p4 ld gg3,gg6 gg3=gg6
p7 ar jk,j
 cr j,l
 bl p8
 bh p10
 ad gg4,gg7 gg4=gg4+gg7
 b p10
p8 ld gg4,gg7  gg4=gg7
 l jl,iky-8(l)
 ar jl,j
 b p11
p6 lr jk,ikyj
 ar jk,k
p10 ar jl,l
p11 adr gg9,gg9
 md gg9,0(p,ij)
 ad gg9,0(a,kl)
 ar il,l
 std gg9,0(a,kl)
 ldr ggg4,gg4
 md gg4,0(p,ik)
 md ggg4,0(q,ik)
 ad gg4,0(a,jl)
 std gg4,0(a,jl)
 ldr ggg3,gg3
 md gg3,0(p,il)
 ad ggg4,0(b,jl)
 std ggg4,0(b,jl)
 md ggg3,0(q,il)
 ad gg3,0(a,jk)
 ldr ggg1,gg1
 std gg3,0(a,jk)
 md gg1,0(p,jl)
 ad ggg3,0(b,jk)
 ld g8,gg8
 std ggg3,0(b,jk)
 md ggg1,0(q,jl)
 ad gg1,0(a,ik)
 adr g8,g8
 ld ggg2,gg2
 std gg1,0(a,ik)
 md g8,0(p,kl)
 ad ggg1,0(b,ik)
 ldr g2,ggg2
 std ggg1,0(b,ik)
 md g2,0(p,jk)
 ad g8,0(a,ij)
 std g8,0(a,ij)
 md ggg2,0(q,jk)
 ad g2,0(a,il)
 ad ggg2,0(b,il)
 std g2,0(a,il)
 std ggg2,0(b,il)
 bct 0,p6000
 lm 2,14,sarea
 br 14
 entry proc21
 using *,15
proc21 stm 2,3,sarea
 lm 0,3,0(1)
 stm 0,3,abas
 lm 0,3,mmmm
q1 st 0,iky(1)
 ar 1,2
 ar 0,1
 bct 3,q1
 lm 2,3,sarea
 br 14
 extrn blkin
gg2 ds d
gg5 ds d
gg6 ds d
gg7 ds d
gg8 ds d
sarea ds 7d
mmmm dc f'-8,0,8,256'
gapb dc a(blkin-8)
abas ds f
pbas ds f
bbas ds f
qbas ds f
zero dc 4f'0'
iky ds 256d
dumg dsect
gspa ds d
g ds 340d
iijj ds 340h
kkll ds 340h
mword ds f
 end
*member name = vecsum
vecsum csect
*
*      function vecsum(a,b,n)
*      REAL  a(1),b(1)
*      vecsum=0.0d0
*      n4=(n/4)*4
*      if(n4.le.0)goto 2
*      do 1 i=1,n4,4
* 1    vecsum=a(i)*b(i)+a(i+1)*b(i+1)
*      >           +a(i+2)*b(i+2)+a(i+3)*b(i+3)+vecsum
* 2    if(n.le.n4)return
*      n4=n4+1
*      do 3 i=n4,n
* 3    vecsum=a(i)*b(i)+vecsum
*      return
*
 using *,15
 stm 2,6,28(13)
 lm 4,6,0(1)
 sdr 0,0
 lm 1,3,const
 a 3,0(6)
 sll 3,3
 bm vecsum2
vecsum1 ld 2,0(1,4)
 ld 4,8(1,4)
 md 2,0(1,5)
 adr 0,2
 md 4,8(1,5)
 adr 0,4
 ld 2,16(1,4)
 ld 4,24(1,4)
 md 2,16(1,5)
 adr 0,2
 md 4,24(1,5)
 adr 0,4
 bxle 1,2,vecsum1
*
vecsum2 a 3,=f'24'
 la 2,8
 sr 1,2
vecsum4 bxh 1,2,vecsum3
 ld 2,0(1,4)
 md 2,0(1,5)
 adr 0,2
 b vecsum4
vecsum3 lm 2,6,28(13)
 br 14
 entry triad
* subroutine triad(n,a,y,x)
* REAL  x(1),y(1),a
* do 1 i=1,n
*1 y(i)=x(i)*a + y(i)
* return
* end
 using *,15
triad stm 2,7,28(13)
 l 7,0(1)
 l 5,4(1)
 l 6,8(1)
 l 4,12(1)
 ld 0,0(5)
 lm 1,3,const
 a 3,0(7)
 sll 3,3
 bm axpy2
axpy1 ld 2,0(1,4)
 ld 4,8(1,4)
 mdr 2,0
 mdr 4,0
 ad 2,0(1,6)
 ad 4,8(1,6)
 std 2,0(1,6)
 std 4,8(1,6)
 ld 2,16(1,4)
 ld 4,24(1,4)
 mdr 2,0
 mdr 4,0
 ad 2,16(1,6)
 ad 4,24(1,6)
 std 2,16(1,6)
 std 4,24(1,6)
 bxle 1,2,axpy1
axpy2 a 3,=f'24'
 la 2,8
 sr 1,2
axpy4 bxh 1,2,axpy3
 ld 2,0(1,4)
 mdr 2,0
 ad 2,0(1,6)
 std 2,0(1,6)
 b axpy4
axpy3 lm 2,7,28(13)
 br 14
_IFN(charmm)
 entry addvec
* subroutine addvec(x,y,z,n)
* REAL  x(1),y(1),z(1)
* do 1 i=1,n
*1 x(i)=y(i)+z(i)
* return
 using *,15
addvec stm 2,7,28(13)
 lm 4,7,0(1)
 l  7,0(0,7)
 ltr 7,7
 bz addvec3
 lm 1,3,const
 ar 3,7
 sll 3,3
 bm addvec2
addvec1 ld 0,0(1,5)
 ld 2,8(1,5)
 ld 4,16(1,5)
 ld 6,24(1,5)
 ad 0,0(1,6)
 ad 2,8(1,6)
 ad 4,16(1,6)
 ad 6,24(1,6)
 std 0,0(1,4)
 std 2,8(1,4)
 std 4,16(1,4)
 std 6,24(1,4)
 bxle 1,2,addvec1
addvec2 a 3,=f'24'
 la 2,8
 sr 1,2
addvec4 bxh 1,2,addvec3
 ld 0,0(1,5)
 ad 0,0(1,6)
 std 0,0(1,4)
 b addvec4
addvec3 lm 2,7,28(13)
 br 14
_ENDIF
 entry add3v
* subroutine add3v(w,x,y,z,n)
* REAL  w(1),x(1),y(1),z(1)
* do 1 i=1,n
*1 w(i)=x(i)+y(i)+z(i)
* return
 using *,15
add3v stm 2,8,28(13)
 lm 4,8,0(1)
 lm 1,3,const
 a 3,0(8)
 sll 3,3
 bm add3v2
add3v1 ld 0,0(1,5)
 ld 2,8(1,5)
 ld 4,16(1,5)
 ld 6,24(1,5)
 ad 0,0(1,6)
 ad 2,8(1,6)
 ad 4,16(1,6)
 ad 6,24(1,6)
 ad 0,0(1,7)
 ad 2,8(1,7)
 ad 4,16(1,7)
 ad 6,24(1,7)
 std 0,0(1,4)
 std 2,8(1,4)
 std 4,16(1,4)
 std 6,24(1,4)
 bxle 1,2,add3v1
add3v2 a 3,=f'24'
 la 2,8
 sr 1,2
add3v4 bxh 1,2,add3v3
 ld 0,0(1,5)
 ad 0,0(1,6)
 ad 0,0(1,7)
 std 0,0(1,4)
 b add3v4
add3v3 lm 2,8,28(13)
 br 14
 entry szero
 entry szero4
* subroutine szero(z,n)
* REAL  z(1)
* do 1 i=1,n
*1 z(i)=0.0d0
* return
 using *,15
szero4 stm 2,5,28(13)
 lm 4,5,0(1)
 l 3,0(5)
 b zero0
 using *,15
szero stm 2,5,28(13)
 lm 4,5,0(1)
 l 3,0(5)
zero0 balr 5,0
 using *,5
 lm 1,2,const
 sdr 0,0
 a 3,const+8
 sll 3,3
 bm zero2
zero1 std 0,0(1,4)
 std 0,8(1,4)
 std 0,16(1,4)
 std 0,24(1,4)
 bxle 1,2,zero1
zero2 a 3,=f'24'
 la 2,8
 sr 1,2
zero5 bxh 1,2,zero3
 std 0,0(1,4)
 b zero5
zero3 lm 2,5,28(13)
 br 14
 entry switch
* subroutine switch (i,j)
* k=i
* i=j
* j=k
* return
switch st 2,28(13)
 l 2,4(1)
 l 1,0(1)
 l 0,0(2)
 mvc 0(4,2),0(1)
 st 0,0(1)
 l 2,28(13)
 br 14
const dc f'0'
 dc f'32'
 dc f'-4'
 end
_ENDIF
_IF(mvs)
*member name = mvsclk
**********************************************************************
*
*             timing routine for mvs/xa operating system
*
*  syntax:    "call timer (cpu,elapse) "  or optionally
*             "call timer (cpu,elapse,cputot,eltot) "
*
*  function:  returns elapsed and cpu time in floating point seconds.
*             the first two arguments give the time from the
*             last call.  the last two arguemnts are optional and
*             give the total times from the first call.
*
*             accuracy at micro second level, but the time for
*             for one call of timer is about 40 micro-seconds on
*             an ibm 3084 running mvs/xa.  therefore the accuracy
*             on the ibm 3084 is at least 50 micro-seconds.
*
*  input:     2 or 4 double words of storage
*
*
*  logic:     gets tod clock via the stck instruction
*             gets cpu clock via the ttimer macro
*             converts clock units to micro seconds
*             converts micro seconds to seconds
*             subtracts current time from time when last called
*             returns result of subtraction
*             optionally returns total times
*
*
*   note:     timer interval set to 19.088 hours by stimer macro.
*             if timer runs for more than 19.088 hours then
*             it is reset to 19.088 hours by another stimer call.
*
**********************************************************************
f0       equ   0
f2       equ   2
f4       equ   4
f6       equ   6
r0       equ   0
r1       equ   1
r2       equ   2
r3       equ   3
r4       equ   4
r5       equ   5
r6       equ   6
r7       equ   7
r8       equ   8
r9       equ   9
r10      equ   10
r11      equ   11
r12      equ   12
r13      equ   13
r14      equ   14
r15      equ   15
timer    csect
         using *,12
         stm   r14,r12,12(r13)     save the regs
         lr    r12,r15             get addressability
         stck  clock               get elapsed time
         mvi   clock,x'4e'         establish the exponent
         ld    f0,clock            clock units to fp representation
         dd    f0,fp4096           convert clock units to micro seconds
         md    f0,millth           convert micro secs to seconds
         ld    f2,savclk           get the last clock value
         std   f0,savclk           save for the next time
         sdr   f0,f2               figure the difference
         l     r2,4(r1)            get address of the parameter
         std   f0,0(r2)            return elapsed seconds
         ld    f6,savcpu           get last value of cputime
         ltdr  f6,f6               check if stimer has been set
         bp    allset              branch if allready set
*        stimer task,micvl=tfhrs   set timer interval to 24 hours
         la    r1,tfhrs
         sr    r15,r15
         la    r0,160(r0,r0)
         sll   r0,24(0)
         svc   47
         mvc   savcpu(8),fptfhrs   set remaining interval to 24 hours
         lm    r14,r12,12(r13)
         mvi   12(r13),x'ff'
         sr    r15,r15
         br    r14
*llset   ttimer ,mic,cpu           get cpu time
allset   la    r0,cpu
         la    r1,2(r0,r0)
         svc   46
         l     r1,24(r13)          restore arguement pointer
         mvi   cpu,x'4e'           establish exponent, positive sign
         ld    f4,cpu              convert clock units to fp
         dd    f4,fp4096           convert clock units to micro seconds
         ld    f6,savcpu           get last cpu call
         std   f4,savcpu           save for next time
         sdr   f6,f4               get difference
         md    f6,millth           convert micro seconds to seconds
         l     r3,0(r1)            get address for cpu argument
         std   f6,0(r3)            return cpu seconds
         ltdr  f2,f2               check if first call (savclk = 0)
         bz    first               if first call do not sum
         ad    f0,eltot            sum total elapsed time
         std   f0,eltot            store new elapsed sum
         ad    f6,cputot           sum total cpu time
         std   f6,cputot           store new cpu sum
first    ltr   r2,r2               check if only two arguments
         bm    finit
         l     r5,8(r1)            get address of cputot
         std   f6,0(r5)            return cpu total to caller
         l     r6,12(r1)           get address of elapsed total
         std   f0,0(r6)            return elapsed total to caller
finit    cd    f4,litleft          compare remaining timer interval
*                                  with a small value (1 hour)
*                                  if timer interval is less than
*                                  1 hour reset interval
         bh    noreset
*        stimer task,micvl=tfhrs   set timer interval to 19 hours
         la    r1,tfhrs
         sr    r15,r15
         la    r0,160(r0,r0)
         sll   r0,24(0)
         svc   47
         mvc   savcpu(8),fptfhrs   set remaining interval to 19 hours
noreset  lm    r14,r12,12(r13)
         mvi   12(r13),x'ff'
         sr    r15,r15
         br    r14
         ds    0d                       align on dword boundry
millth   dc    d'1e-6'                  to convert micro secs to secs
fp4096   dc    x'44',x'1000',xl5'00'  number of clock units/micro sec
clock    dc    x'41',xl7'00'            target of stck
cpu      dc    x'41',xl7'00'            target of stpt
savclk   dc    x'00',xl7'00'            save clock value for next call
savcpu   dc    x'00',xl7'00'            save cpu value for next call
eltot    dc    x'41',xl7'00'            total elapsed time
cputot   dc    x'41',xl7'00'            total cpu time
tfhrs    dc    x'0000ffffffffffff'      timer interval=19.088 hours
fptfhrs  dc    x'4e00000fffffffff'      floating point timer interval
litleft  dc    d'3.6d+9'                one hour
         end
*member name = datim
datim    csect
*
* this module returns the current date, time to the caller in
* printable 8 byte characters. the call format is:
*
*     call datim(date,time)
*
* where date and time are both 8 bytes long. (in fortran, they can
* be declared as REAL  or character*8 variables)
*
         using datim,15
         b     *+84
         dc    x'00',cl7'datim'
save     dc    18f'0'
         stm   14,12,12(13)            save caller's regs
         la    12,12(0,15)             set base register
         using save,12
         drop  15
         st    13,4(0,12)              chain save area back
         st    12,8(0,13)              chain save area forward
         lr    13,12
         lm    8,9,0(1)                save argument addresses
         time  dec                     get time of day
*                                  r0= hh mm ss th  (packed dec digits)
*                                  r1= 00 yy dd ds  (packed dec)
         srl   0,4                     remove fractional seconds
*                                  r0= 0h hm ms st
         st    0,work                  save time for edit
         oi    work+3,x'0f'            or in plus sign bits
         mvc   timework(10),timemask
         ed    timework(10),work       edit hour, min, sec
         mvc   0(8,9),timework+2       move to 2nd argument
         space 1
         lr    0,1                     save year in r0
         n     1,rthalf            r1= 00 00 dd df  (packed dec days)
         srl   0,12                    shift out julian days
         o     0,=f'15'            r0= 00 00 0y yf  (packed dec year)
         st    0,work                  save for leap year calc
         unpk  5(3,8),work+2(2)        unpack year into 1st argument
         mvi   2(8),c'/'                and move in slashes
         mvi   5(8),c'/'
         dp    work(4),fourp(1)        divide year by four, and
         tm    work+3,x'f0'            test remainder eq zero
         la    9,monthtab                assume not leap year
         bnz   *+8                      ~zero, use stand month table
         la    9,leaptab                 zero, use leap year table
         xc    work(4),work            zero work area
         st    1,work+4                convert julian days to binary
         cvb 1,work
         bctr  1,0                      decrement by 1
         sll   1,2                      times four
         la    9,0(1,9)                add to table entry
         mvc   0(2,8),0(9)             move month to 1st argument
         mvc   3(2,8),2(9)             move day   to 1st argument
return   l     13,4(0,13)
         lm    14,12,12(13)
         mvi   12(13),x'ff'
         sr    15,15
         br    14
work     dc    d'0'
rthalf   dc    x'0000ffff'             mask for right half
fourp    dc    x'4f'                   packed decimal four
timework dc    x'00000000000000000000' edit work area
timemask dc    x'402021207a20207a2020' edit mask (hh:mm:ss)
         ltorg
monthtab dc    c'01010102010301040105010601070108' january
         dc    c'01090110011101120113011401150116'
         dc    c'01170118011901200121012201230124'
         dc    c'0125012601270128012901300131'
         dc    c'02010202020302040205020602070208' february
         dc    c'02090210021102120213021402150216'
         dc    c'02170218021902200221022202230224'
         dc    c'0225022602270228'
         dc    c'03010302030303040305030603070308' march
         dc    c'03090310031103120313031403150316'
         dc    c'03170318031903200321032203230324'
         dc    c'0325032603270328032903300331'
         dc    c'04010402040304040405040604070408' april
         dc    c'04090410041104120413041404150416'
         dc    c'04170418041904200421042204230424'
         dc    c'042504260427042804290430'
         dc    c'05010502050305040505050605070508' may
         dc    c'05090510051105120513051405150516'
         dc    c'05170518051905200521052205230524'
         dc    c'0525052605270528052905300531'
         dc    c'06010602060306040605060606070608' june
         dc    c'06090610061106120613061406150616'
         dc    c'06170618061906200621062206230624'
         dc    c'062506260627062806290630'
         dc    c'07010702070307040705070607070708' july
         dc    c'07090710071107120713071407150716'
         dc    c'07170718071907200721072207230724'
         dc    c'0725072607270728072907300731'
         dc    c'08010802080308040805080608070808' august
         dc    c'08090810081108120813081408150816'
         dc    c'08170818081908200821082208230824'
         dc    c'0825082608270828082908300831'
         dc    c'09010902090309040905090609070908' september
         dc    c'09090910091109120913091409150916'
         dc    c'09170918091909200921092209230924'
         dc    c'092509260927092809290930'
         dc    c'10011002100310041005100610071008' october
         dc    c'10091010101110121013101410151016'
         dc    c'10171018101910201021102210231024'
         dc    c'1025102610271028102910301031'
         dc    c'11011102110311041105110611071108' november
         dc    c'11091110111111121113111411151116'
         dc    c'11171118111911201121112211231124'
         dc    c'112511261127112811291130'
         dc    c'12011202120312041205120612071208' december
         dc    c'12091210121112121213121412151216'
         dc    c'12171218121912201221122212231224'
         dc    c'1225122612271228122912301231'
leaptab  dc    c'01010102010301040105010601070108' january
         dc    c'01090110011101120113011401150116'
         dc    c'01170118011901200121012201230124'
         dc    c'0125012601270128012901300131'
         dc    c'02010202020302040205020602070208' february
         dc    c'02090210021102120213021402150216'
         dc    c'02170218021902200221022202230224'
         dc    c'02250226022702280229'
         dc    c'03010302030303040305030603070308' march
         dc    c'03090310031103120313031403150316'
         dc    c'03170318031903200321032203230324'
         dc    c'0325032603270328032903300331'
         dc    c'04010402040304040405040604070408' april
         dc    c'04090410041104120413041404150416'
         dc    c'04170418041904200421042204230424'
         dc    c'042504260427042804290430'
         dc    c'05010502050305040505050605070508' may
         dc    c'05090510051105120513051405150516'
         dc    c'05170518051905200521052205230524'
         dc    c'0525052605270528052905300531'
         dc    c'06010602060306040605060606070608' june
         dc    c'06090610061106120613061406150616'
         dc    c'06170618061906200621062206230624'
         dc    c'062506260627062806290630'
         dc    c'07010702070307040705070607070708' july
         dc    c'07090710071107120713071407150716'
         dc    c'07170718071907200721072207230724'
         dc    c'0725072607270728072907300731'
         dc    c'08010802080308040805080608070808' august
         dc    c'08090810081108120813081408150816'
         dc    c'08170818081908200821082208230824'
         dc    c'0825082608270828082908300831'
         dc    c'09010902090309040905090609070908' september
         dc    c'09090910091109120913091409150916'
         dc    c'09170918091909200921092209230924'
         dc    c'092509260927092809290930'
         dc    c'10011002100310041005100610071008' october
         dc    c'10091010101110121013101410151016'
         dc    c'10171018101910201021102210231024'
         dc    c'1025102610271028102910301031'
         dc    c'11011102110311041105110611071108' november
         dc    c'11091110111111121113111411151116'
         dc    c'11171118111911201121112211231124'
         dc    c'112511261127112811291130'
         dc    c'12011202120312041205120612071208' december
         dc    c'12091210121112121213121412151216'
         dc    c'12171218121912201221122212231224'
         dc    c'1225122612271228122912301231'
         end
_ENDIF
_IF(vmcms)
*member name = elaps
**********************************************************************
* elapse: invoked by "call elapse(t)"
*
*  function: returns elapsed time in floating point seconds
*            from the last call
*
*            accuracy at micro second.
*
*  input: double word of storage
*
*  output: elapsed time since last call in fp seconds.
*          data from first call should be thrown away.
*
*  logic:    -gets tod clock via the stck instruction
*            -converts clock units to micro seconds
*            -converts micro seconds to seconds
*            -subtracts current time from time when last called
*            -returns result of subtraction
*
*   note: 1000 hex clock units equal 1 micro second on the 43xx
*         processors.
*
**********************************************************************
f0       equ   0
f2       equ   2
r0       equ   0
r1       equ   1
r2       equ   2
r3       equ   3
r4       equ   4
r5       equ   5
r6       equ   6
r7       equ   7
r8       equ   8
r9       equ   9
r10      equ   10
r11      equ   11
r12      equ   12
r13      equ   13
r14      equ   14
r15      equ   15
elapse   csect
         using *,12
         stm   r14,r12,12(r13)   save the regs
         lr    r12,r15          get addressability
         stck  clock        get clock
         mvi   clock,x'4e'  establish the exponent
         ld    f0,clock     clock units to fp representation
         dd    f0,fp4096    convert clock units to micro seconds
         md    f0,millth    convert micro secs to seconds
         ld    f2,savclk    get the last clock value
         std   f0,savclk    save for the next time
         sdr   f0,f2        figure the difference
         l     r2,0(r1)     get address of the parameter
         std   f0,0(r2)     return elapsed seconds
         lm    r14,r12,12(r13)
         mvi   12(r13),x'ff'
         sr    r15,r15
         br    r14
         ds    0d                     align on dword boundry
millth   dc    d'1e-6'                to convert micro secs to secs
fp4096   dc    x'44',x'1000',xl5'00'  number of clock units/micro sec
clock    dc    x'41',xl7'00'          target of stck
savclk   dc    x'41',xl7'00'          save clock value for next call
         end
*member name = vmclk
*   subroutine clock(time)
         space
*        definition of register references & assembly parameter
         space
r0       equ   0
r1       equ   1
r2       equ   2
r3       equ   3
r13      equ   13
r14      equ   14
r15      equ   15
f0       equ   0
f2       equ   2
exponent equ   78
         space
         using clock,r15
clock    csect                     enter from calling program
         stm   r0,r3,20(r13)       save registers
         l     r1,0(r1)            get &clock
         lr    r0,r1
         dc    x'8300000c'         diagnose to obtain
*                                  time(0)=date    time(1)=time
*                                  time(2)=cputim  time(3)=tottim
         sr    r2,r2               i=0
         la    r3,8                step=1; imax=1
loop     mvi   16(r1),exponent     convert cputim & tottim to
         ld    f0,16(r1)           d-floating format & scale by 10**-6
         md    f0,scale            to get times in seconds
         std   f0,16(r1)           post cputim & tottim
         ld    f2,lasttime(r2)
         std   f0,lasttime(r2)     update lasttime with cputim & tottim
         ltdr  f2,f2
         bz    loopend             if no prior clock call, skip
         sdr   f0,f2               calculate delcpu & deltot
         std   f0,32(r1)           time(4)=delcpu  time(5)=deltot
loopend  ar    r1,r3               i=i+1
         bxle  r2,r3,loop          if (i.le.imax) recycle loop
         lm    r0,r3,20(r13)       restore registers
         br    r14                 return to calling program
scale    dc    d'1.0d-6'
lasttime dc    d'0,0'              prior values of cputim & tottim
         end
_ENDIF
*member name = alloca
alloc csect
*alloc rmode 24
 bc 15,10(0,15)
 dc x'5'
 dc cl5'alloc'
 stm 14,6,12(13)
 balr 2,0
 using *,2
  mvi  warea,x'00'
 l 4,vblock
 st 13,sarea+4
 la 13,sarea
 l 14,0(0,1)
*lh 15,0(0,14) 15=atmol channel number
 l 15,0(0,14) 15=atmol channel number
 l 6,4(0,1)  6=buffer address
 sll 15,2
 l 3,dcbtab-4(15) 3=dcb address
 using ihadcb,3
 rdjfcb  ((3))
 cli warea,x'00'
 be nofile
 mvc 0(44,4),warea
 open ((3),(output))
 tm dcboflgs,x'10'
 bz nofile
 lh 15,dcbblksi
 sr 5,5 5=block counter
*sth 15,44(0,4)
 st  15,44(0,4)
loop write decb,sf,(3),(6)
 sth 15,warea
 check decb
 cli warea+1,x'0c'
 be loop
 cli warea+1,x'08'
 la 5,1(0,5)
 bne loop
 bct 5,*+4
*mh 5,dcbblksi
 sr 0,0
 lh 1,dcbblksi
 mr 0,5
 srdl 0,12
*sth 5,46(0,4)
 st  1,48(0,4)
 close ((3))
exit l 13,sarea+4
 return (14,6),t
nofile mvi 46(4),x'ff'
 b exit
dcbout ni dcbblksi+1,x'f8' make blksize multiple of 8
 bcr 15,14
warea ds 44f
sarea ds 18f
vblock dc v(disc)
dcbtab dc a(ed0)
 dc a(ed1)
 dc a(ed2)
 dc a(ed3)
 dc a(ed4)
 dc a(ed5)
 dc a(ed6)
 dc  a(ed7)
exls dc x'05'
 dc al3(dcbout)
 dc x'87'
 dc al3(warea)
ed0 dcb ddname=ed0,dsorg=ps,macrf=(wl),devd=da,recfm=f,exlst=exls
ed1 dcb ddname=ed1,dsorg=ps,macrf=(wl),devd=da,recfm=f,exlst=exls
ed2 dcb ddname=ed2,dsorg=ps,macrf=(wl),devd=da,recfm=f,exlst=exls
ed3 dcb ddname=ed3,dsorg=ps,macrf=(wl),devd=da,recfm=f,exlst=exls
ed4 dcb ddname=ed4,dsorg=ps,macrf=(wl),devd=da,recfm=f,exlst=exls
ed5 dcb ddname=ed5,dsorg=ps,macrf=(wl),devd=da,recfm=f,exlst=exls
ed6 dcb ddname=ed6,dsorg=ps,macrf=(wl),devd=da,recfm=f,exlst=exls
ed7 dcb ddname=ed7,dsorg=ps,macrf=(wl),devd=da,recfm=f,exlst=exls
 dcbd dsorg=(ps,da),devd=(da)
 end
_ENDIF
_IF(convex)


; INITIALIZED DATA

	.data
	.align	8
LI:

; UNINITIALIZED DATA

	.bss
	.align	8
LU:
	bs.b	4; +0, function results
	.data
	ds.w	1(0)

; INSTRUCTIONS

	.text
	ds.w	0x3000000
	ds.b	"-O1\0"
	.globl	_leadz_	;ENTRY
_leadz_:
	sub.w	#0x0000008,sp	; #1, 8
	ld.l	#0x000003f,s1	; #4, 63
	ld.l	@0(ap),s0	; #4, INPUT
	lzc	s0,s0		; #4, KNOT
	sub.l	s0,s1		; #4
	mov	s1,s2		; #4
	st.w	s2,LU	; #4, LEADZ
	ld.w	LU,s0	; #5, LEADZ
	rtn	; #5

; INSTRUCTIONS

	.text
	ds.w	0x3000000
	ds.b	"-O1\0"
	.globl	_popcnt_	;ENTRY
_popcnt_:
	sub.w	#0x0000008,sp	; #1, 8
	ld.l	@0(ap),s0	; #5, INPUT
	plc.t	s0,s0		; #5, KNOT
	mov	s0,s1		; #5
	st.w	s1,LU	; #5, POPCNT
	ld.w	LU,s0	; #6, POPCNT
	rtn	; #6
_ENDIF
_IF(alliant,titan,apollo,sun,sgi,hp700)
| integer function leadz(x.r8) ... alliant scalar version
| returns no. of leading zero bits in 64 bit argument
| rjh 7/26/88
.text
     .globl _leadz_
_leadz_:
     link a6,#0          | linkage for traceback is very cheap
     movl a0@,a2         | get address of argument
     vfirst1 a2@+,d0     | examine most sig. 32 bits
     cmpib #-1,d0
     jne gohome
next:
     vfirst1 a2@,d0      | top 32 bits are zero, look at last 32 bits
     cmpib #-1,d0
     jeq allzero
     addb #32,d0
     bras gohome
allzero:
     moveq #64,d0
gohome:
     unlk a6
     rts| integer function popcnt(x.r8) ... alliant version
| return no. of bits set in 64 bit argument
| rjh 7/26/88
.text
     .globl _popcnt_
_popcnt_:
     link a6,#0             | linkage for traceback is very cheap
     movl a0@,a2            | get address
     vcount1 a2@+,d0
     vcount1 a2@,d1
     addb d1,d0
     unlk a6
     rts
_ENDIF
_IF(vax)
          .title    x03mod
;
; mark 7 revised. nag copyright 1979
; written by j.j.du croz, nag central office
;
; this code is written in a generalized form to allow for different
;  lengths of exponent field.  standard vax-11 hardware has 8-bit
;  exponent field.  to allow for exponent field of different length,
;  simply change the symbolic value of explen.  it is assumed however
;  that 8.le.explen.le.15 .
;
          .list
explen    =         8                   ; number of bits in exponent
expbias   =         1@<explen-1>        ; exponent bias
explim    =         1@explen            ; exponent overflow threshold
exppos    =         15-explen           ; starting position of exponent
impbit    =         1@exppos            ; implicit leading bit
fralen    =         64-explen           ; number of bits in fraction
;
; macros for arithmetic shifts on 8-word fields
;
          .macro    ashoctl   sc,base   ; left shift of sc bits
          ashq      #'sc,8+'base,8+'base
          ashl      #-'sc,8+'base,8+'base
          ashq      #'sc,4+'base,4+'base
          ashl      #-'sc,4+'base,4+'base
          ashq      #'sc,base,base
          .endm
          .macro    ashoctr   sc,base   ; right shift of sc bits
          ashq      #-'sc,base,base
          ashl      #'sc,4+'base,4+'base
          ashq      #-'sc,4+'base,4+'base
          ashl      #'sc,8+'base,8+'base
          ashq      #-'sc,8+'base,8+'base
          .endm
;
; macro to unpack fraction of double precision floating-point number
;
          .macro    unpfra    srce,dest
          rotl      #16,4+'srce,dest
          rotl      #16,srce,4+'dest
          insv      #0,#exppos,#explen,6+'dest    ; clear exponent field
          bisw2     #impbit,6+'dest     ; insert leading bit of fraction
          .endm
;
          .psect    $code,long,pic,shr,nowrt
          .entry    x03aaz,^m<r2,r3,r4,r5,r6,r7,r8,r9,r10,r11>
;
; x03aaz accumulates an inner-product in quadruple precision, adding
;  it to a quadruple precision scalar value.
; it follows the pattern of the following code in ibm extended fortran
;
;     subroutine x03aaz(a, b, n, inca, incb, c1, c2, d1, d2)
;     double precision a(*), b(*), c1, c2, d1, d2
;     integer n, inca, incb
;     real*16 sum, term
;     double precision aa, bb
;     integer ia, ib, nn
;     nn = n
;     ia = 1
;     ib = 1
;     sum = qextd(c1)
;     if (c2.eq.0.0) go to 40
;     term = qextd(c2)
;  20 sum = sum + term
;  40 if (nn.le.0) go to 60
;     aa = a(ia)
;     bb = b(ib)
;     ia = ia + inca
;     ib = ib + incb
;     nn = nn - 1
;     if (aa.eq.0.0 .or. bb.eq.0.0) go to 40
;     term = qextd(aa)*qextd(bb)
;     go to 20
;  60 d1 = sum + sum - qextd(dbleq(sum))
;     d2 = sum - qextd(d1)
;     return
;     end
;
; here qextd converts a double precision number to quadruple precision
;  and dbleq truncates a quadruple precision number to double precision.
; the intent of the statement labelled 60 is simply to round sum to
;  double precision.
;
; fetch input arguments
;
          movl      4(ap),r0            ; address of a(1) to r0
          movl      8(ap),r1            ; address of b(1) to r1
          movl      @12(ap),r11         ; value of n to r11
          ashl      #3,@16(ap),r2       ; values of inca and incb to r2
          ashl      #3,@20(ap),r3       ;  and r3 as byte-displacements
          movl      24(ap),r4           ; address of c1 to r4
          movl      28(ap),r5           ; address of c2 to r5
;
; set up pointers to stack
; summary of stack usage
;    sp    points to 4-word workspace
;    sp+8  points to 4-word unpacked fraction of aa
;    sp+16 points to 4-word unpacked fraction of bb
; r7=sp+24 points to 8-word unpacked fraction of term
; r9=sp+40 points to 8-word unpacked fraction of sum
; n.b. r7 and r9 may be interchanged during the quadruple precision
;  addition. they must not be altered in the rest of the code
;
          subl2     #56,sp
          addl3     #24,sp,r7
          addl3     #40,sp,r9
;
; unpack c1 to give initial value of quadruple precision sum.  sum is
;  stored with unbiassed exponent in r8 and with r9 holding address of
;  sign-magnitude 8-word fraction.
;
          extzv     #exppos,#explen,(r4),r8       ; exponent of c1 to r8
          bneq      nzc1
          tstw      (r4)
          bgeq      rm1
          brw       error               ; error if reserved operand
rm1:      clrq      8(r9)               ; clear sum if c1 is zero
          clrq      (r9)
          brb       testc2
nzc1:      unpfra    0(r4),8(r9)          ; unpack fraction of c1
          clrq      (r9)                ; set lower half to zero
          subl2     #expbias,r8         ; remove bias from exponent
;
; fetch c2 and ignore if zero.  otherwise unpack c2 to form quadruple
;  precision value of term,  term is stored with unbiassed exponent in
;  r6 and with r7 holding address of sign-magnitude 8-word fraction.
;
testc2:    extzv     #exppos,#explen,(r5),r6       ; exponent of c2 to r
          bneq      nzc2
          tstw      (r5)
          bgeq      rm2
          brw       error               ; error if reserved operand
rm2:      brw       testn               ; skip if c2 zero
nzc2:      unpfra    0(r5),8(r7)          ; unpack fraction of c2
          clrq      (r7)                ; set lower half to zero
          subl2     #expbias,r6         ; remove bias from exponent
;
; the next section of code performs the addition
;     sum = sum + term
; where sum and term are unpacked quadruple precision numbers
;  and term is known to be non-zero
; the unbiassed exponent of sum is in r8
; the address of the sign-magnitude fraction of sum is in r9
; the unbiassed exponent of term is in r6
; the address of the sign-magnitude fraction of term is in r7
;
add:       tstl      12(r9)              ; test if sum is zero
          bneq      nzsum
          movl      r6,r8               ; if sum is zero then simply
          movq      (r7),(r9)           ;  copy term to sum
          movq      8(r7),8(r9)
          brw       testn
nzsum:     subl3     r6,r8,r6               ; exponent of result is the
          bgeq      difexp              ;  larger of the exponents of
          mnegl     r6,r6               ;  sum and term
          addl2     r6,r8
          movl      r7,r10              ; interchange fractions of sum
          movl      r9,r7               ;  and term to ensure that sum
          movl      r10,r9              ;  has same exponent as result
difexp:    movb      15(r7),r4           ; save sign of term in r4
          movb      15(r9),r5           ; save sign of sum in r5
          clrb      15(r7)              ; clear sign-bits of fractions
          clrb      15(r9)
          tstl      r6
          beql      addfra              ; no shifting if exponents equal
          cmpl      r6,#<fralen+64>
          blss      rm3
          brw       addfin              ; skip if term is insignificant
rm3:      cmpl      r6,#64              ; shift fraction of term to the
          blss      shif1               ;  right using moves to reduce
          movq      8(r7),(r7)          ;  shift count below 32
          clrq      8(r7)
          subl2     #64,r6
shif1:     cmpl      r6,#32
          blss      shif2
          movq      4(r7),(r7)
          movl      12(r7),8(r7)
          subl2     #32,r6
          clrl      12(r7)
shif2:     mnegl     r6,r10              ; now use shift instructions.
          ashq      r10,(r7),(r7)       ;  positive shift count in r6
          ashl      r6,4(r7),4(r7)      ;  negative shift count in r10
          ashq      r10,4(r7),4(r7)
          ashl      r6,8(r7),8(r7)
          ashq      r10,8(r7),8(r7)
addfra:    xorb2     r5,r4               ; compare signs of sum and term
          tstb      r4
          blss      subfra              ; branch if different signs
          addl2     (r7),(r9)           ; add fractions if same signs
          adwc      4(r7),4(r9)
          adwc      8(r7),8(r9)
          adwc      12(r7),12(r9)
          bbs       #<exppos+1>,14(r9),rm4
          brw       addfin
rm4:      ashoctr   1,0(r9)              ; renormalize if fraction
          incl      r8                  ;  overflow
          brw       addfin
subfra:    cmpl      12(r7),12(r9)       ; compare fractions
          blssu      subterm
          bgtru      subsum
          cmpl      8(r7),8(r9)
          blssu      subterm
          bgtru      subsum
          cmpl      4(r7),4(r9)
          blssu      subterm
          bgtru      subsum
          cmpl      (r7),(r9)
          blssu      subterm
          bgtru      subsum
          clrq      (r9)                ; fractions are equal so result
          clrq      8(r9)               ;  is zero
          brb       testn
subsum:    movl      r7,r10              ; interchange fractions so that
          movl      r9,r7               ;  r9 holds address of larger
          movl      r10,r9              ;  fraction
          xorb2     #128,r5               ; also interchange sign
subterm:   subl2     (r7),(r9)
          sbwc      4(r7),4(r9)
          sbwc      8(r7),8(r9)
          sbwc      12(r7),12(r9)
normloop:  bbs       #exppos,14(r9),addfin         ; normalize
          ashoctl   1,0(r9)
          decl      r8
          brb       normloop
addfin:    movb      r5,15(r9)           ; restore sign bit
;
; test if n is negative or zero
;
testn:     tstl      r11
          bgtr      rm6
          brw       packup
;
; fetch exponents of a(ia) (= aa) and b(ib) (= bb)
;
rm6:      extzv     #exppos,#explen,(r0),r4       ; exponent of aa to r4
          bneq      expb
          tstw      (r0)
          bgeq      expb
          brw       error               ; error if reserved operand
expb:      extzv     #exppos,#explen,(r1),r5       ; exponent of bb to r
          bneq      orexp
          tstw      (r1)
          bgeq      orexp
          brw       error               ; error if reserved operand
orexp:    mull3     r4,r5,r10           ; r10=0 if aa*bb=0
;
; store unbiassed exponent of aa*bb in r6 and sign in r10.
; unpack fractions of aa and bb.  store addresses of unsigned fractions
;  in r4 and r5
;
          subl2     #expbias,r4         ; remove bias from exponents
          subl2     #expbias,r5
          addl3     r4,r5,r6            ; exponent of aa*bb to r6
          addl3     #8,sp,r4            ; now set r4 and r5 to point to
          addl3     #16,sp,r5           ;  stack
          unpfra    0(r0),0(r4)           ; unpack fractions of aa and b
          unpfra    0(r1),0(r5)
;
; increment addresses of a(ia) and b(ib) in r0 and r1.  decrement n.
; skip multiplication if aa or bb is zero
;
incadd:    addl2     r2,r0               ; add inca to address of a(ia)
          addl2     r3,r1               ; add incb to address of b(ib)
          decl      r11                 ; decrement n
          tstl      r10
          beql      testn               ; skip if aa or bb is zero
          xorb3     7(r4),7(r5),r10     ; sign of aa*bb to r10
          clrb      7(r4)               ; clear signs of fractions
          clrb      7(r5)
;
; form 8-word sign-magnitude fraction of aa*bb (address in r7)
;
          emul      (r4),(r5),#0,(r7)
          emul      4(r4),(r5),4(r7),4(r7)
          emul      (r4),4(r5),#0,(sp)  ; use 4 words of stack for
          addl2     (sp),4(r7)           ;  workspace
          adwc      4(sp),8(r7)
          emul      4(r4),4(r5),8(r7),8(r7)
          tstl      (r4)                ; test if lower half of fraction
          bgeq      mult1               ;  of aa has sign bit set
          addl2     (r5),4(r7)          ; compensate for negative sign
          adwc      4(r5),8(r7)
          adwc      #0,12(r7)
mult1:     tstl      (r5)                ; test if lower half of fractio
          bgeq      mult2               ;  of bb has sign bit set
          addl2     (r4),4(r7)          ; compensate for negative sign
          adwc      4(r4),8(r7)
          adwc      #0,12(r7)
          tstl      (r4)                ; test if both lower halves have
          bgeq      mult2               ;  sign bits set
          decl      8(r7)               ; compensate for 2 negative
          sbwc      #0,12(r7)           ;  signs
mult2:    ashoctl   explen,0(r7)         ; align fraction
          bbs       #exppos,14(r7),mulfin         ; normalize
          ashoctl   1,0(r7)
          decl      r6
mulfin:    movb      r10,15(r7)          ; insert sign of product
          brw       add
;
; round sum to double precision and pack as double precision number in
;  d1.  set d2 to correct value of sum - d1.
;
packup:   movl      32(ap),r4           ; address of d1 to r4
          movl      36(ap),r5           ; address of d2 to r5
          tstl      12(r9)                ; test for zero sum
          bneq      rm5
          brw       zerod1
rm5:      movb      15(r9),r2           ; sign of d1 in r2
          movb      r2,r3               ; sign of d2 in r3 (may change)
          clrb      15(r9)              ; clear sign of sum
          addl2     #expbias,r8         ; biassed exponent of d1 to r8
          subl3     #fralen,r8,r6       ; biassed exponent of d2 to r6
          ashq      #-explen,(r9),(r9)
          tstb      7(r9)               ; test most significant bit of
          bgeq      packd1              ;  lower half of fraction
          incl      8(r9)               ; add 1 to upper half of
          adwc      #0,12(r9)           ;  fraction to round it
          bbc       #<exppos+1>,14(r9),neglow
          ashq      #-1,8(r9),8(r9)     ; renormalize if necessary
          incl      r8
neglow:    mnegl     4(r9),4(r9)         ; negate lower half of fraction
          mnegl     (r9),(r9)
          sbwc      #0,4(r9)
          xorb2     #128,r3               ; change sign of d2
packd1:    tstl      r8                  ; test for underflow of d1
          bleq      zerod1
          cmpl      r8,#explim          ; test for overflow of d1
          bgeq      error2
          movb      r2,15(r9)           ; insert sign of d1
          insv      r8,#exppos,#explen,14(r9)     ; insert exponent
          rotl      #16,12(r9),(r4)     ; store d1
          rotl      #16,8(r9),4(r4)
          tstl      4(r9)               ; test for zero d2
          bneq      normd2
          tstl      (r9)
          beql      zerod2
normd2:    bbs       #exppos,6(r9),packd2          ; normalize d2
          ashq      #1,(r9),(r9)
          decl      r6
          brb       normd2
packd2:    tstl      r6                  ; test for underflow of d2
          bleq      zerod2
          movb      r3,7(r9)            ; insert sign of d2
          insv      r6,#exppos,#explen,6(r9)      ; insert exponent
          rotl      #16,4(r9),(r5)      ; store d2
          rotl      #16,(r9),4(r5)
          brb       return
zerod1:    clrq      (r4)                ; set d1 to zero
zerod2:    clrq      (r5)                ; set d2 to zero
return:   ret                           ; return
;
; if floating-point error has occurred, set d1 to reserved operand
;
error:     movl      32(ap),r4
          movl      36(ap),r5
error2:    clrq      (r4)
          bisb2      #128,7(r4)
          brb       zerod2
          .end
_ENDIF
_IF(ipsc)
//	asm version of ishft for ipsc-i860
	.text
	.globl		_ishft_
	.align		16
_ishft_:
	ld.l		0(r17), r17	// Load shift value
	ld.l		0(r16), r16	// Load word to be shifted
	subs		0, r17, r18	// Check sign
	bnc		right		// jump to right shift
//
	shl		r17, r16, r16	// do left shift
	bri		r1		// return
	  nop
//
right:	shr		r18, r16, r16	// do right shift
	bri		r1
	  nop
_ENDIF
