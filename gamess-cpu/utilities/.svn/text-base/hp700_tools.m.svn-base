dnl this is the machine dependent file for hp700 with tcgmsg
define(GEN_OPTIONS,`hp700,cio,parallel,tools')dnl
define(GEN_MACHINE,`b')dnl
define(REAL,real*8)dnl
define(COMPLEX,complex*16)dnl
changequote(!,!)dnl
define(INCLUDE,!      include '$1'!)dnl
changequote()dnl
define(setsto,vec_$iinit($3,$1,$2))dnl
define(dgthr,vec_$dgather($2,$4,$1,$3))dnl
define(dsctr,vec_$dscatter($2,$3,$1,$4))dnl
define(addvec,vec_$dadd_vector($2,$3,$4,$1))dnl
define(vclr,`ifelse(TEST_ARGS($2,1),vec_$dinit($1,$3,0.0d0),``vclr($1,$2,$3)'')')dnl
define(dcopy,`ifelse(TEST_ARGS($3,1,$5,1),``vec_$dcopy($2,$4,$1)'',``dcopy($1,$2,$3,$4,$5)'')')dnl
define(daxpy,`ifelse(TEST_ARGS($4,1,$6,1),``vec_$dmult_add($5,$3,$1,$2,$5)'',``daxpy($1,$2,$3,$4,$5,$6)'')')dnl
define(ddot,`ifelse(TEST_ARGS($3,1,$5,1),``vec_$ddot($2,$4,$1)'',``ddot($1,$2,$3,$4,$5)'')')dnl
define(dscal,`ifelse(TEST_ARGS($4,1),``vec_$dmult_constant($3,$1,$2,$3)'',``dscal($1,$2,$3,$4)'')')dnl
define(dsum,`ifelse(TEST_ARGS($3,1),``vec_$dsum($2,$1)'',``dsum($1,$2,$3)'')')dnl
define(dfill,`ifelse(TEST_ARGS($4,1),``vec_$dinit($3,$1,$2)'',``dfill($1,$2,$3,$4)'')')dnl
define(SEND,snd(`$1',`$2',`$3',`$4',`$5'))dnl
define(RECV,rcv(`$1',`$2',`$3',`$4',`$5',`$6',`$7'))dnl
