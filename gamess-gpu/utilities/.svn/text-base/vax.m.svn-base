dnl this is the machine dependent file for the vax-alpha
define(GEN_OPTIONS,`vax,fortio,alpha')dnl
define(GEN_MACHINE,`v')dnl
define(REAL,real*8)dnl
define(COMPLEX,complex*16)dnl
changequote(!,!)dnl
define(INCLUDE,!      include '$1.'!)dnl
changequote()dnl
dnl ----- assembler ------
define(vclr,szero($1,$3))dnl
define(vmul,`ifelse(TEST_ARGS($2,1,$4,1,$6,1),vvtv($7,$5,$1,$3),vmulf($1,$2,$3,$4,$5,$6,$7))')dnl
define(vadd,`ifelse(TEST_ARGS($2,1,$4,1,$6,1),addvec($5,$1,$3,$7),vaddf($1,$2,$3,$4,$5,$6,$7))')dnl
define(vsma,`ifelse(TEST_ARGS($2,1,$5,1,$7,1),gtriad($8,$3,$6,$4,$1),``vsma($1,$2,$3,$4,$5,$6,$7,$8)'')')dnl
define(vsub,`ifelse(TEST_ARGS($2,1,$4,1,$6,1),subvec($5,$1,$3,$7),``vsub($1,$2,$3,$4,$5,$6,$7)'')')dnl
define(vfill,dfill($4,$1,$2,$3))dnl
