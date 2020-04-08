dnl this is the machine dependent file for sun with tcgmsg tools
define(GEN_OPTIONS,`sun,cio,parallel,tools')dnl
define(GEN_MACHINE,`s')dnl
define(REAL,real*8)dnl
define(COMPLEX,complex*16)dnl
changequote(!,!)dnl
define(INCLUDE,!      include '$1'!)dnl
changequote()dnl
