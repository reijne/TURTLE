dnl ******************************************************************* 
dnl * 
dnl *  ***  generate m4 macro  Paul Sherwood (SERC DL 1990)
dnl *
dnl *   implementation of list based block if constructs 
dnl *        
dnl *    This macro definition must be included in all files to be processed
dnl *    A file called machine.m4 must be provided to provide the following
dnl *   
dnl *    1.  definition of MACHINE  to a one-letter key for the machine
dnl *    2.  definition of OPTIONS to a set of option keywords
dnl *    3.  definitions of any function calls which must be substituted
dnl *    4.  definition of INCLUDE to the fortran include statement 
dnl *
dnl *   simple example - 
dnl *
dnl *   define(MACHINE,c)
dnl *   define(OPTIONS,`ibm,3090vf')
dnl *   define(INCLUDE,$insert vp.bigscf.apftn64(sizes))dnl
dnl *   define(szero,`vclr($1,$2,1,$3,$4)')
dnl *
dnl ******************************************************************* 
dnl
dnl
dnl  *** clear index and shift definitions
dnl
define(`G_INDEX',defn(`index'))dnl
define(`G_SHIFT',defn(`shift'))dnl
define(`G_DIVERT',defn(`divert'))dnl
define(`G_INCLUDE',defn(`include'))dnl
define(`G_DEFINE',defn(`define'))dnl
undefine(`divert')dnl
undefine(`index')dnl
undefine(`shift')dnl
undefine(`include')dnl
undefine(`undivert')dnl
undefine(`divnum')dnl
undefine(`incr')dnl
undefine(`eval')dnl
undefine(`substr')dnl
undefine(`defn')dnl
undefine(`decr')dnl
undefine(`define')dnl
undefine(`len')dnl
undefine(`dlen')dnl
undefine(`unix')dnl
undefine(`undefine')dnl
dnl    
dnl  *** generation date stamp
dnl
G_DEFINE(`G_DATEFILE',maketemp(/tmp/m4.XXXXX))dnl
changequote(!,!)dnl
syscmd(echo `date` !dnl! > G_DATEFILE )dnl
changequote()dnl
G_DEFINE(`DATE', `G_INCLUDE(G_DATEFILE)')dnl
dnl
dnl *** one line if and if not tests, 
dnl     argument is a concatenated string of key letters ***
dnl
G_DEFINE(G_CHECKS1,`G_INDEX($1,GEN_MACHINE)')dnl
G_DEFINE(_IF1,   `ifelse(G_CHECKS1($1),-1,`G_DEFINE(`G_ONELINE',0)dnl ',`G_DEFINE(`G_ONELINE',1)')')dnl
G_DEFINE(_IFN1,  `ifelse(G_CHECKS1($1),-1,`G_DEFINE(`G_ONELINE',1)',    `G_DEFINE(`G_ONELINE',0)dnl ')')dnl
G_DEFINE(_ELS1,  `ifelse(G_ONELINE,0,`',`dnl ')')dnl
dnl
dnl *** block if structure, argument(s) are a list of keywords, 
dnl     which are tested against the keywords stored as OPTIONS
dnl
G_DEFINE(G_CHECKSTR, `G_INDEX(G_CHKLIST($@),1)')dnl
G_DEFINE(G_NO,0)G_DEFINE(G_NOTYET,1)G_DEFINE(G_COPY,2)G_DEFINE(G_YES,3)dnl
G_DEFINE(G_LEVEL,G_COPY)dnl
G_DEFINE(G_INIF,G_NO)dnl
G_DEFINE(G_REDEF,`popdef(`$1')pushdef(`$1',`$2')')dnl
G_DEFINE(_IF,     `pushdef(`G_INIF',G_YES)ifelse(G_LEVEL,G_NO,    `pushdef(`G_LEVEL',G_NO)dnl',
                       G_LEVEL,G_NOTYET,`pushdef(`G_LEVEL',G_NO)dnl',
                       G_LEVEL,G_COPY,  `ifelse(G_CHECKSTR($@), -1,
                           `pushdef(`G_LEVEL',G_NOTYET)G_DIVERT(-1)dnl',
                           `pushdef(`G_LEVEL',G_COPY)dnl')')')dnl
G_DEFINE(_ELSEIF, `ifelse(G_INIF,G_NO, `_ERROR(``_ELSEIF'' without ``_IF'')',
                          G_INIF,G_YES, `ifelse(G_LEVEL,G_NO,`dnl',
                       G_LEVEL,G_NOTYET,`ifelse(G_CHECKSTR($@),-1,
                             `dnl',`G_REDEF(`G_LEVEL',G_COPY)G_DIVERT(0)dnl')',
                       G_LEVEL,G_COPY,  `G_REDEF(`G_LEVEL',G_NO)G_DIVERT(-1)dnl')')')dnl
G_DEFINE(_IFN,    `pushdef(`G_INIF',G_YES)ifelse(G_LEVEL,G_NO,    `pushdef(`G_LEVEL',G_NO)dnl',
                       G_LEVEL,G_NOTYET,`pushdef(`G_LEVEL',G_NO)dnl',
                       G_LEVEL,G_COPY,  `ifelse(G_CHECKSTR($@), -1,
                           `pushdef(`G_LEVEL',G_COPY)dnl',
                           `pushdef(`G_LEVEL',G_NOTYET)G_DIVERT(-1)dnl')')')dnl
G_DEFINE(_ELSEIFN, `ifelse(G_INIF,G_NO, `_ERROR(``_ELSEIFN'' without ``_IF'')',
                          G_INIF,G_YES, `ifelse(G_LEVEL,G_NO,`dnl',
                         G_LEVEL,G_NOTYET,`ifelse(G_CHECKSTR($@),-1,
                             `G_REDEF(`G_LEVEL',G_COPY)G_DIVERT(0)dnl',`dnl')',
                         G_LEVEL,G_COPY,  `G_REDEF(`G_LEVEL',G_NO)G_DIVERT(-1)dnl')')')dnl
G_DEFINE(_ELSE,   `ifelse(G_INIF,G_NO, `_ERROR(``_ELSE'' without ``_IF'')',
                          G_INIF,G_YES, `ifelse(G_LEVEL,G_NO,`dnl', 
                       G_LEVEL,G_NOTYET,`G_REDEF(`G_LEVEL',G_COPY)G_DIVERT(0)dnl', 
                       G_LEVEL,G_COPY,  `G_REDEF(`G_LEVEL',G_NO)G_DIVERT(-1)dnl')')')dnl
G_DEFINE(_ENDIF,  `ifelse(G_INIF,G_NO,  `_ERROR(``_ENDIF'' without ``_IF'')',
                          G_INIF,G_YES, `popdef(`G_INIF')popdef(`G_LEVEL')ifelse(G_LEVEL,G_NO,     `dnl',
                                        G_LEVEL,G_NOTYET, `dnl',
	                                G_LEVEL,G_COPY,   `G_DIVERT(0)dnl')')')dnl
G_DEFINE(_MACRO, `ifelse(G_LEVEL,G_NO,`dnl',
                        G_LEVEL,G_NOTYET, `dnl',
                        G_LEVEL,G_COPY,  `G_DEFINE(`$1',`$2')dnl')')dnl
G_DEFINE(_FAIL, `ifelse(G_LEVEL,G_NO,`dnl',
                        G_LEVEL,G_NOTYET, `dnl',
                        G_LEVEL,G_COPY, `_ERROR(`$1')')')dnl
dnl exit due to processing failure
dnl put message on stderr, and on the back of the file 
G_DEFINE(_ERROR, `syscmd(sh -c "echo --- M4 processing failed at this point ---- $1; 1>&2 echo M4 processing failure: $1")m4exit(1)')dnl
dnl ** G_CHECKSTR checks the strings given as options to the IF() or ELSEIF()
dnl **            it returns -1 if no hits (G_CHKLIST evaluates to a string of 0s)
dnl ** G_CHKLIST is a string of n * m characters, each of which is a 0 or 1
dnl               (where n is number of options set in OPTIONS 
dnl                   and m the number of given to IF() or ELSEIF()
dnl ** G_CUR_OPT loops over the options in the source file
G_DEFINE(`G_CHK1',`G_DEFINE(`G_CUR_OPT',$1)G_CHKOPTS(GEN_OPTIONS)')dnl
G_DEFINE(`G_CHK2',`ifelse($1,G_CUR_OPT,1,0)')dnl
G_DEFINE(`G_CHKLIST', `ifelse($2, , `G_CHK1(`$1')', 
                                `G_CHK1(`$1')G_CHKLIST(G_SHIFT($@))')')dnl
G_DEFINE(`G_CHKOPTS', `ifelse($2, , `G_CHK2(`$1')', 
                                 `G_CHK2(`$1')G_CHKOPTS(G_SHIFT($@))')')dnl
dnl
dnl   test macros for argument-conditional substitution - (test is made
dnl   at m4 time, not compile-time)
dnl
G_DEFINE(`G_ARG1',`ifelse(`$1',`$2',1,0)')dnl
G_DEFINE(`G_ARGLIST', `ifelse($3, , `G_ARG1(`$1',`$2')', 
                                `G_ARG1(`$1',`$2')G_ARGLIST(G_SHIFT(G_SHIFT($@)))')')dnl
G_DEFINE(G_CHECKARGSTR, `G_INDEX(G_ARGLIST($@),0)')dnl
dnl combine the -1 test value in the return from the user test macro
G_DEFINE(TEST_ARGS,`G_CHECKARGSTR($@),-1')dnl
ifdef(`GEN_OPTIONS',,`_FAIL(No options have been specified)')dnl
ifdef(`GEN_MACHINE',,`_FAIL(No machine key has been specified)')dnl
ifdef(`INCLUDE',,`_FAIL(No INCLUDE macro has been specified)')dnl
ifdef(`REAL',,`_FAIL(No REAL macro has been specified)')dnl
ifdef(`COMPLEX',,`_FAIL(No COMPLEX macro has been specified)')dnl
dnl no comments recognised
changecom()dnl
changequote(`,`)dnl
m4wrap(`syscmd(rm G_DATEFILE)')dnl
