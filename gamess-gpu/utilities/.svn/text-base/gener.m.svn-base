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
define(`G_IFDEF',defn(`ifdef'))dnl
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
undefine(`ifdef')dnl
undefine(`format')dnl
undefine(`undefine')dnl
dnl
dnl    
dnl  *** generation date stamp
dnl
dnl Suppressed because of problems on Linux
dnl
dnl G_DEFINE(`G_DATEFILE',maketemp(/tmp/m4.XXXXX))dnl
dnl changequote(!,!)dnl
dnl syscmd(echo `date` !dnl! > G_DATEFILE )dnl
dnl changequote(`,')dnl
dnl G_DEFINE(`DATE', `G_INCLUDE(G_DATEFILE)')dnl
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
dnl
G_DEFINE(G_NO,0)G_DEFINE(G_NOTYET,1)G_DEFINE(G_COPY,2)G_DEFINE(G_YES,3)dnl
G_DEFINE(G_LEVEL,G_COPY)dnl
G_DEFINE(G_INIF,G_NO)dnl
G_DEFINE(G_INEX,G_NO)dnl
G_DEFINE(G_REDEF,`popdef(`$1')pushdef(`$1',`$2')')dnl
dnl
dnl Conditional diversion - normally use divert, but not when in extract mode
dnl
G_DEFINE(G_DIVERT2,`ifelse(GEN_EXTRACT,0,`G_DIVERT(`$1')',`ifelse(G_DIVERTFLAG,G_YES,`G_DIVERT(`$1')')')')dnl
dnl
dnl _AND evaluates to true if all its arguments are satisfied
dnl
G_DEFINE(_AND,`ifelse(G_CHKLIST($@),G_CHKLIST2($@),`TRUE',`FALSE')')dnl
dnl
G_DEFINE(_NOT,`ifelse(G_CHECKSTR($@),-1,`TRUE',`FALSE')')dnl
dnl
G_DEFINE(_IF,     `pushdef(`G_INIF',G_YES)ifelse(G_LEVEL,G_NO,    `pushdef(`G_LEVEL',G_NO)dnl',
                       G_LEVEL,G_NOTYET,`pushdef(`G_LEVEL',G_NO)dnl',
                       G_LEVEL,G_COPY,  `ifelse(G_CHECKSTR($@), -1,
                           `pushdef(`G_LEVEL',G_NOTYET)G_DIVERT2(-1)dnl',
                           `pushdef(`G_LEVEL',G_COPY)dnl')')')dnl
G_DEFINE(_ELSEIF, `ifelse(G_INIF,G_NO, `_ERROR(``_ELSEIF'' without ``_IF'')',
                          G_INIF,G_YES, `ifelse(G_LEVEL,G_NO,`dnl',
                       G_LEVEL,G_NOTYET,`ifelse(G_CHECKSTR($@),-1,
                             `dnl',`G_REDEF(`G_LEVEL',G_COPY)G_DIVERT2(0)dnl')',
                       G_LEVEL,G_COPY,  `G_REDEF(`G_LEVEL',G_NO)G_DIVERT2(-1)dnl')')')dnl
G_DEFINE(_IFN,    `pushdef(`G_INIF',G_YES)ifelse(G_LEVEL,G_NO,    `pushdef(`G_LEVEL',G_NO)dnl',
                       G_LEVEL,G_NOTYET,`pushdef(`G_LEVEL',G_NO)dnl',
                       G_LEVEL,G_COPY,  `ifelse(G_CHECKSTR($@), -1,
                           `pushdef(`G_LEVEL',G_COPY)dnl',
                           `pushdef(`G_LEVEL',G_NOTYET)G_DIVERT2(-1)dnl')')')dnl
G_DEFINE(_ELSEIFN, `ifelse(G_INIF,G_NO, `_ERROR(``_ELSEIFN'' without ``_IF'')',
                          G_INIF,G_YES, `ifelse(G_LEVEL,G_NO,`dnl',
                         G_LEVEL,G_NOTYET,`ifelse(G_CHECKSTR($@),-1,
                             `G_REDEF(`G_LEVEL',G_COPY)G_DIVERT2(0)dnl',`dnl')',
                         G_LEVEL,G_COPY,  `G_REDEF(`G_LEVEL',G_NO)G_DIVERT2(-1)dnl')')')dnl
G_DEFINE(_ELSE,   `ifelse(G_INIF,G_NO, `_ERROR(``_ELSE'' without ``_IF'')',
                          G_INIF,G_YES, `ifelse(G_LEVEL,G_NO,`dnl', 
                       G_LEVEL,G_NOTYET,`G_REDEF(`G_LEVEL',G_COPY)G_DIVERT2(0)dnl', 
                       G_LEVEL,G_COPY,  `G_REDEF(`G_LEVEL',G_NO)G_DIVERT2(-1)dnl')')')dnl
G_DEFINE(_ENDIF,  `ifelse(G_INIF,G_NO,  `_ERROR(``_ENDIF'' without ``_IF'')',
                          G_INIF,G_YES, `popdef(`G_INIF')popdef(`G_LEVEL')ifelse(G_LEVEL,G_NO,     `dnl',
                                  G_LEVEL,G_NOTYET, `dnl',
	                          G_LEVEL,G_COPY,   `G_DIVERT2(0)dnl')')')dnl
dnl
dnl Extraction of specific files - acts like _IFN with a single flag, 
dnl unless GEN_EXTRACT is set, in which case the second arg is checked agains GEN_EXTRACTFILE
dnl and if it matches, the section is output. G_DIVERTFLAG should cause the normal 
dnl handling of IF() etc to be resumed while inside the extract block
dnl
G_DEFINE(_EXTRACT,    `pushdef(`G_INEX',G_YES)pushdef(`G_FILE',`$1')ifelse(GEN_EXTRACT,0, `ifelse(G_LEVEL,G_NO,`pushdef(`G_LEVEL',G_NO)dnl',
                                              G_LEVEL,G_NOTYET,`pushdef(`G_LEVEL',G_NO)dnl',
                                              G_LEVEL,G_COPY,  `ifelse(G_CHECKSTR(G_SHIFT($@)), -1,
                                                               `pushdef(`G_LEVEL',G_COPY)dnl',
                                                               `pushdef(`G_LEVEL',G_NOTYET)G_DIVERT2(-1)dnl')')',
                        GEN_EXTRACT,1,`ifelse(G_LEVEL,G_NO,    `pushdef(`G_LEVEL',G_NO)dnl',
                                              G_LEVEL,G_NOTYET,`pushdef(`G_LEVEL',G_NO)dnl',
                                              G_LEVEL,G_COPY,  `ifelse(G_CHECKSTR(G_SHIFT($@)), -1,
                                                                `pushdef(`G_LEVEL',G_NOTYET)',
                                                                `ifelse(G_FILE,GEN_EXTRACTFILE,
                                     `pushdef(`G_LEVEL',G_COPY)pushdef(`G_DIVERTFLAG',G_YES)G_DIVERT(0)dnl',`pushdef(`G_LEVEL',G_NO)')')')')')dnl
dnl
G_DEFINE(_ENDEXTRACT,  `ifelse(G_INEX,G_NO, `_ERROR(``_ENDEXTRACT'' without ``_EXTRACT'')',
                               G_INEX,G_YES, 
                        `popdef(`G_INEX')ifelse(G_LEVEL,G_COPY,
                         ifelse(GEN_EXTRACT,1,G_DIVERT(-1)popdef(`G_DIVERTFLAG')))popdef(`G_LEVEL')ifelse(G_LEVEL,G_NO,`dnl',
                                                                                                          G_LEVEL,G_NOTYET, `dnl',
                               	                                                                          G_LEVEL,G_COPY,`G_DIVERT2(0)dnl')')')dnl
dnl
G_DEFINE(_MACRO, `ifelse(G_LEVEL,G_NO,`dnl',
                        G_LEVEL,G_NOTYET, `dnl',
                        G_LEVEL,G_COPY,  `G_DEFINE(`$1',`$2')dnl')')dnl
G_DEFINE(_FAIL, `ifelse(G_LEVEL,G_NO,`dnl',
                        G_LEVEL,G_NOTYET, `dnl',
                        G_LEVEL,G_COPY, `_ERROR(`$1')')')dnl
dnl
G_DEFINE(_INCLUDE, `ifelse(G_LEVEL,G_COPY,`G_INCLUDE(`$1')dnl')')dnl
dnl
dnl
dnl exit due to processing failure
dnl put message on stderr, and on the back of the file 
G_DEFINE(_ERROR, `syscmd(sh -c "echo --- M4 processing failed at this point ---- $1; 1>&2 echo M4 processing failure: $1")m4exit(1)')dnl
G_DEFINE(_MSG, `syscmd(sh -c "echo $1")')dnl
dnl ** G_CHECKSTR checks the strings given as options to the IF() or ELSEIF()
dnl **            it returns -1 if no hits (G_CHKLIST evaluates to a string of 0s)
dnl
dnl ** G_CHKLIST is a string of characters, 1 for each hit 
dnl (used to be n*m (where n is number of options set in OPTIONS 
dnl                   and m the number of given to IF() or ELSEIF()
dnl but now 0s are left out)
dnl this strategy is flawed as a single field can generate > 1 entry if 
dnl an option is repeated
dnl **  G__CHKLIST2 is the string G_CHLKIST would generate if 
dnl   all options matched (used to implement logical and)
dnl
dnl ** G_CUR_OPT loops over the options in the source file
dnl
dnl  to allow appending of command-line options M4_OPTIONS is
dnl  appended to the GEN_OPTIONS macro.
dnl
dnl  the () are added to protect the embedded , from being interpreted
dnl  by ifelse                 v          v v            v
dnl
G_DEFINE(`GEN_OPTIONSX',ifelse((M4_OPTIONS),(`M4_OPTIONS'),``TRUE,GEN_OPTIONS'',
                                          ``TRUE,GEN_OPTIONS,M4_OPTIONS''))dnl
dnl The keyword TRUE is added so any appearance of TRUE as an option
dnl (e.g. via _AND) will match
dnl
dnl if G_CHK1 returns 11,111 etc it means an option was specified more than once
dnl to fix this G_CHK3 converts these to a single 1 
G_DEFINE(`G_CHK1',`G_DEFINE(`G_CUR_OPT',$1)G_CHKOPTS(GEN_OPTIONSX)')dnl
G_DEFINE(`G_CHK2',`ifelse($1,G_CUR_OPT,1,`')')dnl
G_DEFINE(`G_CHKOPTS', `ifelse($2, , `G_CHK2(`$1')', 
                                 `G_CHK2(`$1')G_CHKOPTS(G_SHIFT($@))')')dnl
G_DEFINE(`G_CHK3',`ifelse(G_CHK1(`$1'),1,1,G_CHK1(`$1'),11,1,G_CHK1(`$1'),111,1,G_CHK1(`$1'),1111,1,G_CHK1(`$1'),11111,1,`')')dnl
G_DEFINE(`G_CHKLIST', `ifelse($2, , `G_CHK3(`$1')', 
                                `G_CHK3(`$1')G_CHKLIST(G_SHIFT($@))')')dnl
dnl
G_DEFINE(`G_CHKLIST2', `ifelse($2, , 1, 1`G_CHKLIST2(G_SHIFT($@))')')dnl
dnl
dnl
dnl File extraction  - if GEN_EXTRACTFILE is set, set flag GEN_EXTRACT to 1 
dnl
G_IFDEF(`GEN_EXTRACTFILE',`G_DEFINE(GEN_EXTRACT,1)',`G_DEFINE(GEN_EXTRACT,0)')dnl
dnl
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
G_IFDEF(`GEN_OPTIONS',,`_FAIL(No options have been specified)')dnl
G_IFDEF(`GEN_MACHINE',,`_FAIL(No machine key has been specified)')dnl
G_IFDEF(`INCLUDE',,`_FAIL(No INCLUDE macro has been specified)')dnl
G_IFDEF(`REAL',,`_FAIL(No REAL macro has been specified)')dnl
G_IFDEF(`COMPLEX',,`_FAIL(No COMPLEX macro has been specified)')dnl
dnl
dnl  Extraction mode - default is to exclude code
dnl
ifelse(GEN_EXTRACT,1,`G_DIVERT(-1)')dnl
dnl
dnl this flag enables the DIVERT2 macro to detect if we are
dnl actively extracting a bit of a file
dnl
G_DEFINE(`G_DIVERTFLAG',G_NO)dnl
dnl
dnl
dnl no comments recognised
changecom()dnl
dnl m4wrap(`syscmd(rm G_DATEFILE)')dnl
