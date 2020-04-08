This is an example of a file containing statements
which appear conditionally when the file is processed with
m4. The file can be processed by passing a machine definition
file (eg sgi.m) the macros (gener.m) and this file through m4.
eg.
     cat sgi.m gener.m example.m | m4 > outfile

1) pre-defined macros from gener.m

  The creation date and time are DATE
  The machine flag is GEN_MACHINE
  The options selected are GEN_OPTIONS
  
1) pre-defined macros from machine file

  REAL for the 64-bit real number representation
  COMPLEX for the 128-bit complex number representation
  INCLUDE(file)  for the include syntax
  

2) Block structured conditional statements. These are controlled by an options
   list (see machine.m4) and effectivley "or" the options together. They
   may be nested. 

_IF(sgi,convex,apollo)
this line is included if any of the options sgi,convex,apollo are selected
_IF(sgi)
this line is contional on sgi
_ELSE
this line is conditional on convex,apollo and not sgi
_ENDIF
_ELSEIF(rs6000)
this line is included if rs6000 is set, without sgi,convex,apollo
_ELSE
this line is included if none of rs6000,sgi,convex,apollo are set
_ENDIF

There is no negation operator, but complements to the if and elseif structures exist
_IFN(convex)
_This line is conditional on option cc not being set
_ELSEIFN(hp700)
This line is conditional on convex being set, but not hp700
_ELSE
This line is included if both convex and hp700 are not set
_ENDIF

3) One-line conditionals. These are controlled by single letter machine
  keys. Only one key may be set at once (see machine.m4). There is no nesting,
  and no else-if construct. 
c
_IF1(abc)This line is included if the machine key is a, b, or c
_ELS1()This line appears when the preceding line does not (ie the machine key is none of a, b, or c)
c
_IFN1(a)This line appears when the machine key is not a
_ELS1()This line apears when the machine key is a
c
