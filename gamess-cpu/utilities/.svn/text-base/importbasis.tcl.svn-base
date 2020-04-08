#!/usr/bin/tcl
#
# Import Basis
# ============
#
# Utility to import basis sets from the EMSL basis set library into GAMESS-UK.
#
# This utility takes 2 file names as arguments. The first file should contain
# a basis set as downloaded from the EMSL basis set library in the GAMESS-UK
# input format. This file is read and processed to produce a file with the 
# second file name which contains the Fortran to implement the basis set in
# the GAMESS-UK internal format.
#
# Limitation: This script cannot handle SP-shells (yet).
#
# Huub van Dam, November 2003

global atoms angmom exp coef nexp oldatom curmom

set oldatom "stom"

set inpfile [open [lindex $argv 0] r]

while {[gets $inpfile line] >= 0} {
   set line [string trim "$line"]
   set line [string tolower "$line"]
   set char [string index "$line" 0]
   if { [string length "$line" ] == 0 } {
   } elseif { [string compare  $char # ] == 0 } {
       if { [string first "basis=" "$line" ] != -1 } {
          set begin [string first \" "$line"]
          incr begin
          set end [string last \" "$line"]
          incr end -1
          set basis [string range "$line" $begin $end]
       }
   } elseif { [string compare $char s ] == 0 || \
              [string compare $char p ] == 0 || \
              [string compare $char d ] == 0 || \
              [string compare $char f ] == 0 || \
              [string compare $char g ] == 0 } {
       set newatom [lindex "$line" 1]
       set curmom  [lindex "$line" 0]
       if { [string compare $newatom $oldatom ] != 0 } {
          set oldatom $newatom
          set nexp($newatom) 0
          lappend atoms $newatom
       }
   } else {
       incr nexp($newatom)
       set coef($newatom,$nexp($newatom)) [lindex "$line" 0]
       set exp($newatom,$nexp($newatom))  [lindex "$line" 1]
       set angmom($newatom,$nexp($newatom))  $curmom
   }
}

close $inpfile



set outfile [open [lindex $argv 1] w]

foreach atom $atoms {

   puts $outfile "c +++ basis ${atom}_$basis"
   for {set i 1} {$i <= $nexp($atom)} {incr i} {
      puts $outfile [format "      e(%4d) =%20.8fd0"  $i $exp($atom,$i)]
      puts $outfile [format "      %s(%4d) =%20.8fd0" $angmom($atom,$i) $i $coef($atom,$i)]
   }
}

close $outfile
