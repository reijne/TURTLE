#!/bin/bash

#   chapter 12 examples using rungamess script
#
export PATH="../../rungamess:../../utilities:${PATH}"

(\cd ../../utilities; make validate)
root=`\cd ../..;pwd`
here=`pwd`
export GAMESS_EXE=$root/bin/gamess
export GAMESS_SCR=$here
export GAMESS_TMP=$here
export GAMESS_LIB=$root/libs
#
# If local disk space is short you may want to 
# edit the next two lines
#
#setenv GAMESS_SCR /scr1/psh
#setenv GAMESS_TMP /scr1/psh
#
time rungamess -t ed2 -t ed3 -n water water_scf
time rungamess -t ed2 -t ed3 -n water water_opt
time rungamess -t ed2 -t ed3 -n water water_restopt
time rungamess -t ed2 -t ed3 -n water water_orhf
time rungamess -t ed2 -t ed3 -n water water_loc
time rungamess -t ed2 -t ed3 -n water water_gvb
(\cd $GAMESS_TMP; /bin/rm -rf water.ed2 water.ed3)
#
time rungamess -t ed3 -t ed2 -n hcl hcl_scf
time rungamess -t ed3 -t ed2 -n hcl hclplus
(\cd $GAMESS_TMP; /bin/rm -rf hcl.ed2 hcl.ed3)
#
time rungamess cubane
time rungamess mg10
time rungamess nitrobenzene
time rungamess tnt
#
time rungamess -t ed2 -t ed3 -n na7mg na7mg_rhf
time rungamess -t ed2 -t ed3 -n na7mg na7mg_orhf
time rungamess -t ed2 -t ed3 -n na7mg na7mg_uhf
#
time rungamess -t ed2 -t ed3 -n na7mg na7mg_ext
time rungamess -t ed2 -t ed3 -n na7mg na7mg_ext_orhf
#
time rungamess -t ed2 -t ed3 -n na7mg na7mg_ecp
time rungamess -t ed2 -t ed3 -n na7mg na7mg_ecp_orhf
(\cd $GAMESS_TMP; /bin/rm -rf na7mg.ed2 na7mg.ed3)
#
time rungamess -t ed3=nico4.ed3 nico4_scf
time rungamess -t ed3=nico4.ed3 nico4_grid
(\cd $GAMESS_TMP; /bin/rm -rf nico4.ed3 nico4_grid.pun)
#
time rungamess -t ed2 -t ed3 -n imino imino_rhf
time rungamess -t ed2 -t ed3 -n imino imino_gvb
(\cd $GAMESS_TMP; /bin/rm -rf imino.ed2 imino.ed3)
#
time rungamess direct_scf
time rungamess berylocene_opt
#
time rungamess hcn_tr
time rungamess hcn_st
time rungamess hcn_js
#
time rungamess -t ed2 -t ed3 -n hsip hsip_scf
time rungamess -t ed2 -t ed3 -n hsip hsip_ts
time rungamess -t ed2 -t ed3 -n hsip hsip_fc
#
time rungamess -t ed2 -t ed3 -n hsip hsip_fcm1
time rungamess -t ed2 -t ed3 -n hsip hsip_tsfcm
time rungamess -t ed2 -t ed3 -n hsip hsip_fcm2
(\cd $GAMESS_TMP; /bin/rm -rf hsip.ed2 hsip.ed3)
#
time rungamess hcn_bf
#
time rungamess -t ed3=ethene.ed3 ethene_opt
time rungamess -t ed3=ethene.ed3 ethene_fcm
#
time rungamess -t ed3=ethene.ed3 ethene_mp2opt
time rungamess -t ed3=ethene.ed3 ethene_mp2fcm
#
time rungamess -t ed3=ethene.ed3 ethene_mp2opt
time rungamess -t ed3=ethene.ed3 ethene_mp2pol
(\cd $GAMESS_TMP; /bin/rm -rf ethene.ed3)
#
time rungamess -t ed3 pyridine
(\cd $GAMESS_TMP; /bin/rm -rf pyridine.ed3)
#
time rungamess -t ed3 -n water water_scf
time rungamess -t ed3 -n water water_cas
#
time rungamess -r casscf -n water water_cas
time rungamess -r casscf -n water water_cas_rest
time rungamess -r casscf -n water water_cas_opt
(\cd $GAMESS_TMP; /bin/rm -rf water.ed1 water.ed10 water.ed11 water.ed2 water.ed3 water.ed4 water.ed6 water.ed9)
#
time rungamess -t ed2 -t ed3 -n beo beo_rhf
time rungamess -t ed2 -t ed3 -n beo beo_casscf_ci
#
time rungamess -t ed2 -t ed3 -n beo beo_rhf
time rungamess -t ed2 -t ed3 -n beo beo_mcscf
#
time rungamess -r mcscf -n beo beo_mcscf
time rungamess -r mcscf -n beo beo_mcscf_rest
time rungamess -t ed3 -n beo beo_mcscf_ci
(\cd $GAMESS_TMP; /bin/rm -rf beo.ed13 beo.ed2 beo.ed3 beo.ed4 beo.ed6)
#
time rungamess -t ed2 -t ed3 -n nh3 nh3_rohf
time rungamess -t ed2 -t ed3 -l table -n nh3 nh3_mrdci
#
time rungamess -r mrdci -l table -n nh3 nh3_sa_tran
time rungamess -r mrdci -l table -n nh3 nh3_select
time rungamess -r mrdci -l table -n nh3 nh3_hamil
time rungamess -r mrdci -l table -n nh3 nh3_diag
time rungamess -r mrdci -l table -n nh3 nh3_analy
(\cd $GAMESS_TMP; /bin/rm -rf nh3.ed2 nh3.ed3 nh3.ftn031 nh3.ftn033 nh3.ftn034 nh3.ftn035 nh3.ftn036)
#
time rungamess -t ed2 -t ed3 -l ed0 -n nicch2 nicch2_rhf
time rungamess -t ed2 -t ed3 -l ed0 -n nicch2 nicch2_swap
time rungamess -t ed2 -t ed3 -l ed0 -n nicch2 nicch2_orhf
time rungamess -t ed2 -t ed3 -t ed4 -t ed6 -l ed0 -n nicch2 nicch2_cas
time rungamess -t ed2 -t ed3 -t ed5 -t ed6 -l ed0 -n nicch2 nicch2_ci
(\cd $GAMESS_TMP; /bin/rm -rf nicch2.ed2 nicch2.ed3 nicch2.ed4 nicch2.ed5 nicch2.ed6)
#
time rungamess -t ed2 -t ed3 -n pyridine2 pyridine2_rhf
time rungamess -r mrdci -l table -n pyridine2 pyridine2_1m1r
time rungamess -r mrdci -l table -n pyridine2 pyridine2_6m1r
time rungamess -r mrdci -l table -n pyridine2 pyridine2_21m10r
time rungamess -r mrdci -l table -n pyridine2 pyridine2_19m10r
(\cd $GAMESS_TMP; /bin/rm -rf pyridine2.ed3 pyridine2.ed2 pyridine2.ftn031 pyridine2.ftn033 pyridine2.ftn034 pyridine2.ftn035 pyridine2.ftn036)
#
time rungamess -r fullci -n h2o_fullci fullci_all 
time rungamess -r fullci -n h2o_fullci fullci_val
#
# cleanup
#
(\cd $GAMESS_TMP; /bin/rm -rf *.ed* *.ftn*)
#
/bin/rm -rf *.pun
/bin/rm -rf beo berylocene_opt cubane direct_scf ethene_fcm ethene_mp2fcm
/bin/rm -rf ethene_mp2opt ethene_mp2pol ethene_opt h2o_fullci hcl hcn_bf
/bin/rm -rf hcn_js hcn_st hcn_tr hsip imino mg10 na7mg nh3 nicch2 nico4_grid
/bin/rm -rf nico4_scf nitrobenzene pyridine pyridine2 tnt water
