V34 :0x2c dlf_store
14 dlf_util_d.F90 S624 0
04/08/2020  10:45:13
use dlf_allocate private
use dlf_global private
use dlf_parameter_module private
enduse
D 61 26 654 4840 653 7
D 214 26 940 200 937 7
D 220 20 297
D 222 23 10 1 306 305 0 1 0 0 1
 300 303 304 300 303 301
D 225 23 7 1 0 71 0 0 0 0 0
 0 71 0 11 71 0
D 228 23 7 1 0 307 0 0 0 0 0
 0 307 0 11 307 0
D 231 23 7 1 0 11 0 0 0 0 0
 0 11 0 11 11 0
D 234 23 7 1 0 307 0 0 0 0 0
 0 307 0 11 307 0
D 237 23 10 1 11 310 0 0 1 0 0
 0 309 11 11 310 310
D 240 23 10 1 11 312 0 0 1 0 0
 0 311 11 11 312 312
D 243 23 10 2 313 319 1 1 0 0 1
 11 314 11 11 314 315
 11 316 317 11 316 318
D 246 23 10 3 320 329 1 1 0 0 1
 11 321 11 11 321 322
 11 323 324 11 323 325
 11 326 327 11 326 328
D 249 23 10 2 330 336 1 1 0 0 1
 11 331 11 11 331 332
 11 333 334 11 333 335
D 252 23 10 3 337 346 1 1 0 0 1
 11 338 11 11 338 339
 11 340 341 11 340 342
 11 343 344 11 343 345
S 624 24 0 0 0 9 1 0 5013 10015 0 A 0 0 0 0 B 0 382 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 382 0 0 0 0 0 0 dlf_store
S 626 23 0 0 0 9 634 624 5044 14 0 A 0 0 0 0 B 0 388 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 rk
S 628 23 0 0 0 9 905 624 5058 14 0 A 0 0 0 0 B 0 389 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 stderr
S 629 23 0 0 0 9 904 624 5065 14 0 A 0 0 0 0 B 0 389 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 stdout
S 631 23 0 0 0 9 918 624 5085 14 0 A 0 0 0 0 B 0 390 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 allocate
S 632 23 0 0 0 9 920 624 5094 14 0 A 0 0 0 0 B 0 390 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 deallocate
S 633 3 0 0 0 7 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
R 634 16 1 dlf_parameter_module rk
S 636 3 0 0 0 19 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 19
S 637 3 0 0 0 19 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 19
S 639 3 0 0 0 7 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 11 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 640 3 0 0 0 7 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 12 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 643 3 0 0 0 7 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 15 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 645 3 0 0 0 7 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 650 3 0 0 0 7 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 16 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
R 653 25 2 dlf_global glob_type
R 654 5 3 dlf_global nvar glob_type
R 655 5 4 dlf_global nat glob_type
R 656 5 5 dlf_global tolerance glob_type
R 657 5 6 dlf_global tolerance_e glob_type
R 658 5 7 dlf_global energy glob_type
R 659 5 8 dlf_global oldenergy glob_type
R 660 5 9 dlf_global toldenergy glob_type
R 661 5 10 dlf_global tinit glob_type
R 662 5 11 dlf_global tatoms glob_type
R 663 5 12 dlf_global maxcycle glob_type
R 664 5 13 dlf_global maxene glob_type
R 665 5 14 dlf_global distort glob_type
R 666 5 15 dlf_global massweight glob_type
R 667 5 16 dlf_global maxdump glob_type
R 668 5 17 dlf_global toldenergy_conv glob_type
R 669 5 18 dlf_global oldenergy_conv glob_type
R 670 5 19 dlf_global task glob_type
R 671 5 20 dlf_global iopt glob_type
R 672 5 21 dlf_global iline glob_type
R 673 5 22 dlf_global imultistate glob_type
R 674 5 23 dlf_global needcoupling glob_type
R 675 5 24 dlf_global imicroiter glob_type
R 676 5 25 dlf_global maxstep glob_type
R 677 5 26 dlf_global scalestep glob_type
R 678 5 27 dlf_global taccepted glob_type
R 679 5 28 dlf_global lbfgs_mem glob_type
R 680 5 29 dlf_global temperature glob_type
R 681 5 30 dlf_global nzero glob_type
R 682 5 31 dlf_global qtsflag glob_type
R 683 5 32 dlf_global timestep glob_type
R 684 5 33 dlf_global fric0 glob_type
R 685 5 34 dlf_global fricfac glob_type
R 686 5 35 dlf_global fricp glob_type
R 687 5 36 dlf_global state_i glob_type
R 688 5 37 dlf_global state_j glob_type
R 689 5 38 dlf_global pf_c1 glob_type
R 690 5 39 dlf_global pf_c2 glob_type
R 691 5 40 dlf_global gp_c3 glob_type
R 692 5 41 dlf_global gp_c4 glob_type
R 693 5 42 dlf_global ln_t1 glob_type
R 694 5 43 dlf_global ln_t2 glob_type
R 695 5 44 dlf_global update glob_type
R 696 5 45 dlf_global maxupd glob_type
R 697 5 46 dlf_global delta glob_type
R 698 5 47 dlf_global soft glob_type
R 699 5 48 dlf_global inithessian glob_type
R 700 5 49 dlf_global carthessian glob_type
R 701 5 50 dlf_global tsrelative glob_type
R 702 5 51 dlf_global minstep glob_type
R 703 5 52 dlf_global icoord glob_type
R 704 5 53 dlf_global icoordinner glob_type
R 705 5 54 dlf_global nivar glob_type
R 706 5 55 dlf_global tcoords2 glob_type
R 707 5 56 dlf_global nimage glob_type
R 708 5 57 dlf_global nebk glob_type
R 709 5 58 dlf_global neb_climb_test glob_type
R 710 5 59 dlf_global neb_freeze_test glob_type
R 711 5 60 dlf_global maxrot glob_type
R 712 5 61 dlf_global tolrot glob_type
R 713 5 62 dlf_global ncons glob_type
R 714 5 63 dlf_global nconn glob_type
R 715 5 64 dlf_global havehessian glob_type
R 716 5 65 dlf_global dump glob_type
R 717 5 66 dlf_global restart glob_type
R 718 5 67 dlf_global cleanup glob_type
R 719 5 68 dlf_global po_pop_size glob_type
R 720 5 69 dlf_global po_radius_base glob_type
R 721 5 70 dlf_global po_contraction glob_type
R 722 5 71 dlf_global po_tolerance_g glob_type
R 723 5 72 dlf_global po_tol_r_base glob_type
R 724 5 73 dlf_global po_distribution glob_type
R 725 5 74 dlf_global po_maxcycle glob_type
R 726 5 75 dlf_global po_init_pop_size glob_type
R 727 5 76 dlf_global po_reset glob_type
R 728 5 77 dlf_global po_mutation_rate glob_type
R 729 5 78 dlf_global po_death_rate glob_type
R 730 5 79 dlf_global po_scalefac glob_type
R 731 5 80 dlf_global po_nsave glob_type
R 732 5 81 dlf_global nprocs glob_type
R 733 5 82 dlf_global ntasks glob_type
R 734 5 83 dlf_global nprocs_per_task glob_type
R 735 5 84 dlf_global iam glob_type
R 736 5 85 dlf_global iam_in_task glob_type
R 737 5 86 dlf_global mytask glob_type
R 738 5 87 dlf_global serial_cycle glob_type
R 739 5 88 dlf_global dotask glob_type
R 740 5 89 dlf_global xcoords glob_type
R 743 5 92 dlf_global xcoords$sd glob_type
R 744 5 93 dlf_global xcoords$p glob_type
R 745 5 94 dlf_global xcoords$o glob_type
R 747 5 96 dlf_global xcoords2 glob_type
R 751 5 100 dlf_global xcoords2$sd glob_type
R 752 5 101 dlf_global xcoords2$p glob_type
R 753 5 102 dlf_global xcoords2$o glob_type
R 755 5 104 dlf_global xgradient glob_type
R 758 5 107 dlf_global xgradient$sd glob_type
R 759 5 108 dlf_global xgradient$p glob_type
R 760 5 109 dlf_global xgradient$o glob_type
R 762 5 111 dlf_global znuc glob_type
R 764 5 113 dlf_global znuc$sd glob_type
R 765 5 114 dlf_global znuc$p glob_type
R 766 5 115 dlf_global znuc$o glob_type
R 768 5 117 dlf_global weight glob_type
R 770 5 119 dlf_global weight$sd glob_type
R 771 5 120 dlf_global weight$p glob_type
R 772 5 121 dlf_global weight$o glob_type
R 774 5 123 dlf_global iweight glob_type
R 776 5 125 dlf_global iweight$sd glob_type
R 777 5 126 dlf_global iweight$p glob_type
R 778 5 127 dlf_global iweight$o glob_type
R 780 5 129 dlf_global mass glob_type
R 782 5 131 dlf_global mass$sd glob_type
R 783 5 132 dlf_global mass$p glob_type
R 784 5 133 dlf_global mass$o glob_type
R 786 5 135 dlf_global icoords glob_type
R 788 5 137 dlf_global icoords$sd glob_type
R 789 5 138 dlf_global icoords$p glob_type
R 790 5 139 dlf_global icoords$o glob_type
R 792 5 141 dlf_global igradient glob_type
R 794 5 143 dlf_global igradient$sd glob_type
R 795 5 144 dlf_global igradient$p glob_type
R 796 5 145 dlf_global igradient$o glob_type
R 798 5 147 dlf_global ihessian glob_type
R 801 5 150 dlf_global ihessian$sd glob_type
R 802 5 151 dlf_global ihessian$p glob_type
R 803 5 152 dlf_global ihessian$o glob_type
R 805 5 154 dlf_global step glob_type
R 807 5 156 dlf_global step$sd glob_type
R 808 5 157 dlf_global step$p glob_type
R 809 5 158 dlf_global step$o glob_type
R 811 5 160 dlf_global spec glob_type
R 813 5 162 dlf_global spec$sd glob_type
R 814 5 163 dlf_global spec$p glob_type
R 815 5 164 dlf_global spec$o glob_type
R 817 5 166 dlf_global icons glob_type
R 820 5 169 dlf_global icons$sd glob_type
R 821 5 170 dlf_global icons$p glob_type
R 822 5 171 dlf_global icons$o glob_type
R 824 5 173 dlf_global iconn glob_type
R 827 5 176 dlf_global iconn$sd glob_type
R 828 5 177 dlf_global iconn$p glob_type
R 829 5 178 dlf_global iconn$o glob_type
R 831 5 180 dlf_global msenergy glob_type
R 833 5 182 dlf_global msenergy$sd glob_type
R 834 5 183 dlf_global msenergy$p glob_type
R 835 5 184 dlf_global msenergy$o glob_type
R 837 5 186 dlf_global msgradient glob_type
R 841 5 190 dlf_global msgradient$sd glob_type
R 842 5 191 dlf_global msgradient$p glob_type
R 843 5 192 dlf_global msgradient$o glob_type
R 845 5 194 dlf_global mscoupling glob_type
R 848 5 197 dlf_global mscoupling$sd glob_type
R 849 5 198 dlf_global mscoupling$p glob_type
R 850 5 199 dlf_global mscoupling$o glob_type
R 852 5 201 dlf_global po_radius glob_type
R 854 5 203 dlf_global po_radius$sd glob_type
R 855 5 204 dlf_global po_radius$p glob_type
R 856 5 205 dlf_global po_radius$o glob_type
R 858 5 207 dlf_global po_tolerance_r glob_type
R 860 5 209 dlf_global po_tolerance_r$sd glob_type
R 861 5 210 dlf_global po_tolerance_r$p glob_type
R 862 5 211 dlf_global po_tolerance_r$o glob_type
R 864 5 213 dlf_global po_radius_scaling glob_type
R 866 5 215 dlf_global po_radius_scaling$sd glob_type
R 867 5 216 dlf_global po_radius_scaling$p glob_type
R 868 5 217 dlf_global po_radius_scaling$o glob_type
R 870 5 219 dlf_global micspec glob_type
R 872 5 221 dlf_global micspec$sd glob_type
R 873 5 222 dlf_global micspec$p glob_type
R 874 5 223 dlf_global micspec$o glob_type
R 876 5 225 dlf_global nicore glob_type
R 877 5 226 dlf_global maxmicrocycle glob_type
R 878 5 227 dlf_global micro_esp_fit glob_type
R 879 5 228 dlf_global macrocoords glob_type
R 883 5 232 dlf_global macrocoords$sd glob_type
R 884 5 233 dlf_global macrocoords$p glob_type
R 885 5 234 dlf_global macrocoords$o glob_type
R 887 5 236 dlf_global g0corr glob_type
R 891 5 240 dlf_global g0corr$sd glob_type
R 892 5 241 dlf_global g0corr$p glob_type
R 893 5 242 dlf_global g0corr$o glob_type
R 895 5 244 dlf_global e0corr glob_type
R 897 5 246 dlf_global e0corr$sd glob_type
R 898 5 247 dlf_global e0corr$p glob_type
R 899 5 248 dlf_global e0corr$o glob_type
R 904 6 253 dlf_global stdout
R 905 6 254 dlf_global stderr
R 918 14 4 dlf_allocate allocate
R 920 14 6 dlf_allocate deallocate
S 921 19 0 0 0 9 1 624 7736 4000 0 A 0 0 0 0 B 0 393 0 0 0 0 0 0 0 927 0 0 0 0 0 0 7 1 0 0 0 0 0 624 0 0 0 0 store_initialise
O 921 1 927
S 922 19 0 0 0 9 1 624 7753 4000 0 A 0 0 0 0 B 0 393 0 0 0 0 0 0 0 928 0 0 0 0 0 0 9 1 0 0 0 0 0 624 0 0 0 0 store_allocate
O 922 1 928
S 923 19 0 0 0 9 1 624 7768 4000 0 A 0 0 0 0 B 0 393 0 0 0 0 0 0 0 929 0 0 0 0 0 0 15 3 0 0 0 0 0 624 0 0 0 0 store_set
O 923 3 931 930 929
S 924 19 0 0 0 9 1 624 7778 4000 0 A 0 0 0 0 B 0 393 0 0 0 0 0 0 0 932 0 0 0 0 0 0 21 3 0 0 0 0 0 624 0 0 0 0 store_get
O 924 3 934 933 932
S 925 19 0 0 0 9 1 624 7788 4000 0 A 0 0 0 0 B 0 393 0 0 0 0 0 0 0 935 0 0 0 0 0 0 23 1 0 0 0 0 0 624 0 0 0 0 store_delete
O 925 1 935
S 926 19 0 0 0 9 1 624 7801 4000 0 A 0 0 0 0 B 0 393 0 0 0 0 0 0 0 936 0 0 0 0 0 0 25 1 0 0 0 0 0 624 0 0 0 0 store_delete_all
O 926 1 936
S 927 27 0 0 0 9 962 624 7736 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 26 0 0 921 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 store_initialise
Q 927 921 0
S 928 27 0 0 0 9 964 624 7753 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 27 0 0 922 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 store_allocate
Q 928 922 0
S 929 27 0 0 0 9 968 624 7768 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 28 0 0 923 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 store_set
Q 929 923 0
S 930 27 0 0 0 9 985 624 7818 10 400000 A 0 0 0 0 B 0 406 0 0 0 0 32 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 store_set_a2
Q 930 923 0
S 931 27 0 0 0 9 997 624 7831 10 400000 A 0 0 0 0 B 0 407 0 0 0 0 33 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 store_set_a3
Q 931 923 0
S 932 27 0 0 0 9 974 624 7778 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 29 0 0 924 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 store_get
Q 932 924 0
S 933 27 0 0 0 9 1012 624 7844 10 400000 A 0 0 0 0 B 0 412 0 0 0 0 34 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 store_get_a2
Q 933 924 0
S 934 27 0 0 0 9 1025 624 7857 10 400000 A 0 0 0 0 B 0 413 0 0 0 0 35 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 store_get_a3
Q 934 924 0
S 935 27 0 0 0 9 980 624 7788 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 30 0 0 925 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 store_delete
Q 935 925 0
S 936 27 0 0 0 9 983 624 7801 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 31 0 0 926 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 store_delete_all
Q 936 926 0
S 937 25 0 0 0 214 1 624 7870 10000014 800014 A 0 0 0 0 B 0 426 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 953 0 0 0 624 0 0 0 0 store_type_r
S 939 3 0 0 0 6 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 40 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 940 5 0 0 0 220 941 624 7883 800004 0 A 0 0 0 0 B 0 427 0 0 0 0 0 0 214 0 0 0 0 0 0 0 0 0 0 0 1 940 0 624 0 0 0 0 tag
S 941 5 0 0 0 7 942 624 2875 800004 0 A 0 0 0 0 B 0 0 0 0 0 40 0 0 214 0 0 0 0 0 0 0 0 0 0 0 940 941 0 624 0 0 0 0 size
S 942 5 6 0 0 222 945 624 7887 10a00004 14 A 0 0 0 0 B 0 429 0 0 0 48 945 0 214 0 947 0 0 0 0 0 0 0 0 944 941 942 946 624 0 0 0 0 array
S 943 6 4 0 0 7 1 624 7893 40800016 0 A 0 0 0 0 B 0 429 0 0 0 0 0 0 0 0 0 0 959 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 z_b_0_1
S 944 5 0 0 0 225 948 624 7901 40822004 1020 A 0 0 0 0 B 0 429 0 0 0 64 0 0 214 0 0 0 0 0 0 0 0 0 0 0 946 944 0 624 0 0 0 0 array$sd
S 945 5 0 0 0 7 946 624 7910 40802001 1020 A 0 0 0 0 B 0 429 0 0 0 48 0 0 214 0 0 0 0 0 0 0 0 0 0 0 942 945 0 624 0 0 0 0 array$p
S 946 5 0 0 0 7 944 624 7918 40802000 1020 A 0 0 0 0 B 0 429 0 0 0 56 0 0 214 0 0 0 0 0 0 0 0 0 0 0 945 946 0 624 0 0 0 0 array$o
S 947 22 1 0 0 9 1 624 7926 40000000 1000 A 0 0 0 0 B 0 429 0 0 0 0 0 942 0 0 0 0 944 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 array$arrdsc
S 948 5 6 0 0 214 951 624 7939 801004 14 A 0 0 0 0 B 0 430 0 0 0 192 951 0 214 0 0 0 0 0 0 0 0 0 0 950 942 948 952 624 0 0 0 0 next
S 949 3 0 0 0 7 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 950 8 1 0 0 228 1 624 7944 40822006 1020 A 0 0 0 0 B 0 430 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 next$sd
S 951 5 0 0 0 7 1 624 7952 40802001 1020 A 0 0 0 0 B 0 430 0 0 0 192 0 0 214 0 0 0 0 0 0 0 0 0 0 0 948 951 0 624 0 0 0 0 next$p
S 952 6 1 0 0 7 1 624 7959 40802010 1020 A 0 0 0 0 B 0 430 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 next$o
S 953 8 5 0 0 231 1 624 7966 40822004 1220 A 0 0 0 0 B 0 431 0 0 0 0 0 214 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 dlf_store$$store_type_r$$$td
S 954 6 6 0 0 214 1 624 7995 34 14 A 0 0 0 0 B 0 433 0 0 0 0 956 0 0 0 0 0 0 0 0 0 0 0 0 955 0 0 957 624 0 0 0 0 first_r
S 955 8 4 0 0 234 1 624 8003 40822016 1020 A 0 0 0 0 B 0 433 0 0 0 0 0 0 0 0 0 0 961 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 first_r$sd
S 956 6 4 0 0 7 957 624 8014 40802011 1020 A 0 0 0 0 B 0 433 0 0 0 0 0 0 0 0 0 0 961 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 first_r$p
S 957 6 4 0 0 7 955 624 8024 40802010 1020 A 0 0 0 0 B 0 433 0 0 0 0 0 0 0 0 0 0 961 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 first_r$o
S 958 6 4 0 0 19 1 624 5177 80003c 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 960 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 tinit
S 959 11 0 0 0 9 913 624 8034 40800010 805000 A 0 0 0 0 B 0 436 0 0 0 8 0 0 943 943 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _dlf_store$6
S 960 11 0 0 0 9 959 624 8047 40800010 805000 A 0 0 0 0 B 0 436 0 0 0 8 0 0 958 958 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _dlf_store$14
S 961 11 0 0 0 9 960 624 8061 40800010 805000 A 0 0 0 0 B 0 436 0 0 0 96 0 0 956 955 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _dlf_store$4
S 962 23 5 0 0 0 963 624 7736 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 921 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 store_initialise
S 963 14 5 0 0 0 1 962 7736 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 5 0 921 0 0 0 0 0 0 0 0 0 0 0 438 0 624 0 0 0 0 store_initialise
F 963 0
S 964 23 5 0 0 0 967 624 7753 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 922 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 store_allocate
S 965 1 3 1 0 30 1 964 7883 14 43000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 tag
S 966 1 3 1 0 7 1 964 2875 14 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 size
S 967 14 5 0 0 0 1 964 7753 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 6 2 922 0 0 0 0 0 0 0 0 0 0 0 449 0 624 0 0 0 0 store_allocate
F 967 2 965 966
S 968 23 5 0 0 0 972 624 7768 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 923 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 store_set
S 969 1 3 1 0 30 1 968 7883 14 43000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 tag
S 970 6 3 1 0 7 1 968 2875 800014 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 size
S 971 7 3 1 0 237 1 968 7887 800214 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 array
S 972 14 5 0 0 0 1 968 7768 200 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 9 3 923 0 0 0 0 0 0 0 0 0 0 0 482 0 624 0 0 0 0 store_set
F 972 3 969 970 971
S 973 6 1 0 0 7 1 968 8074 40800016 3000 A 0 0 0 0 B 0 485 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_309
S 974 23 5 0 0 0 978 624 7778 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 924 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 store_get
S 975 1 3 1 0 30 1 974 7883 14 43000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 tag
S 976 6 3 1 0 7 1 974 2875 800014 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 size
S 977 7 3 2 0 240 1 974 7887 800214 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 array
S 978 14 5 0 0 0 1 974 7778 200 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 13 3 924 0 0 0 0 0 0 0 0 0 0 0 510 0 624 0 0 0 0 store_get
F 978 3 975 976 977
S 979 6 1 0 0 7 1 974 8082 40800016 3000 A 0 0 0 0 B 0 513 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_311
S 980 23 5 0 0 0 982 624 7788 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 925 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 store_delete
S 981 1 3 1 0 30 1 980 7883 14 43000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 tag
S 982 14 5 0 0 0 1 980 7788 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 17 1 925 0 0 0 0 0 0 0 0 0 0 0 537 0 624 0 0 0 0 store_delete
F 982 1 981
S 983 23 5 0 0 0 984 624 7801 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 926 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 store_delete_all
S 984 14 5 0 0 0 1 983 7801 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 19 0 926 0 0 0 0 0 0 0 0 0 0 0 596 0 624 0 0 0 0 store_delete_all
F 984 0
S 985 23 5 0 0 0 989 624 7818 10 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 store_set_a2
S 986 1 3 1 0 30 1 985 7883 14 43000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 tag
S 987 1 3 1 0 7 1 985 2875 14 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 size
S 988 7 3 1 0 243 1 985 7887 20000014 10003000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 array
S 989 14 5 0 0 0 1 985 7818 20000010 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 20 3 0 0 0 0 0 0 0 0 0 0 0 0 624 0 624 0 0 0 0 store_set_a2
F 989 3 986 987 988
S 990 6 1 0 0 7 1 985 8090 40800016 3000 A 0 0 0 0 B 0 628 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_1_1
S 991 6 1 0 0 7 1 985 8098 40800016 3000 A 0 0 0 0 B 0 628 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_2_1
S 992 6 1 0 0 7 1 985 8106 40800016 3000 A 0 0 0 0 B 0 628 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_4_1
S 993 6 1 0 0 7 1 985 8114 40800016 3000 A 0 0 0 0 B 0 628 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_5_1
S 994 6 1 0 0 7 1 985 8122 40800016 3000 A 0 0 0 0 B 0 628 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_6_1
S 995 6 1 0 0 7 1 985 8130 40800016 3000 A 0 0 0 0 B 0 628 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_321
S 996 6 1 0 0 7 1 985 8138 40800016 3000 A 0 0 0 0 B 0 628 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_324
S 997 23 5 0 0 0 1001 624 7831 10 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 store_set_a3
S 998 1 3 1 0 30 1 997 7883 14 43000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 tag
S 999 1 3 1 0 7 1 997 2875 14 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 size
S 1000 7 3 1 0 246 1 997 7887 20000014 10003000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 array
S 1001 14 5 0 0 0 1 997 7831 20000010 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 24 3 0 0 0 0 0 0 0 0 0 0 0 0 632 0 624 0 0 0 0 store_set_a3
F 1001 3 998 999 1000
S 1002 6 1 0 0 7 1 997 8090 40800016 3000 A 0 0 0 0 B 0 636 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_1_1
S 1003 6 1 0 0 7 1 997 8098 40800016 3000 A 0 0 0 0 B 0 636 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_2_1
S 1004 6 1 0 0 7 1 997 8106 40800016 3000 A 0 0 0 0 B 0 636 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_4_1
S 1005 6 1 0 0 7 1 997 8114 40800016 3000 A 0 0 0 0 B 0 636 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_5_1
S 1006 6 1 0 0 7 1 997 8146 40800016 3000 A 0 0 0 0 B 0 636 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_7_1
S 1007 6 1 0 0 7 1 997 8154 40800016 3000 A 0 0 0 0 B 0 636 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_8_1
S 1008 6 1 0 0 7 1 997 8162 40800016 3000 A 0 0 0 0 B 0 636 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_9_1
S 1009 6 1 0 0 7 1 997 8170 40800016 3000 A 0 0 0 0 B 0 636 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_331
S 1010 6 1 0 0 7 1 997 8178 40800016 3000 A 0 0 0 0 B 0 636 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_334
S 1011 6 1 0 0 7 1 997 8186 40800016 3000 A 0 0 0 0 B 0 636 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_337
S 1012 23 5 0 0 0 1017 624 7844 10 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 store_get_a2
S 1013 1 3 1 0 30 1 1012 7883 14 43000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 tag
S 1014 6 3 1 0 7 1 1012 8194 800014 3000 A 0 0 0 0 B 0 640 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 size1
S 1015 6 3 1 0 7 1 1012 8200 800014 3000 A 0 0 0 0 B 0 640 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 size2
S 1016 7 3 2 0 249 1 1012 7887 20000014 10003000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 array
S 1017 14 5 0 0 0 1 1012 7844 20000010 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 28 4 0 0 0 0 0 0 0 0 0 0 0 0 640 0 624 0 0 0 0 store_get_a2
F 1017 4 1013 1014 1015 1016
S 1018 6 1 0 0 7 1 1012 8090 40800016 3000 A 0 0 0 0 B 0 645 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_1_1
S 1019 6 1 0 0 7 1 1012 8098 40800016 3000 A 0 0 0 0 B 0 645 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_2_1
S 1020 6 1 0 0 7 1 1012 8106 40800016 3000 A 0 0 0 0 B 0 645 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_4_1
S 1021 6 1 0 0 7 1 1012 8114 40800016 3000 A 0 0 0 0 B 0 645 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_5_1
S 1022 6 1 0 0 7 1 1012 8122 40800016 3000 A 0 0 0 0 B 0 645 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_6_1
S 1023 6 1 0 0 7 1 1012 8206 40800016 3000 A 0 0 0 0 B 0 645 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_338
S 1024 6 1 0 0 7 1 1012 8214 40800016 3000 A 0 0 0 0 B 0 645 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_341
S 1025 23 5 0 0 0 1031 624 7857 10 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 store_get_a3
S 1026 1 3 1 0 30 1 1025 7883 14 43000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 tag
S 1027 6 3 1 0 7 1 1025 8194 800014 3000 A 0 0 0 0 B 0 651 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 size1
S 1028 6 3 1 0 7 1 1025 8200 800014 3000 A 0 0 0 0 B 0 651 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 size2
S 1029 6 3 1 0 7 1 1025 8222 800014 3000 A 0 0 0 0 B 0 651 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 size3
S 1030 7 3 2 0 252 1 1025 7887 20000014 10003000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 array
S 1031 14 5 0 0 0 1 1025 7857 20000010 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 33 5 0 0 0 0 0 0 0 0 0 0 0 0 651 0 624 0 0 0 0 store_get_a3
F 1031 5 1026 1027 1028 1029 1030
S 1032 6 1 0 0 7 1 1025 8090 40800016 3000 A 0 0 0 0 B 0 657 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_1_1
S 1033 6 1 0 0 7 1 1025 8098 40800016 3000 A 0 0 0 0 B 0 657 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_2_1
S 1034 6 1 0 0 7 1 1025 8106 40800016 3000 A 0 0 0 0 B 0 657 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_4_1
S 1035 6 1 0 0 7 1 1025 8114 40800016 3000 A 0 0 0 0 B 0 657 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_5_1
S 1036 6 1 0 0 7 1 1025 8146 40800016 3000 A 0 0 0 0 B 0 657 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_7_1
S 1037 6 1 0 0 7 1 1025 8154 40800016 3000 A 0 0 0 0 B 0 657 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_8_1
S 1038 6 1 0 0 7 1 1025 8162 40800016 3000 A 0 0 0 0 B 0 657 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_9_1
S 1039 6 1 0 0 7 1 1025 8228 40800016 3000 A 0 0 0 0 B 0 657 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_348
S 1040 6 1 0 0 7 1 1025 8236 40800016 3000 A 0 0 0 0 B 0 657 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_351
S 1041 6 1 0 0 7 1 1025 8244 40800016 3000 A 0 0 0 0 B 0 657 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_354
A 13 2 0 0 0 7 633 0 0 0 13 0 0 0 0 0 0 0 0 0 0 0
A 18 2 0 0 0 7 639 0 0 0 18 0 0 0 0 0 0 0 0 0 0 0
A 20 2 0 0 0 7 640 0 0 0 20 0 0 0 0 0 0 0 0 0 0 0
A 24 2 0 0 0 7 643 0 0 0 24 0 0 0 0 0 0 0 0 0 0 0
A 34 2 0 0 0 7 645 0 0 0 34 0 0 0 0 0 0 0 0 0 0 0
A 71 2 0 0 0 7 650 0 0 0 71 0 0 0 0 0 0 0 0 0 0 0
A 290 2 0 0 0 19 636 0 0 0 290 0 0 0 0 0 0 0 0 0 0 0
A 291 2 0 0 0 19 637 0 0 0 291 0 0 0 0 0 0 0 0 0 0 0
A 297 2 0 0 0 6 939 0 0 0 297 0 0 0 0 0 0 0 0 0 0 0
A 299 1 0 5 0 225 944 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 300 10 0 0 0 7 299 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 18
A 301 10 0 0 300 7 299 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 20
A 302 4 0 0 0 7 301 0 11 0 0 0 0 2 0 0 0 0 0 0 0 0
A 303 4 0 0 0 7 300 0 302 0 0 0 0 1 0 0 0 0 0 0 0 0
A 304 10 0 0 301 7 299 10 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 24
A 305 10 0 0 304 7 299 22 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 34
A 306 10 0 0 305 7 299 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 13
A 307 2 0 0 0 7 949 0 0 0 307 0 0 0 0 0 0 0 0 0 0 0
A 309 1 0 0 0 7 970 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 310 1 0 0 0 7 973 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 311 1 0 0 0 7 976 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 312 1 0 0 0 7 979 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 313 1 0 0 0 7 994 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 314 1 0 0 0 7 990 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 315 1 0 0 0 7 995 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 316 1 0 0 0 7 992 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 317 1 0 0 0 7 991 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 318 1 0 0 0 7 996 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 319 1 0 0 0 7 993 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 320 1 0 0 0 7 1008 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 321 1 0 0 0 7 1002 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 322 1 0 0 0 7 1009 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 323 1 0 0 0 7 1004 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 324 1 0 0 299 7 1003 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 325 1 0 0 0 7 1010 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 326 1 0 0 0 7 1006 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 327 1 0 0 0 7 1005 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 328 1 0 0 0 7 1011 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 329 1 0 0 0 7 1007 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 330 1 0 0 0 7 1022 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 331 1 0 0 0 7 1018 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 332 1 0 0 172 7 1023 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 333 1 0 0 0 7 1020 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 334 1 0 0 0 7 1019 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 335 1 0 0 35 7 1024 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 336 1 0 0 0 7 1021 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 337 1 0 0 0 7 1038 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 338 1 0 0 0 7 1032 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 339 1 0 0 0 7 1039 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 340 1 0 0 0 7 1034 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 341 1 0 0 0 7 1033 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 342 1 0 0 0 7 1040 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 343 1 0 0 0 7 1036 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 344 1 0 0 0 7 1035 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 345 1 0 0 0 7 1041 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 346 1 0 0 0 7 1037 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
Z
T 653 61 0 3 0 0
A 661 19 0 0 1 290 1
A 732 7 0 0 1 11 1
A 733 7 0 0 1 11 1
A 734 7 0 0 1 11 1
A 735 7 0 0 1 10 1
A 736 7 0 0 1 10 1
A 737 7 0 0 1 10 1
A 739 19 0 0 1 291 0
Z
