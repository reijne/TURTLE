V34 :0x2c lbfgs_module
15 dlf_lbfgs_d.F90 S624 0
04/08/2020  10:45:19
use dlf_parameter_module private
enduse
D 61 26 631 1248 630 7
D 67 23 10 1 28 27 0 1 0 0 1
 19 23 25 19 23 21
D 70 23 7 1 0 16 0 0 0 0 0
 0 16 0 11 16 0
D 73 23 10 1 37 36 0 1 0 0 1
 31 34 35 31 34 32
D 76 23 7 1 0 16 0 0 0 0 0
 0 16 0 11 16 0
D 79 23 10 1 46 45 0 1 0 0 1
 40 43 44 40 43 41
D 82 23 7 1 0 16 0 0 0 0 0
 0 16 0 11 16 0
D 85 23 10 1 55 54 0 1 0 0 1
 49 52 53 49 52 50
D 88 23 7 1 0 16 0 0 0 0 0
 0 16 0 11 16 0
D 91 23 10 2 74 73 0 1 0 0 1
 60 63 70 60 63 61
 65 69 72 65 69 67
D 94 23 7 1 0 58 0 0 0 0 0
 0 58 0 11 58 0
D 97 23 10 2 89 88 0 1 0 0 1
 78 81 86 78 81 79
 82 85 87 82 85 83
D 100 23 7 1 0 58 0 0 0 0 0
 0 58 0 11 58 0
D 103 23 10 2 104 103 0 1 0 0 1
 93 96 101 93 96 94
 97 100 102 97 100 98
D 106 23 7 1 0 58 0 0 0 0 0
 0 58 0 11 58 0
D 109 20 106
D 111 23 7 1 0 107 0 0 0 0 0
 0 107 0 11 107 0
D 114 22 7
D 116 22 7
D 118 22 7
D 120 22 7
D 122 22 7
D 124 22 7
D 126 22 7
D 128 23 7 1 0 11 0 0 0 0 0
 0 11 0 11 11 0
D 131 23 7 1 0 107 0 0 0 0 0
 0 107 0 11 107 0
D 134 23 7 1 0 107 0 0 0 0 0
 0 107 0 11 107 0
D 137 20 111
S 624 24 0 0 0 7 1 0 5013 10005 0 A 0 0 0 0 B 0 89 0 0 0 0 0 0 0 0 0 0 116 0 0 0 0 0 0 0 0 89 0 0 0 0 0 0 lbfgs_module
S 626 23 0 0 0 9 628 624 5047 4 0 A 0 0 0 0 B 0 90 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 rk
S 627 3 0 0 0 7 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
R 628 16 1 dlf_parameter_module rk
S 630 25 0 0 0 61 1 624 5053 1000000c 800054 A 0 0 0 0 B 0 91 0 0 0 0 0 0 0 0 0 701 0 0 0 0 0 0 0 700 0 0 0 624 0 0 0 0 lbfgs_type
S 631 5 0 0 0 7 632 624 5064 800004 0 A 0 0 0 0 B 0 92 0 0 0 0 0 0 61 0 0 0 0 0 0 0 0 0 0 0 1 631 0 624 0 0 0 0 n
S 632 5 0 0 0 7 633 624 5066 800004 0 A 0 0 0 0 B 0 93 0 0 0 8 0 0 61 0 0 0 0 0 0 0 0 0 0 0 631 632 0 624 0 0 0 0 m
S 633 5 6 0 0 67 637 624 5068 10a00004 51 A 0 0 0 0 B 0 94 0 0 0 16 637 0 61 0 639 0 0 0 0 0 0 0 0 636 632 633 638 624 0 0 0 0 store
S 634 6 4 0 0 7 646 624 5074 40800006 0 A 0 0 0 0 B 0 94 0 0 0 0 0 0 0 0 0 0 717 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 z_b_0
S 635 3 0 0 0 7 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 16 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 636 5 0 0 0 70 645 624 5080 40822004 1020 A 0 0 0 0 B 0 94 0 0 0 32 0 0 61 0 0 0 0 0 0 0 0 0 0 0 638 636 0 624 0 0 0 0 store$sd
S 637 5 0 0 0 7 638 624 5089 40802001 1020 A 0 0 0 0 B 0 94 0 0 0 16 0 0 61 0 0 0 0 0 0 0 0 0 0 0 633 637 0 624 0 0 0 0 store$p
S 638 5 0 0 0 7 636 624 5097 40802000 1020 A 0 0 0 0 B 0 94 0 0 0 24 0 0 61 0 0 0 0 0 0 0 0 0 0 0 637 638 0 624 0 0 0 0 store$o
S 639 22 1 0 0 9 1 624 5105 40000000 1000 A 0 0 0 0 B 0 94 0 0 0 0 0 633 0 0 0 0 636 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 store$arrdsc
S 640 3 0 0 0 7 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 11 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 641 3 0 0 0 7 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 12 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 643 3 0 0 0 7 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 15 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 644 3 0 0 0 7 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 645 5 6 0 0 73 648 624 5118 10a00004 51 A 0 0 0 0 B 0 95 0 0 0 160 648 0 61 0 650 0 0 0 0 0 0 0 0 647 633 645 649 624 0 0 0 0 store2
S 646 6 4 0 0 7 652 624 5125 40800006 0 A 0 0 0 0 B 0 95 0 0 0 8 0 0 0 0 0 0 717 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 z_b_1
S 647 5 0 0 0 76 651 624 5131 40822004 1020 A 0 0 0 0 B 0 95 0 0 0 176 0 0 61 0 0 0 0 0 0 0 0 0 0 0 649 647 0 624 0 0 0 0 store2$sd
S 648 5 0 0 0 7 649 624 5141 40802001 1020 A 0 0 0 0 B 0 95 0 0 0 160 0 0 61 0 0 0 0 0 0 0 0 0 0 0 645 648 0 624 0 0 0 0 store2$p
S 649 5 0 0 0 7 647 624 5150 40802000 1020 A 0 0 0 0 B 0 95 0 0 0 168 0 0 61 0 0 0 0 0 0 0 0 0 0 0 648 649 0 624 0 0 0 0 store2$o
S 650 22 1 0 0 9 1 624 5159 40000000 1000 A 0 0 0 0 B 0 95 0 0 0 0 0 645 0 0 0 0 647 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 store2$arrdsc
S 651 5 6 0 0 79 654 624 5173 10a00004 51 A 0 0 0 0 B 0 96 0 0 0 304 654 0 61 0 656 0 0 0 0 0 0 0 0 653 645 651 655 624 0 0 0 0 rho
S 652 6 4 0 0 7 658 624 5177 40800006 0 A 0 0 0 0 B 0 96 0 0 0 16 0 0 0 0 0 0 717 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 z_b_2
S 653 5 0 0 0 82 657 624 5183 40822004 1020 A 0 0 0 0 B 0 96 0 0 0 320 0 0 61 0 0 0 0 0 0 0 0 0 0 0 655 653 0 624 0 0 0 0 rho$sd
S 654 5 0 0 0 7 655 624 5190 40802001 1020 A 0 0 0 0 B 0 96 0 0 0 304 0 0 61 0 0 0 0 0 0 0 0 0 0 0 651 654 0 624 0 0 0 0 rho$p
S 655 5 0 0 0 7 653 624 5196 40802000 1020 A 0 0 0 0 B 0 96 0 0 0 312 0 0 61 0 0 0 0 0 0 0 0 0 0 0 654 655 0 624 0 0 0 0 rho$o
S 656 22 1 0 0 9 1 624 5202 40000000 1000 A 0 0 0 0 B 0 96 0 0 0 0 0 651 0 0 0 0 653 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 rho$arrdsc
S 657 5 6 0 0 85 660 624 5213 10a00004 51 A 0 0 0 0 B 0 97 0 0 0 448 660 0 61 0 662 0 0 0 0 0 0 0 0 659 651 657 661 624 0 0 0 0 alpha
S 658 6 4 0 0 7 664 624 5219 40800006 0 A 0 0 0 0 B 0 97 0 0 0 24 0 0 0 0 0 0 717 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 z_b_3
S 659 5 0 0 0 88 663 624 5225 40822004 1020 A 0 0 0 0 B 0 97 0 0 0 464 0 0 61 0 0 0 0 0 0 0 0 0 0 0 661 659 0 624 0 0 0 0 alpha$sd
S 660 5 0 0 0 7 661 624 5234 40802001 1020 A 0 0 0 0 B 0 97 0 0 0 448 0 0 61 0 0 0 0 0 0 0 0 0 0 0 657 660 0 624 0 0 0 0 alpha$p
S 661 5 0 0 0 7 659 624 5242 40802000 1020 A 0 0 0 0 B 0 97 0 0 0 456 0 0 61 0 0 0 0 0 0 0 0 0 0 0 660 661 0 624 0 0 0 0 alpha$o
S 662 22 1 0 0 9 1 624 5250 40000000 1000 A 0 0 0 0 B 0 97 0 0 0 0 0 657 0 0 0 0 659 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 alpha$arrdsc
S 663 5 6 0 0 91 668 624 5263 10a00004 51 A 0 0 0 0 B 0 98 0 0 0 592 668 0 61 0 670 0 0 0 0 0 0 0 0 667 657 663 669 624 0 0 0 0 step
S 664 6 4 0 0 7 665 624 5268 40800006 0 A 0 0 0 0 B 0 98 0 0 0 32 0 0 0 0 0 0 717 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 z_b_4
S 665 6 4 0 0 7 675 624 5274 40800006 0 A 0 0 0 0 B 0 98 0 0 0 40 0 0 0 0 0 0 717 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 z_b_5
S 666 3 0 0 0 7 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 22 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 667 5 0 0 0 94 674 624 5280 40822004 1020 A 0 0 0 0 B 0 98 0 0 0 608 0 0 61 0 0 0 0 0 0 0 0 0 0 0 669 667 0 624 0 0 0 0 step$sd
S 668 5 0 0 0 7 669 624 5288 40802001 1020 A 0 0 0 0 B 0 98 0 0 0 592 0 0 61 0 0 0 0 0 0 0 0 0 0 0 663 668 0 624 0 0 0 0 step$p
S 669 5 0 0 0 7 667 624 5295 40802000 1020 A 0 0 0 0 B 0 98 0 0 0 600 0 0 61 0 0 0 0 0 0 0 0 0 0 0 668 669 0 624 0 0 0 0 step$o
S 670 22 1 0 0 9 1 624 5302 40000000 1000 A 0 0 0 0 B 0 98 0 0 0 0 0 663 0 0 0 0 667 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 step$arrdsc
S 671 3 0 0 0 7 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 17 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 672 3 0 0 0 7 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 18 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 673 3 0 0 0 7 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 21 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 674 5 6 0 0 97 678 624 5314 10a00004 51 A 0 0 0 0 B 0 99 0 0 0 784 678 0 61 0 680 0 0 0 0 0 0 0 0 677 663 674 679 624 0 0 0 0 dgrad
S 675 6 4 0 0 7 676 624 5320 40800006 0 A 0 0 0 0 B 0 99 0 0 0 48 0 0 0 0 0 0 717 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 z_b_6
S 676 6 4 0 0 7 683 624 5326 40800006 0 A 0 0 0 0 B 0 99 0 0 0 56 0 0 0 0 0 0 717 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 z_b_7
S 677 5 0 0 0 100 681 624 5332 40822004 1020 A 0 0 0 0 B 0 99 0 0 0 800 0 0 61 0 0 0 0 0 0 0 0 0 0 0 679 677 0 624 0 0 0 0 dgrad$sd
S 678 5 0 0 0 7 679 624 5341 40802001 1020 A 0 0 0 0 B 0 99 0 0 0 784 0 0 61 0 0 0 0 0 0 0 0 0 0 0 674 678 0 624 0 0 0 0 dgrad$p
S 679 5 0 0 0 7 677 624 5349 40802000 1020 A 0 0 0 0 B 0 99 0 0 0 792 0 0 61 0 0 0 0 0 0 0 0 0 0 0 678 679 0 624 0 0 0 0 dgrad$o
S 680 22 1 0 0 9 1 624 5357 40000000 1000 A 0 0 0 0 B 0 99 0 0 0 0 0 674 0 0 0 0 677 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 dgrad$arrdsc
S 681 5 0 0 0 19 682 624 5370 800004 0 A 0 0 0 0 B 0 100 0 0 0 976 0 0 61 0 0 0 0 0 0 0 0 0 0 0 674 681 0 624 0 0 0 0 tprecon
S 682 5 6 0 0 103 686 624 5378 10a00004 51 A 0 0 0 0 B 0 101 0 0 0 984 686 0 61 0 688 0 0 0 0 0 0 0 0 685 681 682 687 624 0 0 0 0 precon
S 683 6 4 0 0 7 684 624 5385 40800006 0 A 0 0 0 0 B 0 101 0 0 0 64 0 0 0 0 0 0 717 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 z_b_8
S 684 6 4 0 0 7 1 624 5391 40800006 0 A 0 0 0 0 B 0 101 0 0 0 72 0 0 0 0 0 0 717 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 z_b_9
S 685 5 0 0 0 106 689 624 5397 40822004 1020 A 0 0 0 0 B 0 101 0 0 0 1000 0 0 61 0 0 0 0 0 0 0 0 0 0 0 687 685 0 624 0 0 0 0 precon$sd
S 686 5 0 0 0 7 687 624 5407 40802001 1020 A 0 0 0 0 B 0 101 0 0 0 984 0 0 61 0 0 0 0 0 0 0 0 0 0 0 682 686 0 624 0 0 0 0 precon$p
S 687 5 0 0 0 7 685 624 5416 40802000 1020 A 0 0 0 0 B 0 101 0 0 0 992 0 0 61 0 0 0 0 0 0 0 0 0 0 0 686 687 0 624 0 0 0 0 precon$o
S 688 22 1 0 0 9 1 624 5425 40000000 1000 A 0 0 0 0 B 0 101 0 0 0 0 0 682 0 0 0 0 685 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 precon$arrdsc
S 689 5 0 0 0 7 690 624 5439 800004 0 A 0 0 0 0 B 0 105 0 0 0 1176 0 0 61 0 0 0 0 0 0 0 0 0 0 0 682 689 0 624 0 0 0 0 point
S 690 5 0 0 0 7 691 624 5445 800004 0 A 0 0 0 0 B 0 106 0 0 0 1184 0 0 61 0 0 0 0 0 0 0 0 0 0 0 689 690 0 624 0 0 0 0 iter
S 691 5 0 0 0 19 694 624 5450 800004 0 A 0 0 0 0 B 0 107 0 0 0 1192 0 0 61 0 0 0 0 0 0 0 0 0 0 0 690 691 0 624 0 0 0 0 tinit
S 693 3 0 0 0 6 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 40 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 694 5 0 0 0 109 695 624 5456 800004 0 A 0 0 0 0 B 0 108 0 0 0 1200 0 0 61 0 0 0 0 0 0 0 0 0 0 0 691 694 0 624 0 0 0 0 tag
S 695 5 6 0 0 61 698 624 5460 801004 14 A 0 0 0 0 B 0 109 0 0 0 1240 698 0 61 0 0 0 0 0 0 0 0 0 0 697 694 695 699 624 0 0 0 0 next
S 696 3 0 0 0 7 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 697 8 1 0 0 111 1 624 5465 40822006 1020 A 0 0 0 0 B 0 109 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 next$sd
S 698 5 0 0 0 7 1 624 5473 40802001 1020 A 0 0 0 0 B 0 109 0 0 0 1240 0 0 61 0 0 0 0 0 0 0 0 0 0 0 695 698 0 624 0 0 0 0 next$p
S 699 6 1 0 0 7 1 624 5480 40802000 1020 A 0 0 0 0 B 0 109 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 next$o
S 700 8 5 0 0 128 1 624 5487 40822004 1220 A 0 0 0 0 B 0 110 0 0 0 0 0 61 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 lbfgs_module$$lbfgs_type$$td
S 701 6 4 0 0 61 1 624 5516 80004e 0 A 0 0 0 0 B 800 110 0 0 0 0 0 0 0 0 0 0 718 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 ._dtInit0061
S 702 6 6 0 0 61 1 624 5529 24 14 A 0 0 0 0 B 0 111 0 0 0 0 704 0 0 0 0 0 0 0 0 0 0 0 0 703 0 0 705 624 0 0 0 0 lbfgs
S 703 8 4 0 0 131 708 624 5535 40822006 1020 A 0 0 0 0 B 0 111 0 0 0 0 0 0 0 0 0 0 721 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 lbfgs$sd
S 704 6 4 0 0 7 705 624 5544 40802001 1020 A 0 0 0 0 B 0 111 0 0 0 0 0 0 0 0 0 0 721 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 lbfgs$p
S 705 6 4 0 0 7 703 624 5552 40802000 1020 A 0 0 0 0 B 0 111 0 0 0 0 0 0 0 0 0 0 721 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 lbfgs$o
S 706 6 6 0 0 61 1 624 5560 24 14 A 0 0 0 0 B 0 112 0 0 0 0 708 0 0 0 0 0 0 0 0 0 0 0 0 707 0 0 709 624 0 0 0 0 lbfgs_first
S 707 8 4 0 0 134 1 624 5572 40822006 1020 A 0 0 0 0 B 0 112 0 0 0 0 0 0 0 0 0 0 721 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 lbfgs_first$sd
S 708 6 4 0 0 7 709 624 5587 40802001 1020 A 0 0 0 0 B 0 112 0 0 0 0 0 0 0 0 0 0 721 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 lbfgs_first$p
S 709 6 4 0 0 7 707 624 5601 40802000 1020 A 0 0 0 0 B 0 112 0 0 0 0 0 0 0 0 0 0 721 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 lbfgs_first$o
S 710 6 4 0 0 19 1 624 5450 80002c 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 719 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 tinit
S 711 3 0 0 0 19 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 19
S 712 16 0 0 0 19 0 624 5615 800004 400000 A 0 0 0 0 B 0 114 0 0 0 0 0 0 711 108 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 dbg
S 713 6 4 0 0 109 1 624 5619 80002c 0 A 0 0 0 0 B 0 115 0 0 0 0 0 0 0 0 0 0 720 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 newtag
S 714 3 0 0 0 6 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 715 3 0 0 0 137 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 5626 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 20 4 6e 6f 6e 65
S 717 11 0 0 0 9 1 624 5672 40800000 805000 A 0 0 0 0 B 0 116 0 0 0 80 0 0 634 684 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _lbfgs_module$2
S 718 11 0 0 0 9 717 624 5688 40800008 805000 A 0 0 0 0 B 0 116 0 0 0 1248 0 0 701 701 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _lbfgs_module$8
S 719 11 0 0 0 9 718 624 5704 40800008 805000 A 0 0 0 0 B 0 116 0 0 0 8 0 0 710 710 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _lbfgs_module$10
S 720 11 0 0 0 9 719 624 5721 40800008 805000 A 0 0 0 0 B 0 116 0 0 0 40 0 0 713 713 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _lbfgs_module$9
S 721 11 0 0 0 9 720 624 5737 40800000 805000 A 0 0 0 0 B 0 116 0 0 0 192 0 0 704 707 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _lbfgs_module$0
A 13 2 0 0 0 7 627 0 0 0 13 0 0 0 0 0 0 0 0 0 0 0
A 16 2 0 0 0 7 635 0 0 0 16 0 0 0 0 0 0 0 0 0 0 0
A 17 2 0 0 0 7 640 0 0 0 17 0 0 0 0 0 0 0 0 0 0 0
A 18 1 0 1 0 70 636 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 19 10 0 0 0 7 18 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 17
A 20 2 0 0 0 7 641 0 0 0 20 0 0 0 0 0 0 0 0 0 0 0
A 21 10 0 0 19 7 18 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 20
A 22 4 0 0 0 7 21 0 11 0 0 0 0 2 0 0 0 0 0 0 0 0
A 23 4 0 0 7 7 19 0 22 0 0 0 0 1 0 0 0 0 0 0 0 0
A 24 2 0 0 0 7 643 0 0 0 24 0 0 0 0 0 0 0 0 0 0 0
A 25 10 0 0 21 7 18 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 24
A 26 2 0 0 0 7 644 0 0 0 26 0 0 0 0 0 0 0 0 0 0 0
A 27 10 0 0 25 7 18 10 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 26
A 28 10 0 0 27 7 18 13 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 13
A 30 1 0 1 0 76 647 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 31 10 0 0 0 7 30 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 17
A 32 10 0 0 31 7 30 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 20
A 33 4 0 0 0 7 32 0 11 0 0 0 0 2 0 0 0 0 0 0 0 0
A 34 4 0 0 0 7 31 0 33 0 0 0 0 1 0 0 0 0 0 0 0 0
A 35 10 0 0 32 7 30 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 24
A 36 10 0 0 35 7 30 10 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 26
A 37 10 0 0 36 7 30 13 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 13
A 39 1 0 1 0 82 653 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 40 10 0 0 0 7 39 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 17
A 41 10 0 0 40 7 39 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 20
A 42 4 0 0 0 7 41 0 11 0 0 0 0 2 0 0 0 0 0 0 0 0
A 43 4 0 0 0 7 40 0 42 0 0 0 0 1 0 0 0 0 0 0 0 0
A 44 10 0 0 41 7 39 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 24
A 45 10 0 0 44 7 39 10 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 26
A 46 10 0 0 45 7 39 13 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 13
A 48 1 0 1 0 88 659 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 49 10 0 0 0 7 48 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 17
A 50 10 0 0 49 7 48 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 20
A 51 4 0 0 0 7 50 0 11 0 0 0 0 2 0 0 0 0 0 0 0 0
A 52 4 0 0 0 7 49 0 51 0 0 0 0 1 0 0 0 0 0 0 0 0
A 53 10 0 0 50 7 48 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 24
A 54 10 0 0 53 7 48 10 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 26
A 55 10 0 0 54 7 48 13 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 13
A 58 2 0 0 0 7 666 0 0 0 58 0 0 0 0 0 0 0 0 0 0 0
A 59 1 0 3 0 94 667 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 60 10 0 0 0 7 59 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 17
A 61 10 0 0 60 7 59 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 20
A 62 4 0 0 0 7 61 0 11 0 0 0 0 2 0 0 0 0 0 0 0 0
A 63 4 0 0 0 7 60 0 62 0 0 0 0 1 0 0 0 0 0 0 0 0
A 64 2 0 0 0 7 671 0 0 0 64 0 0 0 0 0 0 0 0 0 0 0
A 65 10 0 0 61 7 59 16 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 64
A 66 2 0 0 0 7 672 0 0 0 66 0 0 0 0 0 0 0 0 0 0 0
A 67 10 0 0 65 7 59 19 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 66
A 68 4 0 0 0 7 67 0 11 0 0 0 0 2 0 0 0 0 0 0 0 0
A 69 4 0 0 0 7 65 0 68 0 0 0 0 1 0 0 0 0 0 0 0 0
A 70 10 0 0 67 7 59 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 24
A 71 2 0 0 0 7 673 0 0 0 71 0 0 0 0 0 0 0 0 0 0 0
A 72 10 0 0 70 7 59 22 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 71
A 73 10 0 0 72 7 59 10 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 26
A 74 10 0 0 73 7 59 13 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 13
A 77 1 0 3 0 100 677 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 78 10 0 0 0 7 77 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 17
A 79 10 0 0 78 7 77 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 20
A 80 4 0 0 0 7 79 0 11 0 0 0 0 2 0 0 0 0 0 0 0 0
A 81 4 0 0 0 7 78 0 80 0 0 0 0 1 0 0 0 0 0 0 0 0
A 82 10 0 0 79 7 77 16 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 64
A 83 10 0 0 82 7 77 19 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 66
A 84 4 0 0 0 7 83 0 11 0 0 0 0 2 0 0 0 0 0 0 0 0
A 85 4 0 0 0 7 82 0 84 0 0 0 0 1 0 0 0 0 0 0 0 0
A 86 10 0 0 83 7 77 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 24
A 87 10 0 0 86 7 77 22 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 71
A 88 10 0 0 87 7 77 10 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 26
A 89 10 0 0 88 7 77 13 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 13
A 92 1 0 3 0 106 685 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 93 10 0 0 0 7 92 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 17
A 94 10 0 0 93 7 92 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 20
A 95 4 0 0 0 7 94 0 11 0 0 0 0 2 0 0 0 0 0 0 0 0
A 96 4 0 0 0 7 93 0 95 0 0 0 0 1 0 0 0 0 0 0 0 0
A 97 10 0 0 94 7 92 16 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 64
A 98 10 0 0 97 7 92 19 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 66
A 99 4 0 0 0 7 98 0 11 0 0 0 0 2 0 0 0 0 0 0 0 0
A 100 4 0 0 0 7 97 0 99 0 0 0 0 1 0 0 0 0 0 0 0 0
A 101 10 0 0 98 7 92 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 24
A 102 10 0 0 101 7 92 22 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 71
A 103 10 0 0 102 7 92 10 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 26
A 104 10 0 0 103 7 92 13 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 13
A 106 2 0 0 0 6 693 0 0 0 106 0 0 0 0 0 0 0 0 0 0 0
A 107 2 0 0 0 7 696 0 0 0 107 0 0 0 0 0 0 0 0 0 0 0
A 108 2 0 0 0 19 711 0 0 0 108 0 0 0 0 0 0 0 0 0 0 0
A 109 1 0 0 0 19 710 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 111 2 0 0 0 6 714 0 0 0 111 0 0 0 0 0 0 0 0 0 0 0
A 112 2 0 0 0 137 715 0 0 0 112 0 0 0 0 0 0 0 0 0 0 0
A 113 1 0 0 0 109 713 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
Z
T 630 61 0 0 0 0
A 637 7 114 0 1 2 1
A 636 7 0 16 1 10 1
A 648 7 116 0 1 2 1
A 647 7 0 16 1 10 1
A 654 7 118 0 1 2 1
A 653 7 0 16 1 10 1
A 660 7 120 0 1 2 1
A 659 7 0 16 1 10 1
A 668 7 122 0 1 2 1
A 667 7 0 58 1 10 1
A 678 7 124 0 1 2 1
A 677 7 0 58 1 10 1
A 686 7 126 0 1 2 1
A 685 7 0 58 1 10 0
Z
