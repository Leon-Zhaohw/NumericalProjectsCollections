This is development code. I tried my best to clean it up.

Please check our website for more clean and updated code.

Please also feel free to contact me (Abdalla G. M. Ahmed) by email (abdalla_gafar@hotmail.com) if you have any questions.

The included files should demonstrate the basic functionality of our sampling framework:

= th.cpp:
- compile this to generate strings of the Thue-Morse word.
- There are more efficient ways to do so, but in a range of 192 bits efficiency is irrelevant :)

= optimize.cpp:
- This file includes optimization code to generate self-similar blue-noise tiles, using Lloyed, FPO, and PPO.
- try, for example:
    q=lfc;./optimize -n 300 -p /tmp/test/ $(./th 7 | cut -c1-48) 15 -q $q -d 0.9 -r 0.63 -v 0.03 -T 1 -S4 -c1 -C51
- You will get a table file called "table.dat". A sample table is included in this folder.

= sequence.cpp
- uses the mentioned table to generate sequences of blue noise samples; e.g.
    ./sequence table.dat 1024 > test.txt

= progressive.cpp
- uses the self-similar tile set for stippling; e.g.
    s=188.5; ./progressive density.pgm table.dat -s $s -i 0 > test.eps && epstopdf test.eps
- The input image is pgm format.
- A sample image is included.

= post-optimize.cpp
- Post optimizes the self-similar set for each rank, e.g.
    n=100; P=abl; ./post-optimize -s 0.5 -n $n -d 0.8 -r 0.65 table.dat -n $n -P $P
- You will get a fat table called post-optimized.dat

= stipple.cpp
- Uses the post-optimized tile set for stippling; e.g:
    s=188.5; ./stipple density.pgm post-optimized.dat -s $s -i 0 > test.eps && epstopdf test.eps

