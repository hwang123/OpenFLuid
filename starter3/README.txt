Build Instructions (Athena with Makefile)
==================
$ mkdir build
$ cd build
$ cmake ..
$ make

The sample solution binaries are in sample_solution/

Submission
==========
We only provide official support for developing under the Athena cluster.
Your submission has to build and run over there to receive credit.

Please submit your entire source directory (excluding the build
directory) and compiled binary in inst/.




1. As per the instructions in the handouts:
 	cd /build 
 	./a3 r .01

2. --
3. TAs and office hours
4. My implementation assumes that the cloth has equal dimensions (is a square). Not even sure if this is considered a problem, but for what it's worth I'm pretty positive the issue is a i+1 indexing error in one of my helper functions because the only difference when the cloth isn't a square seems to be a missing spring connected to the fixed points.
5. --
6. Liked this one a lot