This package is gBoost with wildcards with cashing.

need to use:
c++ library, boost
-->install on Ubuntu
     sudo apt-get install libboost-dev

How to make:
$ make
$ cd eval_wild
$ make
$ cd ..

How to use:
$./lpboost -m (minsup) -w (wildcard num) traindata
$./eval_wild  model testdata 
