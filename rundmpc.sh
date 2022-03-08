#!/usr/bin/env bash
#Number of simulations:
r=1
#start Number:
s=6

srce=experiments_small_flock/src1
dest=experiments_small_flock/exp7
f=C\ Code

# rm fitOverHorizon.txt
mkdir $dest
mkdir $srce 
mkdir $dest/result_files
cp C\ Code/gmpc.conf $dest/result_files/
# cp rectangles.txt $dest/result_files/

# g++ -x c "$f"/genconf.c "$f"/conf.c "$f"/common.c "$f"/ziggurat.c -O3 -o genconf
# ./genconf -x -5 -X 5 -y -5 -Y 5 -u -1 -U 1 -v -1 -V 1 -b $r $srce/init_conf_%d.txt

# gcc "$f"/dmpc.c "$f"/reynolds.c "$f"/conf.c "$f"/common.c "$f"/ziggurat.c -O3 -pthread -o dmpc
# ./dmpc -b $r -s $s $srce/init_conf_%d.txt $dest/dmpc_%d_%s.txt -q

g++ -x c "$f"/gmpc.c "$f"/conf.c "$f"/common.c "$f"/ziggurat.c -O3 -pthread -o gmpc
./gmpc -b $r -s $s $srce/init_conf_%d.txt $dest/gmpc_%d_%s.txt

# gcc "$f"/ziggurat.c "$f"/reynolds.c "$f"/conf.c "$f"/common.c -O3 -pthread -o reynolds
# ./reynolds -b $r $srce/init_conf_%d.txt $dest/reynolds_%d_%s.txt

chmod -R 777 experiments_small_flock/*
#-----------------------------------------------------------------
# Start
# dest=experiments_ml/exp29
# mkdir $dest
# mkdir $dest/result_files


# rm rectangles.txt
# cp neurips_results/obstacles/1/rectangles.txt .

# cp C\ Code/gmpc.conf $dest/result_files/
# cp rectangles.txt $dest/result_files/

# ./gmpc -b $r -s $s $srce/init_conf_%d.txt $dest/gmpc_%d_%s.txt

# #-----------------------------------------------------------------
# # Start
# dest=experiments_ml/exp30
# mkdir $dest
# mkdir $dest/result_files


# rm rectangles.txt
# cp neurips_results/obstacles/2/rectangles.txt .

# cp C\ Code/gmpc.conf $dest/result_files/
# cp rectangles.txt $dest/result_files/

# ./gmpc -b $r -s $s $srce/init_conf_%d.txt $dest/gmpc_%d_%s.txt

# #-----------------------------------------------------------------
# # Start
# dest=experiments_ml/exp31
# mkdir $dest
# mkdir $dest/result_files


# rm rectangles.txt
# cp neurips_results/obstacles/3/rectangles.txt .

# cp C\ Code/gmpc.conf $dest/result_files/
# cp rectangles.txt $dest/result_files/

# ./gmpc -b $r -s $s $srce/init_conf_%d.txt $dest/gmpc_%d_%s.txt

# #-----------------------------------------------------------------
# # Start
# dest=experiments_ml/exp32
# mkdir $dest
# mkdir $dest/result_files


# rm rectangles.txt
# cp neurips_results/obstacles/4/rectangles.txt .

# cp C\ Code/gmpc.conf $dest/result_files/
# cp rectangles.txt $dest/result_files/

# ./gmpc -b $r -s $s $srce/init_conf_%d.txt $dest/gmpc_%d_%s.txt

# #-----------------------------------------------------------------

