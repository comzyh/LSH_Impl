# LSH_Impl

实现一个 LSH 查找

# 编译

> g++ test.cpp -std=c++11 -fopenmp -O2 -o test

> ./test --train sift/sift_base.fvecs --test sift/sift_query.fvecs 

# 使用

对于SIFT1M 数据集，建议使用

> ./test -b sift/sift_base.fvecs -q sift/sift_query.fvecs -g sift/sift_groundtruth.ivecs   -w 250 -p 300 -t 64 -f 8

Recall 约为 0.966

如果内存比较受限，可以尝试降低哈希表个数，提升probe

> ./test -b sift/sift_base.fvecs -q sift/sift_query.fvecs -g sift/sift_groundtruth.ivecs   -w 300 -p 500 -t 16 -f 8 -l 100

Recall 约为 0.940