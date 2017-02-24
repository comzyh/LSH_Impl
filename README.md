# LSH_Impl

实现一个 LSH 查找算法

# clone
> git clone git@github.com:comzyh/LSH_Impl.git --recursive

注意 `--recursive` 参数，也不要使用zip 下载，否则会缺失 cmdline 和 flann 库

# 编译

> g++ test.cpp vendor/flann/src/cpp/flann/ext/lz4.c vendor/flann/src/cpp/flann/ext/lz4hc.c -std=c++11 -fopenmp  -O2 -I vendor/flann/src/cpp -o test

# 参数列表

使用[cmdline](https://github.com/tanakh/cmdline/) 进行参数解析
```
usage: ./test --base=string --query=string --ground=string [options] ...
options:
  -b, --base            Base data path. (string)
  -q, --query           Query data path. (string)
  -g, --ground          Ground_Truth data path. (string)
  -t, --table_num       Number of table (int [=16])
  -f, --function_num    Number of function (int [=3])
  -w, --W               W (float [=100])
  -p, --probe           Number of probe (int [=200])
  -l, --limit           imit the query number (int [=-1])
  -k, --kdtree          Using kdtree in flann (int [=0])
  -c, --checks          Checks using by kd-tree (int [=128])
  -?, --help            print this message
```
当参数 `kdtree` 不为 0 时使用 flann 的 KDTreeIndex 进行索引，kd-tree 数量由`kdtree`参数指定。
当使用`kdtree` 参数时，还可以制定 `checks` 参数, 确定使用KDTreeIndex 搜索时访问叶子节点的数量

# 使用

对于SIFT1M 数据集，可以尝试使用

> ./test -b sift/sift_base.fvecs -q sift/sift_query.fvecs -g sift/sift_groundtruth.ivecs   -w 250 -p 300 -t 64 -f 8

Recall 约为 0.966

如果内存比较受限，可以尝试降低哈希表个数，提升probe

> ./test -b sift/sift_base.fvecs -q sift/sift_query.fvecs -g sift/sift_groundtruth.ivecs   -w 250  -p 2500 -t 16 -f 8

Recall 约为 0.960

GIST 数据集同理

> ./test -b gist/gist_base.fvecs -q gist/gist_query.fvecs -g gist/gist_groundtruth.ivecs -w 2 -p 2000 -t 64  -f 14

Recall 约为 0.93

Note: 使用该参数在64bit系统下大约需要消耗 9GB 内存



