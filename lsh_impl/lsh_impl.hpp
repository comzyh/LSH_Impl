#ifndef LSH_IMPL_HPP
#define LSH_IMPL_HPP

#include <cmath>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <functional>
namespace lsh_impl {


//http://www.cnblogs.com/yeahgis/archive/2012/07/13/2590485.html
//Generate std normal distribution
double gaussrand()
{
    static double V1, V2, S;
    static int phase = 0;
    double X;

    if ( phase == 0 ) {
        do {
            double U1 = (double)rand() / RAND_MAX;
            double U2 = (double)rand() / RAND_MAX;

            V1 = 2 * U1 - 1;
            V2 = 2 * U2 - 1;
            S = V1 * V1 + V2 * V2;
        } while (S >= 1 || S == 0);

        X = V1 * sqrt(-2 * log(S) / S);
    } else
        X = V2 * sqrt(-2 * log(S) / S);

    phase = 1 - phase;

    return X;
}

template<typename ElementType>
struct Matrix {
    size_t row, col;
    ElementType *data;
    Matrix(){}
    Matrix(size_t row, size_t col): row(row), col(col) {
        data = new ElementType[col * row];
    }
    ~Matrix() {
        // delete[] data;
    }

    ElementType *operator[](const size_t index) {
        return data + index * col;
    }

    const ElementType *operator[](const size_t index) const {
        return data + index * col;
    }



    static Matrix<ElementType> load_from_file(const char* filename) {
        uint32_t filesize, dimension;
        std::ifstream file(filename, std::ifstream::ate | std::ifstream::binary);
        filesize = file.tellg();
        file.seekg(0);
        file.read(reinterpret_cast<char *>(&dimension), 4);
        lsh_impl::Matrix<ElementType> mat(filesize / (4 + dimension * sizeof(ElementType)), dimension);
        file.seekg(0);
        for (size_t i = 0; i < mat.row; i ++)
        {
            file.ignore(4);
            file.read(reinterpret_cast<char *>(mat.data + mat.col * i), sizeof(ElementType) * mat.col);
        }
        return mat;
    }
};

template<typename ElementType>
struct Vector
{
    size_t d; //dimension
    ElementType *data;
    Vector(){}
    Vector(size_t d): d(d) {
        data = new ElementType[d];
    }
    Vector(size_t d, ElementType *data): d(d), data(data){}

    const ElementType &operator[](const size_t index) const {
        return data[index];
    }

    ElementType &operator[](const size_t index) {
        return data[index];
    }

    float operator *(const ElementType *v) const {
        float result = 0;
        for (size_t i = 0; i < d; i ++) {
            result += (data[i] - v[i]) * (data[i] - v[i]);
        }
        return result;
    }
    float operator*(const Vector<ElementType> &v) const {
        return (*this) * v.data;
    }
};

struct LshIndexParams
{
    size_t table_num; // L
    size_t function_num; // M
    float W; // W
    size_t probe_num; // number proble to use
    LshIndexParams(size_t table_num, size_t function_num, float W, size_t probe_num):
    table_num(table_num), function_num(function_num), W(W), probe_num(probe_num) {}
};

template<typename ElementType>
struct LSH_Table
{
    typedef int32_t hash_element_type; //
    typedef hash_element_type* hash_type; // 
    typedef std::unordered_multimap<
        hash_type, // Key
        size_t, // Value
        std::function<size_t(hash_type)>, // hasher
        std::function<bool(hash_type, hash_type)> //key_equal
        > Table_T;
    Matrix<ElementType> a;
    Vector<ElementType> b;
    Matrix<hash_element_type> hashs;
    size_t feature_size;
    size_t function_num;
    float W;
    Table_T table;
    LSH_Table(size_t feature_size, size_t function_num, float W, size_t row):
    feature_size(feature_size), function_num(function_num), W(W),
    a(function_num, feature_size),b(function_num), hashs(row, function_num) {
        table = Table_T(
            row, //
            [function_num](hash_type hash)-> size_t { // Hasher
                size_t result = 0;
                for (size_t i = 0; i < function_num; i++) {
                    result = result * 47 + *reinterpret_cast<uint32_t*>(hash + i);
                }
                return result;
            },
            [function_num](hash_type u, hash_type v)-> bool {
                for (size_t i = 0; i < function_num; i++) {
                    if (u[i] != v[i]) {
                        return false;
                    }
                }
                return true;
            }
        );
        for (size_t i = 0; i < function_num; i ++) {
            for (size_t j = 0; j < feature_size; j ++) {
                a[i][j] = gaussrand();
            }
            b[i] = rand() / RAND_MAX * W;

        }
    }
    LSH_Table (LSH_Table &&t) {
        std::swap(a, t.a);
        std::swap(b, t.b);
        std::swap(hashs, t.hashs);
        std::swap(feature_size, t.feature_size);
        std::swap(function_num, t.function_num);
        std::swap(W, t.W);
        std::swap(table, t.table);
        t.a.data = t.b.data = nullptr;
        t.hashs.data = nullptr;
    }

    void getKey(const Vector<ElementType> &v, hash_element_type *result) const {
        for (size_t i = 0; i < function_num; i++) {
            result[i] = (v * a[i] + b[i]) / W;
        }
    }
    void add(const Matrix<ElementType> &data) {
        for (size_t i = 0; i < data.row; i++) {
            Vector<hash_element_type> hash(function_num, const_cast<hash_type>(hashs[i]));
            Vector<ElementType> v(feature_size, const_cast<ElementType*>(data[i]));
            getKey(v, hash.data);
            table.insert(std::make_pair(hash.data, i));
        }
    }
    ~LSH_Table() {
        delete[] a.data;
        delete[] b.data;
        delete[] hashs.data;
    }
};

template<typename ElementType>
class LSH_Index
{
public:
    LSH_Index(const Matrix<ElementType> &data, const LshIndexParams &params): data(data), params(params) {
        for (size_t i = 0; i < params.table_num; i++) { 
            //size_t feature_size, size_t function_num, float W, size_t rows
            tables.push_back(LSH_Table<ElementType> (data.col, params.function_num, params.W, data.row));
        }
    }
    void buildIndex() {
        for (size_t i = 0; i < params.table_num; i ++) {
            tables[i].add(data);
            printf("table %4lu initialized\n", i);
        }

    }
    ~LSH_Index(){

    };
private:
    Matrix<ElementType> data;
    LshIndexParams params;
    std::vector<LSH_Table<ElementType> > tables;
};

};
#endif
