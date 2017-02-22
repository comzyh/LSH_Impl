#ifndef LSH_IMPL_HPP
#define LSH_IMPL_HPP

#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <fstream>

const int M_function_num = LSH_IMPL_FUNCTION_NUM;

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
    Matrix(size_t row, size_t col): row(row), col(col) {
        data = new ElementType[col * row];
    }
    ~Matrix() {
        // delete[] data;
    }
    const ElementType *operator[](const size_t index) const{
        return data + index * col;
    }
    ElementType *operator[](const size_t index) {
        return data + index * col;
    }
    static Matrix<ElementType> load_from_file(const char* filename) {
      size_t filesize, dimension;
      std::ifstream file(filename, std::ifstream::ate | std::ifstream::binary);
    	filesize = file.tellg();
    	file.seekg(0);
    	file.read(reinterpret_cast<char *>(&dimension), 4);
    	lsh_impl::Matrix<float> mat(filesize / (4 + dimension * sizeof(ElementType)), dimension);
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

    const ElementType &operator[](const size_t index) const {
        return data[index];
    }

    ElementType &operator[](const size_t index) {
        return data[index];
    }

    float operator *(const ElementType &v) const {
        float result = 0;
        for (size_t i = 0; i < d; i ++) {
            result += (data[i] - v[i]) * (data[i] - v[i]);
        }
        return result;
    }
};

class LSH_Index
{
public:
    LSH_Index();
    ~LSH_Index();

};

}
#endif
