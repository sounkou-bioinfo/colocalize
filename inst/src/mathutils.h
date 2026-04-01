/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef MATHUTILS_H
#define MATHUTILS_H

#define EIGEN_NO_DEBUG

#include "_global.h"
#include <algorithm>
#include "utils.h"
#include "exceptions.h"
#include <Eigen/Eigen>

#define INVERT_MIN_EVAL 1e-10
#define INVERT_EVAL_EST 0.01 

using namespace Eigen;

namespace ProductType {
  enum Type {Raw, Centered, Covariance, Correlation};
}

namespace ProcessType {
  enum Type {Raw, Centered, Normalized};
}

class MathUtils {
public:
  static ProcessType::Type type_convert(ProductType::Type type, bool scale=true) {
    switch (type) {
      case ProductType::Raw: return ProcessType::Raw;
      case ProductType::Centered: return ProcessType::Centered;
      case ProductType::Covariance: return ProcessType::Centered;
      case ProductType::Correlation: return scale ? ProcessType::Normalized : ProcessType::Raw; 
    }
    return ProcessType::Raw;                 
  }

  ///For symmetric matrices only; only uses lower triangle
  ///NOTE: returns improper inverse if not of full rank
  template<typename T>
  static long invert_matrix(Buffer<T>& mat, double cutoff=0, bool root=false) {return invert_matrix(mat, mat, cutoff, root);}

  template<typename T>
  static long invert_matrix(Buffer<T>& mat, Buffer<T>& target, double cutoff=0, bool root=false) {try {
    if (mat.nrow() == 0) throw MathException("matrix inversion", "trying to invert empty matrix");
    if (mat.nrow() != mat.ncol()) throw MathException("matrix inversion", "trying to invert non-square matrix");

    long dim = 0, N = mat.nrow();
    bool est_cutoff = cutoff < 0;
    if (est_cutoff) cutoff = -cutoff;
    if (cutoff < INVERT_MIN_EVAL) cutoff = INVERT_MIN_EVAL;

    if (N == 1) {
      if (mat(0) <= 0) return 0;
      if (&mat != &target) target.set_square(1);
      target(0) = root ? 1/sqrt(mat(0)) : 1/mat(0);
      return 1;
    }

    SelfAdjointEigenSolver<Matrix<T,Dynamic,Dynamic> > eig(Map<Matrix<T,Dynamic,Dynamic> >(mat[0], N, N));
    const T* val = eig.eigenvalues().data();

    if (est_cutoff) {double est = 0;
      for (int i = 0; i < N; i++) est += val[i];
      cutoff *= est / N;
    }

    T* weights = new T[N];    
    for (long i = 0; i < N; i++) {
      if (val[i] > cutoff) {weights[i] = root ? 1/sqrt(val[i]) : 1/val[i]; dim++;}
      else weights[i] = 0;
    }
    if (dim == 0) {target.set_size(0,0); return dim;}

    int diff = N-dim, ncol = root ? dim : N; target.set_size(N, ncol);
    Map<Matrix<T,Dynamic,Dynamic> > out_map(target[0], N, ncol);
    if (root) out_map.noalias() = eig.eigenvectors().block(0, diff, N, dim) * Map<Matrix<T,Dynamic,1> >(weights+diff, dim, 1).asDiagonal();
    else out_map.noalias() = eig.eigenvectors() * Map<Matrix<T,Dynamic,1> >(weights, N, 1).asDiagonal() * eig.eigenvectors().transpose();

    delete[] weights;
    return dim;
    } catch (const MathException& me) {throw;}
      catch (const exception& e) {check_mem_error(e, "inverting matrix"); throw MathException("matrix inversion", "an error occurred when trying to invert matrix", e.what());}
  }

  template<typename T>
  static long rotate(Buffer<T>& data, Buffer<T>& cov, double cutoff=0) {
    Buffer<T> buffer;
    long rank = rotate(data, buffer, cov, cutoff);
    data.swap(buffer);
    return rank;
  }

  ///Assumes data is centered; alters cov
  template<typename T>
  static long rotate(Buffer<T>& data, Buffer<T>& target, Buffer<T>& cov, double cutoff=0) {
    long N = data.nrow(), dim = cov.nrow();
    long rank = invert_matrix(cov, cutoff, true);

    target.set_size(N,rank);
    for (long i_pc = 0; i_pc < rank; i_pc++) {
      for (long i_var = 0; i_var < dim; i_var++) vec_product(data[i_var], target[i_pc], cov(i_var,i_pc), N, i_var>0);
    }

    return rank;
  }

  template<typename T, typename T2>
  static Buffer<T2>& cov_matrix(Buffer<T>& data, Buffer<T2>& mat, bool to_corr=false, int no_obs=0) {int no_var = data.ncol();
    if (no_obs == 0) no_obs = data.nrow();
    
    Buffer<double> stats(no_var, 2); double *means = stats[0], *sds = stats[1];
    mat.set_square(no_var, true);

    for (long i = 0; i < no_var; i++) {
      get_stats(data[i], no_obs, means[i], sds[i]);
      if (sds[i] <= 0) throw MathException("covariance matrix", "variable has variance of zero");
    }
    
    for (long i = 0; i < no_var; i++) {
      mat(i,i) = !to_corr ? sds[i]*sds[i] : 1;
      for (long j = 0; j < i; j++) {double &M1 = means[i], &M2 = means[j]; T2 &coeff = mat(i,j);
        for (T *curr = data[i], *curr2 = data[j], *end = curr+no_obs; curr < end; curr++, curr2++) coeff += (*curr - M1) * (*curr2 - M2);
        coeff /= to_corr ? (no_obs-1) * sds[i] * sds[j] : no_obs-1;
        mat(j,i) = coeff;
      }
    }
    return mat;
  }
  

  template<typename T>
  static Buffer<T>& cov_to_cor(Buffer<T>& mat) {
    if (mat.nrow() != mat.ncol()) throw MathException("matrix normalization", "covariance matrix is not square");

    long dim = mat.nrow();
    if (dim <= 1) {
      if (dim == 1) mat(0) = 1;
      return mat;
    }

    for (int i = 0; i < dim; i++) {
      for (int j = i+1; j < dim; j++) {
        mat(j,i) /= sqrt(mat(i,i)*mat(j,j));
        if (mat(j,i) > 1) mat(j,i) = 1;
        else if (mat(j,i) < -1) mat(j,i) = -1;
        mat(i,j) = mat(j,i);
      }
      mat(i,i) = 1;
    }
    return mat;
  }

  template<typename T>
  static void center(T* buff, long size) {
    double mean = get_mean(buff, size);
    for (T* end = buff + size; buff < end; buff++) *buff -= mean;
  }

  template<typename T>
  static Buffer<T>& center(Buffer<T>& buff) {
    for (int i = 0; i < buff.ncol(); i++) center(buff[i], buff.nrow());
    return buff;
  }

  template<typename T>
  static bool normalize(T* buff, long size, double scale=1) {
    double mean, sd; get_stats(buff, size, mean, sd);
    if (sd == 0) return false;
    if (scale > 0) sd /= scale;
    for (T* curr = buff; curr < buff+size; curr++) *curr = ((*curr) - mean)/sd;
    return true;
  }

  template<typename T>
  static bool normalize(T* buff, T& obs, long size) {
    double mean, sd; get_stats(buff, size, mean, sd);
    if (sd == 0) return false;
    obs = (obs - mean) / sd;
    for (T* curr = buff; curr < buff+size; curr++) *curr = ((*curr) - mean)/sd;
    return true;
  }
  
  template<typename T>
  static bool normalize(Buffer<T>& buff, int no_obs=0) {int failed = 0;
    if (no_obs == 0) no_obs = buff.nrow();
    for (int i = 0; i < buff.ncol(); i++) {
      bool success = normalize(buff[i], no_obs);
      if (!success) {Utils::set_zero(buff[i], no_obs); failed++;} 
    }
    return failed == 0;
  }

  template<typename T> 
  static T clamp(T value, double lower, double upper) {return (value > lower) ? ( (value < upper) ? value : upper ) : lower;}
  
  template<typename T> 
  static void clamp(T* buff, long size, double lower, double upper) {
    for (T* end = buff + size; buff < size; buff++) {
      if (*buff < lower) *buff = lower;
      else if (*buff > upper) *buff = upper;
    }
  }

  template<typename T>   
  static T clamp01(T value) {return clamp(value, 0, 1);}
  template<typename T>   
  static void clamp01(T* buff, long size) {return clamp(buff, size, 0, 1);}  

  template<typename T, typename D>
  static void get_stats(const vector<T>& buff, D& mean, D& sd, bool to_sd=true) {
    mean = 0; sd = 0;
    for (long i = 0; i < buff.size(); i++) {
      mean += buff[i];
      sd += buff[i] * buff[i];
    }
    mean /= buff.size();
    sd = (sd - mean * mean * buff.size()) / (buff.size()-1);    
    if (sd < 0) sd = 0;
    if (to_sd) sd = sqrt(sd);
  }

  template<typename T, typename D>
  static void get_stats(const T* buff, long size, D& mean, D& sd, bool to_sd=true) {
    mean = 0; sd = 0;
    for (const T* curr = buff; curr < buff+size; curr++) {
      mean += *curr;
      sd += *curr * *curr;
    }
    mean /= size;
    sd = (sd - mean * mean * size) / (size-1);
    if (sd < 0) sd = 0;
    if (to_sd) sd = sqrt(sd);
  }
  
  template<typename T, typename D>  
  static void get_stats_skip(const T* buff, long nsteps, int step_size, D& mean, D& sd, bool to_sd=true) {
    if (step_size <= 1) get_stats(buff, nsteps, mean, sd, to_sd);
  
    mean = 0; sd = 0; const T* end = buff + nsteps*step_size;
    for (const T* curr = buff; curr < end; curr+=step_size) {
      mean += *curr;
      sd += *curr * *curr;
    }
    mean /= nsteps;
    sd = (sd - mean * mean * nsteps) / (nsteps-1);
    if (sd < 0) sd = 0;
    if (to_sd) sd = sqrt(sd);
  }

  template<typename T>
  static double get_mean(const vector<T>& buff) {
    double mean = 0;
    for (long i = 0; i < buff.size(); i++) mean += buff[i];
    return mean/double(buff.size());
  }

  template<typename T>
  static double get_mean(const T* buff, long size) {
    double mean = 0;
    for (const T* curr = buff; curr < buff+size; curr++) mean += *curr;
    return mean/double(size);
  }

  template<typename T>
  static double get_sd(const vector<T>& buff) {
    double mean, sd;
    get_stats(buff, mean, sd);
    return sd;
  }
  
  template<typename T>
  static double get_sd(const T* buff, long size) {
    double mean, sd;
    get_stats(buff, size, mean, sd);
    return sd;
  }
  
  template<typename T>
  static double get_var(const T* buff, long size) {
    double mean, var;
    get_stats(buff, size, mean, var, false);
    return var;
  }
  
  template<typename T>
  static double get_sd(const T* buff, long size, double mean) {
    double var = get_var(buff, size, mean);
    return var > 0 ? sqrt(var) : 0;
  }

  template<typename T>
  static double get_var(const T* buff, long size, double mean) {
    double sum = 0, diff; const T* end = buff + size;  
    while (buff < end) {diff = double(*(buff++)) - mean; sum += diff*diff;}
    return sum / (size - 1);
  }

  template<typename T>
  static double get_corr(const T* buff_i, const T* buff_j, long size) {return get_cov(buff_i, buff_j, size, true);}
  
  template<typename T>
  static double get_cov(const T* buff_i, const T* buff_j, long size, bool correlation=false) {
    double mean_i = get_mean(buff_i, size), mean_j = get_mean(buff_j, size), cov = 0;
   
    for (int k = 0; k < size; k++) cov += (buff_i[k] - mean_i) * (buff_j[k] - mean_j);
    if (correlation) {    
      double var_i = 0, var_j = 0;
      for (int k = 0; k < size; k++) {
        var_i += (buff_i[k] - mean_i) * (buff_i[k] - mean_i);    
        var_j += (buff_j[k] - mean_j) * (buff_j[k] - mean_j);            
      }
      if (var_i <= 0 || var_j <= 0) return 0;
      return cov / sqrt(var_i*var_j);
    } else return cov / (size-1);
  }
  
  template<typename T> 
  static double get_mode(T* buff, long size, bool prob=false) {
    Buffer<T> data(buff, size);
    sort(data.begin(), data.end());
    long curr_count = 1, max_count = 1; T max_value = data(0);
    for (long i = 1; i < size; i++) {
      if (data(i) == data(i-1)) curr_count++;
      else {
        if (curr_count > max_count) {max_count = curr_count; max_value = data(i-1);}
        curr_count = 1;
      }                           
    }
    if (curr_count > max_count) {max_count = curr_count; max_value = data(size-1);}

    return prob ? max_count / double(size) : max_value;    
  }

  template<typename T> 
  static T get_median(T* data, long size, bool in_place) {
    if (in_place) {
      sort(data, data+size);
      return (size % 2 == 0) ? (data[size/2]+data[size/2-1])/2.0 : data[size/2];
    } else {
      Buffer<T> buff(data, size); 
      return get_median(buff.begin(), size, true);
    }
  }

  template<typename T> 
  static T get_max(const vector<T>& vec) {
    T curr_max = !vec.empty() ? vec[0] : 0;
    for (long i = 1; i < vec.size(); i++) curr_max = max(curr_max, vec[i]);
    return curr_max;
  }

  template<typename T> 
  static T get_max(T* buff, long size) {
    T max = *(buff++);
    for (T* end = buff + size - 1; buff < end; buff++) {if (*buff > max) max = *buff;}
    return max;
  }

  template<typename T> 
  static T get_min(T* buff, long size) {
    T min = *(buff++);
    for (T* end = buff + size - 1; buff < end; buff++) {if (*buff < min) min = *buff;}
    return min;
  }
  
  template<typename T>
  static double* simple_regression(T* x, T* y, int N) {
    double Sx = 0, Sxx = 0, Sxy = 0, Sy = 0;
    for (T* end = x + N; x < end; x++, y++) {Sx += *x; Sxx += *x * *x; Sxy += *x * *y; Sy += *y;}
 
    double det = N*Sxx - Sx*Sx;
    if (det / (N*N) < 1e-8) return 0;
    
    double* coeff = new double[2];
    coeff[0] = (Sxx*Sy - Sx*Sxy)/det;
    coeff[1] = (N*Sxy - Sx*Sy)/det;
    return coeff;
  }

  template<typename T>
  static long count(T* buff, long len) {long out = 0;
    for (T* end = buff+len; buff < end; buff++) {if (*buff) out++;}
    return out;
  }
 
  template<typename T>
  static long count(vector<T>& data) {long out = 0;
    for (long i = 0; i < data.size(); i++) {if (data[i]) out++;}
    return out;
  }

  template<typename T>
  static long count_value(T* buff, T value, long len) {long out = 0;
    for (T* end = buff+len; buff < end; buff++) {if (*buff == value) out++;}
    return out;
  }

 
  template<typename T>
  static double sum(T* buff, long len) {double out = 0;
    for (T* end = buff+len; buff < end; buff++) out += *buff;
    return out;
  }
  
  template<typename T>
  static double sum(T* buff, double scale, long len) {double out = 0;
    for (T* end = buff+len; buff < end; buff++) out += scale * *buff;
    return out;
  }
  
  template<typename T>
  static long sum(vector<T>& data) {T out = 0;
    for (long i = 0; i < data.size(); i++) out += data[i];
    return out;
  }

  template<typename T>
  static double product(T* buff, long len) {double out = 1;
    for (T* end = buff+len; buff < end; buff++) out *= *buff;
    return out;
  }
   
  template<typename T>
  static double in_product(T* buff, long len) {double out = 0;
    for (T* end = buff+len; buff < end; buff++) out += *buff * *buff;
    return out;
  }

  template<typename T>
  static double in_product(T* buff, long len, int step) {double out = 0;
    for (T* end = buff+len*step; buff < end; buff += step) out += *buff * *buff;
    return out;
  }

  template<typename T1, typename T2>
  static double in_product(T1* buff1, T2* buff2, long len) {double out = 0;
    for (T1* end = buff1+len; buff1 < end; buff1++, buff2++) out += *buff1 * *buff2;
    return out;
  }

  template<typename T1, typename T2>
  static double in_product(T1* buff1, T2* buff2, long len, int step1, int step2) {double out = 0;
    for (T1* end = buff1+len*step1; buff1 < end; buff1 += step1, buff2 += step2) out += *buff1 * *buff2;
    return out;
  }

  template<typename T>
  static void vec_sum(T* in, T* out, long len) {
    T* end = out + len;
    while (out < end) *(out++) += *(in++);
  }

  template<typename T>
  static void vec_sum(T* in1, T* in2, T* out, long len) {
    T* end = out + len;
    while (out < end) *(out++) = *(in1++) + *(in2++);
  }

  template<typename T>
  static void vec_sum(T* in, T* out, double scale, long len) {
    T* end = out + len;
    while (out < end) *(out++) += scale * *(in++);
  }
    
  template<typename T>
  static void vec_sum(T* in1, T* in2, T* out, double sc1, double sc2, long len) {
    T* end = out + len;
    while (out < end) *(out++) = sc1 * *(in1++) + sc2 * *(in2++);
  }

  template<typename T>
  static void vec_diff(T* in, T* out, long len) {
    T* end = out + len;
    while (out < end) *(out++) -= *(in++);
  }

  template<typename T>
  static void vec_diff(T* in1, T* in2, T* out, long len) {
    T* end = out + len;
    while (out < end) *(out++) = *(in1++) - *(in2++);
  }

  template<typename T>
  static void vec_shift(T* in, double shift, long len) {vec_shift(in, in, shift, len);}

  template<typename T>
  static void vec_shift(T* in, T* out, double shift, long len) {
    T* end = out + len;
    while (out < end) *(out++) = *(in++) + shift;
  }

  template<typename T>
  static void vec_scale(T* in, double scale, long len) {vec_scale(in, in, scale, len);}

  template<typename T>
  static void vec_scale(T* in, T* out, double scale, long len) {
    T* end = out + len;
    while (out < end) *(out++) = scale * *(in++);
  }
  
  template<typename T>
  static void vec_transform(T* in, double shift, double scale, long len, bool shift_first) {vec_transform(in, in, shift, scale, len, shift_first);}
  
  template<typename T>
  static void vec_transform(T* in, T* out, double shift, double scale, long len, bool shift_first) {
    if (shift_first) shift *= scale;   
    T* end = out + len;
    while (out < end) *(out++) = scale * *(in++) + shift;
  }
  

  template<typename T>
  static void vec_square(T* in, long len) {vec_square(in, in, len);}

  template<typename T>
  static void vec_square(T* in, T* out, long len) {
    T* end = out + len;
    while (out < end) {*(out++) = *in * *in; in++;}
  }

  template<typename T>
  static void vec_product(T* in, T* out, double scale, long len, bool add) {
    T* end = out + len;
    if (add) {while (out < end) *(out++) += scale * *(in++);}
    else {while (out < end) *(out++) = scale * *(in++);}
  }

  template<typename T>
  static void vec_product(T* in1, T* in2, T* out, long len, bool add) {
    T* end = out + len;
    if (add) {while (out < end) *(out++) += *(in1++) * *(in2++);}
    else {while (out < end) *(out++) = *(in1++) * *(in2++);}
  }

  ///Uses lower triangle only
  template<typename T>
  static T square_trace(Buffer<T>& M) {
    long dim = M.nrow(); T out = 0, add = 0;
    for (long i = 0; i < dim; i++) {
      out += M(i,i)*M(i,i);
      for (T *curr = M[i]+i+1, *end = M[i+1]; curr < end; curr++) add += *curr * *curr;
    }
    return out + 2*add;
  }

  ///Uses lower triangle only
  template<typename T>
  static T corr_trace(Buffer<T>& M) {
    long dim = M.nrow(); T out = 0;
    for (long i = 0; i < dim; i++) {
      T scale = M(i,i), *curr = M[i];
      for (long j = i+1; j < dim; j++) out += curr[j] * curr[j] / scale / M(j,j);
    }
    return 2*out + dim;
  }

  template<typename T>
  static T prod_trace(Buffer<T>& M) {
    return in_product(M[0], M.length());
  }

  ///Compute quadratic form t(V) * M * V, with M symmetric
  template<typename T>
  static T quad_form(Buffer<T>& M, Buffer<T>& V) {return quad_form(M[0], V[0], M.nrow());}
  template<typename T>
  static T quad_form(T* M, T* V, int size) {
    if (size <= 1) return V[0]*V[0]*M[0];

    Map<Matrix<T,Dynamic,1> > V_map(V, size, 1);
    Matrix<T,1,1> out =  V_map.transpose() * Map<Matrix<T,Dynamic,Dynamic> >(M, size, size) * V_map;

    return out(0,0);
  }

  ///Compute cross product t(V1) * M * V1, with M symmetric
  template<typename T>
  static T cross_prod(Buffer<T>& M, Buffer<T>& V1, Buffer<T> V2) {return cross_prod(M[0], V1[0], V2[0], M.nrow());}
  template<typename T>
  static T cross_prod(T* M, T* V1, T* V2, int size) {
    if (size <= 1) return V1[0]*V2[0]*M[0];

    Matrix<T,1,1> out = Map<Matrix<T,1,Dynamic> >(V1, 1, size) * Map<Matrix<T,Dynamic,Dynamic> >(M, size, size) * Map<Matrix<T,1,Dynamic> >(V2, size, 1);
    return out(0,0);
  }


  template<typename T> 
  static Buffer<T>& matrix_scale(Buffer<T>& M, double scale) {
    vec_scale(M[0], scale, M.nrow()*M.ncol());
    return M;  
  }
  
  template<typename T> 
  static Buffer<T>& matrix_square(Buffer<T>& M) {
    vec_square(M[0], M.nrow()*M.ncol());
    return M;  
  }
         
  enum Transpose {TransposeNone=0, TransposeFirst, TransposeSecond};

  template<typename T> ///Assumes out does not overlap M 
  static Buffer<T>& matrix_prod(Buffer<T>& M, Buffer<T>& out, Transpose transpose) {
    if (transpose == TransposeNone && M.nrow() != M.ncol()) {out.set_square(0); return out;}
    int dim = (transpose == TransposeFirst) ? M.ncol() : M.nrow();

    Map<Matrix<T,Dynamic,Dynamic> > in_map = M.matrix();
    out.set_square(dim);

    if (transpose) {
      if (transpose == TransposeFirst) out.matrix().noalias() = in_map.transpose() * in_map;
      else out.matrix().noalias() = in_map * in_map.transpose();
    } else out.matrix().noalias() = in_map * in_map;
    return out;
  }


  template<typename T> ///Assumes out does not overlap M1 or M2
  static Buffer<T>& matrix_prod(Buffer<T>& M1, Buffer<T>& M2, Buffer<T>& out, Transpose transpose=TransposeNone) {
    int nrow = (transpose == TransposeFirst) ? M1.ncol() : M1.nrow();
    int ncol = (transpose == TransposeSecond) ? M2.nrow() : M2.ncol();
    out.set_size(nrow, ncol);

    if (transpose) {
      if (transpose == TransposeFirst) out.matrix().noalias() = M1.matrix().transpose() * M2.matrix();
      else out.matrix().noalias() = M1.matrix() * M2.matrix().transpose();
    } else out.matrix().noalias() = M1.matrix() * M2.matrix();

    return out;
  }
  
  static double pow2(double value) {return value*value;}
  static double pow3(double value) {return value*value*value;}
  static double pow4(double value) {return value*value*value*value;}
  
  static double root3(double value) {static double fact = 1/double(3); return pow(value, fact);}
  static double root4(double value) {return pow(value, 0.25);}
  
};




#endif /**MATHUTILS_H*/
