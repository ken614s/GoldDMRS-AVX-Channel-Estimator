#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <immintrin.h>

#define LEN     144
#define REPEAT 1000

typedef struct {
    float re, im; } cpx;

// Đọc vector phức từ file 
void read_complex(const char *filename, cpx *vec, int len) {
    FILE *f = fopen(filename, "r");
    if (!f) {
        printf("Lỗi mở file %s\n", filename);
        exit(1);
    }
    for (int i = 0; i < len; ++i)
        fscanf(f, "%f %f", &vec[i].re, &vec[i].im);
    fclose(f);
}

// Ghi vector phức ra file
void write_complex(const char *filename, cpx *vec, int len) {
    FILE *f = fopen(filename, "w");
    if (!f) {
        printf("Lỗi tạo file %s\n", filename);
        exit(1);
    }
    for (int i = 0; i < len; ++i) {
        fprintf(f, "%.4f\n", vec[i].re);
        fprintf(f, "%.4f\n", vec[i].im);
    }
    fclose(f);
}

// Hàm tính H_est: H = r_hat * conj(r) 
void compute_H_float(const cpx *r_hat, const cpx *r, cpx *H_out, int len) {
    for (int k = 0; k < len; ++k) {
        H_out[k].re =  r_hat[k].re * r[k].re + r_hat[k].im * r[k].im;
        H_out[k].im =  r_hat[k].im * r[k].re - r_hat[k].re * r[k].im;
    }
}

// Mảng unpack một lần cho AVX, 32-byte aligned
static float _r_re[LEN]   __attribute__((aligned(32)));
static float _r_im[LEN]   __attribute__((aligned(32)));
static float _rh_re[LEN]  __attribute__((aligned(32)));
static float _rh_im[LEN]  __attribute__((aligned(32)));

// Compute H bằng AVX 
void compute_H_avx(const cpx * /*r_hat*/, const cpx * /*r*/, cpx *H_out, int len) {
    static float Hre[LEN] __attribute__((aligned(32)));
    static float Him[LEN] __attribute__((aligned(32)));

    for (int i = 0; i < len; i += 8) {
        __m256 a = _mm256_load_ps(&_rh_re[i]);
        __m256 b = _mm256_load_ps(&_rh_im[i]);
        __m256 c = _mm256_load_ps(&_r_re [i]);
        __m256 d = _mm256_load_ps(&_r_im [i]);

        __m256 re = _mm256_add_ps(_mm256_mul_ps(a,c),
                                  _mm256_mul_ps(b,d));
        __m256 im = _mm256_sub_ps(_mm256_mul_ps(b,c),
                                  _mm256_mul_ps(a,d));

        _mm256_store_ps(&Hre[i], re);
        _mm256_store_ps(&Him[i], im);
    }
    for (int i = 0; i < len; ++i) {
        H_out[i].re = Hre[i];
        H_out[i].im = Him[i];
    }
}

// Timing helper: trả về ms giữa hai timespec
static inline double diff_ms(const struct timespec *t1,
                             const struct timespec *t2) {
    return (t2->tv_sec  - t1->tv_sec ) * 1e3
         + (t2->tv_nsec - t1->tv_nsec) / 1e6;
}

int main() {
    cpx r[LEN], r_hat[LEN];
    cpx H_float[LEN], H_avx[LEN];
    struct timespec t0, t1, t2;
    double min_f = 1e9, max_f = 0, sum_f = 0;
    double min_a = 1e9, max_a = 0, sum_a = 0;
    volatile float sink = 0;

    // 1) Đọc dữ liệu 
    read_complex("dmrs_seq.txt",      r,     LEN);
    read_complex("received_dmrs.txt", r_hat, LEN);

    // 2) Tính H bằng float và ghi ra file
    compute_H_float(r_hat, r, H_float, LEN);
    write_complex("h_est_float.txt", H_float, LEN);

    // 3) Unpack 4 mảng float 1 lần cho AVX
    for (int i = 0; i < LEN; ++i) {
        _r_re [i] = r[i].re;
        _r_im [i] = r[i].im;
        _rh_re[i] = r_hat[i].re;
        _rh_im[i] = r_hat[i].im;
    }

    // 4) Benchmark float 
    for (int i = 0; i < REPEAT; ++i) {
        clock_gettime(CLOCK_MONOTONIC, &t1);
        compute_H_float(r_hat, r, H_float, LEN);
        clock_gettime(CLOCK_MONOTONIC, &t2);
        double dt = diff_ms(&t1, &t2);
        sum_f += dt; 
	if (dt < min_f) min_f = dt; 
	if (dt > max_f) max_f = dt;
        sink += H_float[0].re + H_float[0].im;
    }

    // 5) Benchmark AVX 
    for (int i = 0; i < REPEAT; ++i) {
        clock_gettime(CLOCK_MONOTONIC, &t1);
        compute_H_avx(r_hat, r, H_avx, LEN);
        clock_gettime(CLOCK_MONOTONIC, &t2);
        double dt = diff_ms(&t1, &t2);
        sum_a += dt; if (dt < min_a) min_a = dt; if (dt > max_a) max_a = dt;
        sink += H_avx[0].re + H_avx[0].im;
    }

    // 6) In kết quả benchmark (ms)
    printf("FLOAT: min=%.6f ms, max=%.6f ms, avg=%.6f ms\n",
           min_f, max_f, sum_f/REPEAT);
    printf("AVX  : min=%.6f ms, max=%.6f ms, avg=%.6f ms\n",
           min_a, max_a, sum_a/REPEAT);

    // 7) Ghi kết quả AVX 
    write_complex("h_est_avx.txt", H_avx, LEN);

    return 0;
}
