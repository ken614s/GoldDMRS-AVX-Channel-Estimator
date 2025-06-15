#include <stdio.h>

#define LEN     144
#define REPEAT 1000

typedef struct {
    float re, im; } cpx;

// Đọc vector phức từ file
void read_complex(const char *filename, cpx *vec, int len) {
    FILE *f = fopen(filename, "r");
    if (!f) { 
        printf("Lỗi mở file %s\n", filename);
        return;
    }
    for (int i = 0; i < len; ++i)
        fscanf(f, "%f %f", &vec[i].re, &vec[i].im);
    fclose(f);
}

void write_complex(const char *filename, cpx *vec, int len) {
    FILE *f = fopen(filename, "w");
    if (!f) {
        printf("Lỗi tạo file %s\n", filename);
        return;
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

int main() {
    cpx r[LEN], r_hat[LEN], H_est[LEN], H_float[LEN];

    // 1) Đọc r(k) và r̂(k)
    read_complex("dmrs_seq.txt",  r,     LEN);
    read_complex("received_dmrs.txt",  r_hat, LEN);

    // 2) Tính H_est và ghi ra file
    compute_H_float(r_hat, r, H_est, LEN);
    write_complex("h_est.txt", H_est, LEN);
       return 0;
}

