#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

#define LEN 288
#define NC 1600
#define L 31
#define DMRS_LEN 144

int main() {
    uint32_t c_init = 10062025;
    uint8_t x1[LEN + NC + L] = {0};
    uint8_t x2[LEN + NC + L] = {0};
    uint8_t c_seq[LEN] = {0};

    // Khởi tạo x1: x1(0) = 1, x1(1..30) = 0
    x1[0] = 1;
    for (int i = 1; i < L; ++i)
        x1[i] = 0;

    // Khởi tạo x2 từ c_init (bit thấp nhất là x2[0])
    for (int i = 0; i < L; ++i)
        x2[i] = (c_init >> i) & 1;

    // Sinh x1 theo công thức x1(n+31) = x1(n+3) + x1(n)
    for (int i = L; i < LEN + NC; ++i)
        x1[i] = (x1[i-28] + x1[i-31]) % 2;

    // Sinh x2 theo công thức x2(n+31) = x2(n+3) + x2(n+2) + x2(n+1) + x2(n)
    for (int i = L; i < LEN + NC; ++i)
        x2[i] = (x2[i-28] + x2[i-29] + x2[i-30] + x2[i-31]) % 2;

    // Sinh chuỗi Gold: e(n) = (x1(n+NC) + x2(n+NC)) mod 2
    for (int i = 0; i < LEN; ++i)
        c_seq[i] = (x1[NC + i] + x2[NC + i]) % 2;

    // ---- Sinh DMRS ----
    float r_re[DMRS_LEN], r_im[DMRS_LEN];
    float inv_sqrt2 = 1.0f / sqrtf(2.0f);

    for (int n = 0; n < DMRS_LEN; ++n) {
        r_re[n] = inv_sqrt2 * (1 - 2 * c_seq[2 * n]);
        r_im[n] = inv_sqrt2 * (1 - 2 * c_seq[2 * n + 1]);
    }

    // Xuất ra file dmrs_seq.txt, mỗi dòng: thực ảo (%.4f)
    FILE *f = fopen("dmrs_seq.txt", "w");
    if (f == NULL) {
        printf("Không thể tạo file dmrs_seq.txt!\n");
        return 1;
    }
    for (int n = 0; n < DMRS_LEN; ++n) {
    fprintf(f, "%.4f\n", r_re[n]);  // Ghi phần thực
    fprintf(f, "%.4f\n", r_im[n]);  // Ghi phần ảo
}
    fclose(f);

}
