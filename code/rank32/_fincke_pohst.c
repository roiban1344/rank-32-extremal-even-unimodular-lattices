/* Fincke-Pohst enumeration: the vectors of norm <= bound of an integral Gram matrix G (order n), one of each pair +-v, written to out.
   There is no limit on the number of nodes of the enumeration; if more than cap vectors are found it returns -1.
   Returns the number of vectors and stores the number of nodes in *nodes. Used by rank32/gaussian_keys.py (enum_short).
   Build: gcc -O3 -march=native -shared -fPIC -o <out>.so _fincke_pohst.c -lm (done by gaussian_keys on first use).

   Fincke-Pohst の列挙。整数の Gram 行列 G（n 次）のノルム <= bound のベクトルを +-v の片方ずつ out に書く。
   列挙の節点数の上限は持たず、本数が cap を超えたら -1 を返す。戻り値は本数で、*nodes に節点数を返す。
   rank32/gaussian_keys.py（enum_short）から使う。
   ビルド：gcc -O3 -march=native -shared -fPIC -o <出力>.so _fincke_pohst.c -lm（gaussian_keys が初めて使うときに行う）。 */
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#define NMAX 64
long long fp_enum(int n, const int64_t *G, int bound, int64_t *out, long long cap, long long *nodes) {
    double q[NMAX][NMAX], d[NMAX];
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) q[i][j] = (double)G[i * n + j];
    /* eliminate from below: q(x) = sum_i d_i (x_i + sum_{j > i} q_ij x_j)^2 / 下から消す */
    for (int i = 0; i < n; i++) {
        d[i] = q[i][i];
        if (d[i] <= 0) return -2;
        for (int j = i + 1; j < n; j++) q[i][j] /= d[i];
        for (int k = i + 1; k < n; k++) for (int l = k; l < n; l++) q[k][l] -= q[i][k] * q[i][l] * d[i];
    }
    double rem[NMAX + 1], cen[NMAX];
    int x[NMAX], hi[NMAX];
    long long count = 0, nn = 0;
    const double eps = 1e-6;
    int i = n - 1;
    rem[n] = bound + eps;
    /* sign: the last nonzero coordinate is positive; while all above are 0, start from x_i >= 0 / 符号：最後の 0 でない座標を正にする */
    int top_zero = 1; /* whether all coordinates above i are 0 / i より上がすべて 0 か */
    double s, r;
    cen[i] = 0;
    r = sqrt(rem[n] / d[i]);
    x[i] = 0; hi[i] = (int)floor(r);
    int zero_above[NMAX + 1];
    zero_above[n] = 1;
    x[i] = 0;
    for (;;) {
        if (x[i] > hi[i]) {
            if (++i == n) break;
            x[i]++;
            continue;
        }
        double t = x[i] + cen[i];
        rem[i] = rem[i + 1] - d[i] * t * t;
        nn++;
        if (rem[i] < 0) { x[i]++; continue; }
        zero_above[i] = zero_above[i + 1] && x[i] == 0;
        if (i == 0) {
            if (!zero_above[0]) {
                if (count >= cap) return -1;
                for (int j = 0; j < n; j++) out[count * n + j] = x[j];
                count++;
            }
            x[i]++;
            continue;
        }
        i--;
        s = 0;
        for (int j = i + 1; j < n; j++) s += q[i][j] * x[j];
        cen[i] = s;
        r = sqrt(rem[i + 1] / d[i]);
        hi[i] = (int)floor(r - s);
        x[i] = zero_above[i + 1] ? 0 : (int)ceil(-r - s);
    }
    (void)top_zero;
    *nodes = nn;
    return count;
}
