package com.ndsec.random;

import org.apache.commons.math3.special.Erf;

import java.util.List;

public class util {
    public static final int MATRIX_FORWARD_ELIMINATION = 0;
    public static final int MATRIX_BACKWARD_ELIMINATION = 1;

    public util() {
    }

    public static int pow(int a, int n) {
        int sum = 1;
        for (int i = 0; i < n; i++) {
            sum *= a;
        }
        return sum;
    }

    public static int toInt(List<Boolean> eps, int start, int end) {
        if (end <= start) {
            return 0;
        }
        int len = end - start;
        int sum = 0;
        for (int i = start; i < end; i++) {
            len--;
            if (eps.get(i)) {
                sum += pow(2, len);
            }
        }
        return sum;
    }

    public static int toInt(int[] eps, int start, int end) {
        if (end <= start) return 0;
        int len = end - start;
        int sum = 0;
        for (int i = start; i < end; i++) {
            len--;
            if (eps[i] == 1) {
                sum += pow(2, len);
            }
        }
        return sum;
    }

    public static double psi2p(List<Boolean> eps, int m) {
        int n = eps.size();
        int i, j, k, powLen;
        double sum, numOfBlocks;
        if (m == 0 || m == -1) return 0.0;
        numOfBlocks = n;
        powLen = (int) Math.pow(2, m + 1) - 1;
        int[] P = new int[powLen];
        for (i = 1; i < powLen - 1; i++) {
            P[i] = 0;
        }
        for (i = 0; i < numOfBlocks; i++) {
            k = 1;
            for (j = 0; j < m; j++) {
                if (eps.get((i + j) % n)) {
                    k = 2 * k + 1;
                } else if (!eps.get((i + j) % n)) {
                    k *= 2;
                }
            }
            P[k - 1]++;
        }
        sum = 0.0;
        for (i = (int) Math.pow(2, m) - 1; i < (int) Math.pow(2, m + 1) - 1; i++) {
            sum += Math.pow(P[i], 2);
        }
        sum = (sum * Math.pow(2, m) / n) - n;
        return sum;
    }

    public static double psi2(List<Boolean> eps, int m) {
        int n = eps.size();
        if (m == 0 || m == -1) {
            return 0.0;
        }
        int m_sum = 1 << m;
        int[] eps1 = new int[n + m - 1];
        for (int i = 0; i < n; i++) {
            eps1[i] = eps.get(i) ? 1 : 0;
        }
        for (int i = n, j = 0; i < n + m - 1; i++, j++) {
            eps1[i] = eps.get(j) ? 1 : 0;
        }
        int[] tabs = new int[m_sum];
        for (int i = 0; i < m_sum; i++) {
            tabs[i] = 0;
        }
        int tab;
        for (int i = 0; i < n; i++) {
            tab = toInt(eps1, i, i + m);
            tabs[tab]++;
        }
        double sum = 0.0;
        for (int i = 0; i < m_sum; i++) {
            sum += tabs[i] * tabs[i];
        }
        sum *= m_sum;
        sum /= n;
        sum -= n;
        return sum;
    }

    public static void def_matrix(int M, int Q, int[][] m, int k, List<Boolean> eps) {
        int i, j;
        for (i = 0; i < M; i++)
            for (j = 0; j < Q; j++)
                m[i][j] = eps.get(k * (M * Q) + j + i * M) ? 1 : 0;
    }

    public static void perform_elementary_row_operations(int flag, int i, int M, int Q, int[][] A) {
        int j, k;
        if (flag == MATRIX_FORWARD_ELIMINATION) {
            for (j = i + 1; j < M; j++)
                if (A[j][i] == 1)
                    for (k = i; k < Q; k++)
                        A[j][k] = (A[j][k] + A[i][k]) % 2;
        } else {
            for (j = i - 1; j >= 0; j--)
                if (A[j][i] == 1)
                    for (k = 0; k < Q; k++)
                        A[j][k] = (A[j][k] + A[i][k]) % 2;
        }
    }

    public static int swap_rows(int i, int index, int Q, int[][] A) {
        int p;
        int temp;
        for (p = 0; p < Q; p++) {
            temp = A[i][p];
            A[i][p] = A[index][p];
            A[index][p] = temp;
        }
        return 1;
    }

    public static int find_unit_element_and_swap(int flag, int i, int M, int Q, int[][] A) {
        int index, row_op = 0;
        if (flag == MATRIX_FORWARD_ELIMINATION) {
            index = i + 1;
            while ((index < M) && (A[index][i] == 0))
                index++;
            if (index < M)
                row_op = swap_rows(i, index, Q, A);
        } else {
            index = i - 1;
            while ((index >= 0) && (A[index][i] == 0))
                index--;
            if (index >= 0)
                row_op = swap_rows(i, index, Q, A);
        }

        return row_op;
    }

    public static int determine_rank(int m, int M, int Q, int[][] A) {
        int i, j, rank, allZeroes;
        rank = m;
        for (i = 0; i < M; i++) {
            allZeroes = 1;
            for (j = 0; j < Q; j++) {
                if (A[i][j] == 1) {
                    allZeroes = 0;
                    break;
                }
            }
            if (allZeroes == 1)
                rank--;
        }

        return rank;
    }

    public static int computeRank(int M, int Q, int[][] matrix) {
        int i, rank, m = Math.min(M, Q);
        for (i = 0; i < m - 1; i++) {
            if (matrix[i][i] == 1) {
                perform_elementary_row_operations(MATRIX_FORWARD_ELIMINATION, i, M, Q, matrix);
            } else {
                if (find_unit_element_and_swap(MATRIX_FORWARD_ELIMINATION, i, M, Q, matrix) == 1) {
                    perform_elementary_row_operations(MATRIX_FORWARD_ELIMINATION, i, M, Q, matrix);
                }
            }
        }
        for (i = m - 1; i > 0; i--) {
            if (matrix[i][i] == 1)
                perform_elementary_row_operations(MATRIX_BACKWARD_ELIMINATION, i, M, Q, matrix);
            else {    /* matrix[i][i] = 0 */
                if (find_unit_element_and_swap(MATRIX_BACKWARD_ELIMINATION, i, M, Q, matrix) == 1)
                    perform_elementary_row_operations(MATRIX_BACKWARD_ELIMINATION, i, M, Q, matrix);
            }
        }
        rank = determine_rank(m, M, Q, matrix);
        return rank;
    }

    public static double cephes_normal(double x) {
        double arg, result, sqrt2 = 1.414213562373095048801688724209698078569672;

        if (x > 0) {
            arg = x / sqrt2;
            result = 0.5 * (1 + Erf.erf(arg));
        } else {
            arg = -x / sqrt2;
            result = 0.5 * (1 - Erf.erf(arg));
        }

        return (result);
    }

    public static void calF(double[] qValue, int s, int[] f) {
        for (int i = 0; i < s; ++i) {
            if (qValue[i] >= 0 && qValue[i] < 0.1) {
                f[1]++;
            }
            if (qValue[i] >= 0.1 && qValue[i] < 0.2) {
                f[2]++;
            }
            if (qValue[i] >= 0.2 && qValue[i] < 0.3) {
                f[3]++;
            }
            if (qValue[i] >= 0.3 && qValue[i] < 0.4) {
                f[4]++;
            }
            if (qValue[i] >= 0.4 && qValue[i] < 0.5) {
                f[5]++;
            }
            if (qValue[i] >= 0.5 && qValue[i] < 0.6) {
                f[6]++;
            }
            if (qValue[i] >= 0.6 && qValue[i] < 0.7) {
                f[7]++;
            }
            if (qValue[i] >= 0.7 && qValue[i] < 0.8) {
                f[8]++;
            }
            if (qValue[i] >= 0.8 && qValue[i] < 0.9) {
                f[9]++;
            }
            if (qValue[i] >= 0.9 && qValue[i] < 1.0) {
                f[10]++;
            }
        }
    }

    public static double decisionSample(double[] qValue, int s) {
        int[] f = new int[12];
        calF(qValue, s, f);
        double v = 0.0;
        for (int i = 1; i < 11; ++i) {
            v += Math.pow(f[i] - s / 10.0, 2) / (s / 10.0);
        }
        return cephes.igamc(4.5, v / 2);
    }

    public static double[] getFFT(double[] x) {
        int n = x.length;
        int k = 2;
        while (k < n) {
            k *= 2;
        }
        int m = k - n;
        double[] res;
        if (m != 0) {
            res = new double[k];
            System.arraycopy(x, 0, res, 0, n);
            for (int i = n; i < k; i++) {
                res[i] = 0.0;
            }
        } else {
            res = x;
        }
        return res;
    }

    public static int liner_res(boolean[] var0, int var1) {
        int var2 = 0;
        int var3 = 0;
        int var4 = -1;
        int[] var6 = new int[var1];
        int[] var7 = new int[var1];
        int[] var8 = new int[var1];
        int[] var9 = new int[var1];

        int var10;
        for (var10 = 0; var10 < var1; ++var10) {
            var6[var10] = 0;
            var7[var10] = 0;
            var9[var10] = 0;
            var8[var10] = 0;
        }

        var7[0] = 1;
        for (var6[0] = 1; var2 < var1; ++var2) {
            int var11 = var0[var2] ? 1 : 0;
            for (var10 = 1; var10 <= var3; ++var10) {
                var11 += var7[var10] * ((var0[var2 - var10]) ? 1 : 0);
            }
            var11 %= 2;
            if (var11 == 1) {
                for (var10 = 0; var10 < var1; ++var10) {
                    var9[var10] = var7[var10];
                    var8[var10] = 0;
                }

                for (var10 = 0; var10 < var1; ++var10) {
                    if (var6[var10] == 1) {
                        var8[var10 + var2 - var4] = 1;
                    }
                }

                for (var10 = 0; var10 < var1; ++var10) {
                    var7[var10] = (var7[var10] + var8[var10]) % 2;
                }

                if (var3 <= var2 / 2) {
                    var3 = var2 + 1 - var3;
                    var4 = var2;

                    for (var10 = 0; var10 < var1; ++var10) {
                        var6[var10] = var9[var10];
                    }
                }
            }
        }
        return var3;
    }
}
