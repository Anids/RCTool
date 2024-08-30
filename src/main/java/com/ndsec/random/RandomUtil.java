package com.ndsec.random;

import org.apache.commons.collections4.Bag;
import org.apache.commons.collections4.bag.HashBag;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import static com.ndsec.random.util.cephes_normal;


public class RandomUtil {
    public static final double SQRT_2 = 1.41421356237309504880;

    /**
     * 《GM/T 0005-2021 随机性检测规范》中的附录A的表A.2 1 000 000比特样本检测设置
     *
     * @param eps    检测数据
     * @param result p值和q值存储列表
     */
    public static void sampleCheck(List<Boolean> eps, ArrayList<Double> result) {
        RandomUtil.frequency(eps, result);
        RandomUtil.blockFrequency(eps, 10000, result);
        RandomUtil.pokerDetect(eps, 4, result);
        RandomUtil.pokerDetect(eps, 8, result);
        RandomUtil.serial(eps, 3, result);
        RandomUtil.serial(eps, 5, result);
        RandomUtil.runs(eps, result);
        RandomUtil.runsDistribution(eps, result);
        RandomUtil.longestRunOfBlock(eps, 1, result);
        RandomUtil.longestRunOfBlock(eps, 0, result);
        RandomUtil.binaryDerivate(eps, 3, result);
        RandomUtil.binaryDerivate(eps, 7, result);
        RandomUtil.selfCorrelation(eps, 1, result);
        RandomUtil.selfCorrelation(eps, 2, result);
        RandomUtil.selfCorrelation(eps, 8, result);
        RandomUtil.selfCorrelation(eps, 16, result);
        RandomUtil.rank(eps, result);
        RandomUtil.cumulativeSums(eps, result);
        RandomUtil.approximateEntropy(eps, 2, result);
        RandomUtil.approximateEntropy(eps, 5, result);
        RandomUtil.linearComplexity(eps, 500, result);
        RandomUtil.linearComplexity(eps, 1000, result);
        RandomUtil.universal(eps, result);
        RandomUtil.discreteFourierTransform(eps, result);
    }

    /**
     * A1 单比特频数检测
     *
     * @param eps 检测数据
     * @param res p值和q值
     */
    public static void frequency(List<Boolean> eps, ArrayList<Double> res) {
        int n = eps.size();
        Bag<?> bag = new HashBag<>(eps);
        double Sn;
        Sn = (bag.getCount(true) - bag.getCount(false)) * 1.0;
        double v = Sn / Math.sqrt(n);
        double pValue = Erf.erfc(Math.abs(v) / SQRT_2);
        double qValue = 0.5 * Erf.erfc(v / SQRT_2);
        bag.clear();
        res.add(pValue);
        res.add(qValue);
    }

    /**
     * A2 块内频数检测
     *
     * @param eps 检测数据
     * @param m   子序列长度
     * @param res p值和q值
     */
    public static void blockFrequency(List<Boolean> eps, int m, ArrayList<Double> res) {
        int n = eps.size();
        int N = n / m;
        int blockSum;
        double pi, v = 0.0;
        for (int i = 0; i < N; i++) {
            blockSum = 0;
            for (int j = 0; j < m; j++) {
                blockSum += eps.get(j + i * m) ? 1 : 0;
            }
            pi = (blockSum * 1.0) / (m * 1.0);
            v += (pi - 0.5) * (pi - 0.5);
        }
        double V = 4 * m * v;
        double pValue = cephes.igamc(N / 2.0, V / 2.0);
        res.add(pValue);
        res.add(pValue);
    }

    /**
     * A3 扑克检测
     *
     * @param eps 检测数据
     * @param m   可重叠子序列长度
     * @param res p值和q值
     */
    public static void pokerDetect(List<Boolean> eps, int m, ArrayList<Double> res) {
        int n = eps.size();
        int N = n / m;
        int m_sum = 1 << m;
        int[] tabs = new int[m_sum];
        double sum_ni = 0.0;
        for (int i = 0; i < m_sum; i++) {
            tabs[i] = 0;
        }
        int tab;
        for (int i = 0; i < N; i++) {
            tab = util.toInt(eps, i * m, i * m + m);
            tabs[tab]++;
        }
        for (int i = 0; i < m_sum; i++) {
            sum_ni += tabs[i] * tabs[i];
        }
        sum_ni *= m_sum;
        sum_ni /= N;
        sum_ni -= N;
        double pValue = cephes.igamc((m_sum - 1) / 2.0, sum_ni / 2.0);
        res.add(pValue);
        res.add(pValue);
    }

    /**
     * A4 重叠子序列检测
     *
     * @param eps 检测数据
     * @param m   可重叠子序列长度
     * @param res p值和q值（四值）
     */
    public static void serial(List<Boolean> eps, int m, ArrayList<Double> res) {
        if (m < 2) {
            res.add(0.0);
            res.add(0.0);
            res.add(0.0);
            res.add(0.0);
            return;
        }
        double psiM0 = util.psi2p(eps, m);
        double psiM1 = util.psi2p(eps, m - 1);
        double psiM2 = util.psi2p(eps, m - 2);
        double delM1 = psiM0 - psiM1;
        double delM2 = psiM0 - 2.0 * psiM1 + psiM2;
        double pValue1 = cephes.igamc(Math.pow(2, m - 2), delM1 / 2.0);
        double pValue2 = cephes.igamc(Math.pow(2, m - 3), delM2 / 2.0);
        res.add(pValue1);
        res.add(pValue1);
        res.add(pValue2);
        res.add(pValue2);
    }

    /**
     * A5 游程总数检测
     *
     * @param eps 检测数据
     * @param res p值和q值
     */
    public static void runs(List<Boolean> eps, ArrayList<Double> res) {
        int S, k;
        double pi, V, erfc_arg, pValue, qValue;
        int n = eps.size();
        S = 0;
        for (k = 0; k < n; k++) {
            if (eps.get(k)) S++;
        }
        pi = S / (n * 1.0);
        if (Math.abs(pi - 0.5) > (2.0 / Math.sqrt(n))) {
            pValue = 0.0;
            qValue = 0.0;
        } else {
            V = 1;
            for (k = 1; k < n; k++) {
                if (eps.get(k) != eps.get(k - 1)) {
                    V++;
                }
            }
            erfc_arg = (V - 2.0 * n * pi * (1 - pi)) / (2.0 * pi * (1 - pi) * Math.sqrt(2 * n));
            pValue = Erf.erfc(Math.abs(erfc_arg));
            qValue = 0.5 * Erf.erfc(erfc_arg);
        }
        res.add(pValue);
        res.add(qValue);
    }

    /**
     * A6 游程分布检测
     *
     * @param eps 检测数据
     * @param res p值和q值
     */
    public static void runsDistribution(List<Boolean> eps, ArrayList<Double> res) {
        int n = eps.size();
        double[] ei = new double[n + 1];
        int k = 0;
        for (int i = 1; i <= n; i++) {
            ei[i - 1] = (n - i + 3) / Math.pow(2, i + 2);
            if (ei[i - 1] >= 5) {
                k = i;
            }
        }
        double[] bi = new double[n + 1];
        double[] gi = new double[n + 1];
        eps.add(!eps.get(n - 1));
        boolean flag = eps.get(0);
        int j = 1;
        for (int i = 1; i <= n; i++) {
            if (eps.get(i) != flag) {
                if (j <= k) {
                    if (!flag) {
                        gi[j]++;
                    } else {
                        bi[j]++;
                    }
                } else {
                    if (!flag) {
                        gi[k]++;
                    } else {
                        bi[k]++;
                    }
                }
                flag = eps.get(i);
                j = 1;
            } else {
                ++j;
            }
        }
        eps.remove(eps.size() - 1);
        double T = 0.0;
        for (int i = 1; i <= k; i++) {
            T += bi[i] + gi[i];
        }
        double[] e1 = new double[k + 1];
        for (int i = 1; i <= k; i++) {
            if (i == k) {
                e1[i] = T / Math.pow(2, i);
            } else {
                e1[i] = T / Math.pow(2, i + 1);
            }
        }
        double V = 0.0;
        for (int i = 1; i <= k; i++) {
            V += Math.pow(bi[i] - e1[i], 2) / e1[i] + Math.pow(gi[i] - e1[i], 2) / e1[i];
        }
        double pValue = cephes.igamc(k - 1, V / 2.0);
        res.add(pValue);
        res.add(pValue);
    }

    /**
     * A7 块内最大游程检测
     *
     * @param eps 检测数据
     * @param fg  游程值：1、最大“1”游程；0、最大“0”游程
     * @param res p值和q值
     */
    public static void longestRunOfBlock(List<Boolean> eps, int fg, ArrayList<Double> res) {
        int n = eps.size();
        int run, v_n_obs, N, i, j, K, M;
        int[] V = new int[7];
        int[] nu = new int[7];
        double pValue, chi2;
        double[] pi = new double[7];
        if (n < 128) {
            res.add(0.0);
            res.add(0.0);
            return;
        }
        if (n < 6272) {
            K = 3;
            M = 8;
            V[0] = 1;
            V[1] = 2;
            V[2] = 3;
            V[3] = 4;
            pi[0] = 0.2148;
            pi[1] = 0.3672;
            pi[2] = 0.2305;
            pi[3] = 0.1875;
        } else if (n < 750000) {
            K = 5;
            M = 128;
            V[0] = 4;
            V[1] = 5;
            V[2] = 6;
            V[3] = 7;
            V[4] = 8;
            V[5] = 9;
            pi[0] = 0.1174;
            pi[1] = 0.2430;
            pi[2] = 0.2494;
            pi[3] = 0.1752;
            pi[4] = 0.1027;
            pi[5] = 0.1124;
        } else {
            K = 6;
            M = 10000;
            V[0] = 10;
            V[1] = 11;
            V[2] = 12;
            V[3] = 13;
            V[4] = 14;
            V[5] = 15;
            V[6] = 16;
            pi[0] = 0.086632;
            pi[1] = 0.208201;
            pi[2] = 0.248419;
            pi[3] = 0.193913;
            pi[4] = 0.121458;
            pi[5] = 0.068011;
            pi[6] = 0.073366;
        }

        N = n / M;
        for (i = 0; i < N; i++) {
            v_n_obs = 0;
            run = 0;
            for (j = 0; j < M; j++) {
                int ax = eps.get(i * M + j) ? 1 : 0;
                if (ax == fg) {
                    run++;
                    if (run > v_n_obs)
                        v_n_obs = run;
                } else
                    run = 0;
            }
            if (v_n_obs < V[0])
                nu[0]++;
            for (j = 0; j <= K; j++) {
                if (v_n_obs == V[j])
                    nu[j]++;
            }
            if (v_n_obs > V[K])
                nu[K]++;
        }
        chi2 = 0.0;
        for (i = 0; i <= K; i++) {
            double v = nu[i] - N * pi[i];
            chi2 += (v * v) / (N * pi[i]);
        }
        pValue = cephes.igamc(K / 2.0, chi2 / 2.0);
        res.add(pValue);
        res.add(pValue);

    }

    /**
     * A8 二元推导检测
     *
     * @param eps 检测数据
     * @param k   二元推导次数
     * @param res p值和q值
     */
    public static void binaryDerivate(List<Boolean> eps, int k, ArrayList<Double> res) {
        int n = eps.size();
        int Sn_k = 0;
        int n_k = n - k;
        double V, pValue = 0.0, sqrt2 = 1.41421356237309504880;
        int i, j;
        int[] epsExt = new int[n];
        for (i = 0; i < n; i++) {
            epsExt[i] = eps.get(i) ? 1 : 0;
        }
        for (i = 0; i < k; i++) {
            for (j = 0; j < n - 1; j++) {
                epsExt[j] = epsExt[j] ^ epsExt[j + 1];
            }
        }
        for (i = 0; i < n_k; i++) {
            Sn_k += (2 * epsExt[i]) - 1;
        }
        V = Sn_k / Math.sqrt(n_k);
        pValue = Erf.erfc(Math.abs(V) / sqrt2);
        double qValue = Erf.erfc(V / sqrt2) * 0.5;
        res.add(pValue);
        res.add(qValue);
    }

    /**
     * A9 自相关检测
     *
     * @param eps 检测数据
     * @param d   逻辑左移位数
     * @param res p值和q值
     */
    public static void selfCorrelation(List<Boolean> eps, int d, ArrayList<Double> res) {
        int n = eps.size();
        int i;
        int n_d = n - d;
        int Ad = 0;
        double V, pValue, sqrt2 = 1.41421356237309504880;
        for (i = 0; i < n_d; ++i) {
            Ad += (eps.get(i) ^ eps.get(i + d)) ? 1 : 0;
        }
        V = 2 * (Ad - (n_d / 2.0)) / Math.sqrt(n_d);
        pValue = Erf.erfc(Math.abs(V) / sqrt2);
        double qValue = 0.5 * Erf.erfc(V / sqrt2);
        res.add(pValue);
        res.add(qValue);
    }

    /**
     * A10 矩阵秩检测
     * 按标准设置常量计算
     *
     * @param eps 检测数据
     * @param res p值和q值
     */
    public static void rank(List<Boolean> eps, ArrayList<Double> res) {
        int n = eps.size();
        int ret = 0;
        int N, i, k, r;
        double pValue, product, chi_squared, arg1, p_32, p_31, p_30, R, F_32, F_31, F_30;
        int[][] matrix = new int[32][32];
        N = n / (32 * 32);
        if (N == 0.e0) {
            res.add(0.0);
            res.add(0.0);
            return;
        } else {
            p_32 = 0.2888;
            p_31 = 0.5776;
            p_30 = 0.1336;
            F_32 = 0;
            F_31 = 0;
            for (k = 0; k < N; k++) {
                util.def_matrix(32, 32, matrix, k, eps);
                R = util.computeRank(32, 32, matrix);
                if (R == 32)
                    F_32++;
                if (R == 31)
                    F_31++;
            }
            F_30 = N - (F_32 + F_31);
            chi_squared = Math.pow(F_32 - N * p_32, 2) / (N * p_32) + Math.pow(F_31 - N * p_31, 2) / (N * p_31)
                    + Math.pow(F_30 - N * p_30, 2) / (N * p_30);
            pValue = cephes.igamc(1, chi_squared / 2.0);
            res.add(pValue);
            res.add(pValue);
        }
    }

    /**
     * A10 矩阵秩检测
     * 按原理计算
     *
     * @param eps 检测数据
     * @param res p值和q值
     */
    public static void rank_bak(List<Boolean> eps, ArrayList<Double> res) {
        int n = eps.size();
        int ret = 0;
        int N, i, k, r;
        double pValue, product, chi_squared, arg1, p_32, p_31, p_30, R, F_32, F_31, F_30;
        int[][] matrix = new int[32][32];
        N = n / (32 * 32);
        if (N == 0.e0) {
            res.add(0.0);
            res.add(0.0);
            return;
        } else {
            r = 32;
            product = 1;
            for (i = 0; i < r; i++) {
                product *= ((1.e0 - Math.pow(2, i - 32)) * (1.e0 - Math.pow(2, i - 32))) / (1.e0 - Math.pow(2, i - r));
            }
            p_32 = product;
            r = 31;
            product = 1;
            for (i = 0; i < r; i++) {
                product *= ((1.e0 - Math.pow(2, i - 32)) * (1.e0 - Math.pow(2, i - 32))) / (1.e0 - Math.pow(2, i - r));
            }
            p_31 = Math.pow(2, r * (32 + 32 - r) - 32 * 32) * product;
            p_30 = 1 - (p_32 + p_31);
            F_32 = 0;
            F_31 = 0;
            for (k = 0; k < N; k++) {
                util.def_matrix(32, 32, matrix, k, eps);
                R = util.computeRank(32, 32, matrix);
                if (R == 32)
                    F_32++;
                if (R == 31)
                    F_31++;
            }
            F_30 = N - (F_32 + F_31);
            chi_squared = Math.pow(F_32 - N * p_32, 2) / (N * p_32) + Math.pow(F_31 - N * p_31, 2) / (N * p_31)
                    + Math.pow(F_30 - N * p_30, 2) / (N * p_30);
            pValue = cephes.igamc(1, chi_squared / 2.0);
            res.add(pValue);
            res.add(pValue);
        }
    }

    /**
     * A11 累加和检测
     *
     * @param eps 检测数据
     * @param res p值和q值
     */
    public static void cumulativeSums(List<Boolean> eps, ArrayList<Double> res) {
        int n = eps.size();
        int S, sup, inf, z = 0, zrev = 0, k;
        double sum1, sum2, p1Value, p2Value;
        S = 0;
        sup = 0;
        inf = 0;
        for (k = 0; k < n; k++) {
            if (eps.get(k)) S++;
            else S--;
            if (S > sup)
                sup++;
            if (S < inf)
                inf--;
            z = Math.max(sup, -inf);
            zrev = Math.max(sup - S, S - inf);
        }
        // forward
        sum1 = 0.0;

        for (k = (-n / z + 1) / 4; k <= (n / z - 1) / 4; k++) {
            sum1 += cephes_normal(((4 * k + 1) * z) / Math.sqrt(n));
            sum1 -= cephes_normal(((4 * k - 1) * z) / Math.sqrt(n));
        }
        sum2 = 0.0;
        for (k = (-n / z - 3) / 4; k <= (n / z - 1) / 4; k++) {
            sum2 += cephes_normal(((4 * k + 3) * z) / Math.sqrt(n));
            sum2 -= cephes_normal(((4 * k + 1) * z) / Math.sqrt(n));
        }
        p1Value = 1.0 - sum1 + sum2;
        // backward
        sum1 = 0.0;
        for (k = (-n / zrev + 1) / 4; k <= (n / zrev - 1) / 4; k++) {
            sum1 += cephes_normal(((4 * k + 1) * zrev) / Math.sqrt(n));
            sum1 -= cephes_normal(((4 * k - 1) * zrev) / Math.sqrt(n));
        }
        sum2 = 0.0;
        for (k = (-n / zrev - 3) / 4; k <= (n / zrev - 1) / 4; k++) {
            sum2 += cephes_normal(((4 * k + 3) * zrev) / Math.sqrt(n));
            sum2 -= cephes_normal(((4 * k + 1) * zrev) / Math.sqrt(n));
        }
        p2Value = 1.0 - sum1 + sum2;
        res.add(p1Value);
        res.add(p1Value);
        res.add(p2Value);
        res.add(p2Value);
    }

    /**
     * A12 近似熵检测
     *
     * @param eps 检测数据
     * @param m   可重叠子序列长度
     * @param res p值和q值
     */
    public static void approximateEntropy(List<Boolean> eps, int m, ArrayList<Double> res) {
        int n = eps.size();
        int i, j, k, r, blockSize, seqLength, powLen, index;
        double sum, numOfBlocks, apen, chi_squared, p_value;
        double[] ApEn = new double[2];
        int[] P;
        seqLength = n;
        r = 0;
        for (blockSize = m; blockSize <= m + 1; blockSize++) {
            if (blockSize == 0) {
                ApEn[0] = 0.00;
                r++;
            } else {
                numOfBlocks = (double) seqLength;
                powLen = (int) Math.pow(2, blockSize + 1) - 1;
                P = new int[powLen];
                for (i = 1; i < powLen - 1; i++)
                    P[i] = 0;
                for (i = 0; i < numOfBlocks; i++) { /* COMPUTE FREQUENCY */
                    k = 1;
                    for (j = 0; j < blockSize; j++) {
                        k <<= 1;
                        if (eps.get((i + j) % seqLength))
                            k++;
                    }
                    P[k - 1]++;
                }
                /* DISPLAY FREQUENCY */
                sum = 0.0;
                index = (int) Math.pow(2, blockSize) - 1;
                for (i = 0; i < (int) Math.pow(2, blockSize); i++) {
                    if (P[index] > 0)
                        sum += P[index] * Math.log(P[index] / numOfBlocks);
                    index++;
                }
                sum /= numOfBlocks;
                ApEn[r] = sum;
                r++;
            }
        }
        apen = ApEn[0] - ApEn[1];
        chi_squared = 2.0 * seqLength * (Math.log(2) - apen);
        p_value = cephes.igamc(Math.pow(2, m - 1), chi_squared / 2.0);
        res.add(p_value);
        res.add(p_value);
    }

    /**
     * A13 线性复杂度检测
     * 按标准计算
     *
     * @param eps 检测数据
     * @param M   子序列长度
     * @param res p值和q值
     */
    public static void linearComplexity(List<Boolean> eps, int M, ArrayList<Double> res) {
        int n = eps.size();
        int N = n / M;
        double[] nu = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double[] pi = new double[]{0.010417, 0.031250, 0.125, 0.500, 0.250, 0.062500, 0.020833};
        double chi2 = 0.0;
        boolean[] v11 = new boolean[M];
        double v15 = (double) M / 2.0 + (9.0 + (M % 2 == 0 ? -1.0 : 1.0)) / 36.0 - ((double) M / 3.0 + 0.2222222222222222) / Math.pow(2.0, M);
        Iterator<Boolean> v17 = eps.iterator();
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                v11[j] = v17.next();
            }
            int v12 = util.liner_res(v11, M);
            double v13 = (M % 2 == 0 ? 1.0 : -1.0) * ((double) v12 - v15) + 0.2222222222222222;
            if (v13 <= -2.5) {
                nu[0]++;
            } else if (v13 <= -1.5) {
                nu[1]++;
            } else if (v13 <= -0.5) {
                nu[2]++;
            } else if (v13 <= 0.5) {
                nu[3]++;
            } else if (v13 <= 1.5) {
                nu[4]++;
            } else if (v13 <= 2.5) {
                nu[5]++;
            } else {
                nu[6]++;
            }
        }
        for (int i = 0; i < 7; i++) {
            chi2 += Math.pow(nu[i] - N * pi[i], 2) / (N * pi[i]);
        }
        double p_value = cephes.igamc(3.0, chi2 / 2.0);
        res.add(p_value);
        res.add(p_value);
    }

    /**
     * A13 线性复杂度检测
     * 按原理计算
     *
     * @param eps 检测数据
     * @param M   子序列长度
     * @param res p值和q值
     */
    public static void linearComplexity_bak(List<Boolean> eps, int M, ArrayList<Double> res) {
        int n = eps.size();
        int i, ii, j, d, N, L, m, N_, parity, sign, K = 6;
        double p_value, T_, mean, chi2;
        double[] nu = new double[7];
        double[] pi = {0.01047, 0.031250, 0.125, 0.500, 0.250, 0.062500, 0.020833};
        int[] T = new int[M];
        int[] P = new int[M];
        int[] B_ = new int[M];
        int[] C = new int[M];
        N = n / M;
        for (i = 0; i < K + 1; i++)
            nu[i] = 0.00;
        for (ii = 0; ii < N; ii++) {
            for (i = 0; i < M; i++) {
                B_[i] = 0;
                C[i] = 0;
                T[i] = 0;
                P[i] = 0;
            }
            L = 0;
            m = -1;
            d = 0;
            C[0] = 1;
            B_[0] = 1;

            /* DETERMINE LINEAR COMPLEXITY */
            N_ = 0;
            while (N_ < M) {
                d = eps.get(ii * M + N_) ? 1 : 0;
                for (i = 1; i <= L; i++)
                    d += C[i] * (eps.get(ii * M + N_ - i) ? 1 : 0);
                d = d % 2;
                if (d == 1) {
                    for (i = 0; i < M; i++) {
                        T[i] = C[i];
                        P[i] = 0;
                    }
                    for (j = 0; j < M; j++)
                        if (B_[j] == 1)
                            P[j + N_ - m] = 1;
                    for (i = 0; i < M; i++)
                        C[i] = (C[i] + P[i]) % 2;
                    if (L <= N_ / 2) {
                        L = N_ + 1 - L;
                        m = N_;
                        for (i = 0; i < M; i++)
                            B_[i] = T[i];
                    }
                }
                N_++;
            }
            if ((parity = (M + 1) % 2) == 0)
                sign = -1;
            else
                sign = 1;
            mean = M / 2.0 + (9.0 + sign) / 36.0 - 1.0 / Math.pow(2, M) * (M / 3.0 + 2.0 / 9.0);
            if ((parity = M % 2) == 0)
                sign = 1;
            else
                sign = -1;
            T_ = sign * (L - mean) + 2.0 / 9.0;

            if (T_ <= -2.5)
                nu[0]++;
            else if (T_ > -2.5 && T_ <= -1.5)
                nu[1]++;
            else if (T_ > -1.5 && T_ <= -0.5)
                nu[2]++;
            else if (T_ > -0.5 && T_ <= 0.5)
                nu[3]++;
            else if (T_ > 0.5 && T_ <= 1.5)
                nu[4]++;
            else if (T_ > 1.5 && T_ <= 2.5)
                nu[5]++;
            else if (T_ > 2.5)
                nu[6]++;
        }
        chi2 = 0.00;
        for (i = 0; i < K + 1; i++) {
            chi2 += Math.pow(nu[i] - N * pi[i], 2) / (N * pi[i]);
        }
        p_value = cephes.igamc(3.0, chi2 / 2.0);
        res.add(p_value);
        res.add(p_value);
    }

    /**
     * A14 通用统计检测
     *
     * @param eps 检测数据
     * @param res p值和q值
     */
    public static void universal(List<Boolean> eps, ArrayList<Double> res) {
        int n = eps.size();
        int ret = 0;
        int i, j, p, L, Q, K;
        double arg, sqrt2, sigma, phi, sum, p_value, c;
        long decRep;
        double[] expected_value = {0, 0, 0, 0, 0, 0, 5.2177052, 6.1962507, 7.1836656,
                8.1764248, 9.1723243, 10.170032, 11.168765,
                12.168070, 13.167693, 14.167488, 15.167379};
        double[] variance = {0, 0, 0, 0, 0, 0, 2.954, 3.125, 3.238, 3.311, 3.356, 3.384,
                3.401, 3.410, 3.416, 3.419, 3.421};
        L = 5;
        if (n >= 387840) L = 6;
        if (n >= 904960) L = 7;
        if (n >= 2068480) L = 8;
        if (n >= 4654080) L = 9;
        if (n >= 10342400) L = 10;
        if (n >= 22753280) L = 11;
        if (n >= 49643520) L = 12;
        if (n >= 107560960) L = 13;
        if (n >= 231669760) L = 14;
        if (n >= 496435200) L = 15;
        if (n >= 1059061760) L = 16;

        Q = (int) (Math.pow(2, L)) * 10;
        K = (int) (Math.floor(n / L) - (double) Q);

        p = (int) Math.pow(2, L);
        long[] T = new long[p];
        if (L < 6 || (double) Q < 10 * Math.pow(2, L)) {
            res.add(0.0);
            res.add(0.0);
            return;
        }
        c = 0.7 - 0.8 / (double) L + (4 + 32 / (double) L) * Math.pow(K, -3 / (double) L) / 15;
        sigma = c * Math.sqrt(variance[L] / (double) K);
        sqrt2 = Math.sqrt(2);
        sum = 0.0;
        for (i = 0; i < p; i++)
            T[i] = 0;
        for (i = 1; i <= Q; i++) {        /* INITIALIZE TABLE */
            decRep = 0;
            for (j = 0; j < L; j++)
                decRep += (eps.get((i - 1) * L + j) ? 1 : 0) * (long) Math.pow(2, L - 1 - j);
            T[(int) decRep] = i;
        }
        for (i = Q + 1; i <= Q + K; i++) {    /* PROCESS BLOCKS */
            decRep = 0;
            for (j = 0; j < L; j++)
                decRep += (eps.get((i - 1) * L + j) ? 1 : 0) * (long) Math.pow(2, L - 1 - j);
            sum += Math.log(i - T[(int) decRep]) / Math.log(2);
            T[(int) decRep] = i;
        }
        phi = (sum / (double) K);
        arg = (phi - expected_value[L]) / (sqrt2 * sigma);
        p_value = Erf.erfc(Math.abs(arg));
        double q_value = Erf.erfc(arg) * 0.5;
        res.add(p_value);
        res.add(q_value);
    }

    /**
     * A15 离散傅里叶检测
     *
     * @param eps 检测数据
     * @param res p值和q值
     */
    public static void discreteFourierTransform(List<Boolean> eps, ArrayList<Double> res) {
        int n = eps.size();
        double[] X = new double[n];
        double T = Math.sqrt(2.995732274 * n);
        double N0 = 0.95 * n / 2.0;
        int N1 = 0;
        for (int i = 0; i < n; i++) {
            if (eps.get(i)) {
                X[i] = 1.0;
            } else {
                X[i] = -1.0;
            }
        }
        X = util.getFFT(X);
        FastFourierTransformer F0 = new FastFourierTransformer(DftNormalization.STANDARD);
        Complex[] v10 = F0.transform(X, TransformType.FORWARD);
        for (int i = 0; i < n / 2 - 1; i++) {
            if (v10[i].abs() < T) {
                N1++;
            }
        }
        double V = (N1 - N0) / Math.sqrt(0.95 * 0.05 * n / 3.8);
        double pValue = Erf.erfc(Math.abs(V) / Math.sqrt(2.0));
        double qValue = Erf.erfc(V / Math.sqrt(2.0)) / 2.0;
        res.add(pValue);
        res.add(qValue);
    }


}
