package com.ndsec.random;

import org.apache.commons.math3.special.Gamma;

public class cephes {
    public static final double MAX_LOG = 709.782712893384;
    public static final double MAC_HEP = 1.1102230246251565E-16;
    public static final double BIG_INV = 2.220446049250313E-16;
    public static final double BIG = 4.503599627370496E15;

    public cephes() {
    }


    public static double igam(double a, double x) {
        if (!(x <= 0.0) && !(a <= 0.0)) {
            if (x > 1.0 && x > a) {
                return 1.0 - igamc(a, x);
            } else {
                double ax = a * Math.log(x) - x - Gamma.logGamma(a);
                if (ax < -MAX_LOG) {
                    return 0.0;
                } else {
                    ax = Math.exp(ax);
                    double r = a;
                    double c = 1.0;
                    double ans = 1.0;
                    do {
                        r += 1.0;
                        c *= x / r;
                        ans += c;
                    } while (c / ans > MAC_HEP);
                    return ans * ax / a;
                }
            }
        } else {
            return 0.0;
        }
    }

    public static double igamc(double a, double x) {
        if (!(x <= 0.0) && !(a <= 0.0)) {
            if (!(x < 1.0) && !(x < a)) {
                double ax = a * Math.log(x) - x - Gamma.logGamma(a);
                if (ax < -MAX_LOG) {
                    return 0.0;
                } else {
                    ax = Math.exp(ax);
                    double y = 1.0 - a;
                    double z = x + y + 1.0;
                    double c = 0.0;
                    double pkm2 = 1.0;
                    double qkm2 = x;
                    double pkm1 = x + 1.0;
                    double qkm1 = z * x;
                    double ans = pkm1 / qkm1;

                    double t;
                    do {
                        c++;
                        y++;
                        z += 2.0;
                        double yc = y * c;
                        double pk = pkm1 * z - pkm2 * yc;
                        double qk = qkm1 * z - qkm2 * yc;
                        if (qk != 0) {
                            double r = pk / qk;
                            t = Math.abs((ans - r) / r);
                            ans = r;
                        } else {
                            t = 1.0;
                        }
                        pkm2 = pkm1;
                        pkm1 = pk;
                        qkm2 = qkm1;
                        qkm1 = qk;
                        if (Math.abs(pk) > BIG) {
                            pkm2 *= BIG_INV;
                            pkm1 *= BIG_INV;
                            qkm2 *= BIG_INV;
                            qkm1 *= BIG_INV;
                        }
                    } while (t > MAC_HEP);
                    return ans * ax;
                }
            } else {
                return 1.0 - igam(a, x);
            }
        } else {
            return 1.0;
        }
    }
}
