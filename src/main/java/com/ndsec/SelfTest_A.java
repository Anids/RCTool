package com.ndsec;

import com.ndsec.common.Utils;
import com.ndsec.random.RandomUtil;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Objects;

public class SelfTest_A {
    public static void main(String[] args) throws IOException {
        ArrayList<Double> totalRes = new ArrayList<>();
        File[] files = Utils.getFiles("random_check_data_128K");
        for (File file : Objects.requireNonNull(files)) {
            ArrayList<Boolean> eps = new ArrayList<>();
            Utils.readBinaryFile(file, eps);
            RandomUtil.frequency(eps, totalRes);
//            RandomUtil.blockFrequency(eps, 10000, totalRes);
            eps.clear();
        }

        double[] pValueList = new double[files.length];
        double[] qValueList = new double[files.length];
        int count1 = 0, count2 = 0;
        for (int i = 0; i < totalRes.size(); i++) {
            if (i % 2 == 0) {
                pValueList[count1] = totalRes.get(i);
                count1++;
            } else {
                qValueList[count2] = totalRes.get(i);
                count2++;
            }
        }
        int pCount = Utils.getPassOfP(pValueList);
        boolean pRes = Utils.statPValue(pValueList);
        double qStat = Utils.getStatQValue(qValueList);
        boolean qRes = Utils.statQValue(qValueList);
        System.out.println("P值： " + pCount);
        System.out.println("Q值： " + Double.parseDouble(String.format("%.6f", qStat)));
    }
}
