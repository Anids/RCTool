package com.ndsec;

import com.ndsec.random.RandomUtil;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

public class SelfTest {
    // 《GM/T 0005-2021 随机性检测规范》上附录C检测数据（2组）
    public static final String TEST_DATA_1 = "11001100000101010110110001001100111000000000001001001101010100010001001111010110100000001101011111001100111001101101100010110010";
    public static final String TEST_DATA_2 = "1100100100001111110110101010001000100001011010001100001000110100110001001100011001100010100010111000";

    public static final double[] STAND_DATA_LIST = {0.215925, 0.892038, 0.706438, 0.706438, 0.213734, 0.213734, 0.436868, 0.436868, 0.723674, 0.723674, 0.620729, 0.310364, 0.970152, 0.970152, 0.839299, 0.839299, 0.180598, 0.180598, 0.039669, 0.980166, 0.790080, 0.395040, 0.307543, 0.307543, 0.219194, 0.219194, 0.114866, 0.114866, 0.235301, 0.235301, 0.844721, 0.844721, 0.282568, 0.141284, 0.654721, 0.327360};

    public static final String[] ITEM_CHECK_NAME = {"A1 单比特频数检测", "A2 块内频数检测", "A3 扑克检测", "A4.1 重叠子序列检测", "A4.2 重叠子序列检测", "A5 游程总数检测", "A6 游程分布检测", "A7 块内最大0游程检测", "A7 块内最大1游程检测", "A8 二元推导检测", "A9 自相关检测", "A10 矩阵秩检测", "A11.1 累加和检测", "A11.2 累加和检测", "A12 近似熵检测", "A13 线性复杂度检测", "A14 Maurer通用统计检测", "A15 离散傅里叶检测"};


    public static ArrayList<Double> result = new ArrayList<>();

    public static void StringToBoolean(String a, List<Boolean> b) {
        for (int i = 0; i < a.length(); i++) {
            b.add(a.charAt(i) == '1');
        }
    }

    public static void frequencyTest() {
        List<Boolean> res = new ArrayList<>();
        StringToBoolean(TEST_DATA_1, res);
        RandomUtil.frequency(res, result);
    }

    public static void blockFrequencyTest() {
        List<Boolean> res = new ArrayList<>();
        StringToBoolean(TEST_DATA_2, res);
        RandomUtil.blockFrequency(res, 10, result);

    }

    public static void PokerDetectTest() {
        List<Boolean> res = new ArrayList<>();
        StringToBoolean(TEST_DATA_1, res);
        RandomUtil.pokerDetect(res, 4, result);
    }

    public static void serialTest() {
        List<Boolean> res = new ArrayList<>();
        StringToBoolean(TEST_DATA_1, res);
        RandomUtil.serial(res, 2, result);
    }

    public static void runsTest() {
        List<Boolean> res = new ArrayList<>();
        StringToBoolean(TEST_DATA_1, res);
        RandomUtil.runs(res, result);
    }

    public static void runsDistributionTest() {
        List<Boolean> res = new ArrayList<>();
        StringToBoolean(TEST_DATA_1, res);
        RandomUtil.runsDistribution(res, result);
    }

    public static void LongestRunTest(int x) {
        List<Boolean> res = new ArrayList<>();
        StringToBoolean(TEST_DATA_1, res);
        RandomUtil.longestRunOfBlock(res, x, result);
    }

    public static void binaryDerivateTest() {
        List<Boolean> res = new ArrayList<>();
        StringToBoolean(TEST_DATA_1, res);
        RandomUtil.binaryDerivate(res, 3, result);
    }

    public static void selfCorrelationTest() {
        List<Boolean> res = new ArrayList<>();
        StringToBoolean(TEST_DATA_1, res);
        RandomUtil.selfCorrelation(res, 1, result);
    }

    public static void rankTest() {
        List<Boolean> res = new ArrayList<>();
        readFile(new File("src/main/resources/e/EData.txt"), res);
        List<Boolean> eps_f = new ArrayList<>();
//        System.out.println(res.size());
        int count = 1000000;
        for (int i = 0; i < count; i++) {
            eps_f.add(res.get(i));
        }
        RandomUtil.rank(eps_f, result);
    }

    public static void cumulativeSumsTest() {
        List<Boolean> res = new ArrayList<>();
        StringToBoolean(TEST_DATA_2, res);
        RandomUtil.cumulativeSums(res, result);
    }

    public static void approximateEntropyTest() {
        List<Boolean> res = new ArrayList<>();
        StringToBoolean(TEST_DATA_2, res);
        RandomUtil.approximateEntropy(res, 2, result);
    }

    public static void linearComplexityTest() {
        List<Boolean> res = new ArrayList<>();
        readFile(new File("src/main/resources/e/EData.txt"), res);
        List<Boolean> eps = new ArrayList<>();
        int count = 1000000;
        for (int i = 0; i < count; i++) {
            eps.add(res.get(i));
        }
        RandomUtil.linearComplexity(eps, 1000, result);
    }

    public static void universalTest() {
        List<Boolean> res = new ArrayList<>();
        readFile(new File("src/main/resources/e/EData.txt"), res);
        List<Boolean> eps = new ArrayList<>();
        int count = 1000000;
        for (int i = 0; i < count; i++) {
            eps.add(res.get(i));
        }
        RandomUtil.universal(eps, result);
    }

    public static void discreteFourierTransformTest() {
        List<Boolean> res = new ArrayList<>();
        StringToBoolean(TEST_DATA_2, res);
        RandomUtil.discreteFourierTransform(res, result);
    }

    public static void readFile(File file, List<Boolean> eps) {
        if (file.isFile() && file.exists()) {
            try {
                InputStreamReader inputStreamReader = new InputStreamReader(Files.newInputStream(file.toPath()));
                BufferedReader bufferedReader = new BufferedReader(inputStreamReader);
                String lineStr;
                try {
                    while ((lineStr = bufferedReader.readLine()) != null) {
                        for (int i = 0; i < lineStr.length(); i++) {
                            if (lineStr.charAt(i) == '1') {
                                eps.add(true);
                            }
                            if (lineStr.charAt(i) == '0') {
                                eps.add(false);
                            }
                        }
                    }
                    bufferedReader.close();
                    inputStreamReader.close();
                } catch (Exception e) {
                    e.printStackTrace();
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    public static boolean testRandom() {
        frequencyTest();
        blockFrequencyTest();
        PokerDetectTest();
        serialTest();
        runsTest();
        runsDistributionTest();
        LongestRunTest(0);
        LongestRunTest(1);
        binaryDerivateTest();
        selfCorrelationTest();
        rankTest();
        cumulativeSumsTest();
        approximateEntropyTest();
        linearComplexityTest();
        universalTest();
        discreteFourierTransformTest();
        boolean flag = true;
        for (int i = 0; i < result.size(); i += 2) {
            String s1 = String.format("%.6f", result.get(i));
            double s11 = Double.parseDouble(s1);
            String s2 = String.format("%.6f", result.get(i + 1));
            double s22 = Double.parseDouble(s2);
            if ((s11 == STAND_DATA_LIST[i]) && (s22 == STAND_DATA_LIST[i + 1])) {
                System.out.println(ITEM_CHECK_NAME[i / 2] + " [PASS]");
            } else {
                System.out.println(ITEM_CHECK_NAME[i / 2] + " [FAIL]");
                System.out.println("\tExpect - PValue : " + STAND_DATA_LIST[i]);
                System.out.println("\tActual - PValue : " + result.get(i));
                System.out.println("\tExpect - QValue : " + STAND_DATA_LIST[i + 1]);
                System.out.println("\tActual - QValue : " + result.get(i + 1));
                flag = false;
            }
        }
        return flag;
    }


    public static void main(String[] args) {
        if (testRandom()) {
            System.out.println("Self test pass!");
        } else {
            System.out.println("Self test fail!");
        }
    }
}
