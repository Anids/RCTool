package com.ndsec;

import com.ndsec.common.Utils;
import com.ndsec.task.sampleTask;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

public class RandomTest_A {
    public static final int RANDOM_SUM = 1000;
    private static final int THREAD_SIZE = Runtime.getRuntime().availableProcessors();

    public static final String[] sampleSchemeStr = {"A1 单比特频数检测", "A2 块内频数检测", "A3.1 扑克检测", "A3.2 扑克检测", "A4.1 重叠子序列检测",
            "A4.2 重叠子序列检测", "A4.3 重叠子序列检测", "A4.4 重叠子序列检测", "A5 游程总数检测", "A6 游程分布检测", "A7.1 块内最大游程检测", "A7.2 块内最大游程检测", "A8.1 二元推导检测",
            "A8.2 二元推导检测", "A9.1 自相关检测", "A9.2 自相关检测", "A9.3 自相关检测", "A9.4 自相关检测", "A10 矩阵秩检测",
            "A11.1 累加和检测", "A11.2 累加和检测", "A12.1 近似熵检测", "A12.2 近似熵检测", "A13.1 线性复杂度检测", "A13.2 线性复杂度检测", "A14 通用统计检测", "A15 离散傅里叶检测"};

    public static void threadSampleTest(String filePath) throws IOException {
        File[] files = Utils.getFiles(filePath);
        ArrayList<Double> totalRes = new ArrayList<>();
        long startTime = System.currentTimeMillis();
        long[] end = new long[1];
        ExecutorService threadPool = Executors.newWorkStealingPool(THREAD_SIZE);
        AtomicInteger sum = new AtomicInteger(0);
        Future<ArrayList<Double>>[] futures = new Future[RANDOM_SUM];
        if (Objects.requireNonNull(files).length == RANDOM_SUM) {
            for (int i = 0; i < RANDOM_SUM; i++) {
                List<Boolean> eps = new ArrayList<>();
                Utils.readBinaryFile(files[i], eps);
                futures[i] = threadPool.submit(new sampleTask(eps, files[i].getName(), (String rt) -> {
                    int num = sum.incrementAndGet();
                    end[0] = System.currentTimeMillis();
                    System.out.println("file/" + rt + ", finished! id=" + num + ", cost time = " + (end[0] - startTime) + "ms");
                }));
            }

            for (int i = 0; i < RANDOM_SUM; i++) {
                try {
                    ArrayList<Double> rs = futures[i].get();
                    totalRes.addAll(rs);
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }

            // 关闭线程池
            try {
                threadPool.shutdown();
                if (!threadPool.awaitTermination(100, TimeUnit.SECONDS)) {
                    threadPool.shutdownNow();
                    if (!threadPool.awaitTermination(100, TimeUnit.SECONDS)) {
                        System.err.println("线程池无法及时关闭");
                    }
                }
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
            System.out.println("完成" + RANDOM_SUM + "组数据的计算，花费总时间：" + (end[0] - startTime) + "ms");
        }

        double[] pValueList = new double[RANDOM_SUM];
        double[] qValueList = new double[RANDOM_SUM];
        boolean re = true;
        double[][] totalResult = new double[RANDOM_SUM][totalRes.size() / RANDOM_SUM];
        for (int i = 0; i < RANDOM_SUM; i++) {
            for (int j = 0; j < (totalRes.size() / RANDOM_SUM); j++) {
                totalResult[i][j] = totalRes.get(i * (totalRes.size() / RANDOM_SUM) + j);
            }
        }
        for (int i = 0; i < (totalRes.size() / RANDOM_SUM); i += 2) {
            for (int j = 0; j < RANDOM_SUM; j++) {
                pValueList[j] = totalResult[j][i];
                qValueList[j] = totalResult[j][i + 1];
            }
            int pCount = Utils.getPassOfP(pValueList);
            boolean pRes = Utils.statPValue(pValueList);
            double qStat = Double.parseDouble(String.format("%.6f", Utils.getStatQValue(qValueList)));
            boolean qRes = Utils.statQValue(qValueList);
            int num = (i / 2) + 1;
            System.out.println("(" + num + ")、" + sampleSchemeStr[i / 2] + "\tresult: P(" + pCount + "," + pRes + ") , Q(" + qStat + "," + qRes + ")");
            if (!pRes || !qRes) {
                re = false;
            }
        }
        System.out.println("Spend time : " + (System.currentTimeMillis() - startTime) + "ms");
        if (re) {
            System.out.println("All test Passed!");
        } else {
            System.out.println("All test Failed!");
        }
    }

    public static void main(String[] args) {
        System.out.println("RCTools V1.0.1");
        System.out.println("Usage:");
        System.out.println("java -cp RCTool.jar com.ndsec.RandomTest_A <path includes 1000 [125000B , 131072B] files>");
        if (args.length < 1) {
            System.err.println("请命令行输入待检测文件所在文件夹的路径");
            System.exit(-1);
        }
        File file = new File(args[0]);
        if (!file.isDirectory()) {
            System.err.println("请命令行输入待检测文件所在文件夹的路径");
            System.exit(-2);
        }
        File[] files = Utils.getFiles(args[0]);
        if (Objects.requireNonNull(files).length < 1000) {
            System.err.println("指定路径下文件的数量少于1000");
            System.exit(-3);
        }
        for (File file1 : files) {
            if (file1.length() < 125000 || file1.length() > 131072) {
                System.err.println("指定文件夹下存在文件：" + file1.getName() + "， 其大小为：" + file1.length() + "，不再125000B到131072B之间");
            }
        }
        try {
            threadSampleTest(args[0]);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
}
