package com.ndsec.task;

import com.ndsec.random.RandomUtil;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.function.Consumer;

public class sampleTask implements Callable<ArrayList<Double>> {
    List<Boolean> eps;

    String fileName;
    Consumer<String> callback;

    public sampleTask(List<Boolean> eps, String fileName, Consumer<String> callback) {
        this.eps = eps;
        this.fileName = fileName;
        this.callback = callback;
    }

    @Override
    public ArrayList<Double> call() throws Exception {
        ArrayList<Double> res = new ArrayList<>();
        RandomUtil.sampleCheck(eps, res);
        callback.accept(fileName);
        return res;
    }
}
