package com.ndsec.common;

import com.ndsec.random.cephes;

import java.io.*;
import java.nio.file.Files;
import java.util.List;

public class Utils {
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

    public static void readBinaryFile(File file, List<Boolean> eps) throws IOException {
        if (file.isFile() && file.exists()) {
            BufferedInputStream in = new BufferedInputStream(new FileInputStream(file));
            ByteArrayOutputStream out = new ByteArrayOutputStream(1024);
            byte[] temp = new byte[1024];
            int size = 0;
            while ((size = in.read(temp)) != -1) {
                out.write(temp, 0, size);
            }
            in.close();
            byte[] content = out.toByteArray();
            for (byte b : content) {
                if (eps.size() == 1000000) break;
                for (int j = 7; j >= 0; j--) {
                    if ((byte) ((b >> j) & 0x1) == 0x1) {
                        eps.add(true);
                    } else {
                        eps.add(false);
                    }
                }
            }
        }
    }

    public static boolean statPValue(double[] p) {
        int count = 0;
        for (double v : p) {
            if (v >= 0.01) count++;
        }
        return count / (p.length * 1.0) >= 0.981;
    }

    public static int getPassOfP(double[] p) {
        int count = 0;
        for (double v : p) {
            if (v >= 0.01) count++;
        }
        return count;
    }

    public static double getStatQValue(double[] q) {
        int[] f = new int[12];
        for (double value : q) {
            if (value >= 0 && value < 0.1) {
                f[1]++;
            }
            if (value >= 0.1 && value < 0.2) {
                f[2]++;
            }
            if (value >= 0.2 && value < 0.3) {
                f[3]++;
            }
            if (value >= 0.3 && value < 0.4) {
                f[4]++;
            }
            if (value >= 0.4 && value < 0.5) {
                f[5]++;
            }
            if (value >= 0.5 && value < 0.6) {
                f[6]++;
            }
            if (value >= 0.6 && value < 0.7) {
                f[7]++;
            }
            if (value >= 0.7 && value < 0.8) {
                f[8]++;
            }
            if (value >= 0.8 && value < 0.9) {
                f[9]++;
            }
            if (value >= 0.9 && value <= 1.0) {
                f[10]++;
            }
        }
        double v = 0.0;
        for (int i = 1; i < 11; i++) {
            v += Math.pow(f[i] - q.length / 10.0, 2) / (q.length / 10.0);
        }
        return cephes.igamc(4.5, v / 2);
    }

    public static boolean statQValue(double[] q) {
        return getStatQValue(q) >= 0.0001;
    }

    public static File[] getFiles(String dirPath) {
        File file = new File(dirPath);
        if (file.isDirectory()) {
            return file.listFiles();
        }
        return null;
    }

    public static boolean read_file_bits(String filePath, long pos1, long pos2, List<Boolean> eps) {
        if (pos1 >= pos2) {
            System.err.println("pos1 err!");
            return false;
        }
        if (pos2 > new File(filePath).length()) {
            System.err.println("pos2 err!");
            return false;
        }
        try (RandomAccessFile file = new RandomAccessFile(filePath, "r")) {
            file.seek(pos1);
            for (long pos = pos1; pos < pos2; pos++) {
                byte b = file.readByte();
                for (int i = 7; i >= 0; i--) {
                    boolean bit = (b & (1 << i)) != 0;
                    if (eps.size() < 1000000) eps.add(bit);
                }
            }
            if (eps.size() == 1000000) return true;
        } catch (Exception e) {
            e.printStackTrace();
        }
        return false;
    }
}
