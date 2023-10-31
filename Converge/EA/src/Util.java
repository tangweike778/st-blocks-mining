import javafx.util.Pair;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class Util {

    public static double sigmoid_base = 0;

    public static void setSigmoid_base(double sigmoid_base) {
        Util.sigmoid_base = sigmoid_base;
    }

    /**
     * 一维数组相减操作函数
     * @param arr1
     * @param arr2
     * @return
     */
    public static int[] arrSubstract(int[] arr1, int[] arr2, boolean[] bools) throws Exception {
        // 检查两个数组的维度是否相同
        if (arr1.length != arr2.length) {
            System.out.println("两个数组的维度不同，无法进行减法运算。");
            throw new Exception();
        }

        int[] result = arr1.clone();
        for (int i = 0; i < arr1.length; i++) {
            if(bools[i]){
                result[i] -= arr2[i];
            }
        }

        return result;
    }


    /**
     * 统计每一个维度的大小
     * @param file_path：输入文件路径
     * @param splitStr：数据文件分隔符
     * @return 三个维度的大小和记录的条数
     */
    public static int[] computerSize(String file_path, String splitStr) throws Exception {
        BufferedReader reader = null;
        int dimensions;
        //提前读取一行，判断其维度
        try {
            reader = new BufferedReader(new FileReader(file_path));
            dimensions = reader.readLine().split(",").length;
        } catch (Exception e) {
            throw new Exception(e);
        } finally {
            if(reader != null){
                reader.close();
            }
        }

        //统计三个维度的大小和记录的条数
        int[] size = new int[dimensions];
        String[] attributes;
        List<Set<String>> dimensionList = new ArrayList<>();
        for(int i = 0; i < dimensions - 1; i++){
            dimensionList.add(new HashSet<>());
        }
        String line;
        try {
            reader = new BufferedReader(new FileReader(file_path));
            while((line = reader.readLine()) != null){
                //统计记录条数
                size[dimensions - 1]++;
                attributes = line.split(splitStr);
                for(int i = 0; i < dimensions - 1; i++){
                    dimensionList.get(i).add(attributes[i]);
                    //统计最大最小时间
                    Converge.globalMinTime = Math.min(Integer.parseInt(attributes[dimensions - 4]), Converge.globalMinTime);
                    Converge.globalMaxTime = Math.max(Integer.parseInt(attributes[dimensions - 4]), Converge.globalMaxTime);
                    Converge.globalMinLatitude = Math.min(Double.parseDouble(attributes[dimensions - 3]), Converge.globalMinLatitude);
                    Converge.globalMaxLatitude = Math.max(Double.parseDouble(attributes[dimensions - 3]), Converge.globalMaxLatitude);
                    Converge.globalMinlongitude = Math.min(Double.parseDouble(attributes[dimensions - 2]), Converge.globalMinlongitude);
                    Converge.globalMaxlongitude = Math.max(Double.parseDouble(attributes[dimensions - 2]), Converge.globalMaxlongitude);
                }
            }
            for(int i = 0; i < dimensions - 1; i++){
                size[i] = dimensionList.get(i).size();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }finally {
            try {
                if(reader != null){
                    reader.close();
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        return size;

    }

    /**
     * 通过经纬度计算空间块索引
     * 输入的是原始的经纬度
     * @param latitude
     * @param longitude
     * @return
     */
    public static int computeSpaceIndex(double latitude, double longitude){
        return computeDirectSpaceIndex(transformLatitude(latitude), transformLongitude(longitude));
    }

    /**
     * 通过时间和经纬度计算时空块索引
     * @param time
     * @param latitude
     * @param longitude
     * @return
     */
    public static int computeTimeSpaceIndex(int time, double latitude, double longitude){
        return transformTime(time) * Converge.PARTITION_NUMBER * Converge.PARTITION_NUMBER + computeSpaceIndex(latitude, longitude);
    }

    /**
     * 通过时间和经纬度计算时空块索引
     * @param time
     * @param latitude
     * @param longitude
     * @return
     */
    public static int computeTimeSpaceIndexDirectly(int time, int latitude, int longitude){
        return time * Converge.PARTITION_NUMBER * Converge.PARTITION_NUMBER + computeDirectSpaceIndex(latitude, longitude);
    }

    /**
     * 通过经纬度计算空间块索引
     * 输入的是处理过的经纬度， 是整数
     * @param latitude
     * @param longitude
     * @return
     */
    public static int computeDirectSpaceIndex(int latitude, int longitude){
        return latitude * Converge.PARTITION_NUMBER + longitude;
    }



    public static double computeMeasurement(double density, double temporal, double spatial){
        double sigmoid = 1.0 / (1 + Math.exp(-density + sigmoid_base));
        return sigmoid * Converge.DENSITY_PROPORTION + (1 - temporal) * Converge.TEMPORAL_PROPORTION + (1 - spatial) * Converge.SPATIAL_PROPORTION;
    }


    /**
     * 找到某一个维度中存在状态为true的属性
     * @param arr
     * @return
     */
    public static List<Integer> findTrueIndexes(boolean[] arr) {
        // 创建存储 true 索引下标的数组
        List<Integer> attrList = new ArrayList<>();

        // 遍历数组，将 true 的索引下标存入 trueIndexes 数组中
        for (int i = 0; i < arr.length; i++) {
            if (arr[i]) {
                attrList.add(i);
            }
        }

        return attrList;
    }

    /**
     * 计算子张量密度
     *
     * @return
     */
    public static double computeDensity(long mass, long size){
        // return size == 0 ? 0 : mass / (size * 1.0) * (Converge.CORE_COLUMN_LENGTH);
        //TODO：修改了度量
        return size == 0 ? 0 : mass / (size * 1.0) * (Converge.globalSize.length - 1);
        // return size == 0 ? 0 : mass ;
    }

    public static double computeDensity2(long mass, long size){
        // return size == 0 ? 0 : mass / (size * 1.0) * (Converge.CORE_COLUMN_LENGTH);
        //TODO：修改了度量
        return size == 0 ? 0 : mass / (size * 1.0) * (Converge.CORE_COLUMN_LENGTH);
        // return size == 0 ? 0 : mass ;
    }

    public static double computeDensity(SubTensor subTensor){
        // return subTensor.getMass() * (float)(Converge.CORE_COLUMN_LENGTH) / computeVolume(subTensor);
        //TODO：修改了度量
        return subTensor.getMass() * (float)(Converge.globalSize.length - 1) / computeVolume(subTensor);
    }

    public static double computeDensity2(SubTensor subTensor){
        // return subTensor.getMass() * (float)(Converge.CORE_COLUMN_LENGTH) / computeVolume(subTensor);
        //TODO：修改了度量
        return subTensor.getMass() * (float)(Converge.CORE_COLUMN_LENGTH) / computeVolume2(subTensor);
    }

    public static double computeTotalDensity(SubTensor subTensor){
        // double density = subTensor.getMass() * (float)(Converge.CORE_COLUMN_LENGTH) / computeVolume(subTensor);
        //TODO：修改了度量
        double density = subTensor.getMass() * (float)(Converge.globalSize.length - 1) / computeVolume(subTensor);
        // double density = subTensor.getMass();
        double temporalRange = computeTemporalRange(subTensor.getMaxTime(), subTensor.getMinTime());
        double spatialRange = computeOfSpatialRange(subTensor.getMinLatitude(), subTensor.getMaxLatitude(), subTensor.getMinLongitude(), subTensor.getMaxLongitude());
        return computeMeasurement(density, temporalRange, spatialRange);
    }

/*    public static double computeWeightDensity(SubTensor subTensor){
        int[] size = subTensor.getSize();
        double coreMass = 0;
        for(int i = 0; i < Converge.CORE_COLUMN_LENGTH; i++){
            coreMass += size[i];
        }
        double core = coreMass * Converge.CORE_PROPORTION;
        double time = (subTensor.getMaxTime() - subTensor.getMinTime() + 1) * Converge.timeStep * Converge.TIME_PROPORTION;
        double space = (((subTensor.getMaxLatitude() - subTensor.getMinLatitude() + 1) * Converge.latitudeStep * 0.5)
                + ((subTensor.getMaxLongitude() - subTensor.getMinLongitude() + 1) * Converge.longitudeStep * 0.5)) * Converge.SPACE_PROPORTION;
        //TODO:更改为mass/core + mass / core + mass / space
        return Util.computeWeightDensity(subTensor.getMass(), core, time, space);
    }*/

    /*public static double computeWeightDensity(long mass, double core, double time, double space){
        return mass * 1.0 / core + mass * 1.0 / time + mass * 1.0 / space;
        // return mass * 1.0 / (core + time + space);
    }*/

    /**
     * 根据要划分的份数h，将time转换为对应的序号
     * @param time
     * @return
     */
    public static int transformTime(int time){
        if(time == Converge.globalMinTime)return 1;
        return (int) Math.ceil((time - Converge.globalMinTime) * 1.0 / Converge.timeStep);
    }

    public static int transformLatitude(double latitude){
        if(latitude == Converge.globalMinLatitude)return 1;
        return (int) Math.ceil((latitude - Converge.globalMinLatitude) / Converge.latitudeStep);
    }

    public static int transformLongitude(double longitude){
        if(longitude == Converge.globalMinlongitude)return 1;
        return (int) Math.ceil((longitude - Converge.globalMinlongitude) / Converge.longitudeStep);
    }


    public static double computeTemporalRange(int maxTime, int minTime) {
        return  (maxTime - minTime + 1) * 1.0 / (Converge.globalMaxTime - Converge.globalMinTime + 1);
    }

    //TODO:计算范围需要修改
    public static double computeOfSpatialRange(double minLatitude, double maxLatitude, double minLongitude, double maxLongitude) {
        double localArea = (maxLatitude - minLatitude + 1) * (maxLongitude - minLongitude + 1);
        double globalArea = (Converge.globalMaxLatitude - Converge.globalMinLatitude + 1)
                * (Converge.globalMaxlongitude - Converge.globalMinlongitude + 1);
        return localArea / globalArea;
    }

    /**
     * 计算子张量的体积，也可称为size
     * @param subTensor 子张量
     * @return 子张量体积
     */
    public static long computeVolume(SubTensor subTensor){
        int[] size = subTensor.getSize();
        int coreSize = 0;
        for(int i = 0; i < Converge.globalSize.length - 1; i++){
            coreSize += size[i];
        }
        //TODO:修改了度量
        /*for(int i = 0; i < Converge.CORE_COLUMN_LENGTH; i++){
            coreSize += size[i];
        }*/
        return coreSize;
    }

    public static long computeVolume2(SubTensor subTensor){
        int[] size = subTensor.getSize();
        int coreSize = 0;
        for(int i = 0; i < Converge.CORE_COLUMN_LENGTH; i++){
            coreSize += size[i];
        }
        //TODO:修改了度量
        /*for(int i = 0; i < Converge.CORE_COLUMN_LENGTH; i++){
            coreSize += size[i];
        }*/
        return coreSize;
    }
}
