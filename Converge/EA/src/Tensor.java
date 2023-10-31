import java.util.*;

public class Tensor {


    int mass; //Tensor的总质量

    int[] size; //Tensor各个维度的大小


    int[][] attVals; //保存每一条记录每个维度的信息

    double[][] spaceVals; //保存每一条维度的空间信息

    int[] value; //保存每一条记录的count

    List<int[]> coreDimensionMass;

    int[] timeMass;

    int[] timeSpaceMass;

    int[] latitudeMass;

    int[] longitudeMass;


    boolean[] timeSpaceExisted;    //存储存在的时空块序号

    List<Map<Integer, List<Integer>>> selectDataByCoreDimension;

    Map<Integer, List<Integer>> selectDataByTime; //通过时间维度找到count不为0的数据项

    Map<Integer, List<Integer>> selectDataByTimeSpace; //通过时空块索引找到属于该时空块的数据

    public Tensor() {
    }

    public Tensor(int mass, int[] size, int[][] attVals, double[][] spaceVals, int[] value, List<int[]> coreDimensionMass, int[] timeMass, int[] timeSpaceMass, int[] latitudeMass, int[] longitudeMass,  boolean[] timeSpaceExisted, List<Map<Integer, List<Integer>>> selectDataByCoreDimension, Map<Integer, List<Integer>> selectDataByTime, Map<Integer, List<Integer>> selectDataByTimeSpace) {
        this.mass = mass;
        this.size = size;
        this.attVals = attVals;
        this.spaceVals = spaceVals;
        this.value = value;
        this.coreDimensionMass = coreDimensionMass;
        this.timeMass = timeMass;
        this.timeSpaceMass = timeSpaceMass;
        this.latitudeMass = latitudeMass;
        this.longitudeMass = longitudeMass;
        this.timeSpaceExisted = timeSpaceExisted;
        this.selectDataByCoreDimension = selectDataByCoreDimension;
        this.selectDataByTime = selectDataByTime;
        this.selectDataByTimeSpace = selectDataByTimeSpace;
    }

    public List<int[]> getCoreDimensionMass() {
        return coreDimensionMass;
    }

    public int[] getTimeSpaceMass() {
        return timeSpaceMass;
    }

    public void setTimeSpaceMass(int[] timeSpaceMass) {
        this.timeSpaceMass = timeSpaceMass;
    }

    public boolean[] getTimeSpaceExisted() {
        return timeSpaceExisted;
    }

    public Map<Integer, List<Integer>> getSelectDataByTimeSpace() {
        return selectDataByTimeSpace;
    }

    public int[] getLatitudeMass() {
        return latitudeMass;
    }

    public void setLatitudeMass(int[] latitudeMass) {
        this.latitudeMass = latitudeMass;
    }

    public int[] getlongitudeMass() {
        return longitudeMass;
    }

    public void setlongitudeMass(int[] longitudeMass) {
        this.longitudeMass = longitudeMass;
    }

    public double[][] getSpaceVals() {
        return spaceVals;
    }

    public int[] getTimeMass() {
        return timeMass;
    }

    public void setTimeMass(int[] timeMass) {
        this.timeMass = timeMass;
    }


    public int getMass() {
        return mass;
    }

    public void setMass(int mass) {
        this.mass = mass;
    }

    public int[] getSize() {
        return size;
    }

    public void setSize(int[] size) {
        this.size = size;
    }

    public int[][] getAttVals() {
        return attVals;
    }

    public int[] getValue() {
        return value;
    }

}
