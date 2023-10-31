import java.io.Serializable;
import java.util.*;

public class SubTensor implements Serializable {

    private int mass; //质量信息

    private int minTime; //最小的时间戳

    private int maxTime; // 最大的时间戳

    private double minLongitude; //最小的经度

    private double maxLongitude; //最大的经度

    private double minLatitude; //最小的纬度

    private double maxLatitude; //最大的纬度

    private int[] size; //每一个维度的大小

    private boolean[] timeExisted; //时间维度索引

    private List<boolean[]> coreDimensionExisted;

    boolean[] timeSpaceExisted;    //存储存在的时空块序号

    private boolean[] latitudeExisted; //纬度维度索引

    private boolean[] longitudeExisted; //经度维度索引

    List<int[]> coreDimensionMass; //所有核心维度的质量

    int[] timeMass; //时间维度质量大小

    int[] timeSpaceMass;

    int[] latitudeMass; //纬度维度质量大小

    int[] longitudeMass; //经度维度质量大小

    private final int[][] attVals; //保存每条记录每一个维度的信息

    private double[][] spaceVals; //保存每条记录的空间信息

    private final int[] value; //保存每一条记录的count

    private final List<Map<Integer, List<Integer>>> selectDataByCoreDimension;

    private Map<Integer, List<Integer>> selectDataByTime; //通过时间维度找到count不为0的数据在value和attVals的索引

    private Map<Integer, List<Integer>> selectDataByLatitude; //通过时间维度找到count不为0的数据在value和attVals的索引

    private Map<Integer, List<Integer>> selectDataByLongitude; //通过时间维度找到count不为0的数据在value和attVals的索引

    public SubTensor(SubTensor subTensor){
        this.mass = subTensor.getMass();
        this.size = subTensor.getSize().clone();
        this.attVals = subTensor.getAttVals().clone();
        this.selectDataByCoreDimension = subTensor.getSelectDataByCoreDimension();
        this.coreDimensionMass = new ArrayList<>();
        List<int[]> subTensorCoreDimensionMass = subTensor.getCoreDimensionMass();
        for (int[] masses : subTensorCoreDimensionMass) {
            coreDimensionMass.add(masses.clone());
        }
        this.value = subTensor.getValue().clone();
        this.minTime = subTensor.getMinTime();
        this.maxTime = subTensor.getMaxTime();
    }

    public SubTensor(int mass, int minTime, int maxTime, double minLongitude, double maxLongitude, double minLatitude, double maxLatitude, int[] size, boolean[] timeExisted, List<boolean[]> coreDimensionExisted, boolean[] timeSpaceExisted, boolean[] latitudeExisted, boolean[] longitudeExisted, List<int[]> coreDimensionMass, int[] timeMass, int[] timeSpaceMass, int[] latitudeMass, int[] longitudeMass, int[][] attVals, double[][] spaceVals, int[] value, List<Map<Integer, List<Integer>>> selectDataByCoreDimension, Map<Integer, List<Integer>> selectDataByTime, Map<Integer, List<Integer>> selectDataByLatitude, Map<Integer, List<Integer>> selectDataByLongitude) {
        this.mass = mass;
        this.minTime = minTime;
        this.maxTime = maxTime;
        this.minLongitude = minLongitude;
        this.maxLongitude = maxLongitude;
        this.minLatitude = minLatitude;
        this.maxLatitude = maxLatitude;
        this.size = size;
        this.timeExisted = timeExisted;
        this.coreDimensionExisted = coreDimensionExisted;
        this.timeSpaceExisted = timeSpaceExisted;
        this.latitudeExisted = latitudeExisted;
        this.longitudeExisted = longitudeExisted;
        this.coreDimensionMass = coreDimensionMass;
        this.timeMass = timeMass;
        this.timeSpaceMass = timeSpaceMass;
        this.latitudeMass = latitudeMass;
        this.longitudeMass = longitudeMass;
        this.attVals = attVals;
        this.spaceVals = spaceVals;
        this.value = value;
        this.selectDataByCoreDimension = selectDataByCoreDimension;
        this.selectDataByTime = selectDataByTime;
        this.selectDataByLatitude = selectDataByLatitude;
        this.selectDataByLongitude = selectDataByLongitude;
    }

    public int[] getLatitudeMass() {
        return latitudeMass;
    }


    public int[] getlongitudeMass() {
        return longitudeMass;
    }

    public int getMass() {
        return mass;
    }

    public double getMaxLongitude() {
        return maxLongitude;
    }

    public void setMass(int mass) {
        this.mass = mass;
    }

    public int getMinTime() {
        return minTime;
    }

    public int getMaxTime() {
        return maxTime;
    }


    public double getMinLongitude() {
        return minLongitude;
    }

    public double getMaxlongitude() {
        return maxLongitude;
    }

    public int[] getTimeSpaceMass() {
        return timeSpaceMass;
    }

    public double getMinLatitude() {
        return minLatitude;
    }

    public double getMaxLatitude() {
        return maxLatitude;
    }

    public int[] getSize() {
        return size;
    }

    public void setSize(int[] size) {
        this.size = size;
    }

    public boolean[] getTimeExisted() {
        return timeExisted;
    }

    public boolean[] getLatitudeExisted() {
        return latitudeExisted;
    }

    public boolean[] getLongitudeExisted() {
        return longitudeExisted;
    }

    public List<boolean[]> getCoreDimensionExisted() {
        return coreDimensionExisted;
    }
    public List<int[]> getCoreDimensionMass() {
        return coreDimensionMass;
    }

    public int[] getTimeMass() {
        return timeMass;
    }

    public int[][] getAttVals() {
        return attVals;
    }

    public double[][] getSpaceVals() {
        return spaceVals;
    }

    public Map<Integer, List<Integer>> getSelectDataByLatitude() {
        return selectDataByLatitude;
    }

    public Map<Integer, List<Integer>> getSelectDataByLongitude() {
        return selectDataByLongitude;
    }

    public int[] getValue() {
        return value;
    }

    public List<Map<Integer, List<Integer>>> getSelectDataByCoreDimension() {
        return selectDataByCoreDimension;
    }

    public Map<Integer, List<Integer>> getSelectDataByTime() {
        return selectDataByTime;
    }



    /**
     * 模拟将传入的数据集全部添加后的密度
     * 添加的时候只模拟time、user、object三个维度的size的增加，空间维度不管
     * @param tensor : 全局张量
     * @param dataList : 需要添加到Subtensor中的数据集
     */
    public double imitateAppendDatas(Tensor tensor, List<Integer> dataList) {
        //模拟需要的数据：mass、size、所有的Existed、
        int mass = getMass();
        int[] size = getSize().clone();
        boolean[] timeExisted = getTimeExisted().clone();
        List<boolean[]> coreExistedCopy = new ArrayList<>();
        for (boolean[] coreExisted : coreDimensionExisted) {
            coreExistedCopy.add(coreExisted.clone());
        }

        int tempMinTime = minTime;
        int tempMaxTime = maxTime;
        double tempMinLatitude = minLatitude;
        double tempMaxLatitude = maxLatitude;
        double tempMinLongitude = minLongitude;
        double tempMaxLongitude = maxLongitude;

        int[] value = tensor.getValue();
        int[][] attVals = tensor.getAttVals();
        for(Integer data : dataList){
            if(value[data] == 0)continue;
            mass += value[data];
            for(int i = 0; i < Converge.globalSize.length - 3; i++){
                int attr = attVals[i][data];
                if(i < Converge.CORE_COLUMN_LENGTH){
                    boolean[] existed = coreExistedCopy.get(i);
                    if(!existed[attr]){
                        size[i]++;
                        existed[attr] = true;
                    }
                }else if(i == Converge.TIME_COLUMN){
                    if(!timeExisted[Util.transformTime(attr)]){
                        size[i]++;
                        timeExisted[Util.transformTime(attr)] = true;
                    }
                }


            }
            tempMinTime = Math.min(tempMinTime, attVals[Converge.TIME_COLUMN][data]);
            tempMaxTime = Math.max(tempMaxTime, attVals[Converge.TIME_COLUMN][data]);
            tempMinLatitude = Math.min(tempMinLatitude, spaceVals[0][data]);
            tempMaxLatitude = Math.max(tempMaxLatitude, spaceVals[0][data]);
            tempMinLongitude = Math.min(tempMinLongitude, spaceVals[1][data]);
            tempMaxLongitude = Math.max(tempMaxLongitude, spaceVals[1][data]);

        }

        double temporalRange = Util.computeTemporalRange(tempMaxTime, tempMinTime);
        double spatialRange = Util.computeOfSpatialRange(tempMinLatitude, tempMaxLatitude, tempMinLongitude, tempMaxLongitude);

        /*long coreSize = 0;
        for(int i = 0; i < Converge.CORE_COLUMN_LENGTH; i++){
            coreSize += size[i];
        }*/
        long totalSize = 0;
        for(int i = 0; i < Converge.globalSize.length - 1; i++){
            totalSize += size[i];
        }

        // return Util.computeDensity(mass, coreSize);
        return Util.computeMeasurement(Util.computeDensity(mass, totalSize), temporalRange, spatialRange);
    }

    /**
     * 向子张量中添加数据
     * @param tensor
     * @param dataList 要添加的数据
     */
    public void appendDatas(Tensor tensor, List<Integer> dataList) {
        int[] dataValues = tensor.getValue();
        for(Integer data : dataList){
            if(dataValues[data] == 0)continue;
            int currentCount = dataValues[data];
            //更新mass
            mass += currentCount;
            value[data] = currentCount;
            for (int j = 0; j < Converge.globalSize.length - 3; j++){
                int attribute = attVals[j][data];
                if (j < Converge.CORE_COLUMN_LENGTH){
                    int[] masses = coreDimensionMass.get(j);
                    Map<Integer, List<Integer>> selectDataByCore = selectDataByCoreDimension.get(j);
                    masses[attribute] += currentCount;
                    if(!selectDataByCore.containsKey(attribute)){
                        selectDataByCore.put(attribute, new ArrayList<>());
                    }
                    selectDataByCore.get(attribute).add(data);
                    boolean[] coreExisted = coreDimensionExisted.get(j);
                    if(!coreExisted[attribute]){
                        coreExisted[attribute] = true;
                        size[j]++;
                    }
                } else if(j == Converge.TIME_COLUMN) {
                    int processedTime = Util.transformTime(attribute);
                    timeMass[processedTime] += currentCount;
                    minTime = Math.min(attribute, minTime);
                    maxTime = Math.max(attribute, maxTime);
                    if(!selectDataByTime.containsKey(processedTime)){
                        selectDataByTime.put(processedTime, new ArrayList<>());
                    }
                    selectDataByTime.get(processedTime).add(data);
                    if(!timeExisted[processedTime]){
                        timeExisted[processedTime] = true;
                        size[j]++;
                    }
                }
            }

            int latitude = Util.transformLatitude(spaceVals[0][data]);
            int longitude = Util.transformLongitude(spaceVals[1][data]);
            minLatitude = Math.min(spaceVals[0][data], minLatitude);
            maxLatitude = Math.max(spaceVals[0][data], maxLatitude);
            minLongitude = Math.min(spaceVals[1][data], minLongitude);
            maxLongitude = Math.max(spaceVals[1][data], maxLongitude);
            latitudeMass[latitude] += currentCount;
            longitudeMass[longitude] += currentCount;

            if(!selectDataByLatitude.containsKey(latitude)){
                selectDataByLatitude.put(latitude, new ArrayList<>());
            }
            selectDataByLatitude.get(latitude).add(data);
            if(!selectDataByLongitude.containsKey(longitude)){
                selectDataByLongitude.put(longitude, new ArrayList<>());
            }
            selectDataByLongitude.get(longitude).add(data);

        }

        //更新空间维度的size和存在状态
        for(int k = Util.transformLatitude(minLatitude); k <= Util.transformLatitude(maxLatitude); k++){
            if(!latitudeExisted[k]){
                latitudeExisted[k] = true;
                size[Converge.LATITUDE_COLUMN]++;
            }
        }
        for(int k = Util.transformLongitude(minLongitude); k <= Util.transformLongitude(maxLongitude); k++){
            if(!longitudeExisted[k]){
                longitudeExisted[k] = true;
                size[Converge.LONGITUDE_COLUMN]++;
            }
        }

    }

    /**
     * 更新极值信息
     */
    public void updatePeakInfomation() {
        maxTime = Integer.MIN_VALUE;
        minTime = Integer.MAX_VALUE;
        maxLatitude = -Double.MAX_VALUE;
        minLatitude = Double.MAX_VALUE;
        maxLongitude = -Double.MAX_VALUE;
        minLongitude = Double.MAX_VALUE;
        for (int i = 0; i < value.length; i++) {
            if(value[i] != 0){
                maxTime = Math.max(attVals[Converge.TIME_COLUMN][i], maxTime);
                minTime = Math.min(attVals[Converge.TIME_COLUMN][i], minTime);
                maxLatitude = Math.max(spaceVals[0][i], maxLatitude);
                minLatitude = Math.min(spaceVals[0][i], minLatitude);
                maxLongitude = Math.max(spaceVals[1][i], maxLongitude);
                minLongitude = Math.min(spaceVals[1][i], minLongitude);
            }
        }
    }
}
