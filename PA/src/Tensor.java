import com.sun.jndi.toolkit.url.Uri;

import java.io.Serializable;
import java.util.*;

public class Tensor implements Cloneable{

    private int minTime;
    private int maxTime;
    private int minLatitude;
    private int maxLatitude;
    private int minLongitude;
    private int maxLongitude;
    private int mass;  //tuple数量
    private int[] size; //tensor各个维度的大小
    private final int[][] attVals; //保存每一条记录每个维度的信息
    private final int[][] spaceVals; //保存每一条记录空间经纬度信息
    private int[] value; //保存每一条记录的count

    //注：每个维度k号索引对应的是k+1属性值
    private List<List<Slice>> slices; //保存每一个维度的slice
    private List<Map<Integer, Integer>> dataOfAttr; //保存每个维度中每个属性的记录数量


    @Override
    public Tensor clone() throws CloneNotSupportedException {
        Tensor tensor = (Tensor) super.clone();

        tensor.setSize(size.clone());
        tensor.setValue(value.clone());

        //attVals
        int[][] attVals1 = tensor.getAttVals();
        for (int i = 0; i < attVals.length; i++) {
            attVals1[i] = attVals[i].clone();
        }

        //spaceVals
        int[][] spaceVals1 = tensor.getSpaceVals();
        for (int i = 0; i < spaceVals.length; i++) {
            spaceVals1[i] = spaceVals[i].clone();
        }

        //slices
        List<List<Slice>> slicesCopy = new ArrayList<>();
        for (List<Slice> sliceList : slices) {
            List<Slice> sliceListCopy = new ArrayList<>();
            for (Slice s : sliceList) {
                sliceListCopy.add(s.clone());
            }
            slicesCopy.add(sliceListCopy);
        }
        tensor.setSlices(slicesCopy);

        //dataOfAttr
        List<Map<Integer, Integer>> dataOfAttrCopy = new ArrayList<>();
        for (Map<Integer, Integer> map : dataOfAttr) {
            Map<Integer, Integer> mapCopy = new HashMap<>();
            mapCopy.putAll(map);
            dataOfAttrCopy.add(mapCopy);
        }
        tensor.setDataOfAttr(dataOfAttrCopy);

        tensor.updatePeakInfomation();


        return tensor;
    }

    /**
     * 更新极值信息
     */
    public void updatePeakInfomation() {
        maxTime = Integer.MIN_VALUE;
        minTime = Integer.MAX_VALUE;
        maxLatitude = Integer.MIN_VALUE;
        minLatitude = Integer.MAX_VALUE;
        maxLongitude = Integer.MIN_VALUE;
        minLongitude = Integer.MAX_VALUE;

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

    @Override
    public String toString() {
        return "Tensor{" +
                "minTime=" + minTime +
                ", maxTime=" + maxTime +
                ", minLatitude=" + minLatitude +
                ", maxLatitude=" + maxLatitude +
                ", minLongitude=" + minLongitude +
                ", maxLongitude=" + maxLongitude +
                ", mass=" + mass +
                ", size=" + Arrays.toString(size) +
                ", attVals=" + Arrays.toString(attVals) +
                ", spaceVals=" + Arrays.toString(spaceVals) +
                ", value=" + Arrays.toString(value) +
                ", slices=" + slices +
                ", dataOfAttr=" + dataOfAttr +
                '}';
    }

    public Tensor(Tensor tensor) {
        this.mass = tensor.getMass();
        // this.size = Arrays.copyOfRange(Converge.globalSize, 0, Converge.globalSize.length - 1);
        this.size = tensor.getSize().clone();
        this.attVals = tensor.getAttVals().clone();
        this.spaceVals = tensor.getSpaceVals().clone();
        this.value = tensor.getValue().clone();
        this.slices = tensor.getSlices();
        this.dataOfAttr = new ArrayList<>();
        for(int i = 0; i < Converge.globalSize.length - 1; i++){
            Map<Integer, Integer> parentDataMap = tensor.getDataOfAttr().get(i);
            HashMap<Integer, Integer> dataMap = new HashMap<>();
            for (Integer key : parentDataMap.keySet()) {
                dataMap.put(key, parentDataMap.get(key));
            }
            this.dataOfAttr.add(dataMap);
        }
    }

    public Tensor(int minTime, int maxTime, int minLatitude, int maxLatitude, int minLongitude, int maxLongitude, int mass, int[] size, int[][] attVals, int[][] spaceVals, int[] value, List<List<Slice>> slices, List<Map<Integer, Integer>> dataOfAttr) {
        this.minTime = minTime;
        this.maxTime = maxTime;
        this.minLatitude = minLatitude;
        this.maxLatitude = maxLatitude;
        this.minLongitude = minLongitude;
        this.maxLongitude = maxLongitude;
        this.mass = mass;
        this.size = size;
        this.attVals = attVals;
        this.spaceVals = spaceVals;
        this.value = value;
        this.slices = slices;
        this.dataOfAttr = dataOfAttr;
    }

    public Tensor createSubTensorBydataList(List<Integer> dataList){
        int minTime = Integer.MAX_VALUE;
        int maxTime = Integer.MIN_VALUE;
        int minLatitude = Integer.MAX_VALUE;
        int maxLatitude = Integer.MIN_VALUE;
        int minLongitude = Integer.MAX_VALUE;
        int maxLongitude = Integer.MIN_VALUE;
        int newMass = 0;
        int[] newSize = new int[Converge.globalSize.length - 1];
        int[][] newAttVals = this.getAttVals();
        int[][] newSpaceVals = this.getSpaceVals();
        //创建一个和全局tensor的value属性一样大小的value，初始化所有值为0
        int[] value = this.getValue();
        int[] newValue = new int[value.length];
        List<Map<Integer, Integer>> newDataOfAttr = new ArrayList<>();
        //初始化dataOfAttr
        for (int i = 0; i < Converge.globalSize.length - 1; i++) {
            newDataOfAttr.add(new HashMap<>());
        }

        //初始化slices
        List<List<Slice>> newSlices = new ArrayList<>();
        //初始化user、item维度slices
        for (int j = 0; j < Converge.CORE_COLUMN_LENGTH; j++) {
            List<Slice> tempSlices = new ArrayList<>();
            for (int k = 0; k < Converge.globalSize[j]; k++) {
                tempSlices.add(new Slice(k + 1, false));
            }
            newSlices.add(tempSlices);
        }

        //初始化时间、经纬度维度的slices
        List<Slice> tempSlices1 = new ArrayList<>();
        for (int k = 0; k < Converge.globalMaxTime; k++) {
            tempSlices1.add(new Slice(k + 1,  false));
        }
        newSlices.add(tempSlices1);
        List<Slice> tempSlices2 = new ArrayList<>();
        for (int k = 0; k < Converge.globalMaxLatitude; k++) {
            tempSlices2.add(new Slice(k + 1,  false));
        }
        newSlices.add(tempSlices2);
        List<Slice> tempSlices3 = new ArrayList<>();
        for (int k = 0; k < Converge.globalMaxlongitude; k++) {
            tempSlices3.add(new Slice(k + 1,  false));
        }
        newSlices.add(tempSlices3);

        for (Integer data : dataList) {
            //修改newValue
            newValue[data] = value[data];
            //统计总质量
            newMass += newValue[data];
            for (int i = 0; i < Converge.CORE_COLUMN_LENGTH; i++) {
                //修改dataOfAttr
                Map<Integer,Integer> datas = newDataOfAttr.get(i);
                Integer tempValue = datas.get(newAttVals[i][data]);
                if(tempValue == null){
                    datas.put(newAttVals[i][data], 0);
                    tempValue = 0;
                    newSize[i]++;
                }
                datas.replace(newAttVals[i][data], tempValue + value[data]);

            }

            //处理时间的dataOfAttr
            Map<Integer, Integer> datas3 = newDataOfAttr.get(Converge.TIME_COLUMN);
            Integer tempValue = datas3.get(newAttVals[Converge.TIME_COLUMN][data]);
            if(tempValue == null){
                datas3.put(newAttVals[Converge.TIME_COLUMN][data], 0);
                tempValue = 0;
                newSize[Converge.TIME_COLUMN]++;
            }
            datas3.replace(newAttVals[Converge.TIME_COLUMN][data], tempValue + value[data]);

            //修改经纬度的dataOfAttr
            Map<Integer, Integer> datas1 = newDataOfAttr.get(Converge.LATITUDE_COLUMN);
            Map<Integer, Integer> datas2 = newDataOfAttr.get(Converge.LONGITUDE_COLUMN);
            //使经度对应的keyValue自增1
            int latitude = (int) newSpaceVals[0][data];
            int longitude = (int) newSpaceVals[1][data];
            Integer latitudeValue = datas1.get(latitude);
            if (latitudeValue == null) {
                datas1.put(latitude, 0);
                latitudeValue = 0;
                newSize[Converge.LATITUDE_COLUMN]++;
            }
            datas1.replace(latitude, latitudeValue + value[data]);

            //使维度对应的keyValue自增1
            Integer longitudeValue = datas2.get(longitude);
            if (longitudeValue == null) {
                datas2.put(longitude, 0);
                longitudeValue = 0;
                newSize[Converge.LONGITUDE_COLUMN]++;
            }
            datas2.replace(longitude, longitudeValue + value[data]);

            //更新user、item、time纬度的slices信息
            for (int mode = 0; mode < Converge.globalSize.length - 1; mode++) {
                int attr = 0;
                //取出当前维度的slices
                List<Slice> tempSlices = newSlices.get(mode);
                if(mode < Converge.CORE_COLUMN_LENGTH){
                    attr = newAttVals[mode][data];
                } else if (mode == Converge.TIME_COLUMN) {
                    attr = newAttVals[mode][data];
                } else if(mode == Converge.LATITUDE_COLUMN){
                    attr = latitude;
                }else if(mode == Converge.LONGITUDE_COLUMN){
                    attr = longitude;
                }
                Slice slice = tempSlices.get(attr - 1);
                slice.setExists(true);

                if(mode == Converge.TIME_COLUMN){
                    minTime = Math.min(minTime, newAttVals[mode][data]);
                    maxTime = Math.max(maxTime, newAttVals[mode][data]);
                }else if(mode == Converge.LATITUDE_COLUMN){
                    minLatitude = Math.min(minLatitude, newSpaceVals[0][data]);
                    maxLatitude = Math.max(maxLatitude, newSpaceVals[0][data]);
                }else if(mode == Converge.LONGITUDE_COLUMN){
                    minLongitude = Math.min(minLongitude, newSpaceVals[1][data]);
                    maxLongitude = Math.max(maxLongitude, newSpaceVals[1][data]);
                }

                //修改切片包含的数据索引
                slice.dataIndexs.add(data);
            }

        }

        if(newMass == 0)return null;
        return new Tensor(minTime, maxTime, minLatitude, maxLatitude, minLongitude, maxLongitude,newMass, newSize, newAttVals, newSpaceVals, newValue, newSlices, newDataOfAttr);
    }

    /**
     * 模拟删除指定切片
     * @param sliceInformation
     * @param index
     * @return
     */
    public boolean preDeleteSlice(int[][] sliceInformation, int index) {
        //根据参数找到预删除的切片
        int attr = sliceInformation[0][index];
        int mode = sliceInformation[1][index];
        Slice slice = slices.get(mode).get(attr - 1);
        int[] valueClone = value.clone();

        //将sllice的存在状态转换为二维数组
        boolean[][] sliceExisted = new boolean[Converge.globalSize.length - 1][];
        for (int i = 0; i < slices.size(); i++) {
            if(i < Converge.CORE_COLUMN_LENGTH){
                sliceExisted[i] = new boolean[Converge.globalSize[i] + 1];
            }else if(i == Converge.TIME_COLUMN){
                sliceExisted[i] = new boolean[Converge.globalMaxTime + 1];
            }else if(i == Converge.LATITUDE_COLUMN){
                sliceExisted[i] = new boolean[Converge.globalMaxLatitude + 1];
            }else if(i == Converge.LONGITUDE_COLUMN){
                sliceExisted[i] = new boolean[Converge.globalMaxlongitude + 1];
            }
            List<Slice> sliceList = slices.get(i);
            for (Slice tempSlice : sliceList) {
                sliceExisted[i][tempSlice.getIndex()] = tempSlice.isExists();
            }
        }

        List<Map<Integer, Integer>> tempDataOfAttr = new ArrayList<>();
        for (Map<Integer, Integer> map : dataOfAttr) {
            Map<Integer, Integer> mapCopy = new HashMap<>();
            mapCopy.putAll(map);
            tempDataOfAttr.add(mapCopy);
        }

        int tempMinTime = minTime;
        int tempMaxTime = maxTime;
        int tempMinLatitude = minLatitude;
        int tempMaxLatitude = maxLatitude;
        int tempMinLongitude = minLongitude;
        int tempMaxLongitude = maxLongitude;

        //拷贝size属性
        int[] tempSize = size.clone();
        //取出切片包含的数据,遍历删除。删除过程中更新size属性和mass属性。
        Vector<Integer> dataIndexs = slice.getDataIndexs();
        int sliceMass = 0;
        for (Integer data : dataIndexs) {
            if(value[data] == 0)continue;
            sliceMass += value[data];
            //根据data查找到每个维度的属性，根据属性更新tempDataOfAttr并更新size。
            for(int i = 0; i < Converge.globalSize.length - 1; i++){
                int tempAttr = 0;
                if(i < Converge.CORE_COLUMN_LENGTH){
                    tempAttr = attVals[i][data];
                } else if (i == Converge.TIME_COLUMN) {
                    tempAttr = attVals[i][data];
                } else if(i == Converge.LATITUDE_COLUMN){
                    tempAttr = spaceVals[0][data];
                }else if(i == Converge.LONGITUDE_COLUMN){
                    tempAttr = spaceVals[1][data];
                }
                int attrValue = tempDataOfAttr.get(i).get(tempAttr);
                tempDataOfAttr.get(i).replace(tempAttr, attrValue - valueClone[data]);
                if(tempDataOfAttr.get(i).get(tempAttr) == 0){
                    //通过切片判断该维度的该属性是否存在，只有存在才修改其size
                    if(sliceExisted[i][tempAttr]){
                        //对于经纬度和时间不需要更改其size
                        tempSize[i]--;
                    }
                    sliceExisted[i][tempAttr] = false;
                }

                //更新最大最小时间和最大最小经纬度和时间
                if(i == Converge.TIME_COLUMN){
                    if(attVals[i][data] == tempMinTime){
                        for(int time = tempMinTime; time <= tempMaxTime; time++){
                            if(sliceExisted[i][time]){
                                tempMinTime = time;
                                break;
                            }
                        }
                    }else if(attVals[i][data] == tempMaxTime){
                        for(int time = tempMaxTime; time >= tempMinTime; time--){
                            if(sliceExisted[i][time]){
                                tempMaxTime = time;
                                break;
                            }
                        }
                    }
                }else if(i == Converge.LATITUDE_COLUMN){
                    if(spaceVals[0][data] == tempMinLatitude){
                        for(int latitude = minLatitude; latitude <= maxLatitude; latitude++){
                            if(sliceExisted[i][latitude]){
                                tempMinLatitude = latitude;
                                break;
                            }
                        }
                    }else if(spaceVals[0][data] == tempMaxLatitude){
                        for(int latitude = maxLatitude; latitude >= minLatitude; latitude--){
                            if(sliceExisted[i][latitude]){
                                tempMaxLatitude = latitude;
                                break;
                            }
                        }
                    }
                }else if(i == Converge.LONGITUDE_COLUMN){
                    if(spaceVals[1][data] == tempMinLongitude){
                        for(int longitude = minLongitude; longitude < maxLongitude; longitude++){
                            if(sliceExisted[i][longitude]){
                                tempMinLongitude = longitude;
                                break;
                            }
                        }
                    }else if(spaceVals[1][data] == tempMaxLongitude){
                        for(int longitude = maxLongitude; longitude >= minLongitude; longitude--){
                            if(sliceExisted[i][longitude]){
                                tempMaxLongitude = longitude;
                                break;
                            }
                        }
                    }
                }
            }
            //将当前data的value设置为0，防止更新极值时更新到该值。
            valueClone[data] = 0;
        }

        double densityBefore = Util.computeDensity(this);
        int coreSize = 0;
        for (int i = 0; i < Converge.CORE_COLUMN_LENGTH; i++) {
            coreSize += tempSize[i];
        }
        double densityAfter = (mass - sliceMass) * (float)Converge.CORE_COLUMN_LENGTH / coreSize;

        double temporalRangeBefore = Util.computeTemporalRange(minTime, maxTime);
        double temporalRangeAfter = Util.computeTemporalRange(tempMinTime, tempMaxTime);

        double spatialRangeBefore = Util.computeSpatialRange(minLatitude, maxLatitude, minLongitude, maxLongitude);
        double spatialRangeAfter = Util.computeSpatialRange(tempMinLatitude, tempMaxLatitude, tempMinLongitude, tempMaxLongitude);

        double measurementBefore = Util.computeMeasurement(densityBefore, temporalRangeBefore, spatialRangeBefore);
        double measurementAfter = Util.computeMeasurement(densityAfter, temporalRangeAfter, spatialRangeAfter);

        return measurementAfter > measurementBefore;
    }

    /**
     * 模拟删除指定切片
     * @return
     */
    public boolean preDeleteSlice(int attr, int mode) {
        if(attr <= 0)return false;
        //根据参数找到预删除的切片
        Slice slice = slices.get(mode).get(attr - 1);
        int[] valueClone = value.clone();

        //将sllice的存在状态转换为二维数组
        boolean[][] sliceExisted = new boolean[Converge.globalSize.length - 1][];
        for (int i = 0; i < slices.size(); i++) {
            if(i < Converge.CORE_COLUMN_LENGTH){
                sliceExisted[i] = new boolean[Converge.globalSize[i] + 1];
            }else if(i == Converge.TIME_COLUMN){
                sliceExisted[i] = new boolean[Converge.globalMaxTime + 1];
            }else if(i == Converge.LATITUDE_COLUMN){
                sliceExisted[i] = new boolean[Converge.globalMaxLatitude + 1];
            }else if(i == Converge.LONGITUDE_COLUMN){
                sliceExisted[i] = new boolean[Converge.globalMaxlongitude + 1];
            }
            List<Slice> sliceList = slices.get(i);
            for (Slice tempSlice : sliceList) {
                sliceExisted[i][tempSlice.getIndex()] = tempSlice.isExists();
            }
        }

        List<Map<Integer, Integer>> tempDataOfAttr = new ArrayList<>();
        for (Map<Integer, Integer> map : dataOfAttr) {
            Map<Integer, Integer> mapCopy = new HashMap<>();
            mapCopy.putAll(map);
            tempDataOfAttr.add(mapCopy);
        }

        int tempMinTime = minTime;
        int tempMaxTime = maxTime;
        int tempMinLatitude = minLatitude;
        int tempMaxLatitude = maxLatitude;
        int tempMinLongitude = minLongitude;
        int tempMaxLongitude = maxLongitude;

        //拷贝size属性
        int[] tempSize = size.clone();
        //取出切片包含的数据,遍历删除。删除过程中更新size属性和mass属性。
        Vector<Integer> dataIndexs = slice.getDataIndexs();
        int sliceMass = 0;
        for (Integer data : dataIndexs) {
            if(value[data] == 0)continue;
            sliceMass += value[data];
            //根据data查找到每个维度的属性，根据属性更新tempDataOfAttr并更新size。
            for(int i = 0; i < Converge.globalSize.length - 1; i++){
                int tempAttr = 0;
                if(i < Converge.CORE_COLUMN_LENGTH){
                    tempAttr = attVals[i][data];
                } else if (i == Converge.TIME_COLUMN) {
                    tempAttr = attVals[i][data];
                } else if(i == Converge.LATITUDE_COLUMN){
                    tempAttr = spaceVals[0][data];
                }else if(i == Converge.LONGITUDE_COLUMN){
                    tempAttr = spaceVals[1][data];
                }
                int attrValue = tempDataOfAttr.get(i).get(tempAttr);
                tempDataOfAttr.get(i).replace(tempAttr, attrValue - valueClone[data]);
                if(tempDataOfAttr.get(i).get(tempAttr) == 0){
                    //通过切片判断该维度的该属性是否存在，只有存在才修改其size
                    if(sliceExisted[i][tempAttr]){
                        tempSize[i]--;
                    }
                    sliceExisted[i][tempAttr] = false;
                }

                //更新最大最小时间和最大最小经纬度和时间
                if(i == Converge.TIME_COLUMN){
                    if(attVals[i][data] == tempMinTime){
                        for(int time = tempMinTime; time <= tempMaxTime; time++){
                            if(sliceExisted[i][time]){
                                tempMinTime = time;
                                break;
                            }
                        }
                    }else if(attVals[i][data] == tempMaxTime){
                        for(int time = tempMaxTime; time >= tempMinTime; time--){
                            if(sliceExisted[i][time]){
                                tempMaxTime = time;
                                break;
                            }
                        }
                    }
                }else if(i == Converge.LATITUDE_COLUMN){
                    if(spaceVals[0][data] == tempMinLatitude){
                        for(int latitude = minLatitude; latitude <= maxLatitude; latitude++){
                            if(sliceExisted[i][latitude]){
                                tempMinLatitude = latitude;
                                break;
                            }
                        }
                    }else if(spaceVals[0][data] == tempMaxLatitude){
                        for(int latitude = maxLatitude; latitude >= minLatitude; latitude--){
                            if(sliceExisted[i][latitude]){
                                tempMaxLatitude = latitude;
                                break;
                            }
                        }
                    }
                }else if(i == Converge.LONGITUDE_COLUMN){
                    if(spaceVals[1][data] == tempMinLongitude){
                        for(int longitude = minLongitude; longitude < maxLongitude; longitude++){
                            if(sliceExisted[i][longitude]){
                                tempMinLongitude = longitude;
                                break;
                            }
                        }
                    }else if(spaceVals[1][data] == tempMaxLongitude){
                        for(int longitude = maxLongitude; longitude >= minLongitude; longitude--){
                            if(sliceExisted[i][longitude]){
                                tempMaxLongitude = longitude;
                                break;
                            }
                        }
                    }
                }
            }
            //将当前data的value设置为0，防止更新极值时更新到该值。
            valueClone[data] = 0;
        }

        double densityBefore = Util.computeDensity(this);
        int coreSize = 0;
        for (int i = 0; i < Converge.globalSize.length - 1; i++) {
            coreSize += tempSize[i];
        }
        double densityAfter = (mass - sliceMass) * (float)Converge.CORE_COLUMN_LENGTH / coreSize;

        double temporalRangeBefore = Util.computeTemporalRange(minTime, maxTime);
        double temporalRangeAfter = Util.computeTemporalRange(tempMinTime, tempMaxTime);

        double spatialRangeBefore = Util.computeSpatialRange(minLatitude, maxLatitude, minLongitude, maxLongitude);
        double spatialRangeAfter = Util.computeSpatialRange(tempMinLatitude, tempMaxLatitude, tempMinLongitude, tempMaxLongitude);

        double measurementBefore = Util.computeMeasurement(densityBefore, temporalRangeBefore, spatialRangeBefore);
        double measurementAfter = Util.computeMeasurement(densityAfter, temporalRangeAfter, spatialRangeAfter);

        return measurementAfter > measurementBefore;
        // return densityAfter > densityBefore;
    }

    /**
     * 用于更新最大最小值的比较器
     */
    class AttrComparator implements Comparator<Integer>{

        @Override
        public int compare(Integer o1, Integer o2) {
            if(value[o1] == 0 && value[o2] != 0){
                return -1;
            }else if(value[o1] != 0 && value[o2] == 0){
                return 1;
            }else if(value[o1] == 0 && value[o2] == 0){
                return 0;
            }else{
                return Double.compare(spaceVals[0][o1], spaceVals[0][o2]);
            }
        }
    }

    /**
     * 删除指定切片
     * @param sliceInformation 0维度为数据的属性值，1维度数据为所属维度
     * @param index
     */
    public void deleteSlice(int[][] sliceInformation, int index) {
        //1、根据参数找到预删除的切片
        int attr = sliceInformation[0][index];
        int mode = sliceInformation[1][index];
        Slice slice = slices.get(mode).get(attr - 1);

        int deleteMass = 0;
        //找到该切片包含的数据
        Vector<Integer> dataIndexs = slice.getDataIndexs();
        for (int k = 0; k < dataIndexs.size(); k++) {
            Integer data = dataIndexs.get(k);
            if(value[data] == 0)continue;
            deleteMass += value[data];

            //更新slices
            for(int i = 0; i < Converge.globalSize.length - 1; i++){
                int tempAttr = 0;
                if(i < Converge.CORE_COLUMN_LENGTH){
                    tempAttr = attVals[i][data];
                } else if (i == Converge.TIME_COLUMN) {
                    tempAttr = attVals[i][data];
                } else if(i == Converge.LATITUDE_COLUMN){
                    tempAttr = (int) spaceVals[0][data];
                }else if(i == Converge.LONGITUDE_COLUMN){
                    tempAttr = (int) spaceVals[1][data];
                }

                //更新Slice
                Slice tempSlice = slices.get(i).get(tempAttr - 1);

                //更新dataOfAttr
                Map<Integer, Integer> attrMap = dataOfAttr.get(i);
                int attrValue = attrMap.get(tempAttr);
                attrMap.replace(tempAttr, attrValue - value[data]);
                if(attrMap.get(tempAttr) == 0){
                    attrMap.remove(tempAttr);
                    //通过切片判断该维度的该属性是否存在，只有存在才修改其size
                    if(tempSlice.isExists()){
                        //对于经纬度和时间不需要更改其size
                        size[i]--;
                    }
                    tempSlice.setExists(false);
                }


                //更新最大最小时间信息
                if(i == Converge.TIME_COLUMN){
                    if(attVals[i][data] == minTime){
                        for(int time = minTime; time <= maxTime; time++){
                            if(slices.get(i).get(time - 1).isExists()){
                                minTime = time;
                                break;
                            }
                        }
                    }else if(attVals[i][data] == maxTime){
                        for(int time = maxTime; time >= minTime; time--){
                            if(slices.get(i).get(time - 1).isExists()){
                                maxTime = time;
                                break;
                            }
                        }
                    }
                }else if(i == Converge.LATITUDE_COLUMN){
                    //更新最小最大经纬度信息
                    if(spaceVals[0][data] == minLatitude){
                        for(int latitude = minLatitude; latitude <= maxLatitude; latitude++){
                            if(slices.get(i).get(latitude - 1).isExists()){
                                minLatitude = latitude;
                                break;
                            }
                        }
                    }else if(spaceVals[0][data] == maxLatitude){
                        for(int latitude = maxLatitude; latitude >= minLatitude; latitude--){
                            if(slices.get(i).get(latitude - 1).isExists()){
                                maxLatitude = latitude;
                                break;
                            }
                        }
                    }
                }else if(i == Converge.LONGITUDE_COLUMN){
                    if(spaceVals[1][data] == minLongitude){
                        for(int longitude = minLongitude; longitude < maxLongitude; longitude++){
                            if(slices.get(i).get(longitude - 1).isExists()){
                                minLongitude = longitude;
                                break;
                            }
                        }
                    }else if(spaceVals[1][data] == maxLongitude){
                        for(int longitude = maxLongitude; longitude >= minLongitude; longitude--){
                            if(slices.get(i).get(longitude - 1).isExists()){
                                maxLongitude = longitude;
                                break;
                            }
                        }
                    }
                }
            }
            //将对应数据value值置为0
            value[data] = 0;
        }
        //更新mass
        mass-=deleteMass;

    }

    /**
     * 删除指定切片
     */
    public void deleteSlice(int attr, int mode) {
        //1、根据参数找到预删除的切片
        Slice slice = slices.get(mode).get(attr - 1);

        int deleteMass = 0;
        //找到该切片包含的数据
        Vector<Integer> dataIndexs = slice.getDataIndexs();
        for (int k = 0; k < dataIndexs.size(); k++) {
            Integer data = dataIndexs.get(k);
            if(value[data] == 0)continue;
            deleteMass += value[data];

            //更新slices
            for(int i = 0; i < Converge.globalSize.length - 1; i++){
                int tempAttr = 0;
                if(i < Converge.CORE_COLUMN_LENGTH){
                    tempAttr = attVals[i][data];
                } else if (i == Converge.TIME_COLUMN) {
                    tempAttr = attVals[i][data];
                } else if(i == Converge.LATITUDE_COLUMN){
                    tempAttr = spaceVals[0][data];
                }else if(i == Converge.LONGITUDE_COLUMN){
                    tempAttr = spaceVals[1][data];
                }

                //更新Slice
                Slice tempSlice = slices.get(i).get(tempAttr - 1);

                //更新dataOfAttr
                Map<Integer, Integer> attrMap = dataOfAttr.get(i);
                int attrValue = attrMap.get(tempAttr);
                attrMap.replace(tempAttr, attrValue - value[data]);
                if(attrMap.get(tempAttr) == 0){
                    attrMap.remove(tempAttr);
                    //通过切片判断该维度的该属性是否存在，只有存在才修改其size
                    if(tempSlice.isExists()){
                        //对于经纬度和时间不需要更改其size
                        size[i]--;
                    }
                    tempSlice.setExists(false);
                }

                //更新最大最小时间信息
                if(i == Converge.TIME_COLUMN){
                    if(attVals[i][data] == minTime){
                        for(int time = minTime; time <= maxTime; time++){
                            if(slices.get(i).get(time - 1).isExists()){
                                minTime = time;
                                break;
                            }
                        }
                    }else if(attVals[i][data] == maxTime){
                        for(int time = maxTime; time >= minTime; time--){
                            if(slices.get(i).get(time - 1).isExists()){
                                maxTime = time;
                                break;
                            }
                        }
                    }
                }else if(i == Converge.LATITUDE_COLUMN){
                    //更新最小最大经纬度信息
                    if(spaceVals[0][data] == minLatitude){
                        for(int latitude = minLatitude; latitude <= maxLatitude; latitude++){
                            if(slices.get(i).get(latitude - 1).isExists()){
                                minLatitude = latitude;
                                break;
                            }
                        }
                    }else if(spaceVals[0][data] == maxLatitude){
                        for(int latitude = maxLatitude; latitude >= minLatitude; latitude--){
                            if(slices.get(i).get(latitude - 1).isExists()){
                                maxLatitude = latitude;
                                break;
                            }
                        }
                    }
                }else if(i == Converge.LONGITUDE_COLUMN){
                    if(spaceVals[1][data] == minLongitude){
                        for(int longitude = minLongitude; longitude < maxLongitude; longitude++){
                            if(slices.get(i).get(longitude - 1).isExists()){
                                minLongitude = longitude;
                                break;
                            }
                        }
                    }else if(spaceVals[1][data] == maxLongitude){
                        for(int longitude = maxLongitude; longitude >= minLongitude; longitude--){
                            if(slices.get(i).get(longitude - 1).isExists()){
                                maxLongitude = longitude;
                                break;
                            }
                        }
                    }
                }
            }
            //将对应数据value值置为0
            value[data] = 0;
        }
        //更新mass
        mass-=deleteMass;

    }

    /**
     * 模拟划分，并判断划分之后子张量是否大于划分前的子张量
     * @param dataList1
     * @param dataList2
     * @return
     */
    public double[] preDivideSlice(List<Integer> dataList1, List<Integer> dataList2, int mode) {
        double[] result = new double[2];

        if(dataList1.size() == 0 || dataList2.size() == 0)return result;

        int[] tempMinTime = new int[2];
        int[] tempMaxTime = new int[2];
        double[] tempMinLatitude = new double[2];
        double[] tempMaxLatitude = new double[2];
        double[] tempMinLongitude = new double[2];
        double[] tempMaxLongitude = new double[2];

        //通过mass和dataOfAttr属性计算每个切片集的密度
        Integer totalMass1 = 0;
        Integer totalMass2 = 0;
        List<Set<Integer>> totalSizeList1 = new ArrayList<>();
        List<Set<Integer>> totalSizeList2 = new ArrayList<>();
        for (int i = 0; i < Converge.globalSize.length - 1; i++) {
            totalSizeList1.add(new HashSet<>());
            totalSizeList2.add(new HashSet<>());
        }
        for (int i = 0; i < dataList1.size(); i++) {
            totalMass1 += value[dataList1.get(i)];
            for (int j = 0; j < Converge.globalSize.length - 1; j++) {
                if(j < Converge.globalSize.length - 3){
                    totalSizeList1.get(j).add(attVals[j][dataList1.get(i)]);
                }else if(j == Converge.LATITUDE_COLUMN){
                    totalSizeList1.get(j).add(spaceVals[0][dataList1.get(i)]);
                } else if (j == Converge.LONGITUDE_COLUMN) {
                    totalSizeList1.get(j).add(spaceVals[1][dataList1.get(i)]);
                }
            }
        }
        for (int i = 0; i < dataList2.size(); i++) {
            totalMass2 += value[dataList2.get(i)];
            for (int j = 0; j < Converge.globalSize.length - 1; j++) {
                if(j < Converge.globalSize.length - 3){
                    totalSizeList2.get(j).add(attVals[j][dataList2.get(i)]);
                }else if(j == Converge.LATITUDE_COLUMN){
                    totalSizeList2.get(j).add(spaceVals[0][dataList2.get(i)]);
                } else if (j == Converge.LONGITUDE_COLUMN) {
                    totalSizeList2.get(j).add(spaceVals[1][dataList2.get(i)]);
                }
            }
        }

        //找到两个数据集的最大最小时空
        tempMinTime[0] = totalSizeList1.get(Converge.TIME_COLUMN).stream().min(Integer::compareTo).get();
        tempMaxTime[0] = totalSizeList1.get(Converge.TIME_COLUMN).stream().max(Integer::compareTo).get();
        tempMinLatitude[0] = totalSizeList1.get(Converge.LATITUDE_COLUMN).stream().min(Integer::compareTo).get();
        tempMaxLatitude[0] = totalSizeList1.get(Converge.LATITUDE_COLUMN).stream().max(Integer::compareTo).get();
        tempMinLongitude[0] = totalSizeList1.get(Converge.LONGITUDE_COLUMN).stream().min(Integer::compareTo).get();
        tempMaxLongitude[0] = totalSizeList1.get(Converge.LONGITUDE_COLUMN).stream().max(Integer::compareTo).get();

        tempMinTime[1] = totalSizeList2.get(Converge.TIME_COLUMN).stream().min(Integer::compareTo).get();
        tempMaxTime[1] = totalSizeList2.get(Converge.TIME_COLUMN).stream().max(Integer::compareTo).get();
        tempMinLatitude[1] = totalSizeList2.get(Converge.LATITUDE_COLUMN).stream().min(Integer::compareTo).get();
        tempMaxLatitude[1] = totalSizeList2.get(Converge.LATITUDE_COLUMN).stream().max(Integer::compareTo).get();
        tempMinLongitude[1] = totalSizeList2.get(Converge.LONGITUDE_COLUMN).stream().min(Integer::compareTo).get();
        tempMaxLongitude[1] = totalSizeList2.get(Converge.LONGITUDE_COLUMN).stream().max(Integer::compareTo).get();


        //密度条件
        int size1 = 0;
        int size2 = 0;
        int otherSize1 = 0;
        int otherSize2 = 0;
        for (int i = 0; i < Converge.globalSize.length - 1; i++) {
            // if(i < Converge.CORE_COLUMN_LENGTH){
                size1 += totalSizeList1.get(i).size();
                size2 += totalSizeList2.get(i).size();
            // }else if(i == mode){
            //     otherSize1 += totalSizeList1.get(i).size();
            //     otherSize2 += totalSizeList2.get(i).size();
            // }
        }


        double densityAfter1 = Util.computeDensity(totalMass1, size1);
        double densityAfter2 = Util.computeDensity(totalMass2, size2);
        // double densityAfter1 = Util.computeDensity(totalMass1, size1);
        // double densityAfter2 = Util.computeDensity(totalMass2, size2);

        //时间范围条件
        double timeRangeAfter1 = Util.computeTemporalRange(tempMinTime[0],tempMaxTime[0]);
        double timeRangeAfter2 = Util.computeTemporalRange(tempMinTime[1],tempMaxTime[1]);
        //空间范围条件
        double spatialRangeAfter1 = Util.computeSpatialRange(tempMinLatitude[0], tempMaxLatitude[0], tempMinLongitude[0], tempMaxLongitude[0]);
        double spatialRangeAfter2 = Util.computeSpatialRange(tempMinLatitude[1], tempMaxLatitude[1], tempMinLongitude[1], tempMaxLongitude[1]);

        result[0] = Util.computeMeasurement(densityAfter1, timeRangeAfter1, spatialRangeAfter1);
        result[1] = Util.computeMeasurement(densityAfter2, timeRangeAfter2, spatialRangeAfter2);

        // result[0] = densityAfter1;
        // result[1] = densityAfter2;

        return result;
    }

    /**
     * 删除子张量中mode维度的i属性
     * @param i 删除的属性
     * @param mode 删除的用户维度
     */
    public void deleteAttr(int i, int mode){
        Map<Integer, Integer> attrMasses = dataOfAttr.get(mode);
        int attrMass = attrMasses.get(i);
        //更新mass
        mass = mass - attrMass;

        //更新mode维度属性存在状态
        //更新modeMass
        List<Integer> dataList = slices.get(mode).get(i - 1).getDataIndexs();
        if(dataList != null){
            for (Integer data : dataList) {
                if(value[data] == 0)continue;
                int tempValue = value[data];
                //将value设置为0，代表该纪录被删除
                value[data] = 0;
                for(int k = 0; k < Converge.globalSize.length - 1; k++){
                    int attr = 0;
                    if(k < Converge.CORE_COLUMN_LENGTH){
                        attr = attVals[k][data];
                    } else if (k == Converge.TIME_COLUMN) {
                        attr = attVals[k][data];
                    } else if(k == Converge.LATITUDE_COLUMN){
                        attr = (int) spaceVals[0][data];
                    }else if(k == Converge.LONGITUDE_COLUMN){
                        attr = (int) spaceVals[1][data];
                    }

                    int attrValue = dataOfAttr.get(k).get(attr);
                    dataOfAttr.get(k).replace(attr, attrValue - tempValue);
                    if(dataOfAttr.get(k).get(attr) == 0){
                        Slice slice = slices.get(k).get(attr - 1);
                        //通过切片判断该维度的该属性是否存在，只有存在才修改其size
                        // if(k < Converge.CORE_COLUMN_LENGTH){
                        //对于经纬度和时间不需要更改其size
                        size[k]--;
                        // }
                        slice.setExists(false);
                    }

                    /*AttrComparator comparator = new AttrComparator();

                    //更新最大最小时间信息
                    if(k == Converge.TIME_COLUMN){
                        if(attVals[k][data] == minTime){
                            for(int time = minTime; time <= maxTime; time++){
                                if(slices.get(k).get(time - 1).isExists()){
                                    Vector<Integer> dataIndexsInTime = slices.get(k).get(time - 1).getDataIndexs();
                                    Integer minTimeIndex = dataIndexsInTime.stream().min(comparator).get();
                                    minTime = attVals[k][minTimeIndex];
                                    break;
                                }
                            }
                        }else if(attVals[k][data] == maxTime){
                            for(int time = maxTime; time >= minTime; time--){
                                if(slices.get(k).get(time - 1).isExists()){
                                    Vector<Integer> dataIndexsInTime = slices.get(k).get(time - 1).getDataIndexs();
                                    Integer maxTimeIndex = dataIndexsInTime.stream().max(comparator).get();
                                    maxTime = attVals[k][maxTimeIndex];
                                    break;
                                }
                            }
                        }
                    }else if(k == Converge.LATITUDE_COLUMN){
                        //更新最小最大经纬度信息
                        if(spaceVals[0][data] == minLatitude){
                            for(int latitude = minLatitude; latitude <= maxLatitude; latitude++){
                                if(slices.get(k).get(latitude - 1).isExists()){
                                    Vector<Integer> dataIndexInLatitude = slices.get(k).get(latitude - 1).getDataIndexs();
                                    Integer minLatitudeIndex = dataIndexInLatitude.stream().min(comparator).get();
                                    minLatitude = spaceVals[0][minLatitudeIndex];
                                    break;
                                }
                            }
                        }else if(spaceVals[0][data] == maxLatitude){
                            for(int latitude = maxLatitude; latitude >= minLatitude; latitude--){
                                if(slices.get(k).get(latitude - 1).isExists()){
                                    Vector<Integer> dataIndexInLatitude = slices.get(k).get(latitude - 1).getDataIndexs();
                                    Integer maxLatitudeIndex = dataIndexInLatitude.stream().max(comparator).get();
                                    maxLatitude = spaceVals[0][maxLatitudeIndex];
                                    break;
                                }
                            }
                        }
                    }else if(k == Converge.LONGITUDE_COLUMN){
                        if(spaceVals[1][data] == minLongitude){
                            for(int longitude = minLongitude; longitude < maxLongitude; longitude++){
                                if(slices.get(k).get(longitude - 1).isExists()){
                                    Vector<Integer> dataIndexInLongitude = slices.get(k).get(longitude - 1).getDataIndexs();
                                    Integer minLongitudeIndex = dataIndexInLongitude.stream().min(comparator).get();
                                    minLongitude = spaceVals[1][minLongitudeIndex];
                                    break;
                                }
                            }
                        }else if(spaceVals[1][data] == maxLongitude){
                            for(int longitude = maxLongitude; longitude >= minLongitude; longitude--){
                                if(slices.get(k).get(longitude - 1).isExists()){
                                    Vector<Integer> dataIndexInLongitude = slices.get(k).get(longitude - 1).getDataIndexs();
                                    Integer maxLongitudeIndex = dataIndexInLongitude.stream().max(comparator).get();
                                    maxLongitude = spaceVals[1][maxLongitudeIndex];
                                    break;
                                }
                            }
                        }
                    }*/

                }

            }
        }
    }

    public void setSlices(List<List<Slice>> slices) {
        this.slices = slices;
    }

    public void setDataOfAttr(List<Map<Integer, Integer>> dataOfAttr) {
        this.dataOfAttr = dataOfAttr;
    }

    public void setValue(int[] value) {
        this.value = value;
    }

    public List<Map<Integer, Integer>> getDataOfAttr() {
        return dataOfAttr;
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

    public int[][] getSpaceVals() {
        return spaceVals;
    }

    public int[] getValue() {
        return value;
    }

    public List<List<Slice>> getSlices() {
        return slices;
    }

    public int getMinTime() {
        return minTime;
    }

    public int getMaxTime() {
        return maxTime;
    }

    public int getMinLatitude() {
        return minLatitude;
    }

    public int getMaxLatitude() {
        return maxLatitude;
    }

    public int getMinLongitude() {
        return minLongitude;
    }

    public int getMaxLongitude() {
        return maxLongitude;
    }

}
