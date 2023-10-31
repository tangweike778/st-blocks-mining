import java.io.*;
import java.util.*;

public class Converge{

    //全局属性：用于转换属性值
    static List<Map<String, Integer>> convert = new ArrayList<>();
    static List<Map<Integer, String>> recoverDicts = new ArrayList<>();

    public static double globalMaxLatitude = Double.MIN_VALUE;
    public static double globalMinLatitude = Double.MAX_VALUE;
    public static double globalMaxlongitude = Double.MIN_VALUE;
    public static double globalMinlongitude = Double.MAX_VALUE;

    public static int globalMaxTime = Integer.MIN_VALUE;

    public static int globalMinTime = Integer.MAX_VALUE;

    public static int[] globalSize;

    static double timeStep;
    static double latitudeStep;
    static double longitudeStep;

    static int CORE_COLUMN_LENGTH;
    static int TIME_COLUMN;
    static int LATITUDE_COLUMN;
    static int LONGITUDE_COLUMN;

    static int PARTITION_NUMBER = 50;

    static int DENSITY_PROPORTION = 1;
    static int TEMPORAL_PROPORTION = 1;
    static int SPATIAL_PROPORTION = 1;

    /**
     *
     * @param args
     * input_path: 文件输入路径
     * output_path: 文件输出路径
     * splitStr：分隔符
     * k: 需要找寻的子张量个数k，默认k=0（表示只执行一次findSubTensor函数）
     */
    public static void main(String[] args) throws Exception {
        //利用缓冲流处理读取并处理文件
        if(args.length < 4){
            System.err.println("please enter parameters in the following format: ");
            System.err.println("input_path output_path splitStr k(Number of sub-tensor to search) Partition_number Density_proportion Temporal_proportion Spatial_proportion.");
            return;
        }
        final String input_path = args[0];
        final String output_path = args[1];
        final String splitStr = args[2];
        int k = Integer.parseInt(args[3]);
        if(args.length >= 5){
            PARTITION_NUMBER = Integer.parseInt(args[4]);
        }
        if (args.length >= 8) {
            DENSITY_PROPORTION = Integer.parseInt(args[5]);
            TEMPORAL_PROPORTION = Integer.parseInt(args[6]);
            SPATIAL_PROPORTION = Integer.parseInt(args[7]);
        }
        System.out.println("input_path: " + input_path);
        System.out.println("output_path: " + output_path);
        System.out.println("splitStr: \"" + splitStr + "\"");
        System.out.println("k (Number of sub-tensor to search): " + k + " (Default: 0)");

        long start = System.currentTimeMillis();

        System.out.println("Converting data into tensors...");

        long startTime = System.currentTimeMillis();
        preComputeParameter(input_path, splitStr);
        long endTime = System.currentTimeMillis();
        System.out.println("查找各个维度大小信息的时间为："+ (endTime - startTime) + "ms");


        startTime = System.currentTimeMillis();
        Tensor tensor = createTensor(input_path, splitStr);
        endTime = System.currentTimeMillis();
        System.out.println("构造父张量的时间为："+ (endTime - startTime) + "ms");

        //计算时空步长
        // computeTemporalSpatialStep();

        //调用findSubTensor()函数获取到k个子张量信息，该信息使用SubTensor对象数组保存。
        //将子张量放入subTensors的过程是从小到大，从后到前放入的。
        List<SubTensor> subTensors = findSubTensors(tensor, k);

        System.out.println("Running time: " + (System.currentTimeMillis() - start + 0.0) / 1000 + " seconds");

        //根据输出文件路径创建输出流，为每一个子张量信息创建相应的文件夹，并且存入相应维度的数据信息。
        //得到k个子张量信息之后，使用对应的字典将数据还原。
        System.out.println("A total of " + subTensors.size() + " dense sub-tensors were found.");
        System.out.println("Writing outputs...\n");
        outputTensorList(tensor,subTensors, output_path, subTensors.size());
        System.out.println("Outputs were written: " + (System.currentTimeMillis() - start + 0.0) / 1000 + " seconds was taken.");
    }

    /**
     * 计算全局参数
     * @param input_path
     * @param splitStr
     * @throws Exception
     */
    private static void preComputeParameter(String input_path, String splitStr) throws Exception {
        globalSize = Util.computerSize(input_path,splitStr);
        CORE_COLUMN_LENGTH = globalSize.length - 4;
        TIME_COLUMN = globalSize.length - 4;
        LATITUDE_COLUMN = globalSize.length - 3;
        LONGITUDE_COLUMN = globalSize.length - 2;
        timeStep = (globalMaxTime - globalMinTime) * 1.0 / PARTITION_NUMBER;
        latitudeStep = (globalMaxLatitude - globalMinLatitude) / PARTITION_NUMBER;
        longitudeStep = (globalMaxlongitude - globalMinlongitude) / PARTITION_NUMBER;
    }

    /**
     * 计算时空步长并存储在全局变量中
     */
/*    private static void computeTemporalSpatialStep(){
        int coreSize = 0;
        for(int i = 0; i < CORE_COLUMN_LENGTH; i++){
            coreSize += globalSize[i];
        }
        //计算时空步长
        timeStep = (coreSize) * 1.0 / (globalMaxTime - globalMinTime + 1);
        latitudeStep = (coreSize) * 1.0 / (Util.ceilSpaceByLatitude(globalMaxLatitude) - Util.ceilSpaceByLatitude(globalMinLatitude) + 1);
        longitudeStep = (coreSize) * 1.0 / (Util.ceilSpaceBylongitude(globalMaxlongitude) - Util.ceilSpaceBylongitude(globalMinlongitude) + 1);
    }*/

    /**
     * 根据子张量数组和输出文件路径，将子张量输出到磁盘中
     *
     * @param tensor
     * @param tensorList  子张量数组
     * @param output_path 输出路径
     */
    public static void outputTensorList(Tensor tensor, List<SubTensor> tensorList, String output_path, int totalNum) {
        File directory = null;
        try {
            //1、先判断文件加是否存在？不存在则创建对应文件夹
            directory = new File(output_path);
            if(!directory.exists() || !directory.isDirectory()){
                //文件夹不存在或不是文件夹，新建对应文件夹
                directory.mkdirs();
            }
            //2、遍历子张量数组，每一个张量存在一个文件中
            for (int i = 0; i < totalNum; i++) {
                SubTensor subTensor = tensorList.get(i);
                int[] size = subTensor.getSize();
                System.out.println("Tensor " + (i + 1) + ":");
                System.out.println("Mass: " + subTensor.getMass());
                System.out.println("Valume: " + Util.computeVolume(subTensor));
                System.out.println("Density: " + Util.computeDensity(subTensor));
                StringBuilder s = new StringBuilder();
                for (int j = 0; j < globalSize.length - 4; j++) {
                    s.append(size[j]).append(" x ");
                }
                s.replace(s.length() - 3, s.length(), "");
                System.out.println("The size of the " + (i + 1) + "th sub-tensor is " + s);

                //将指定子张量属性写到磁盘中
                outputProperties(tensor, subTensor, i + 1, output_path);

                //将指定子张量数据写到磁盘中
                ourputData(subTensor, i + 1, output_path);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void ourputData(SubTensor subTensor, int i, String output_path) throws IOException {
        String subTensorDataFilePath = null;
        File subTensorDataFile = null;
        FileWriter wt2 = null;
        subTensorDataFilePath = output_path + "/" + "block_" + i + ".tuples";
        subTensorDataFile = new File(subTensorDataFilePath);

        //输出张量数据文件
        if(subTensorDataFile.exists()){
            subTensorDataFile.delete();
        }else{
            subTensorDataFile.createNewFile();
        }
        wt2 = new FileWriter(subTensorDataFile);
        StringBuilder sb2 = new StringBuilder();
        //获取恢复数据的字典
        Map<Integer, String> timeRecoverDict = recoverDicts.get(TIME_COLUMN);
        int[][] attVals = subTensor.getAttVals();
        double[][] spaceVals = subTensor.getSpaceVals();
        int[] values = subTensor.getValue();
        //遍历子张量selectDataByTime，获得记录位置信息
        //希望输出的格式信息：
        //通过记录位置信息结合attVals、value输出该记录
        for (int k = 0; k < values.length; k++) {
            if(values[k] == 0)continue;
            for(int mode = 0; mode < globalSize.length - 4; mode++ ){
                sb2.append(recoverDicts.get(mode).get(attVals[mode][k])).append(",");
            }
            sb2.append(attVals[TIME_COLUMN][k]  + "," + (int)spaceVals[0][k] + "," + (int)spaceVals[1][k] + "," + values[k] + "\n");
        }
        wt2.append(sb2);
        wt2.flush();
        wt2.close();
    }

    /**
     * 输出对应子张量的各维度的属性和整体信息
     * @param subTensor
     * @param i
     * @param output_path
     * @throws IOException
     */
    public static void outputProperties(Tensor tensor, SubTensor subTensor, int i, String output_path) throws IOException {
        String subTensorProFilePath = null;
        File subTensorProFile = null;
        FileWriter wt1 = null;
        int[] size = subTensor.getSize();
        subTensorProFilePath = output_path + "/" + "block_" + i + ".properties";
        subTensorProFile = new File(subTensorProFilePath);

        //输出张量属性文件
        if(subTensorProFile.exists()){
            subTensorProFile.delete();
        }else{
            subTensorProFile.createNewFile();
        }
        wt1 = new FileWriter(subTensorProFilePath);

        //找出每一个维度存在的属性
        List<List<Integer>> modeAttrs = new ArrayList<>();
        List<boolean[]> coreDimensionExisted = subTensor.getCoreDimensionExisted();
        for (int mode = 0; mode < CORE_COLUMN_LENGTH; mode++){
            boolean[] coreExisted = coreDimensionExisted.get(mode);
            modeAttrs.add(Util.findTrueIndexes(coreExisted));
        }
        // modeAttrs.add(Util.findTrueIndexes(subTensor.getUserExisted()));
        // modeAttrs.add(Util.findTrueIndexes(subTensor.getObjectExisted()));
        modeAttrs.add(Util.findTrueIndexes(subTensor.getTimeExisted()));

        //希望的文件格式：
        StringBuilder sb1 = new StringBuilder();
        sb1.append("Tensor " + i + ":\n");
        sb1.append("Mass: "+ subTensor.getMass() +"\n");
        sb1.append("Valume: "+ Util.computeVolume(subTensor) +"\n");
        sb1.append("Mass: "+ Util.computeDensity(subTensor.getMass(), Util.computeVolume(subTensor)) +"\n");
        sb1.append("The size of the " + i + "th sub-tensor is ");
        for (int j = 0; j < CORE_COLUMN_LENGTH; j++) {
            sb1.append(size[j] + " x ");
        }
        sb1.append(size[TIME_COLUMN] + " x " + size[LATITUDE_COLUMN] + " x " + size[LONGITUDE_COLUMN] + "\n");
        for (int j = 0; j < CORE_COLUMN_LENGTH; j++) {
            sb1.append(j + ":\n");
            Map<Integer, String> recoverDict = recoverDicts.get(j);
            List<Integer> modeAttr = modeAttrs.get(j);
            for (Integer attr : modeAttr) {
                sb1.append(recoverDict.get(attr) + ", ");
            }
            sb1.replace(sb1.length() - 2, sb1.length(), "\n");
        }

        //查询子张量的最大最小时间和经纬度
        subTensor.updatePeakInfomation();

        double globalArea = (globalMaxlongitude - globalMinlongitude + 1) * (globalMaxLatitude - globalMinLatitude + 1);
        double localArea = (subTensor.getMaxLatitude() - subTensor.getMinLatitude() + 1) * (subTensor.getMaxlongitude() - subTensor.getMinLongitude() + 1);
        System.out.println("空间范围的百分比为：" + ((localArea * 100.0) / globalArea) + "%");
        System.out.println("时间范围的百分比为：" + ((subTensor.getMaxTime() - subTensor.getMinTime() + 1) * 100.0) / (globalMaxTime * 1.0) + "%");
        System.out.println();
        wt1.append(sb1);
        wt1.flush();
        wt1.close();
    }

    /**
     * 找到指定k个密集子张量
     *
     * @param tensor 全局张量
     * @param k      要求返回的子张量个数
     * @return 子张量数组
     */
    public static List<SubTensor> findSubTensors(Tensor tensor, int k) {
        PriorityQueue<Pair<SubTensor, Double>> subTensorPriorityQueue = new PriorityQueue<>(k, new Comparator<Pair<SubTensor, Double>>() {

            //TODO:修改度量时需要修改
            @Override
            public int compare(Pair<SubTensor, Double> o1, Pair<SubTensor, Double> o2) {
                double a = Util.computeDensity(o1.getKey());
                double b = Util.computeDensity(o2.getKey());
                return Double.compare(a, b);
            }
        });

        findSubTensorsInOneRound(tensor, subTensorPriorityQueue, k, 0);
        System.gc();

        //将top-k个子张量转换为列表并返回
        List<SubTensor> subTensors = new ArrayList<>();
        for (int i = 0; i < k; i++){
            if(subTensorPriorityQueue.peek() == null)break;
            SubTensor subTensor = subTensorPriorityQueue.poll().getKey();
            subTensors.add(subTensor);
        }
        Collections.reverse(subTensors);

        return subTensors;
    }


    /**
     * 返回一轮密集子张量的筛选列表d
     * @param tensor 全局张量
     * @return 子张量列表
     */
    public static void findSubTensorsInOneRound(Tensor tensor, PriorityQueue<Pair<SubTensor, Double>> subTensorPriorityQueue, int k, int circleNum) {
        //计算空间块的平均质量,遍历空间块,将空间块中小于平均空间块密度的空间块删除。
        double timeSpaceAvgMass = tensor.getMass() * 1.0 / tensor.getSelectDataByTimeSpace().size();
        int[] timeSpaceMass = tensor.getTimeSpaceMass();
        boolean[] timeSpaceExisted = tensor.getTimeSpaceExisted().clone();
        for (int i = 0; i < timeSpaceExisted.length; i++) {
            if(timeSpaceMass[i] < timeSpaceAvgMass){
                //逻辑删除
                timeSpaceExisted[i] = false;
            }
        }

        //将剩下的时空块按照其质量进行排序。
        List<Integer> timespaceIndexList = new ArrayList<>();
        for (int i = 0; i < timeSpaceExisted.length; i++) {
            if(timeSpaceExisted[i])timespaceIndexList.add(i);
        }
        Integer[] timeSpaceIndexsBySortMass = new Integer[timespaceIndexList.size()];
        timespaceIndexList.toArray(timeSpaceIndexsBySortMass);
        sort(timeSpaceIndexsBySortMass, timeSpaceMass, 0, timespaceIndexList.size() - 1);

        //选取密度最大的时空块开始向外扩展,根据timeSpaceIndexArray找到质量最大的时空块。
        for (int i = 0; i < timeSpaceIndexsBySortMass.length; i++) {
            if(timeSpaceIndexsBySortMass[i] == null || !timeSpaceExisted[timeSpaceIndexsBySortMass[i]]){
                continue;
            }
            int timeSpaceBlockIndex = timeSpaceIndexsBySortMass[i];
            //从tensor中的selectDataByTimeSpace属性中取出数据，将取出的最大密度时空块删除。
            List<Integer> dataList = tensor.getSelectDataByTimeSpace().get(timeSpaceBlockIndex);
            tensor.getSelectDataByTimeSpace().remove(timeSpaceBlockIndex);
            timeSpaceExisted[timeSpaceBlockIndex] = false;

            //创建用于记录下一次可以扩散的时空块序号的数据结构
            //其中timeSpaceState：1表示该块属于可拼接的时空块，0表示未判断的时空块，-1表示以及存在当前subTensor中的时空块
            int[] timeSpaceState = new int[timeSpaceMass.length];
            Set<Integer> expansionTimeSpaceBlocks = new HashSet<>();
            timeSpaceState[timeSpaceBlockIndex] = 1;
            expansionTimeSpaceBlocks.add(timeSpaceBlockIndex);
            //将当前block空间中邻接六个方向的block块查找出来加入expansionTimeSpaceBlocks中
            findAdjacencyBlock(timeSpaceState, expansionTimeSpaceBlocks, timeSpaceBlockIndex);

            //可能会出现dataList有数据，但是数据都为0造成subTensor为空，此时应该遍历下一个块。
            SubTensor subTensor = createSubTensorByDatas(tensor, dataList);
            if (subTensor == null)continue;
            //计算后续sigmoid函数需要平移的距离
            Util.setSigmoid_base(Util.computeDensity(subTensor));
            while(true){
                Set<Integer> addBlock = new HashSet<>();
                double maxDensity = Util.computeTotalDensity(subTensor);
                for (Integer blockIndex : expansionTimeSpaceBlocks) {
                    double density = preComputeDensityAfterExtension(subTensor, tensor, blockIndex);
                    //下降不超过百分之25
                    if(density > maxDensity){
                        addBlock.add(blockIndex);
                    }
                }

                //检验全部加入后密度变化，若加入后密度大于加入前密度则执行扩展操作。
                double densityBatchAdd = 0;
                try {
                    densityBatchAdd = preComputeDensityAfterBatchExtension(subTensor, tensor, addBlock);
                } catch (Exception e) {
                    System.out.println(maxDensity);
                    System.out.println(subTensor.getMass() * (float)(Converge.CORE_COLUMN_LENGTH) / Util.computeVolume(subTensor));
                    System.out.println(Util.computeTemporalRange(subTensor.getMaxTime(), subTensor.getMinTime()));
                    System.out.println(Util.computeOfSpatialRange(subTensor.getMinLatitude(), subTensor.getMaxLatitude(), subTensor.getMinLongitude(), subTensor.getMaxLongitude()));
                    throw new RuntimeException(e);
                }
                if(densityBatchAdd > maxDensity){
                    extendSubtensorInDirect(subTensor, tensor, addBlock);
                    for (Integer blockIndex : expansionTimeSpaceBlocks) {
                        if(addBlock.contains(blockIndex)){
                            findAdjacencyBlock(timeSpaceState, addBlock, blockIndex);
                        }
                    }
                    expansionTimeSpaceBlocks = addBlock;
                }else{
                    processCore(tensor, subTensor, k, subTensorPriorityQueue, circleNum);
                    break;
                }

            }
        }

    }

    private static double preComputeDensityAfterBatchExtension(SubTensor subTensor, Tensor tensor, Set<Integer> addBlock) {
        List<Integer> dataInTimeSpaceBlock = new ArrayList<>();
        for (Integer blockIndex : addBlock) {
            dataInTimeSpaceBlock.addAll(tensor.getSelectDataByTimeSpace().get(blockIndex));
        }
        if(dataInTimeSpaceBlock.size() == 0){
            return 0;
        }

        //3、将数据模拟加入到subtensor中，并计算密度。
        double densityAfterAppend = subTensor.imitateAppendDatas(tensor, dataInTimeSpaceBlock);

        //4、将密度返回。
        return densityAfterAppend;
    }

    /**
     * 根据timeSpaceBlockIndex，将其六个方向的block加入到expansionTimeSpaceBlocks中
     * @param timeSpaceState
     * @param expansionTimeSpaceBlocks
     * @param timeSpaceBlockIndex
     */
    private static void findAdjacencyBlock(int[] timeSpaceState, Set<Integer> expansionTimeSpaceBlocks, int timeSpaceBlockIndex) {
        int[] expansionBlocks = new int[6];
        expansionBlocks[0] = timeSpaceBlockIndex + 1;
        expansionBlocks[1] = timeSpaceBlockIndex - 1;
        expansionBlocks[2] = timeSpaceBlockIndex + PARTITION_NUMBER;
        expansionBlocks[3] = timeSpaceBlockIndex - PARTITION_NUMBER;
        expansionBlocks[4] = timeSpaceBlockIndex + (PARTITION_NUMBER * PARTITION_NUMBER);
        expansionBlocks[5] = timeSpaceBlockIndex - (PARTITION_NUMBER * PARTITION_NUMBER);

        expansionTimeSpaceBlocks.remove(timeSpaceBlockIndex);
        timeSpaceState[timeSpaceBlockIndex] = -1;

        for (int blockIndex : expansionBlocks) {
            if(blockIndex > 0 && blockIndex < timeSpaceState.length && timeSpaceState[blockIndex] == 0){
                timeSpaceState[blockIndex] = 1;
                expansionTimeSpaceBlocks.add(blockIndex);
            }
        }

    }

    public static SubTensor createSubTensorByDatas(Tensor tensor, List<Integer> dataList) {
        int newMass = 0;
        int newMinTime = Integer.MAX_VALUE;
        int newMaxTime = Integer.MIN_VALUE;
        double newMinLongitude = Double.MAX_VALUE;
        double newMaxLongitude = -Double.MAX_VALUE;
        double newMinLatitude = Double.MAX_VALUE;
        double newMaxLatitude = -Double.MAX_VALUE;
        int[] newSize = new int[globalSize.length - 1];
        int timeSpaceMaxIndex = PARTITION_NUMBER * PARTITION_NUMBER * PARTITION_NUMBER + PARTITION_NUMBER * PARTITION_NUMBER + PARTITION_NUMBER;
        //获取记录数量
        //初始时用户维度、商品维度和原始张量一样，时间维度只有当前时间属性为true
        boolean[] newTimeExisted = new boolean[PARTITION_NUMBER + 1];
        List<boolean[]> newCoreDimensionExisted = new ArrayList<>();
        for(int i = 0; i < CORE_COLUMN_LENGTH; i++){
            newCoreDimensionExisted.add(new boolean[globalSize[i] + 1]);
        }
        boolean[] newTimeSpaceExisted = new boolean[timeSpaceMaxIndex + 1];
        boolean[] newLatitudeExisted = new boolean[PARTITION_NUMBER + 1];
        boolean[] newLongitudeExisted = new boolean[PARTITION_NUMBER + 1];

        List<int[]> newCoreDimensionMass = new ArrayList<>();
        for(int i = 0; i < CORE_COLUMN_LENGTH; i++){
            newCoreDimensionMass.add(new int[globalSize[i] + 1]);
        }
        // int[] newUserMass = new int[globalSize[0] + 1];
        // int[] newObjectMass = new int[globalSize[1] + 1];
        int[] newTimeMass = new int[PARTITION_NUMBER + 1];
        int[] newTimeSpaceMass = new int[timeSpaceMaxIndex + 1];
        int[] newLatitudeMass = new int[PARTITION_NUMBER + 1];
        int[] newLongitudeMass = new int[PARTITION_NUMBER + 1];

        int[][] newAttVals = tensor.getAttVals();
        double[][] newSpaceVals = tensor.getSpaceVals();
        int[] newValue = new int[tensor.getValue().length];
        // Map<Integer, List<Integer>> newSelectDataByUser = new HashMap<>();
        // Map<Integer, List<Integer>> newSelectDataByObject = new HashMap<>();
        List<Map<Integer, List<Integer>>> newSelectDataByCoreDimension = new ArrayList<>();
        for(int i = 0; i < CORE_COLUMN_LENGTH; i++){
            newSelectDataByCoreDimension.add(new HashMap<>());
        }
        Map<Integer, List<Integer>> newSelectDataByTime = new HashMap<>();
        Map<Integer, List<Integer>> newSelectDataByLatitude = new HashMap<>();
        Map<Integer, List<Integer>> newSelectDataByLongitude = new HashMap<>();

        //当前子张量记录索引
        int[] value = tensor.getValue();
        if(dataList == null)return null;
        for(Integer data : dataList){
            //每条数据的质量
            int dataMass = value[data];
            //数据已经被输出
            if(dataMass == 0)continue;
            //修改质量
            newMass += value[data];
            newValue[data] = dataMass;
            for(int j = 0; j < globalSize.length - 3; j++){
                //拿到对应的维度的属性值
                int attribute = newAttVals[j][data];
                if (j < CORE_COLUMN_LENGTH){
                    int[] masses = newCoreDimensionMass.get(j);
                    masses[attribute] += dataMass;
                    Map<Integer, List<Integer>> newSelectDataByCore = newSelectDataByCoreDimension.get(j);
                    if(!newSelectDataByCore.containsKey(attribute)){
                        newSelectDataByCore.put(attribute, new ArrayList<>());
                    }
                    newSelectDataByCore.get(attribute).add(data);
                    boolean[] newCoreExisted = newCoreDimensionExisted.get(j);
                    if(!newCoreExisted[attribute]){
                        newCoreExisted[attribute] = true;
                        newSize[j]++;
                    }
                }else if (j == TIME_COLUMN) {
                    newTimeMass[Util.transformTime(attribute)] += dataMass;
                    newMinTime = Math.min(attribute, newMinTime);
                    newMaxTime = Math.max(attribute, newMaxTime);
                    if(!newSelectDataByTime.containsKey(Util.transformTime(attribute))){
                        newSelectDataByTime.put(Util.transformTime(attribute), new ArrayList<>());
                    }
                    newSelectDataByTime.get(Util.transformTime(attribute)).add(data);
                    if(!newTimeExisted[Util.transformTime(attribute)]){
                        newTimeExisted[Util.transformTime(attribute)] = true;
                        newSize[j]++;
                    }
                }
            }
            //拿到空间数据
            double latitude = newSpaceVals[0][data];
            double longitude = newSpaceVals[1][data];
            int processedLatitude = Util.transformLatitude(latitude);
            int processedLongitude = Util.transformLongitude(longitude);

            newMinLatitude = Math.min(latitude, newMinLatitude);
            newMaxLatitude = Math.max(latitude, newMaxLatitude);
            newMinLongitude = Math.min(longitude, newMinLongitude);
            newMaxLongitude = Math.max(longitude, newMaxLongitude);

            //计算时空块索引地址
            int timeSpaceBlockIndex = Util.computeTimeSpaceIndexDirectly(Util.transformTime(newAttVals[TIME_COLUMN][data]), processedLatitude, processedLongitude);
            //处理经纬度质量是向上取整
            newLatitudeMass[processedLatitude] += dataMass;
            newLongitudeMass[processedLongitude] += dataMass;
            newTimeSpaceMass[timeSpaceBlockIndex] += value[data];
            if(!newTimeSpaceExisted[timeSpaceBlockIndex]){
                newTimeSpaceExisted[timeSpaceBlockIndex] = true;
            }

            //更新Latitude和Longitude的selectData
            if(!newSelectDataByLatitude.containsKey(processedLatitude)){
                newSelectDataByLatitude.put(processedLatitude, new ArrayList<>());
            }
            newSelectDataByLatitude.get(processedLatitude).add(data);
            if(!newSelectDataByLongitude.containsKey(processedLongitude)){
                newSelectDataByLongitude.put(processedLongitude, new ArrayList<>());
            }
            newSelectDataByLongitude.get(processedLongitude).add(data);

            //初始化latitudeExisted和longitudeExisted
            int latitudeSize = 0;
            for(int k = Util.transformLatitude(newMinLatitude); k <= Util.transformLatitude(newMaxLatitude); k++){
                newLatitudeExisted[k] = true;
                latitudeSize++;
            }
            int longitudeSize = 0;
            for(int k = Util.transformLongitude(newMinLongitude); k <= Util.transformLongitude(newMaxLongitude); k++){
                newLongitudeExisted[k] = true;
                longitudeSize++;
            }
            newSize[LATITUDE_COLUMN] = latitudeSize;
            newSize[LONGITUDE_COLUMN] = longitudeSize;
        }
        if(newMass == 0)return null;
        return new SubTensor(newMass, newMinTime, newMaxTime, newMinLongitude, newMaxLongitude, newMinLatitude, newMaxLatitude, newSize, newTimeExisted, newCoreDimensionExisted, newTimeSpaceExisted, newLatitudeExisted, newLongitudeExisted, newCoreDimensionMass, newTimeMass,  newTimeSpaceMass, newLatitudeMass, newLongitudeMass, newAttVals, newSpaceVals, newValue, newSelectDataByCoreDimension, newSelectDataByTime, newSelectDataByLatitude, newSelectDataByLongitude);
    }

    /**
     * 按照maxDensityDir执行真正的时空块扩展操作
     * @param subTensor
     * @param tensor
     * @param maxDensityDir
     */
    /*private static void extendSubtensorInDirect(SubTensor subTensor, Tensor tensor, int maxDensityDir) {
        //1、根据当前子张量保存的时空块序号计算出指定方向的时空块序号。
        List<Integer> appendTimeSpaceBlocks = findTimeSpaceBlockForExtension( subTensor, maxDensityDir);

        //2、使用tensor中的selectDataByTimeSpace取出指定方向时空块序号的数据。
        List<Integer> dataList = new ArrayList<>();
        for (Integer blockIndex : appendTimeSpaceBlocks) {
            if(!tensor.getTimeSpaceExisted()[blockIndex])continue;
            List<Integer> dataInTimeSpaceBlock = tensor.getSelectDataByTimeSpace().get(blockIndex);
            if(dataInTimeSpaceBlock != null){
                dataList.addAll(dataInTimeSpaceBlock);
                //4、从tensor中删除拼接的block
                tensor.getTimeSpaceExisted()[blockIndex] = false;
            }
        }

        //3、将数据加入到subtensor中
        subTensor.appendDatas(tensor, dataList);

    }*/

    /**
     * 按照maxDensityDir执行真正的时空块扩展操作
     * @param subTensor
     * @param tensor
     */
    private static void extendSubtensorInDirect(SubTensor subTensor, Tensor tensor, Set<Integer> addBlock) {
        List<Integer> dataInTimeSpaceBlock = new ArrayList<>();
        for (Integer blockIndex : addBlock) {
            dataInTimeSpaceBlock.addAll(tensor.getSelectDataByTimeSpace().get(blockIndex));

            //4、从tensor中删除拼接的block
            tensor.getTimeSpaceExisted()[blockIndex] = false;
            tensor.getSelectDataByTimeSpace().remove(blockIndex);
        }

        //3、将数据加入到subtensor中
        subTensor.appendDatas(tensor, dataInTimeSpaceBlock);
    }

    /**
     * 模拟拼接完毕后的子张量密度
     * @param subTensor
     * @param tensor
     * @param dir
     * @return
     */
    /*private static double preComputeDensityAfterExtension(SubTensor subTensor, Tensor tensor, int dir) {
        //1、根据当前子张量保存的时空块序号计算出指定方向的时空块序号。
        //调用findTimeSpaceBlockForExtension方法找到要扩展的时空块
        List<Integer> appendTimeSpaceBlocks = findTimeSpaceBlockForExtension( subTensor, dir);

        //2、使用tensor中的selectDataByTimeSpace取出指定方向时空块序号的数据。
        List<Integer> dataList = new ArrayList<>();
        for (Integer blockIndex : appendTimeSpaceBlocks) {
            if(!tensor.getTimeSpaceExisted()[blockIndex])continue;
            List<Integer> dataInTimeSpaceBlock = tensor.getSelectDataByTimeSpace().get(blockIndex);
            if(dataInTimeSpaceBlock != null){
                dataList.addAll(dataInTimeSpaceBlock);
            }
        }

        //3、将数据模拟加入到subtensor中，并计算密度。
        double densityAfterAppend = 0;
        if(dataList.size() != 0){
            densityAfterAppend = subTensor.imitateAppendDatas(tensor, dataList, dir);
        }

        //4、将密度返回。
        return densityAfterAppend;
    }*/

    /**
     * 模拟拼接完毕后的子张量密度
     * @param subTensor
     * @param tensor
     * @return
     */
    private static double preComputeDensityAfterExtension(SubTensor subTensor, Tensor tensor, int blockIndex) {
        List<Integer> dataInTimeSpaceBlock = tensor.getSelectDataByTimeSpace().get(blockIndex);
        if(dataInTimeSpaceBlock == null || dataInTimeSpaceBlock.size() == 0){
            return 0;
        }

        //3、将数据模拟加入到subtensor中，并计算密度。
        double densityAfterAppend = subTensor.imitateAppendDatas(tensor, dataInTimeSpaceBlock);

        //4、将密度返回。
        return densityAfterAppend;
    }


    public static void updateTensorByRemoveSubTensor(Tensor tensor, SubTensor subTensor) {
        //更新mass
        tensor.setMass(tensor.getMass() - subTensor.getMass());

        //更新value，删除某一个tuple就是把子张量中对应的value设置为0
        int[] value = tensor.getValue();
        int[] subValue = subTensor.getValue();
        for (int i = 0; i < value.length; i++) {
            if(subValue[i] != 0){
                value[i] = 0;
            }
        }

        try {
            //更新Mass每个维度质量信息
            for(int i = 0; i < CORE_COLUMN_LENGTH; i++){
                int[] coreDimensionMass = tensor.getCoreDimensionMass().get(i);
                int[] subTensorCoreDimensionMass = subTensor.getCoreDimensionMass().get(i);
                boolean[] subTensorCoreDimensionExisted = subTensor.getCoreDimensionExisted().get(i);
                tensor.getCoreDimensionMass().set(i, Util.arrSubstract(coreDimensionMass, subTensorCoreDimensionMass, subTensorCoreDimensionExisted));
            }
            tensor.setTimeMass(Util.arrSubstract(tensor.getTimeMass(), subTensor.getTimeMass(), subTensor.getTimeExisted()));
            tensor.setLatitudeMass(Util.arrSubstract(tensor.getLatitudeMass(), subTensor.getLatitudeMass(), subTensor.getLatitudeExisted()));
            tensor.setlongitudeMass(Util.arrSubstract(tensor.getlongitudeMass(), subTensor.getlongitudeMass(), subTensor.getLongitudeExisted()));
            tensor.setTimeSpaceMass(Util.arrSubstract(tensor.getTimeSpaceMass(), subTensor.getTimeSpaceMass(), tensor.getTimeSpaceExisted()));
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private static void updateCoreLikeMzoom(SubTensor subTensor){
        //创建一个用于遍历删除的压缩版subTensor
        SubTensor tempTensor = new SubTensor(subTensor);
        final IMinHeap[] heaps = createHeaps(subTensor);
        int mass = tempTensor.getMass();
        int[] cardinality = tempTensor.getSize().clone();
        final int[][] modeToAttVals = createModeToAttVals(cardinality, subTensor);
        int sumOfCardinality = 0;
        for(int i = 0; i < CORE_COLUMN_LENGTH; i++){
            sumOfCardinality += cardinality[i];
        }

        //copy coreDimensionMass
        List<int[]> coreDimensionMass = new ArrayList<>();
        List<int[]> subTensorCoreDimensionMass = subTensor.getCoreDimensionMass();
        for (int[] masses : subTensorCoreDimensionMass) {
            coreDimensionMass.add(masses.clone());
        }

        int delIndex = 0;
        int[] deleteModes = new int[sumOfCardinality];
        int[] deleteAttrs = new int[sumOfCardinality];

        int maxIters = 0;
        double maxDensityAmongIters = Util.computeDensity2(tempTensor);

        int i = 0;

        while(i < sumOfCardinality){
            //1、先对所有core维度进行按质量从小到大的排序
            //2、每次遍历三个维度中质量最小的属性并将其假删除并记录删除的属性和密度
            //3、记录这个过程中密度达到的最大值和maxItem，最后恢复到maxItem。
            int maxMode = 0;
            double maxDensityAmongModes = -Double.MIN_VALUE;
            for(int mode = 0; mode < CORE_COLUMN_LENGTH; mode++){
                final Pair<Integer, Integer> pair = heaps[mode].peek();
                if(pair != null){
                    double tempDensity = Util.computeDensity2(mass - pair.getValue(), Util.computeVolume2(tempTensor) - 1);
                    if(tempDensity >= maxDensityAmongModes){
                        maxMode = mode;
                        maxDensityAmongModes = tempDensity;
                    }
                }
            }
            Pair<Integer, Integer> pair = heaps[maxMode].poll();

            mass -= pair.getValue();
            double density = Util.computeDensity2(tempTensor.getMass() - pair.getValue(), Util.computeVolume2(tempTensor) - 1);
            //后续需要更新tempTensor中的Mass和Size
            tempTensor.setMass(tempTensor.getMass() - pair.getValue());
            int[] tempSize = tempTensor.getSize();
            tempSize[maxMode]--;
            tempTensor.setSize(tempSize);

            if(density > maxDensityAmongIters){
                maxDensityAmongIters = density;
                maxIters = i + 1;
            }
            deleteModes[delIndex] = maxMode;
            deleteAttrs[delIndex++] = pair.getKey();
            i++;

            //该删除只处理coreDimensions，并且不实时更新其它维度的size
            removeAndUpdateAttValMassesLikeMzoom(tempTensor, maxMode, pair.getKey(), coreDimensionMass, heaps);

        }

        //真正的子张量删除阶段
        for(int index = 0; index < maxIters; index++){
            deleteAttrInCoreUpdate(subTensor, deleteAttrs[index],deleteModes[index]);
        }
    }

    /**
     * 针对每一个子张量，通过减少用户或商品维度提升其密度。
     * @param subTensor 时间维度上密集的子张量
     * @return
     */
    private static void updateCoreLikeDcube(SubTensor subTensor){
        //创建一个用于遍历删除的压缩版subTensor
        SubTensor tempTensor = new SubTensor(subTensor);
        int mass = tempTensor.getMass();
        int[] cardinality = tempTensor.getSize().clone();
        final int[][] modeToAttVals = createModeToAttVals(cardinality, subTensor);
        int sumOfCardinality = 0;
        for(int i = 0; i < CORE_COLUMN_LENGTH; i++){
            sumOfCardinality += cardinality[i];
        }
        List<int[]> coreDimensionMass = new ArrayList<>();
        List<int[]> subTensorCoreDimensionMass = subTensor.getCoreDimensionMass();
        for (int[] masses : subTensorCoreDimensionMass) {
            coreDimensionMass.add(masses.clone());
        }
        int[] modeToAliveNum = tempTensor.getSize().clone();
        int[] modeToRemovedNum = new int[CORE_COLUMN_LENGTH];

        int delIndex = 0;
        int[] deleteModes = new int[sumOfCardinality];
        int[] deleteAttrs = new int[sumOfCardinality];

        int maxIters = 0;
        double maxDensityAmongIters = Util.computeDensity2(tempTensor);

        int i = 0;

        while(i < sumOfCardinality){
            int maxMode = 0;
            double maxDensityAmongModes = -Double.MIN_VALUE;
            for(int mode = 0; mode < CORE_COLUMN_LENGTH; mode++){
                if(modeToAliveNum[mode] > 0){
                    double threshold = mass * 1.0/ modeToAliveNum[mode];
                    int numToRemove = 0;
                    long removedMassSum = 0;
                    int[] attValToMass;
                    int[] attVals = modeToAttVals[mode];
                    // if(mode == 0){
                    //     attValToMass = userMass;
                    // }else{
                    //     attValToMass = objectMass;
                    // }
                    attValToMass = coreDimensionMass.get(mode);

                    for(int j = modeToRemovedNum[mode]; j < cardinality[mode]; j++){
                        int attVal = attVals[j];
                        if(attValToMass[attVal] <= threshold){
                            numToRemove++;
                            removedMassSum += attValToMass[attVal];
                        }
                    }

                    if(numToRemove >= 1){
                        double tempDensity = Util.computeDensity2(mass - removedMassSum, Util.computeVolume2(tempTensor) - numToRemove);
                        if(tempDensity >= maxDensityAmongModes){
                            maxMode = mode;
                            maxDensityAmongModes = tempDensity;
                        }
                    }else{
                        System.out.println("Sanity Check!");
                    }
                }
            }

            double threshold = mass * 1.0 / modeToAliveNum[maxMode];
            final int[] attValToMass;
            final boolean[] attValsToRemove = new boolean[cardinality[maxMode]];
            attValToMass = coreDimensionMass.get(maxMode);
            sort(modeToAttVals[maxMode], attValToMass, modeToRemovedNum[maxMode], cardinality[maxMode] - 1);
            int[] attVals = modeToAttVals[maxMode];

            for(int j = modeToRemovedNum[maxMode]; j < cardinality[maxMode]; j++){
                int attVal = attVals[j];
                if(attValToMass[attVal] <= threshold){
                    mass -= attValToMass[attVal];
                    double density = Util.computeDensity2(tempTensor.getMass() - attValToMass[attVal], Util.computeVolume2(tempTensor) - 1);
                    //后续需要更新tempTensor中的Mass和Size
                    tempTensor.setMass(tempTensor.getMass() - attValToMass[attVal]);
                    int[] tempSize = tempTensor.getSize();
                    tempSize[maxMode]--;
                    tempTensor.setSize(tempSize);

                    if(density > maxDensityAmongIters){
                        maxDensityAmongIters = density;
                        maxIters = i + 1;
                    }
                    modeToRemovedNum[maxMode]++;
                    modeToAliveNum[maxMode]--;
                    deleteModes[delIndex] = maxMode;
                    deleteAttrs[delIndex++] = attVal;
                    i++;
                    attValsToRemove[j] = true;
                }else {
                    break;
                }
            }

            removeAndUpdateAttValMassesLikeDcube(tempTensor, maxMode, attVals, attValsToRemove, coreDimensionMass);

        }

        //真正的子张量删除阶段
        for(int index = 0; index < maxIters; index++){
            deleteAttrInCoreUpdate(subTensor, deleteAttrs[index],deleteModes[index]);
        }
    }

    static class AttrComparator implements Comparator<Integer>{

        int[] value;

        double[][] spaceVals;

        public AttrComparator(int[] value, double[][] spaceVals) {
            this.value = value;
            this.spaceVals = spaceVals;
        }

        @Override
        public int compare(Integer o1, Integer o2) {
            if(value[o1] == 0 && value[o2] != 0){
                return -1;
            }else if(value[o1] != 0 && value[o2] == 0){
                return 1;
            }else if(value[o1] == 0 && value[o2] == 0){
                return 0;
            }else{
                return -Double.compare(spaceVals[1][o1], spaceVals[1][o2]);
            }
        }
    }

    private static void updateCoreLikeMbiz(SubTensor subTensor) {
        updateCoreLikeMzoom(subTensor);

        int[] globalValue = subTensor.getValue().clone();
        List<Map<Integer, List<Integer>>> initSelectDataByCoreDimension = subTensor.getSelectDataByCoreDimension();
        List<Map<Integer, List<Integer>>> globalAttributeToValuesToTuples = new ArrayList<>();
        for (Map<Integer, List<Integer>> map : initSelectDataByCoreDimension) {
            Map<Integer, List<Integer>> copyMap = new HashMap<>();
            for (Integer key : map.keySet()) {
                List<Integer> copyList = new ArrayList<>();
                for (Integer value: map.get(key)){
                    copyList.add(value);
                }
                copyMap.put(key, copyList);
            }
            globalAttributeToValuesToTuples.add(copyMap);
        }

        final byte[] nonmemberCounts = new byte[globalValue.length];
        int[][] attributes = subTensor.getAttVals();
        List<boolean[]> coreDimensionExisted = subTensor.getCoreDimensionExisted();

        int[][] attributeToValueToMassChange = new int[CORE_COLUMN_LENGTH][];
        for (int mode = 0; mode < CORE_COLUMN_LENGTH; mode++) {
            attributeToValueToMassChange[mode] = new int[globalSize[mode] + 1];
        }

        for(int i=0; i< globalValue.length; i++){
            int measureValue = globalValue[i];
            if(measureValue == 0)continue;

            //记录有多少tuple已经被删除了
            byte nonmemberCount = 0;
            int nonmemberMode = 0;
            for(int mode=0; mode<CORE_COLUMN_LENGTH; mode++){
                if(!coreDimensionExisted.get(mode)[attributes[mode][i]]) {
                    nonmemberCount += 1;
                    nonmemberMode = mode;
                }
            }

            nonmemberCounts[i] = nonmemberCount;

            if(nonmemberCount==0) { // in the block
                for(int mode=0; mode<CORE_COLUMN_LENGTH; mode++){
                    attributeToValueToMassChange[mode][attributes[mode][i]] += measureValue;
                }
            }
            else if(nonmemberCount == 1){
                attributeToValueToMassChange[nonmemberMode][attributes[nonmemberMode][i]] += measureValue;
            }
        }

        final IMinHeap[] inHeaps = new IMinHeap[CORE_COLUMN_LENGTH];
        final IMaxHeap[] outHeaps = new IMaxHeap[CORE_COLUMN_LENGTH];
        final int[] attributeToCardinalities = new int[CORE_COLUMN_LENGTH];
        int sumOfCardinalities = 0;

        //初始化大小顶堆
        for(int mode = 0; mode < CORE_COLUMN_LENGTH; mode++) {
            boolean[] attValToBeIncluded = coreDimensionExisted.get(mode);
            int[] attValToMassChange = attributeToValueToMassChange[mode];
            IMinHeap inHeap = new HashIndexedMinHeap(globalSize[mode] + 1);
            IMaxHeap outHeap = new HashIndexedMaxHeap(globalSize[mode] + 1);
            for(int index = 1; index < globalSize[mode] + 1; index++) {
                if(attValToBeIncluded[index]) {
                    inHeap.insert(index, attValToMassChange[index]);
                    attributeToCardinalities[mode] += 1;
                    sumOfCardinalities ++;
                }
                else {
                    outHeap.insert(index, attValToMassChange[index]);
                }
            }
            inHeaps[mode] = inHeap;
            outHeaps[mode] = outHeap;
        }

        double currentScore = Util.computeDensity2(subTensor);
        while(true) {

            double previousScore = currentScore;
            int maxMode = 0;
            boolean action = false; //false: remove, true: insert
            // if(sumOfCardinalities > lower) {
            if(sumOfCardinalities > 0) {
                for (int mode = 0; mode < CORE_COLUMN_LENGTH; mode++) {
                    final Pair<Integer, Integer> pair = inHeaps[mode].peek();
                    if (pair != null) {
                        double tempScore = Util.computeDensity2(subTensor.getMass() - pair.getValue(), Util.computeVolume2(subTensor) - 1);
                        // double tempScore = measure.ifRemoved(attribute, 1, pair.getValue());
                        if (tempScore > currentScore) {
                            maxMode = mode;
                            action = false;
                            currentScore = tempScore;
                        }
                    }
                }
            }

            if(sumOfCardinalities < Integer.MAX_VALUE) {
                for (int mode = 0; mode < CORE_COLUMN_LENGTH; mode++) {
                    final Pair<Integer, Integer> pair = outHeaps[mode].peek();
                    if (pair != null) {
                        double tempScore = Util.computeDensity2(subTensor.getMass() + pair.getValue(), Util.computeVolume2(subTensor) + 1);
                        // double tempScore = measure.ifInserted(attribute, 1, pair.getValue());
                        if (tempScore > currentScore) {
                            maxMode = mode;
                            action = true;
                            currentScore = tempScore;
                        }
                    }
                }
            }

            if(currentScore == previousScore) { // terminates
                break;
            }
            int[] values = subTensor.getValue();
            // int maxTime = subTensor.getMaxTime();
            // int minTime = subTensor.getMinTime();
            // double maxLatitude = subTensor.getMaxLatitude();
            // double minLatitude = subTensor.getMinLatitude();
            // double maxLongitude = subTensor.getMaxLongitude();
            // double minLongitude = subTensor.getMinLongitude();
            List<int[]> coreDimensionMass = subTensor.getCoreDimensionMass();

            //更新用户维度属性存在状态
            boolean[] timeExisted = subTensor.getTimeExisted();
            boolean[] latitudeExisted = subTensor.getLatitudeExisted();
            boolean[] longitudeExisted = subTensor.getLongitudeExisted();

            double[][] spaceVals = subTensor.getSpaceVals();
            int[] size = subTensor.getSize();
            Map<Integer, List<Integer>> selectDataByTime = subTensor.getSelectDataByTime();
            Map<Integer, List<Integer>> selectDataByLatitude = subTensor.getSelectDataByLatitude();
            Map<Integer, List<Integer>> selectDataByLongitude = subTensor.getSelectDataByLongitude();
            int[][] attVals = subTensor.getAttVals();
            int[] timeMass = subTensor.getTimeMass();
            int[] latitudeMasses = subTensor.getLatitudeMass();
            int[] longitudeMasses = subTensor.getlongitudeMass();

            if(action == false) { //remove
                Pair<Integer, Integer> pair = inHeaps[maxMode].poll();
                int attValToRemove = pair.getKey();
                //更新mass
                subTensor.setMass(subTensor.getMass() - pair.getValue());

                sumOfCardinalities--;
                attributeToCardinalities[maxMode]--;

                //update degree in
                int massSumOut = 0;
                List<Integer> tuples = globalAttributeToValuesToTuples.get(maxMode).get(attValToRemove);
                for (int tuple : tuples) {
                    if(globalValue[tuple] == 0)continue;
                    byte nonmemberCount = nonmemberCounts[tuple];

                    if(nonmemberCount > 1) {
                        nonmemberCounts[tuple]++;
                    }
                    else if(nonmemberCount == 1){
                        int measureValue = values[tuple];
                        int nonMemberMode = 0;
                        for(int mode = 0; mode < CORE_COLUMN_LENGTH; mode++) {
                            if(mode != maxMode) {
                                if(!coreDimensionExisted.get(mode)[attributes[mode][tuple]]) {
                                    nonMemberMode = mode;
                                    break;
                                }
                            }
                        }
                        nonmemberCounts[tuple]++;
                        int attVal = attributes[nonMemberMode][tuple];
                        outHeaps[nonMemberMode].updatePriority(attVal, outHeaps[nonMemberMode].getPriority(attVal) - measureValue);
                    }
                    else if(nonmemberCount==0){ //nonmember count == 0;
                        int measureValue = globalValue[tuple];
                        massSumOut += measureValue;
                        nonmemberCounts[tuple]++;
                        for (int mode = 0; mode < globalSize.length - 1; mode++) {
                            int attr = 0;
                            if(mode < Converge.CORE_COLUMN_LENGTH){
                                if(mode != maxMode) {
                                    int attVal = attributes[mode][tuple];
                                    inHeaps[mode].updatePriority(attVal, inHeaps[mode].getPriority(attVal) - measureValue);
                                }
                                attr = attVals[mode][tuple];
                                int[] coreMasses = coreDimensionMass.get(mode);
                                coreMasses[attr] -= globalValue[tuple];
                            }else if(mode == Converge.TIME_COLUMN){
                                attr = Util.transformTime(attVals[mode][tuple]);
                                timeMass[attr] -= globalValue[tuple];
                            }else if(mode == Converge.LATITUDE_COLUMN){
                                attr = Util.transformLatitude(spaceVals[0][tuple]);
                                latitudeMasses[attr] -= globalValue[tuple];
                            }else if(mode == Converge.LONGITUDE_COLUMN){
                                attr = Util.transformLongitude(spaceVals[1][tuple]);
                                longitudeMasses[attr] -= globalValue[tuple];
                            }

                            // AttrComparator comparator = new AttrComparator(values, spaceVals);

                            if(mode >= Converge.CORE_COLUMN_LENGTH){
                                if(mode == Converge.TIME_COLUMN){
                                    if(timeMass[attr] == 0){
                                        timeExisted[attr] = false;
                                        size[TIME_COLUMN]--;
                                        selectDataByTime.remove(attr);
                                    }

                                }else if(mode == Converge.LATITUDE_COLUMN){
                                    if(latitudeMasses[attr] == 0){
                                        latitudeExisted[attr] = false;
                                        size[LATITUDE_COLUMN]--;
                                        selectDataByLatitude.remove(attr);
                                    }
                                }else if(mode == Converge.LONGITUDE_COLUMN){
                                    if(longitudeMasses[attr] == 0){
                                        longitudeExisted[attr] = false;
                                        size[LONGITUDE_COLUMN]--;
                                        selectDataByLongitude.remove(attr);
                                    }

                                }
                            }
                        }

                        values[tuple] = 0;
                    }
                    else {
                        System.out.println("error-2!");
                        System.exit(0);
                    }
                }
                outHeaps[maxMode].insert(attValToRemove, massSumOut);
                coreDimensionExisted.get(maxMode)[attValToRemove] = false;
                size[maxMode]--;

            }
            else { //insert
                Pair<Integer, Integer> pair = outHeaps[maxMode].poll();
                int attValToInsert = pair.getKey();
                //更新mass
                subTensor.setMass(subTensor.getMass() + pair.getValue());


                sumOfCardinalities++;
                attributeToCardinalities[maxMode]++;

                //update degree in
                int massSumIn = 0;
                List<Integer> tuples = globalAttributeToValuesToTuples.get(maxMode).get(attValToInsert);
                for (int tuple : tuples) {
                    if(globalValue[tuple] == 0)continue;
                    byte nonmemberCount = nonmemberCounts[tuple];

                    if(nonmemberCount > 2) {
                        nonmemberCounts[tuple]--;
                        continue;
                    }
                    else if(nonmemberCount == 2){
                        int measureValue = globalValue[tuple];
                        int nonMemberMode = 0;
                        for(int mode = 0; mode < CORE_COLUMN_LENGTH; mode++) {
                            if(mode != maxMode) {
                                if(!coreDimensionExisted.get(mode)[attributes[mode][tuple]]) {
                                    nonMemberMode = mode;
                                    break;
                                }
                            }
                        }
                        nonmemberCounts[tuple]--;
                        int index = attributes[nonMemberMode][tuple];
                        outHeaps[nonMemberMode].updatePriority(index, outHeaps[nonMemberMode].getPriority(index) + measureValue);
                    }
                    else if(nonmemberCount == 1){ //nonmember count == 1;
                        int measurevALUE = globalValue[tuple];
                        massSumIn += measurevALUE;
                        nonmemberCounts[tuple]--;
                        for (int mode = 0; mode < globalSize.length - 1; mode++) {
                            int attr = 0;
                            if(mode < Converge.CORE_COLUMN_LENGTH){
                                if(mode != maxMode) {
                                    int index = attributes[mode][tuple];
                                    inHeaps[mode].updatePriority(index, inHeaps[mode].getPriority(index) + measurevALUE);
                                }
                                attr = attVals[mode][tuple];
                                int[] coreMasses = coreDimensionMass.get(mode);
                                coreMasses[attr] += globalValue[tuple];
                            }else if(mode == Converge.TIME_COLUMN){
                                attr = Util.transformTime(attVals[mode][tuple]);
                                timeMass[attr] += globalValue[tuple];
                                // minTime = Math.min(minTime, attVals[mode][tuple]);
                                // maxTime = Math.max(maxTime, attVals[mode][tuple]);
                                if (!timeExisted[attr]){
                                    size[TIME_COLUMN]++;
                                    timeExisted[attr] = true;
                                }
                                selectDataByTime.putIfAbsent(attr, new ArrayList<>());
                                selectDataByTime.get(attr).add(tuple);
                            }else if(mode == Converge.LATITUDE_COLUMN){
                                attr = Util.transformLatitude(spaceVals[0][tuple]);
                                latitudeMasses[attr] += globalValue[tuple];
                                // maxLatitude = Math.max(spaceVals[0][tuple], maxLatitude);
                                // minLatitude = Math.min(spaceVals[0][tuple], minLatitude);
                                if (!latitudeExisted[attr]) {
                                    size[LATITUDE_COLUMN]++;
                                    latitudeExisted[attr] = true;
                                }
                                selectDataByTime.putIfAbsent(attr, new ArrayList<>());
                                selectDataByTime.get(attr).add(tuple);
                            }else if(mode == Converge.LONGITUDE_COLUMN){
                                attr = Util.transformLongitude(spaceVals[1][tuple]);
                                longitudeMasses[attr] += globalValue[tuple];
                                // minLongitude = Math.min(spaceVals[1][tuple], minLongitude);
                                // maxLongitude = Math.max(spaceVals[1][tuple], maxLongitude);
                                if (!longitudeExisted[attr]) {
                                    size[LONGITUDE_COLUMN]++;
                                    longitudeExisted[attr] = true;
                                }
                                selectDataByTime.putIfAbsent(attr, new ArrayList<>());
                                selectDataByTime.get(attr).add(tuple);
                            }
                        }
                        values[tuple] = globalValue[tuple];

                    }
                    else { //nonmemberCount
                        System.out.println(coreDimensionExisted.get(maxMode)[attributes[maxMode][tuple]]);
                        System.out.println("error-1!");
                        System.exit(0);
                    }
                }
                inHeaps[maxMode].insert(attValToInsert, massSumIn);
                coreDimensionExisted.get(maxMode)[attValToInsert] = true;
                size[maxMode]++;
            }

            //更新最小最大经纬度信息
            // subTensor.setMinTime(minTime);
            // subTensor.setMaxTime(maxTime);
            // subTensor.setMinLatitude(minLatitude);
            // subTensor.setMaxLatitude(maxLatitude);
            // subTensor.setMinLongitude(minLongitude);
            // subTensor.setMaxLongitude(maxLongitude);
            currentScore = Util.computeDensity2(subTensor);

        }

    }

    private static void removeAndUpdateAttValMassesLikeDcube(SubTensor tensor, int maxMode, int[] attVals, boolean[] attValsToRemove, List<int[]> coreDimensionMass) {
        int len = attValsToRemove.length;
        List<Integer> dataList;
        int[][] attValsInMode = tensor.getAttVals();
        List<int[]> originalCoreDimensionMass = tensor.getCoreDimensionMass();
        int[] originalMainMasses = originalCoreDimensionMass.get(maxMode);
        int[] mainMass = coreDimensionMass.get(maxMode);
        // int[] originalObjectMass = tensor.getObjectMass();
        // int[] originalUserMass = tensor.getUserMass();
        int[] value = tensor.getValue();
        //1、找到要删除的所有tuple
        //2、根据tuple找到要质量减小的userMass或objectMass
        //3、将tuple的value置为0，将对于删除的userMass或ObjectMass置为0
        for(int i = 0; i < len; i++){
            if(attValsToRemove[i]){
                int index = attVals[i];
                mainMass[index] = 0;
                originalMainMasses[index] = 0;
                Map<Integer, List<Integer>> selectDataByCore = tensor.getSelectDataByCoreDimension().get(maxMode);
                dataList = selectDataByCore.get(index);
                if(dataList == null)continue;
                for(int data : dataList){
                    if(value[data] == 0)continue;
                    for(int otherMode = 0; otherMode < coreDimensionMass.size(); otherMode++){
                        if(otherMode == maxMode)continue;
                        int attval = attValsInMode[otherMode][data];
                        int[] otherMasses = coreDimensionMass.get(otherMode);
                        otherMasses[attval] -= value[data];
                        int[] originalOtherMasses = originalCoreDimensionMass.get(otherMode);
                        originalOtherMasses[attval] -= value[data];
                    }
                    // objectMass[attval] -= value[data];
                    // originalObjectMass[attval] -= value[data];
                    value[data] = 0;
                }
            }
        }
    }

    private static IMinHeap[] createHeaps(SubTensor subTensor) {
        List<int[]> coreDimensionMass = subTensor.getCoreDimensionMass();
        IMinHeap[] heaps = new IMinHeap[CORE_COLUMN_LENGTH];
        for (int i = 0; i < CORE_COLUMN_LENGTH; i++) {
            boolean[] coreDimensionExisted = subTensor.getCoreDimensionExisted().get(i);
            IMinHeap heap = new HashIndexedMinHeap(globalSize[i] + 1);
            int[] masses = coreDimensionMass.get(i);
            for(int index = 0; index < coreDimensionExisted.length; index++) {
                if(coreDimensionExisted[index]){
                    heap.insert(index, masses[index]);
                }
            }
            heaps[i]  = heap;
        }
        return heaps;
    }

    private static void removeAndUpdateAttValMassesLikeMzoom(SubTensor tensor, int maxMode, int attVal, List<int[]> coreDimensionMass, IMinHeap[] heaps) {
        List<Integer> dataList;
        int[][] attValsInMode = tensor.getAttVals();
        List<int[]> originalCoreDimensionMass = tensor.getCoreDimensionMass();
        int[] originalMainMasses = originalCoreDimensionMass.get(maxMode);
        int[] mainMass = coreDimensionMass.get(maxMode);
        // int[] originalObjectMass = tensor.getObjectMass();
        // int[] originalUserMass = tensor.getUserMass();
        int[] value = tensor.getValue();
        //1、找到要删除的所有tuple
        //2、根据tuple找到要质量减小的userMass或objectMass
        //3、将tuple的value置为0，将对于删除的userMass或ObjectMass置为0
        mainMass[attVal] = 0;
        originalMainMasses[attVal] = 0;
        Map<Integer, List<Integer>> selectDataByCore = tensor.getSelectDataByCoreDimension().get(maxMode);
        dataList = selectDataByCore.get(attVal);
        if(dataList == null)return;
        for(int data : dataList){
            if(value[data] == 0)continue;
            for(int otherMode = 0; otherMode < coreDimensionMass.size(); otherMode++){
                if(otherMode == maxMode)continue;
                int att = attValsInMode[otherMode][data];
                int[] otherMasses = coreDimensionMass.get(otherMode);
                if(otherMode < CORE_COLUMN_LENGTH){
                    heaps[otherMode].updatePriority(att, otherMasses[att] - value[data]);
                }
                otherMasses[att] -= value[data];
                int[] originalOtherMasses = originalCoreDimensionMass.get(otherMode);
                originalOtherMasses[att] -= value[data];
            }
            value[data] = 0;
        }
    }

    private static void removeAndUpdateAttValMassesLikeMbiz(SubTensor tensor, int maxMode, int attVal, List<int[]> coreDimensionMass, IMinHeap[] heaps) {
        List<Integer> dataList;
        int[][] attValsInMode = tensor.getAttVals();
        List<int[]> originalCoreDimensionMass = tensor.getCoreDimensionMass();
        int[] originalMainMasses = originalCoreDimensionMass.get(maxMode);
        int[] mainMass = coreDimensionMass.get(maxMode);
        // int[] originalObjectMass = tensor.getObjectMass();
        // int[] originalUserMass = tensor.getUserMass();
        int[] value = tensor.getValue();
        //1、找到要删除的所有tuple
        //2、根据tuple找到要质量减小的userMass或objectMass
        //3、将tuple的value置为0，将对于删除的userMass或ObjectMass置为0
        mainMass[attVal] = 0;
        originalMainMasses[attVal] = 0;
        Map<Integer, List<Integer>> selectDataByCore = tensor.getSelectDataByCoreDimension().get(maxMode);
        dataList = selectDataByCore.get(attVal);
        if(dataList == null)return;
        for(int data : dataList){
            if(value[data] == 0)continue;
            for(int otherMode = 0; otherMode < coreDimensionMass.size(); otherMode++){
                if(otherMode == maxMode)continue;
                int att = attValsInMode[otherMode][data];
                int[] otherMasses = coreDimensionMass.get(otherMode);
                if(otherMode < CORE_COLUMN_LENGTH){
                    heaps[otherMode].updatePriority(att, otherMasses[att] - value[data]);
                }
                otherMasses[att] -= value[data];
                int[] originalOtherMasses = originalCoreDimensionMass.get(otherMode);
                originalOtherMasses[att] -= value[data];
            }
            value[data] = 0;
        }
    }

    /**
     * 将attributes数组中的属性值，按照其切片质量排序
     * @param attributes
     * @param masses
     * @param left
     * @param right
     */
    private static void sort(int[] attributes, int[] masses, int left, int right) {

        if (attributes == null || attributes.length == 0)
            return;

        if (left >= right)
            return;

        int middle = left + (right - left) / 2;
        int pivot = masses[attributes[middle]];

        int i = left, j = right;
        while (i <= j) {
            while (masses[attributes[i]] < pivot) {
                i++;
            }

            while (masses[attributes[j]] > pivot) {
                j--;
            }

            if (i <= j) {
                int temp = attributes[i];
                attributes[i] = attributes[j];
                attributes[j] = temp;
                i++;
                j--;
            }
        }

        if (left < j)
            sort(attributes, masses, left, j);

        if (right > i)
            sort(attributes, masses, i, right);
    }

    private static void sort(Integer[] attributes, int[] masses, int left, int right) {
        if (attributes == null || attributes.length == 0)
            return;

        if (left >= right)
            return;

        int middle = left + (right - left) / 2;
        int pivot = masses[attributes[middle]];

        int i = left, j = right;
        while (i <= j) {
            while (masses[attributes[i]] > pivot) {
                i++;
            }

            while (masses[attributes[j]] < pivot) {
                j--;
            }

            if (i <= j) {
                int temp = attributes[i];
                attributes[i] = attributes[j];
                attributes[j] = temp;
                i++;
                j--;
            }
        }

        if (left < j)
            sort(attributes, masses, left, j);

        if (right > i)
            sort(attributes, masses, i, right);
    }

    /**
     * 创建一个二维数组，按照顺序存放质量不为0的属性索引
     * @return
     */
    private static int[][] createModeToAttVals(int[] cardinality, SubTensor subTensor) {
        int[][] modeToIndices = new int[CORE_COLUMN_LENGTH][];
        // int[] userMass = subTensor.getUserMass();
        // int[] objectMass = subTensor.getObjectMass();
        for(int mode = 0; mode < CORE_COLUMN_LENGTH; mode++){
            int[] masses = subTensor.getCoreDimensionMass().get(mode);
            int[] indices = new int[cardinality[mode]];
            int num = 0;
            for(int index = 1; index < masses.length; index++){
                if(masses[index] != 0){
                    indices[num++] = index;
                }
            }
            modeToIndices[mode] = indices;
        }
        return modeToIndices;
    }

    /**
     * 删除子张量中i商品维度
     * @param subTensor 子张量信息
     * @param i 删除的用户维度
     */
    public static void deleteAttrInCore(SubTensor subTensor, int i, int mode){
        // int[] objectMass = subTensor.getObjectMass();
        List<int[]> coreDimensionMass = subTensor.getCoreDimensionMass();
        int[] mainMasses = coreDimensionMass.get(mode);
        //删除子张量用户维度该属性,值得一提的是删除操作不修改attVals和value属性，只需要把value属性中的对应记录设置为0即可。
        //更新mass
        subTensor.setMass(subTensor.getMass() - mainMasses[i]);
        //更新用户维度规模
        int[] size = subTensor.getSize();
        size[mode]--;

        //更新用户维度属性存在状态
        boolean[] coreExisted = subTensor.getCoreDimensionExisted().get(mode);
        // boolean[] objectExisted = subTensor.getObjectExisted();
        coreExisted[i] = false;

        //更新Mass
        Map<Integer, List<Integer>> selectDataByCore = subTensor.getSelectDataByCoreDimension().get(mode);
        List<Integer> dataList = selectDataByCore.get(i);
        int[][] attVals = subTensor.getAttVals();
        // int[] userMass = subTensor.getUserMass();
        // int[] objectMass = subTensor.getObjectMass();
        int[] timeMass = subTensor.getTimeMass();
        int[] value = subTensor.getValue();
        if(dataList != null){
            for (Integer data : dataList) {
                if(value[data] == 0)continue;
                for(int coreMode = 0; coreMode < coreDimensionMass.size(); coreMode++){
                    int coreAttr = attVals[coreMode][data];
                    int[] coreMasses = coreDimensionMass.get(coreMode);
                    coreMasses[coreAttr] -= value[data];
                }
                // int userAtt = attVals[0][data];
                // userMass[userAtt] -= value[data];
                int timeAtt = attVals[TIME_COLUMN][data];
                timeMass[timeAtt] -= value[data];
                value[data] = 0;

                //TODO:需要检查最大最小经纬度和时间是否发生变化
            }
        }
        selectDataByCore.remove(i);
    }

    /**
     * 删除子张量中i商品维度
     * @param subTensor 子张量信息
     * @param i 删除的用户维度
     */
    public static void deleteAttrInCoreUpdate(SubTensor subTensor, int i, int mode){
        int maxTime = subTensor.getMaxTime();
        int minTime = subTensor.getMinTime();
        double maxLatitude = subTensor.getMaxLatitude();
        double minLatitude = subTensor.getMinLatitude();
        double maxLongitude = subTensor.getMaxLongitude();
        double minLongitude = subTensor.getMinLongitude();

        List<int[]> coreDimensionMass = subTensor.getCoreDimensionMass();
        int[] mainMasses = coreDimensionMass.get(mode);
        double[][] spaceVals = subTensor.getSpaceVals();
        //删除子张量用户维度该属性,值得一提的是删除操作不修改attVals和value属性，只需要把value属性中的对应记录设置为0即可。
        //更新mass
        subTensor.setMass(subTensor.getMass() - mainMasses[i]);
        //更新用户维度规模
        int[] size = subTensor.getSize();
        // size[mode]--;

        //更新用户维度属性存在状态
        List<boolean[]> coreDimensionExisted = subTensor.getCoreDimensionExisted();
        boolean[] timeExisted = subTensor.getTimeExisted();
        boolean[] latitudeExisted = subTensor.getLatitudeExisted();
        boolean[] longitudeExisted = subTensor.getLongitudeExisted();
        // coreExisted[i] = false;

        //更新Mass
        Map<Integer, List<Integer>> selectDataByTime = subTensor.getSelectDataByTime();
        Map<Integer, List<Integer>> selectDataByLatitude = subTensor.getSelectDataByLatitude();
        Map<Integer, List<Integer>> selectDataByLongitude = subTensor.getSelectDataByLongitude();
        Map<Integer, List<Integer>> selectDataByCore = subTensor.getSelectDataByCoreDimension().get(mode);
        List<Integer> dataList = selectDataByCore.get(i);
        int[][] attVals = subTensor.getAttVals();
        // int[] userMass = subTensor.getUserMass();
        // int[] objectMass = subTensor.getObjectMass();
        int[] timeMass = subTensor.getTimeMass();
        int[] latitudeMasses = subTensor.getLatitudeMass();
        int[] longitudeMasses = subTensor.getlongitudeMass();
        int[] value = subTensor.getValue();
        if(dataList != null){
            for (Integer data : dataList) {
                if(value[data] == 0)continue;
                for(int m = 0; m < globalSize.length - 1; m++){
                    int attr = 0;
                    if(m < Converge.CORE_COLUMN_LENGTH){
                        attr = attVals[m][data];
                        int[] coreMasses = coreDimensionMass.get(m);
                        coreMasses[attr] -= value[data];
                    }else if(m == Converge.TIME_COLUMN){
                        attr = Util.transformTime(attVals[m][data]);
                        timeMass[attr] -= value[data];
                    }else if(m == Converge.LATITUDE_COLUMN){
                        attr = Util.transformLatitude(spaceVals[0][data]);
                        latitudeMasses[attr] -= value[data];
                    }else if(m == Converge.LONGITUDE_COLUMN){
                        attr = Util.transformLongitude(spaceVals[1][data]);
                        longitudeMasses[attr] -= value[data];
                    }

                    if(m < Converge.CORE_COLUMN_LENGTH && coreDimensionMass.get(m)[attr] == 0){
                        //通过切片判断该维度的该属性是否存在，只有存在才修改其size
                        if(coreDimensionExisted.get(m)[attr]){
                            //对于经纬度和时间不需要更改其size
                            size[m]--;
                        }
                        coreDimensionExisted.get(m)[attr] = false;
                    }else if(m == Converge.TIME_COLUMN && timeMass[attr] == 0){
                        timeExisted[attr] = false;
                    }else if(m == Converge.LATITUDE_COLUMN && latitudeMasses[attr] == 0){
                        latitudeExisted[attr] = false;
                    }else if(m == Converge.LONGITUDE_COLUMN && longitudeMasses[attr] == 0){
                        longitudeExisted[attr] = false;
                    }

                    // AttrComparator comparator = new AttrComparator(value, spaceVals);

                    //更新最小最大经纬度信息
                    /*if(m == Converge.TIME_COLUMN){
                        if(attVals[mode][data] == minTime){
                            for(int time = Util.transformTime(minTime); time <= Util.transformTime(maxTime); time++){
                                if(timeExisted[time]){
                                    List<Integer> dataIndexsInTime = selectDataByTime.get(time);
                                    Integer minTimeIndex = dataIndexsInTime.stream().min(comparator).get();
                                    minTime = attVals[mode][minTimeIndex];
                                    break;
                                }
                            }
                        }else if(attVals[mode][data] == maxTime){
                            for(int time = Util.transformTime(maxTime); time >= Util.transformTime(minTime); time--){
                                if(timeExisted[time]){
                                    List<Integer> dataIndexsInTime = selectDataByTime.get(time);
                                    Integer maxTimeIndex = dataIndexsInTime.stream().max(comparator).get();
                                    maxTime = attVals[mode][maxTimeIndex];
                                    break;
                                }
                            }
                        }
                    }else if(m == Converge.LATITUDE_COLUMN){
                        if(spaceVals[0][data] == minLatitude){
                            for(int latitude = Util.transformLatitude(minLatitude); latitude <= Util.transformLatitude(maxLatitude); latitude++){
                                if(latitudeExisted[latitude]){
                                    List<Integer> dataIndexsInLatitude = selectDataByLatitude.get(latitude);
                                    Integer minLatitudeIndex = dataIndexsInLatitude.stream().min(comparator).get();
                                    minLatitude = attVals[mode][minLatitudeIndex];
                                    break;
                                }
                            }
                        }else if(spaceVals[0][data] == maxLatitude){
                            for(int latitude = Util.transformLatitude(maxLatitude); latitude >= Util.transformLatitude(minLatitude); latitude--){
                                if(latitudeExisted[latitude]){
                                    List<Integer> dataIndexsInLatitude = selectDataByLatitude.get(latitude);
                                    Integer maxLatitudeIndex = dataIndexsInLatitude.stream().max(comparator).get();
                                    maxLatitude = attVals[mode][maxLatitudeIndex];
                                    break;
                                }
                            }
                        }
                    }else if(m == Converge.LONGITUDE_COLUMN){
                        if(spaceVals[1][data] == minLongitude){
                            for(int longitude = Util.transformLongitude(minLongitude); longitude <= Util.transformLongitude(maxLongitude); longitude++){
                                if(longitudeExisted[longitude]){
                                    List<Integer> dataIndexsInLongitude = selectDataByLongitude.get(longitude);
                                    Integer minLongitudeIndex = dataIndexsInLongitude.stream().min(comparator).get();
                                    minLongitude = attVals[mode][minLongitudeIndex];
                                    break;
                                }
                            }
                        }else if(spaceVals[1][data] == maxLongitude){
                            for(int longitude = Util.transformLongitude(maxLongitude); longitude >= Util.transformLongitude(minLongitude); longitude--){
                                if(longitudeExisted[longitude]){
                                    List<Integer> dataIndexsInLongitude = selectDataByLongitude.get(longitude);
                                    Integer maxLongitudeIndex = dataIndexsInLongitude.stream().max(comparator).get();
                                    maxLongitude = attVals[mode][maxLongitudeIndex];
                                    break;
                                }
                            }
                        }
                    }*/


                }
                value[data] = 0;

            }
        }
        selectDataByCore.remove(i);
    }

    /**
     * 整合处理子张量的操作，包括对其它两个维度的贪心，判断是否为前k个密集子张量，以及从tensor中去除
     * @param tensor
     * @param subTensor
     * @param k
     * @param subTensorPriorityQueue
     */
    public static void processCore(Tensor tensor, SubTensor subTensor, int k, PriorityQueue<Pair<SubTensor, Double>> subTensorPriorityQueue, int circleNum){
        Double density1 = Util.computeDensity(subTensor);
        // Double density2 = subTensorPriorityQueue.peek() == null ? 0.0 : subTensorPriorityQueue.peek().getValue();
        // //第二轮后初始重量小于队列中最小子张量初始mass的不进行update操作，降低时间复杂度
        // //第一轮中初始mass不到队列中最小张量一半的不进行update操作
        // if(subTensorPriorityQueue.size() == k && (circleNum > 0 && density1 <= density2 || circleNum == 0 && density1 <= density2 / 2)){
        //     updateTensorByRemoveSubTensor(tensor, subTensor);
        //     return;
        // }

        // updateCoreLikeDcube(subTensor);
        // updateCoreLikeMzoom(subTensor);
        updateCoreLikeMbiz(subTensor);
        //TODO:添加了isSingle函数忽略1X1X1的子张量
        /*if(!isSingle(subTensor)){
            if(subTensorPriorityQueue.size() >= k){
                SubTensor minSubTensor = subTensorPriorityQueue.peek().getKey();
                if(Util.computeDensity(minSubTensor) < Util.computeDensity(subTensor)){
                    subTensorPriorityQueue.poll();
                    subTensorPriorityQueue.offer(new Pair<>(subTensor,density1));
                }
            }else{
                subTensorPriorityQueue.offer(new Pair<>(subTensor,density1));
            }
        }*/

        //TODO：改度量时此处要改
        if(subTensorPriorityQueue.size() >= k){
            SubTensor minSubTensor = subTensorPriorityQueue.peek().getKey();
            if(Util.computeDensity(minSubTensor) < Util.computeDensity(subTensor)){
                subTensorPriorityQueue.poll();
                subTensorPriorityQueue.offer(new Pair<>(subTensor,density1));
            }
        }else{
            subTensorPriorityQueue.offer(new Pair<>(subTensor,density1));
        }
        //将子张量数组中所有子数组从全局张量tensor中去除
        updateTensorByRemoveSubTensor(tensor, subTensor);
    }

    /**
     * 判断coreDimension是否都为1，是则丢弃该子张量
     */
    static boolean isSingle(SubTensor subTensor){
        boolean flag = true;
        int[] size = subTensor.getSize();
        for(int i = 0; i < CORE_COLUMN_LENGTH; i++){
            if(size[i] > 1)flag = false;
        }
        return flag;
    }

    /**
     * 根据输入数据，将数据转换为张量进行存储
     * @param file_path：文件输入路径
     * @param splitStr：文件数据分隔符
     * @return 原始数据张量
     */
    public static Tensor createTensor(String file_path, String splitStr){
        //初始化需要创建的属性
        String line;
        String[] attributes;
        int[] size = Arrays.copyOfRange(globalSize, 0, globalSize.length - 1);
        int[][] attVals = new int[globalSize.length - 3][globalSize[globalSize.length - 1]];
        double[][] spaceVals = new double[2][globalSize[globalSize.length - 1]];
        int[] value = new int[globalSize[globalSize.length - 1]];
        int mass = 0;
        List<int[]> coreDimensionMass = new ArrayList<>();
        for(int i = 0; i < CORE_COLUMN_LENGTH; i++){
            coreDimensionMass.add(new int[globalSize[i] + 1]);
        }
        int timeSpaceMaxIndex = PARTITION_NUMBER * PARTITION_NUMBER * PARTITION_NUMBER + PARTITION_NUMBER * PARTITION_NUMBER + PARTITION_NUMBER;
        int[] timeSpaceMass = new int[timeSpaceMaxIndex + 1];
        int[] timeMass = new int[PARTITION_NUMBER + 1];
        int[] latitudeMass = new int[PARTITION_NUMBER + 1];
        int[] longitudeMass = new int[PARTITION_NUMBER + 1];
        //由于各个维度的属性都是从1开始，因此需要加1
        boolean[] timeExisted = new boolean[PARTITION_NUMBER + 1];
        boolean[] timeSpaceExisted = new boolean[timeSpaceMaxIndex + 1];
        List<Map<Integer, List<Integer>>> newSelectDataByCoreDimension = new ArrayList<>();
        for(int i = 0; i < CORE_COLUMN_LENGTH; i++){
            newSelectDataByCoreDimension.add(new HashMap<>());
        }
        Map<Integer, List<Integer>> selectDataByTime = new HashMap<>();
        Map<Integer, List<Integer>> selectDataByTimeSpace = new HashMap<>();
        int i = 0; //显示当前记录为第几条记录
        //modeId用于给属性值进行转换
        //初始化modeId,全部填充为1
        Integer[] modeId = new Integer[globalSize.length - 3];
        Arrays.fill(modeId, 1);
        //初始化timeExisted数组，对于0索引应为false, 防止后续遍历到0认为其存在
        Arrays.fill(timeExisted, 1, timeExisted.length - 1, true);
        for (int j = 0; j < globalSize.length - 3; j++) {
            convert.add(new HashMap<>());
            recoverDicts.add(new HashMap<>());
        }

        //开始读取数据
        BufferedReader reader;
        try {
            reader = new BufferedReader(new FileReader(file_path));
            while ((line = reader.readLine()) != null){
                attributes = line.split(splitStr);
                value[i] = Integer.parseInt(attributes[globalSize.length - 1]);
                for(int mode = 0; mode < globalSize.length - 4; mode++){
                    //转换数据格式并存储属性映射关系,导入数据到attVals、value,存储每一个维度的质量
                    Map<String ,Integer> currentMap = convert.get(mode);
                    Map<Integer, String> recoverMap = recoverDicts.get(mode);
                    String key = attributes[mode];
                    if(!currentMap.containsKey(key)){
                        currentMap.put(key, modeId[mode]);
                        recoverMap.put(modeId[mode]++, key);
                    }
                    attVals[mode][i] = currentMap.get(key);
                    //统计时间和核心维度中属性的质量信息
                    int[] masses = coreDimensionMass.get(mode);
                    masses[attVals[mode][i]] += value[i];
                    Map<Integer, List<Integer>> newSelectDataByCore = newSelectDataByCoreDimension.get(mode);
                    if(!newSelectDataByCore.containsKey(attVals[mode][i])){
                        newSelectDataByCore.put(attVals[mode][i], new ArrayList<>());
                    }
                    newSelectDataByCore.get(attVals[mode][i]).add(i);
                }

                //处理时间维度信息
                attVals[TIME_COLUMN][i] = Integer.parseInt(attributes[TIME_COLUMN]);
                timeMass[Util.transformTime(attVals[TIME_COLUMN][i])] += value[i];
                //记录selectDataByTime
                if(!selectDataByTime.containsKey(Util.transformTime(attVals[TIME_COLUMN][i]))){
                    selectDataByTime.put(Util.transformTime(attVals[TIME_COLUMN][i]), new ArrayList<>());
                }
                selectDataByTime.get(Util.transformTime(attVals[TIME_COLUMN][i])).add(i);

                //统计总质量
                mass += value[i];
                //处理经纬度信息，将其转换为h * h * h的时空块
                double latitude = Double.parseDouble(attributes[LATITUDE_COLUMN]);
                double longitude = Double.parseDouble(attributes[LONGITUDE_COLUMN]);
                spaceVals[0][i] = latitude;
                spaceVals[1][i] = longitude;
                //计算时空块索引地址
                int timeSpaceIndex = Util.computeTimeSpaceIndex(attVals[TIME_COLUMN][i], latitude, longitude);
                if(!selectDataByTimeSpace.containsKey(timeSpaceIndex)){
                    selectDataByTimeSpace.put(timeSpaceIndex, new ArrayList<>());
                    timeSpaceExisted[timeSpaceIndex] = true;
                }
                selectDataByTimeSpace.get(timeSpaceIndex).add(i);
                timeSpaceMass[timeSpaceIndex] += value[i];
                latitudeMass[Util.transformLatitude(latitude)] += value[i];
                longitudeMass[Util.transformLongitude(longitude)] += value[i];

                i++;
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return new Tensor(mass, size, attVals, spaceVals, value, coreDimensionMass, timeMass, timeSpaceMass , latitudeMass, longitudeMass, timeSpaceExisted ,newSelectDataByCoreDimension, selectDataByTime, selectDataByTimeSpace);
    }

}
