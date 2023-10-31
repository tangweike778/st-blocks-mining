import java.io.Serializable;
import java.util.*;

public class Slice implements Cloneable {

    private int index; //切片对应的属性值（转换之后的）

    //是使用数据表示所有切片是否存在还是使用vector存储存在的数据
    //两者的差异在于，当数据量比较大时前者过于浪费空间，后者在于删除时需要遍历。
    public Vector<Integer> dataIndexs; //包含数据的索引

    private boolean isExists; //该切片是否存在

    public Slice(int index, boolean isExists) {
        this.index = index;
        dataIndexs = new Vector<>();
        this.isExists = isExists;
    }


    public Vector<Integer> getDataIndexs() {
        return dataIndexs;
    }

    public boolean isExists() {
        return isExists;
    }

    public void setExists(boolean exists) {
        isExists = exists;
    }

    public int getIndex() {
        return index;
    }

    public void setDataIndexs(Vector<Integer> dataIndexs) {
        this.dataIndexs = dataIndexs;
    }

    public Slice clone() throws CloneNotSupportedException {
        Slice s = (Slice) super.clone();
        Vector<Integer> dataIndexsCopy = new Vector<>();
        for (Integer dataIndex : dataIndexs) {
            dataIndexsCopy.add(dataIndex);
        }
        s.setDataIndexs(dataIndexsCopy);
        return s;
    }

}
