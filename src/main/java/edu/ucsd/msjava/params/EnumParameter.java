package edu.ucsd.msjava.params;

import java.util.ArrayList;

public class EnumParameter extends IntParameter {

    private int defaultValue = Integer.MIN_VALUE;
    private ArrayList<String> descriptions = new ArrayList<String>();

    public EnumParameter(String key) {
        super(key, null, null);
        super.minValue(0);
    }

    public EnumParameter(ParamManager.ParamNameEnum paramInfo) {
        super(paramInfo);
        super.minValue(0);
    }

    public EnumParameter setMinIndex(int minIndex) {
        super.minValue(minIndex);
        return this;
    }

    public EnumParameter registerEntry(String description) {
        descriptions.add(description);
        return this;
    }

    public EnumParameter setDefault() {
        this.defaultValue = getMinValue() + descriptions.size() - 1;
        super.defaultValue(defaultValue);
        return this;
    }

    protected int getCurIndex() {
        return getMinValue() + descriptions.size();
    }

    protected int getMinValue() {
        if (super.minValue == null)
            return 0;
        else
            return super.minValue;
    }

    @Override
    public String getName() {
        if (super.getName() != null)
            return super.getName();
        StringBuffer buf = new StringBuffer();
        for (int i = super.minValue; i < getMinValue() + descriptions.size(); i++) {
            if (i > getMinValue())
                buf.append("/");
            buf.append(i);
        }
        return buf.toString();
    }

    @Override
    public String getDescription() {
        StringBuffer buf = new StringBuffer();
        if (super.getDescription() != null) {
            buf.append(super.getDescription() + ", ");
            buf.append("Default: " + this.defaultValue);
            return buf.toString();
        }
        for (int i = super.minValue; i < getMinValue() + descriptions.size(); i++) {
            if (i > getMinValue())
                buf.append(", ");
            buf.append(i + ": " + descriptions.get(i - getMinValue()));
            if (i == defaultValue)
                buf.append(" (Default)");
        }
        return buf.toString();
    }

    @Override
    public String parse(String value) {
        super.maxValue(getMinValue() + descriptions.size());
        return super.parse(value);
    }
}
