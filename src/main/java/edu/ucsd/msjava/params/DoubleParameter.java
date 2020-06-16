package edu.ucsd.msjava.params;

public class DoubleParameter extends NumberParameter<Double> {

    public DoubleParameter(ParamManager.ParamNameEnum paramInfo) {
        super(paramInfo.getKey(), paramInfo.getName(), paramInfo.getDescription());
        setAdditionalDescription(paramInfo.getAdditionalDescription());
    }

    public DoubleParameter(String key, String name, String description) {
        super(key, name, description);
        super.minValue = Double.NEGATIVE_INFINITY;
        super.maxValue = Double.POSITIVE_INFINITY;
    }

    @Override
    public String parse(String value) {
        try {
            super.value = Double.valueOf(value);

            if (minValue == null)
                minValue = Double.NEGATIVE_INFINITY;

            if (maxValue == null)
                maxValue = Double.POSITIVE_INFINITY;

            String range = (super.isMinInclusive ? "[" : "(") + minValue + "," + maxValue + (super.isMaxInclusive ? "]" : ")");
            if (this.value < minValue ||
                    this.value > maxValue
                    || !super.isMinInclusive && this.value.equals(minValue)
                    || !super.isMaxInclusive && this.value.equals(maxValue))
                return "must be in the range " + range;
        } catch (NumberFormatException e) {
            return "must be a double";
        }
        return null;
    }
}
