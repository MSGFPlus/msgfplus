package edu.ucsd.msjava.params;

public class IntParameter extends NumberParameter<Integer> {

    public IntParameter(ParamManager.ParamNameEnum paramInfo) {
        super(paramInfo.getCommandlineName(), paramInfo.getName(), paramInfo.getDescription());
    }

    public IntParameter(String key, String name, String description) {
        super(key, name, description);
        super.minValue = 0;
        super.maxValue = Integer.MAX_VALUE;
    }

    @Override
    public String parse(String value) {
        try {
            super.value = Integer.valueOf(value);
            if (this.value == null) {
                return "Value cannot be null";
            }

            if (minValue == null && maxValue == null) {
                // Skip the range check
                return null;
            }

            if (minValue == null) {
                minValue = Integer.MIN_VALUE;
            }

            if (maxValue == null) {
                maxValue = Integer.MAX_VALUE;
            }

            String range = (super.isMinInclusive ? "[" : "(") + minValue + "," + maxValue + (super.isMaxInclusive ? "]" : ")");
            if (this.value < minValue || this.value > maxValue
                    || !super.isMinInclusive && this.value.equals(minValue)
                    || !super.isMaxInclusive && this.value.equals(maxValue)) {

                if (super.isMinInclusive && super.isMaxInclusive)
                    return "must be in the range " + minValue + " to " + maxValue;
                else if (super.isMinInclusive)
                    return "must be in the range " + minValue + " to " + (maxValue - 1);
                else if (!super.isMinInclusive && super.isMaxInclusive)
                    return "must be in the range " + (minValue + 1) + " to " + maxValue;
                else
                    return "must be in the range " + range;
            }
        } catch (NumberFormatException e) {
            return "must be an integer";
        }
        return null;
    }
}
