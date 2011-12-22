package params;

import java.util.ArrayList;

public class EnumParameter extends IntParameter {

	private int defaultValue = Integer.MIN_VALUE;
	private ArrayList<String> descriptions = new ArrayList<String>();
	public EnumParameter(String key) {
		super(key, null, null);
		super.minValue(0);
	}

	public EnumParameter(String key, String name) {
		super(key, name, null);
		super.minValue(0);
	}

	public EnumParameter(String key, String name, String description) {
		super(key, name, description);
		super.minValue(0);
	}
	
	public EnumParameter setMinIndex(int minIndex)
	{
		super.minValue(minIndex);
		return this;
	}
	
	public EnumParameter registerEntry(String description)
	{
		descriptions.add(description);
		return this;
	}
	
	public EnumParameter setDefault()
	{
		this.defaultValue = super.minValue+descriptions.size()-1; 
		super.defaultValue(defaultValue);
		return this;
	}

	protected int getCurIndex()
	{
		return super.minValue + descriptions.size();
	}
	
	@Override
	public String getName()
	{
		if(super.getName() != null)
			return super.getName();
		StringBuffer buf = new StringBuffer();
		for(int i=super.minValue; i<super.minValue+descriptions.size(); i++)
		{
			if(i > super.minValue)
				buf.append("/");
			buf.append(i);
		}
		return buf.toString();
	}
	
	@Override 
	public String getDescription()
	{
		StringBuffer buf = new StringBuffer();
		if(super.getDescription() != null)
		{
			buf.append(super.getDescription()+", ");
			buf.append("Default: " + this.defaultValue);
			return buf.toString();
		}
		for(int i=super.minValue; i<super.minValue+descriptions.size(); i++)
		{
			if(i > super.minValue)
				buf.append(", ");
			buf.append(i+": " + descriptions.get(i-super.minValue));
			if(i == defaultValue)
				buf.append(" (Default)");
		}
		return buf.toString();
	}
	
	@Override
	public String parse(String value) {
		super.maxValue(super.minValue+descriptions.size());
		return super.parse(value);
	}	
}
