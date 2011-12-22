package params;

import java.util.ArrayList;

public class ObjectEnumParameter<T extends ParamObject> extends EnumParameter {

	private ArrayList<T> objectList = new ArrayList<T>();
	public ObjectEnumParameter(String key, String name) 
	{
		super(key, name);
	}

	public ObjectEnumParameter<T> registerObject(T obj)
	{
		super.registerEntry(obj.getParamDescription());
		objectList.add(obj);
		return this;
	}
	
	public T getObject()
	{
		int value = getValue();
		return objectList.get(value-minValue);
	}
	
	@Override
	public String getValueAsString()
	{
		return getObject().getParamDescription();
	}
}
