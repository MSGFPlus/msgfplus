package params;

import java.lang.reflect.ParameterizedType;
import java.util.HashMap;

public class ParamManager {
	private HashMap<String,Parameter> params;
	
	public ParamManager()
	{
		params = new HashMap<String,Parameter>();
	}

	public boolean addParameter(Parameter param)
	{
		if(params.containsKey(param.getKey()))
		{
			System.err.println("ParamManager: duplicate keys (" + param.getKey() + ")");
			return false;
		}
		params.put(param.getKey(), param);
		return true;
	}
	
	public Parameter getParameter(String key)
	{
		return params.get(key);
	}
	
	public <T extends Parameter> T getFileParameter(String key)
	{
//		Parameter param = params.get(key);
//		return (T)param;
//		Parameter instance = (Parameter) ((ParameterizedType)getClass().getGenericSuperclass()).getActualTypeArguments()[0];
//		if(param.getClass() != instance.getClass())
//			return (T)param;
		return null;
	}
}
