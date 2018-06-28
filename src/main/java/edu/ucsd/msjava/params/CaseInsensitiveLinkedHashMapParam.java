package edu.ucsd.msjava.params;

import java.util.LinkedHashMap;

/**
 * Case insensitive LinkedHashMap (Key:String, Value:Parameter)
 * from https://stackoverflow.com/a/8237007/1179467
 */
public class CaseInsensitiveLinkedHashMapParam extends LinkedHashMap<String, Parameter> {

    @Override
    public Parameter put(String key, Parameter value) {
        return super.put(key.toLowerCase(), value);
    }

    // not @Override because that would require the key parameter to be of type Object
    public Parameter get(String key) {
        return super.get(key.toLowerCase());
    }

    public boolean containsKey(String key) {
        return super.containsKey(key.toLowerCase());
    }
}
