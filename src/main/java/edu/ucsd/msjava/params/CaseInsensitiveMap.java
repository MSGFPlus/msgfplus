package edu.ucsd.msjava.params;

import java.util.HashMap;

/**
 * Case insensitive HashMap (Key:String, Value:String)
 * from https://stackoverflow.com/a/8237007/1179467
 */
public class CaseInsensitiveMap extends HashMap<String, String> {

    @Override
    public String put(String key, String value) {
        return super.put(key.toLowerCase(), value);
    }

    // not @Override because that would require the key parameter to be of type Object
    public String get(String key) {
        return super.get(key.toLowerCase());
    }

    public boolean containsKey(String key) {
        return super.containsKey(key.toLowerCase());
    }
}
