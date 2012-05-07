package msutil;

//import java.util.ArrayList;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Hashtable;
import java.util.LinkedList;
import java.util.StringTokenizer;

public abstract class IonType {
    // IonType.InternalIon
    public static class InternalIon extends IonType {
        public InternalIon(String name, int charge, float offset) {
            super(name, charge, offset);
        }
        public InternalIon(int charge, float offset) {
        	super("I_"+charge+"_"+Math.round(offset), charge, offset);
        }
    }

    // added by kyowon
    public static class CyclicIon extends IonType {
        public CyclicIon(String name, int charge, float offset) {
            super(name, charge, offset);
        }
        public CyclicIon(int charge, float offset) {
        	super("C_"+charge+"_"+Math.round(offset), charge, offset);
        }
    }
    
 // added by kyowon
    public static class PrecursorIon extends IonType {
        public PrecursorIon(String name, int charge, float offset) {
            super(name, charge, offset);
        }
        public PrecursorIon(int charge, float offset) {
        	super("R_"+charge+"_"+Math.round(offset), charge, offset);
        }
    }
    
    // IonType.PrefixIon
    public static class PrefixIon extends IonType {
        public PrefixIon(String name, int charge, float offset) {
            super(name, charge, offset);
        }
        public PrefixIon(int charge, float offset) {
        	super("P_"+charge+"_"+Math.round(offset), charge, offset);
        }
    }

    // IonType.SuffixIon
    public static class SuffixIon extends IonType {
        public SuffixIon(String name, int charge, float offset) {
            super(name, charge, offset);
        }
        public SuffixIon(int charge, float offset) {
        	super("S_"+charge+"_"+Math.round(offset), charge, offset);
        }
    }

    
    public String toString()
    {
    	return name+"("+charge+","+offset+")";
    }
    
    public boolean equals(Object o){
    	if(this == o) return true;
    	if ( !(o instanceof IonType) ) return false;
    	IonType io = (IonType)o;
    	return io.name.equals(this.name) && io.charge == this.charge && io.offset == this.offset;
    }
    
    public int hashCode(){
    	return this.name.hashCode() * this.charge * new Float(this.offset).hashCode();
    }
    
    private String name;
    private int charge;
    private float offset;

    // kyowon added it
  
    
    protected IonType(String name, int charge, float offset) // Only to be used by child classes
    {
        this.name = name;
        this.charge = charge;
        this.offset = offset;
    }

    public int getCharge() { return charge; }
    public String getName() { return name; }
    public float getOffset() { return offset; }
    public boolean isPrefixIon()	{ return this instanceof PrefixIon; }
    public boolean isSuffixIon()	{ return this instanceof SuffixIon; }
    public float getMz(float mass)
    {
    	return mass/charge + offset;
    }

    public float getMass(float mz)
    {
    	return (mz - offset) * charge;
    }
    
    
    /**
     * Return ion type from string.
     * Ion name format: a/b/c  a=[sp], (s: suffixIon, p: prefixIon, i: internalIon, r: precursorIon),  b=charge,  c=offset  
     *  or
     * Ion name format: [abcxyz][+-]c,  (c=['H''H2''H2O''NH3'NH'] or c=offset)
     * Examples: y2-12.02   a+1.002-H2O  i/2/+1.23  s/1/-22.11  b-H2O-NH3
     * Returns null if format is not valid or ion does not exist
     * @param name
     * @return
     */
    public static IonType getIonType(String name)
    {
        if (name==null || name.length()==0) return null;
        // Ion name format: a/b/c  a=[spi]  b=charge  c=offset
        if (name.startsWith("s/") || name.startsWith("p/") || name.startsWith("i/") || name.startsWith("r/")) {
            StringTokenizer s = new StringTokenizer(name, "/", false);
            s.nextToken();
            if (!s.hasMoreTokens()) return null;
            String t=s.nextToken();
            try {
                int charge=Integer.parseInt(t.replace("+", ""));
                if (!s.hasMoreTokens()) return null;
                t=s.nextToken();
                float offset=Float.parseFloat(t);
                IonType it;
                if (name.startsWith("s"))
                    it=new IonType.SuffixIon(name, charge, offset);
                else if (name.startsWith("p"))
                    it= new IonType.PrefixIon(name, charge, offset);
                else if (name.startsWith("i"))
                    it= new IonType.InternalIon(name, charge, offset);
                else
                    it= new IonType.PrecursorIon(name, charge, offset);
                
                return it;
            } catch (NumberFormatException e) {
                return null;
            }
        }

        // Ion name format: [abcxyz][+-]c  c=['H''H2''H2O''NH3'NH'] or c=offset
        StringTokenizer s = new StringTokenizer(name, "+-", true);
        String token=s.nextToken();
        IonType base=ionTable.get(token);
        if (base==null) return null;
        float offset=0;
        // Add og subtract H2O, NH3, H, H2, ...
        while(s.hasMoreTokens()) {
            token = s.nextToken(); // + or -
            int sign;
            if (token.equals("+")) sign=1;
            else sign=-1;
            if (!s.hasMoreTokens()) throw new Error();
            token = s.nextToken();
            Float offs=compositionOffsetTable.get(token);
            if (offs==null) {
                try {
                    offs=Float.parseFloat(token);
                } catch (NumberFormatException e) {
                    return null;
                }
            }
            offset+=sign*offs;
        }
        IonType it;
        if (base instanceof PrefixIon)
            it=new PrefixIon(name, base.charge, base.offset+offset/base.charge);
        else if (base instanceof SuffixIon)
            it=new SuffixIon(name, base.charge, base.offset+offset/base.charge);
        else if (base instanceof InternalIon)
            it=new InternalIon(name, base.charge, base.offset+offset/base.charge);
        else it=null;
        return it;
    }

    public static ArrayList<IonType> getAllKnownIonTypes(int maxCharge, boolean removeRedundancy)
    {
    	return getAllKnownIonTypes(maxCharge, removeRedundancy, false);
    }
    
    public static ArrayList<IonType> getAllKnownIonTypes(int maxCharge, boolean removeRedundancy, boolean addPhosphoNL)
    {
    	return getAllKnownIonTypes(maxCharge, removeRedundancy, "H3PO4");
    }
    
    public static ArrayList<IonType> getAllKnownIonTypes(int maxCharge, boolean removeRedundancy, String nlString)
    {
        String[] base ={
                "x","y","z","a","b","c", //"x2","y2","z2","a2","b2","c2",
            };
        String[] extension={
          "", "-H2O", "-H2O-H2O", "-NH3", "-NH3-NH3", "-NH3-H2O","+n", "+n2", "-H"
        	};
        
        String[] nlExt;
        if(nlString != null)
        {
        	String[] token = nlString.split(",");
        	nlExt = new String[token.length];
        	for(int i=0; i<token.length; i++)
        	{
        		String nl = token[i].trim();
        		
        	}
        	nlExt = new String[] {"", "-H3PO4"};
        }
        else
        	nlExt = new String[] {""};
        
        ArrayList<IonType> ionList = new ArrayList<IonType>();
        for(int charge=1; charge<=maxCharge; charge++)
        {
            for (int i = 0; i < base.length; i++) {
                for (int j = 0; j < extension.length; j++) {
                	if(i==5 && j == 3)// c-NH3
                		continue;
               		for(int k=0; k<nlExt.length; k++)
                	{
                    	IonType ion = IonType.getIonType(base[i]+(charge > 1 ? charge : "")+extension[j]+nlExt[k]);
                    	assert(ion != null): base[i]+extension[j]+nlExt[k];
                        ionList.add(ion);
                	}
                }
            }
        }

        Collections.sort(ionList, new Comparator<IonType>() {
            public int compare(IonType i1, IonType i2) {
            	if(i1.getCharge() < i2.getCharge())
            		return -1;
            	else if (i1.getOffset()<i2.getOffset()) return -1;
                else if (i1.getOffset()>i2.getOffset()) return 1;
                return 0;
            }
        });
        
        if(!removeRedundancy)
        	return ionList;
        else
        {
            LinkedList<IonType> newIonList = new LinkedList<IonType>();
            for(int i=1; i<ionList.size(); i++)
            {
            	IonType prevIon = ionList.get(i-1);
            	IonType curIon = ionList.get(i);
            	if(curIon.getOffset()-prevIon.getOffset() < 0.1f && 
            			curIon.getCharge() == prevIon.getCharge() &&
            			curIon.isPrefixIon() == prevIon.isPrefixIon())
            	{
            		if(curIon.getName().length() < prevIon.getName().length())
            		{
            			newIonList.removeLast();
            			newIonList.add(curIon);
            		}
            	}
            	else
            	{
            		newIonList.add(curIon);
            	}
            }
            return new ArrayList<IonType>(newIonList);
        }
    }

    protected static Hashtable<String, IonType> ionTable;
    protected static Hashtable<String, Float> compositionOffsetTable;
    public final static IonType Y = new SuffixIon("y", 1, (float)Composition.OFFSET_Y);
    public final static IonType Z = new SuffixIon("z", 1, (float)(Y.offset-(Composition.NH2)));
    public final static IonType X = new SuffixIon("x", 1, (float)(Y.offset+Composition.CO));
    public final static IonType B = new PrefixIon("b", 1, (float)Composition.OFFSET_B);
    public final static IonType A = new PrefixIon("a", 1, (float)(B.offset-Composition.CO));
    public final static IonType C = new PrefixIon("c", 1, (float)(B.offset+Composition.NH3));
    public final static IonType NOISE = new PrefixIon("noise", 0, 0);
    
    // Composition (int C,   int H,      int N,      int O,      int S)
    // Mass         12.0f,   1.0078250f, 14.003074f, 15.994915f, 31.9720718f
    static {
        ionTable = new Hashtable<String, IonType>();
        ionTable.put("x", X); //+63.03697
        ionTable.put("y", Y); //+19.01839
        ionTable.put("z", Z); //+4.012321 => +3
        ionTable.put("a", A); //-27.00246
        ionTable.put("b", B); //+1.00794
        ionTable.put("c", C); //+16.0188
        
        for(int charge=2; charge<=4; charge++)
        {
            ionTable.put("x"+charge, new SuffixIon("x"+charge, charge, (float)((X.offset+Composition.PROTON*(charge-1))/charge)));
            ionTable.put("y"+charge, new SuffixIon("y"+charge, charge, (float)((Y.offset+Composition.PROTON*(charge-1))/charge)));
            ionTable.put("z"+charge, new SuffixIon("z"+charge, charge, (float)((Z.offset+Composition.PROTON*(charge-1))/charge)));
            ionTable.put("a"+charge, new PrefixIon("a"+charge, charge, (float)((A.offset+Composition.PROTON*(charge-1))/charge)));
            ionTable.put("b"+charge, new PrefixIon("b"+charge, charge, (float)((B.offset+Composition.PROTON*(charge-1))/charge)));
            ionTable.put("c"+charge, new PrefixIon("c"+charge, charge, (float)((C.offset+Composition.PROTON*(charge-1))/charge)));
        	
        }
//        ionTable.put("x2", new SuffixIon("x2", 2, (float)((X.offset+Composition.PROTON)/2)));
//        ionTable.put("y2", new SuffixIon("y2", 2, (float)((Y.offset+Composition.PROTON)/2)));
//        ionTable.put("z2", new SuffixIon("z2", 2, (float)((Z.offset+Composition.PROTON)/2)));
//        ionTable.put("a2", new PrefixIon("a2", 2, (float)((A.offset+Composition.PROTON)/2)));
//        ionTable.put("b2", new PrefixIon("b2", 2, (float)((B.offset+Composition.PROTON)/2)));
//        ionTable.put("c2", new PrefixIon("c2", 2, (float)((C.offset+Composition.PROTON)/2)));
//        // Internal ions "i_.."
//        ionTable.put("i_a",  new InternalIon("i_a",  1, (float)A.offset));
//        ionTable.put("i_a2", new InternalIon("i_a2", 2, (float)((A.offset+Composition.PROTON)/2)));
//        ionTable.put("i_b",  new InternalIon("i_b",  1, (float)(B.offset)));
//        ionTable.put("i_b2", new InternalIon("i_b2", 2, (float)((B.offset+Composition.PROTON)/2)));
//        ionTable.put("i_c",  new InternalIon("i_c",  1, (float)(C.offset)));
//        ionTable.put("i_c2", new InternalIon("i_c2", 2, (float)((C.offset+Composition.PROTON)/2)));
//        ionTable.put("i_x",  new InternalIon("i_x",  1, (float)(X.offset)));
//        ionTable.put("i_x2", new InternalIon("i_x2", 2, (float)((X.offset+Composition.PROTON)/2)));
//        ionTable.put("i_y",  new InternalIon("i_y",  1, (float)(Y.offset)));
//        ionTable.put("i_y2", new InternalIon("i_y2", 2, (float)((Y.offset+Composition.PROTON)/2)));
//        ionTable.put("i_z",  new InternalIon("i_z",  1, (float)(Z.offset)));
//        ionTable.put("i_z2", new InternalIon("i_z2", 2, (float)((Z.offset+Composition.PROTON)/2)));

        compositionOffsetTable = new Hashtable<String, Float>();
        compositionOffsetTable.put("H2O", (float)Composition.H2O);
        compositionOffsetTable.put("NH3", (float)Composition.NH3);
        compositionOffsetTable.put("NH", (float)Composition.NH);
        compositionOffsetTable.put("n", (float)Composition.ISOTOPE);
        compositionOffsetTable.put("n2", (float)Composition.ISOTOPE2);
        compositionOffsetTable.put("H", (float)Composition.H);
        compositionOffsetTable.put("H3PO4", (float)(Composition.H*3+Composition.P+Composition.O*4));
    }
    public static void main(String[] args) {
    	ArrayList<IonType> allIons = IonType.getAllKnownIonTypes(3, true, true);
    	for(IonType ion : allIons)
    		System.out.println(ion.getName()+"\t"+ion.getOffset());
    }
}
