package edu.ucsd.msjava.msdbsearch;

import java.io.*;
import java.util.*;
import java.util.Map.Entry;

import edu.ucsd.msjava.sequences.Constants;
import edu.ucsd.msjava.sequences.Sequence;



//import suffixarray.ByteSequence;


/**
 * An implementation of the Sequence class allowing a fasta file to be used as 
 * the database.
 * @author sangtae
 *
 */
public class CompactFastaSequence implements Sequence {

	public static final String SEQ_FILE_EXTENSION = ".cseq";
	public static final String ANNOTATION_FILE_EXTENSION = ".canno";
	
	//this is the file in which the sequence was generated
	private String baseFilepath;

	private TreeMap<Integer,String> annotations;

	// the contents of the sequence concatenated into a long string
	private byte[] sequence;

	// the number of characters in the buffer
	private int size;

	// the alphabet map
	private HashMap<Character, Byte> alpha2byte;

	// the reverse translation map
	private HashMap<Byte, Character> byte2alpha;

	// the string representation of the alphabet
	private String alphabetString;

	// the identifier for this sequence
	private int id;

	private boolean truncateAnnotation = false;	// if true, store annotations only before first blank
	
	/***** CONSTRUCTORS *****/
	/**
	 * Constructor. The alphabet will be created dynamically according from the 
	 * fasta file.
	 * @param filepath the path to the fasta file.
	 */
	public CompactFastaSequence(String filepath) {
		this(filepath, Constants.AMINO_ACIDS_20);
	}

	/**
	 * Constructor using the specified alphabet set. If there is a letter not in
	 * the alphabet, it will be encoded as the TERMINATOR byte.
	 * @param filepath the path to the fasta file.
	 * @param alphabet the specifications alphabet string. This could take the 
	 *        predefined AminoAcid strings defined in this class or customized strings.
	 * @param seqExtension the extension to use for the sequence file.
	 */
	private CompactFastaSequence(String filepath, String alphabet) {

		String[] tokens = filepath.split("\\.");
		String extension = tokens[tokens.length-1];
		String basepath = filepath.substring(0, filepath.length()-extension.length()-1);

		this.baseFilepath = basepath;
		if(!extension.equalsIgnoreCase("fasta") && !extension.equalsIgnoreCase("fa")) {
			System.err.println("Input error: not a fasta file");
			System.exit(-1);
		}

		String metaFile = basepath + ANNOTATION_FILE_EXTENSION;
		String sequenceFile = basepath + SEQ_FILE_EXTENSION;
		if(!new File(metaFile).exists() || !new File(sequenceFile).exists()) {      
			createObjectFromRawFile(filepath, alphabet);
		}

		int metaId = readMetaInfo();
		int seqId = readSequence();
		
		if(metaId == seqId && metaId != 0) {
			initializeAlphabet(this.alphabetString);
			this.id = metaId;
		}
		else {
			System.err.println("The sequence file and meta file have different ids: " + basepath);
			System.err.println("The problem can be solved by recreating the files");
			System.exit(-1);
		}
	}

	public CompactFastaSequence truncateAnnotation()
	{
		truncateAnnotation = true;
		return this;
	}
	
	/***** CLASS METHODS *****/  
	public Set<Byte> getAlphabetAsBytes() {
		return this.byte2alpha.keySet();
	}

	public Collection<Character> getAlphabet() {
		ArrayList<Character> results = new ArrayList<Character>();
		for (char c : this.byte2alpha.values()) 
			if (c!='_') results.add(c);
		return results;
	}

	public boolean isTerminator(long position) {
		return getByteAt(position)==Constants.TERMINATOR;
	}

	public char toChar(byte b) {
		if (byte2alpha.containsKey(b)) return byte2alpha.get(b);
		return '?';
	}

	public int getAlphabetSize() {
		return this.byte2alpha.size();
	}

	public long getSize() {
		return this.size;
	}

	public byte getByteAt(long position) {
		// forget boundary check for faster access
//		if(position >= this.size) return Constants.TERMINATOR;
//		return this.sequence.get((int)position);
		return this.sequence[(int)position];
	}

	public String getSubsequence(long start, long end) {
		if(start >= end || end > this.size)     return null;
		char[] seq = new char[(int)(end-start)];
		for(long i=start; i<end; i++) {
			seq[(int)(i-start)] = toChar(this.sequence[(int)i]);
		}
		return new String(seq);
	}

	public char getCharAt(long position) {
		return toChar(this.sequence[(int)position]);
	}

	public String toString(byte[] sequence) {
		String retVal = "";
		for(byte item : sequence) {
			Character c = byte2alpha.get(item);
			if(c != null)     retVal += c;
			else              retVal += '?';
		}
		return retVal;
	}

	public byte toByte(char c) {
		return alpha2byte.get(c);
	}

	public byte[] getBytes(int start, int end) {
		byte[] result = new byte[end-start];
		for (int i = start; i <end; i++) {
			result[i-start] = getByteAt(i);
		}
		return result;
	}

	public boolean isInAlphabet(char c) {
		return alpha2byte.containsKey(c);
	}

	public boolean isValid(long position) {
		if (isTerminator(position)) return false;
		if (isInAlphabet(getCharAt(position))) return true;
		return false;
	}

	public int getId() {
		return this.id;
	}

	public String getAnnotation(long position) {
		Entry<Integer, String> entry = annotations.higherEntry((int)position);
		if(entry != null)
			return entry.getValue();
		else
			return null;
	}

	public long getStartPosition(long position) {
		Integer startPos = annotations.floorKey((int)position);
		if (startPos==null) {
			return 0;
		}
		return startPos;
	}

	public String getMatchingEntry(long position) {
		Integer start = annotations.floorKey((int)position);	 // always "_" at start
		Integer end = annotations.higherKey((int)position);	   // exclusive
		if (start==null)    start = 0;
		if (end==null)      end = (int)this.getSize();
		while (!isValid(end-1)) end--;     // ensure that the last character is valid (exclusive)
		return this.getSubsequence(start+1, end);
	}

	public String getMatchingEntry(String name) {
		return null;
	}

	/**
	 * Setter method.
	 * @param baseFilepath set the baseFilepath for this object. The baseFilepath
	 *        has no extension.
	 */
	public void setBaseFilepath(String baseFilepath) { this.baseFilepath = baseFilepath; }

	/**
	 * Getter method.
	 * @return the baseFilename with properties described in the setter method.
	 */
	public String getBaseFilepath() { return this.baseFilepath; }
	
	/***** HELPER METHODS *****/
	// helper method, initialize the alphabet with given colon separated string
	private void initializeAlphabet(String s) {
		String[] tokens = s.split(":");
		this.alpha2byte = new HashMap<Character, Byte>();
		this.byte2alpha = new HashMap<Byte, Character>();
		this.byte2alpha.put(Constants.TERMINATOR, Constants.TERMINATOR_CHAR);
		this.byte2alpha.put(Constants.INVALID_CHAR_CODE, Constants.INVALID_CHAR);
		byte value = 2;
		for(byte i = 0; i < tokens.length; i++, value++) {
			for(int j = 0; j < tokens[i].length(); j++) {
				alpha2byte.put(tokens[i].charAt(j), value);
			}
			byte2alpha.put(value, tokens[i].charAt(0));
		}
	}

	// helper method to read and write the processed files given the alphabet
	private void createObjectFromRawFile(String filepath, String alphabet) {
		initializeAlphabet(alphabet);
		int size = 0;
		int id = UUID.randomUUID().hashCode();

		String seqFilepath = this.baseFilepath + SEQ_FILE_EXTENSION;
		String metaFilepath = this.baseFilepath + ANNOTATION_FILE_EXTENSION;
		
		// read the fasta file
		try {
			BufferedReader in = new BufferedReader(new FileReader(filepath));

			DataOutputStream seqOut = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(seqFilepath)));
			seqOut.writeInt(size);
			seqOut.writeInt(id);
			
			PrintStream metaOut = new PrintStream(new BufferedOutputStream(new FileOutputStream(metaFilepath)));
			metaOut.println(id);
			metaOut.println(alphabet);
			
			Integer offset = 0;
			String annotation = null;
			String s;
			
			// write protein sequences 
			while((s = in.readLine()) != null) {

				// this is a regular fasta line
				if(!s.startsWith(">")) {
					for(int index = 0; index < s.length(); index++) {
						Byte encoded = alpha2byte.get(s.charAt(index));
						if(encoded != null) {
							seqOut.writeByte(encoded);
						}
						else {
							seqOut.writeByte(Constants.INVALID_CHAR_CODE);
						}
					}
					offset += s.length();
				}

				// annotation line
				else {
					seqOut.writeByte(Constants.TERMINATOR);
					if(annotation != null)
						metaOut.println(offset+":"+annotation);
					// remember for the next annotation
					offset++;
					if(this.truncateAnnotation)
						annotation = s.substring(1).split("\\s+")[0];
					else
						annotation = s.substring(1);
				}
			}
			
			seqOut.writeByte(Constants.TERMINATOR);
			offset++;
			// the offset always points to the terminator of this sequence
			
			metaOut.println(offset+":"+annotation);
			size = offset;
			in.close();

			metaOut.flush();
			metaOut.close();
			
			seqOut.close();
			seqOut.close();
			
			// replace size
			RandomAccessFile raf = new RandomAccessFile(seqFilepath, "rw");
			raf.seek(0);
			raf.writeInt(size);
			raf.close();
		} 
		catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		} 
	}

	// read the metainformation file
	private int readMetaInfo() {
		String filepath = this.baseFilepath + ANNOTATION_FILE_EXTENSION;
		try {
			BufferedReader in = new BufferedReader(new FileReader(filepath));
			int id = Integer.parseInt(in.readLine());
			this.alphabetString = in.readLine().trim();
//			this.boundaries = new TreeSet<Long>();
//			for(String line = in.readLine(); line != null; line = in.readLine()) {
//				String[] tokens = line.split(":", 2);
//				this.boundaries.add(Long.parseLong(tokens[0]));
//			}
			this.annotations = new TreeMap<Integer, String>();
			for(String line = in.readLine(); line != null; line = in.readLine()) {
				String[] tokens = line.split(":", 2);
				this.annotations.put(Integer.parseInt(tokens[0]), tokens[1]);
			}
			in.close();
			return id;
		}
		catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return 0; 
	}

	// read the sequence in binary
	private int readSequence() {
		String filepath = this.baseFilepath + SEQ_FILE_EXTENSION;
		try {
			// read the first integer which encodes for the size of the file
			DataInputStream in = new DataInputStream(new BufferedInputStream(new FileInputStream(filepath)));
			int size = in.readInt();
			this.size = size;
			int id = in.readInt(); 

			sequence = new byte[size];
			in.read(sequence);

			in.close();
			return id;
		}
		catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return 0;
	}
	
	public int getNumProteins()
	{
		return annotations.keySet().size();
	}
	
	public float getRatioUniqueProteins()
	{
		int numProteins = 0;
		ArrayList<Integer> proteinLastIndexList = new ArrayList<Integer>(annotations.keySet());
		HashMap<Integer,ArrayList<Integer>> lengthProtIndexMap = new HashMap<Integer,ArrayList<Integer>>();
		int fromIndex = 0;
		for(int i=0; i<proteinLastIndexList.size(); i++)
		{
			int toIndex = proteinLastIndexList.get(i);
			int length = toIndex-fromIndex;
			ArrayList<Integer> list = lengthProtIndexMap.get(length);
			if(list == null)
			{
				list = new ArrayList<Integer>();
				lengthProtIndexMap.put(length, list);
			}
			list.add(i);
			fromIndex = toIndex;
		}
		
		int numUniqueProteins = 0;
		for(int length : lengthProtIndexMap.keySet())
		{
			ArrayList<Integer> protIndexList = lengthProtIndexMap.get(length);
			if(protIndexList.size() > 500)
				continue;
			numProteins += protIndexList.size();
			boolean[] isRedundant = new boolean[protIndexList.size()];
			for(int i=0; i<protIndexList.size(); i++)
			{
				if(isRedundant[i])
					continue;
				int toIndex1 = proteinLastIndexList.get(protIndexList.get(i));
				for(int j=i+1; j<protIndexList.size(); j++)
				{
					if(isRedundant[j])
						continue;
					int toIndex2 = proteinLastIndexList.get(protIndexList.get(j));
					boolean isIdentical = true;
					for(int l=0; l<length; l++)
					{
						if(sequence[toIndex1-1-l] != sequence[toIndex2-1-l])
						{
							isIdentical = false;
							break;
						}
					}
					if(isIdentical)
					{
						isRedundant[i] = isRedundant[j] = true;
//						System.out.println(annotations.get(toIndex1) + " = " + annotations.get(toIndex2));
						break;
					}
				}
				if(!isRedundant[i])
					numUniqueProteins++;
			}
		}
		return numUniqueProteins/(float)numProteins;
	}
	
}