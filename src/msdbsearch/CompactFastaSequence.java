package msdbsearch;

import java.io.*;
import java.nio.ByteBuffer;
import java.util.*;
import java.util.Map.Entry;


import sequences.Constants;
import sequences.Sequence;

//import suffixarray.ByteSequence;


/**
 * An implementation of the Sequence class allowing a fasta file to be used as 
 * the database.
 * @author sangtae
 *
 */
public class CompactFastaSequence implements Sequence {

	//this is the file in which the sequence was generated
	private String baseFilepath;

	// used for writing the encoded binary sequence.
	private final String seqExtension;

	// records the starting position of each protein
	private TreeSet<Long> boundaries;

	// maps the header strings of the fasta entries to the position of the terminators
//	private TreeMap<String,Integer> header2ends;

	// the contents of the sequence concatenated into a long string
	private ByteBuffer sequence;

	// the original serialized fasta file
//	private ByteBuffer original;

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

	/***** HELPER METHODS *****/
	// helper method, initialize the alphabet with given colon separated string
	private void initializeAlphabet(String s) {
		String[] tokens = s.split(":");
		this.alpha2byte = new HashMap<Character, Byte>();
		this.byte2alpha = new HashMap<Byte, Character>();
		this.byte2alpha.put(Constants.TERMINATOR, Constants.TERMINATOR_CHAR);
		for(byte i = 0, value = 1; i < tokens.length; i++, value++) {
			for(int j = 0; j < tokens[i].length(); j++) {
				alpha2byte.put(tokens[i].charAt(j), value);
			}
			byte2alpha.put(value, tokens[i].charAt(0));
		}
	}

	// the other helper method when the hashmap is not known before hand
	// Modified by Sangtae to reduce memory usage
	private void createObjectFromRawFile(String filepath) {  

		HashMap<Character, Byte> alpha2byte = new HashMap<Character, Byte>();
		String alphabet = "";
		byte alphabetSize = 1;
		int size = 0;
		int id = UUID.randomUUID().hashCode();

		String seqFilepath = this.baseFilepath + this.seqExtension;
		String metaFilepath = this.baseFilepath + this.seqExtension + "anno";
		
		// read the fasta file
		try {
			BufferedReader in = new BufferedReader(new FileReader(filepath));

			DataOutputStream seqOut = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(seqFilepath)));
			seqOut.writeInt(size);
			seqOut.writeInt(id);
			
			File tempMetaOutFile = File.createTempFile(this.baseFilepath, ".msgfdbmeta");
			tempMetaOutFile.deleteOnExit();
			
			File tempOriginalOutFile = File.createTempFile(this.baseFilepath, "msgfdbseq");
			tempOriginalOutFile.deleteOnExit();
			DataOutputStream originalOut = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(tempOriginalOutFile)));;
			
			PrintStream tempMetaOut = new PrintStream(new BufferedOutputStream(new FileOutputStream(tempMetaOutFile)));
			
			Integer offset = 0;
			String annotation = null;
			String s;              // 
			while((s = in.readLine()) != null) {

				// this is a regular fasta line
				if(!s.startsWith(">")) {
					for(int index = 0; index < s.length(); index++) {
						Byte encoded = alpha2byte.get(s.charAt(index));
						if(encoded != null) {
							seqOut.writeByte(encoded);
						}
						else {
							seqOut.writeByte(alphabetSize);
							alpha2byte.put(s.charAt(index), alphabetSize++);
							alphabet += ":" + s.charAt(index);
						}
						originalOut.writeChar(s.charAt(index));
					}
					offset += s.length();
				}

				// annotation line
				else {
					seqOut.writeByte(Constants.TERMINATOR);
					originalOut.writeChar(Constants.TERMINATOR_CHAR);
					if(annotation != null)
						tempMetaOut.println(offset+":"+annotation);
					// remember for the next annotation
					offset++;
					annotation = s.substring(1);
				}
			}
			seqOut.writeByte(Constants.TERMINATOR);
			originalOut.writeChar(Constants.TERMINATOR_CHAR);
			offset++;
			// the offset always points to the terminator of this sequence
			tempMetaOut.println(offset+":"+annotation);
			size = offset;
			in.close();

			originalOut.flush();
			originalOut.close();
			tempMetaOut.flush();
			tempMetaOut.close();
			
			// append originalOut to seqOut
			DataInputStream originalIn = new DataInputStream(new BufferedInputStream(new FileInputStream(tempOriginalOutFile)));
			for(int i=0; i<size; i++)
				seqOut.writeByte(originalIn.readChar());
			originalIn.close();
			seqOut.close();
			tempOriginalOutFile.delete();
			
			// replace size
			RandomAccessFile raf = new RandomAccessFile(seqFilepath, "rw");
			raf.seek(0);
			raf.writeInt(size);
			raf.close();
			
			// write meta info file
			PrintStream metaOut = new PrintStream(new BufferedOutputStream(new FileOutputStream(metaFilepath)));
			metaOut.println(size);
			metaOut.println(id);
			metaOut.println(alphabet.substring(1));
			
			BufferedReader metaIn = new BufferedReader(new FileReader(tempMetaOutFile));
			while((s=metaIn.readLine()) != null)
				metaOut.println(s);
			metaIn.close();
			tempMetaOutFile.delete();
			metaOut.close();
		} 
		catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		} 
	}

	// helper method to read and write the processed files given the alphabet
	private void createObjectFromRawFile(String filepath, String alphabet) {
		HashMap<Character, Byte> alpha2byte = new HashMap<Character, Byte>();
		initializeAlphabet(alphabet);
		int size = 0;
		int id = UUID.randomUUID().hashCode();

		String seqFilepath = this.baseFilepath + this.seqExtension;
		String metaFilepath = this.baseFilepath + this.seqExtension + "anno";
		
		// read the fasta file
		try {
			BufferedReader in = new BufferedReader(new FileReader(filepath));

			DataOutputStream seqOut = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(seqFilepath)));
			seqOut.writeInt(size);
			seqOut.writeInt(id);
			
			File tempMetaOutFile = File.createTempFile(this.baseFilepath, ".msgfdbmeta");
			tempMetaOutFile.deleteOnExit();
			
			File tempOriginalOutFile = File.createTempFile(this.baseFilepath, "msgfdbseq");
			tempOriginalOutFile.deleteOnExit();
			DataOutputStream originalOut = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(tempOriginalOutFile)));;
			
			PrintStream tempMetaOut = new PrintStream(new BufferedOutputStream(new FileOutputStream(tempMetaOutFile)));
			
			Integer offset = 0;
			String annotation = null;
			String s;              // 
			while((s = in.readLine()) != null) {

				// this is a regular fasta line
				if(!s.startsWith(">")) {
					for(int index = 0; index < s.length(); index++) {
						Byte encoded = alpha2byte.get(s.charAt(index));
						if(encoded != null) {
							seqOut.writeByte(encoded);
						}
						else {
							seqOut.writeByte(Constants.TERMINATOR);
						}
						originalOut.writeChar(s.charAt(index));
					}
					offset += s.length();
				}

				// annotation line
				else {
					seqOut.writeByte(Constants.TERMINATOR);
					originalOut.writeChar(Constants.TERMINATOR_CHAR);
					if(annotation != null)
						tempMetaOut.println(offset+":"+annotation);
					// remember for the next annotation
					offset++;
					annotation = s.substring(1);
				}
			}
			seqOut.writeByte(Constants.TERMINATOR);
			originalOut.writeChar(Constants.TERMINATOR_CHAR);
			offset++;
			// the offset always points to the terminator of this sequence
			tempMetaOut.println(offset+":"+annotation);
			size = offset;
			in.close();

			originalOut.flush();
			originalOut.close();
			tempMetaOut.flush();
			tempMetaOut.close();
			
			// append originalOut to seqOut
			DataInputStream originalIn = new DataInputStream(new BufferedInputStream(new FileInputStream(tempOriginalOutFile)));
			for(int i=0; i<size; i++)
				seqOut.writeByte(originalIn.readChar());
			originalIn.close();
			seqOut.close();
			tempOriginalOutFile.delete();
			
			// replace size
			RandomAccessFile raf = new RandomAccessFile(seqFilepath, "rw");
			raf.seek(0);
			raf.writeInt(size);
			raf.close();
			
			// write meta info file
			PrintStream metaOut = new PrintStream(new BufferedOutputStream(new FileOutputStream(metaFilepath)));
			metaOut.println(size);
			metaOut.println(id);
			metaOut.println(alphabet);
			
			BufferedReader metaIn = new BufferedReader(new FileReader(tempMetaOutFile));
			while((s=metaIn.readLine()) != null)
				metaOut.println(s);
			metaIn.close();
			tempMetaOutFile.delete();
			metaOut.close();
		} 
		catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		} 
	}

	// read the metainformation file
	private int readMetaInfo() {
		String filepath = this.baseFilepath + this.seqExtension + "anno";
		try {
			BufferedReader in = new BufferedReader(new FileReader(filepath));
			this.size = Integer.parseInt(in.readLine());     
			int id = Integer.parseInt(in.readLine());
			this.alphabetString = in.readLine().trim();
			this.boundaries = new TreeSet<Long>();
			for(String line = in.readLine(); line != null; line = in.readLine()) {
				String[] tokens = line.split(":", 2);
				this.boundaries.add(Long.parseLong(tokens[0]));
			}
//			this.annotations = new TreeMap<Integer, String>();
//			for(String line = in.readLine(); line != null; line = in.readLine()) {
//				String[] tokens = line.split(":", 2);
//				this.annotations.put(Integer.parseInt(tokens[0]), tokens[1]);
//			}
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
		String filepath = this.baseFilepath + this.seqExtension;
		try {
			// read the first integer which encodes for the size of the file
			DataInputStream in = new DataInputStream(new BufferedInputStream(new FileInputStream(filepath)));
			int size = in.readInt();
			int id = in.readInt(); 

			// Modified by Sangtae
//			System.out.println("MaxMem: " + Runtime.getRuntime().maxMemory());
//			System.out.println("FreeMem1: " + Runtime.getRuntime().freeMemory());
//			System.out.println("TotalMem1: " + Runtime.getRuntime().totalMemory());
			byte[] sequenceArr = new byte[size];
//			System.out.println("FreeMem2: " + Runtime.getRuntime().freeMemory());
//			System.out.println("TotalMem2: " + Runtime.getRuntime().totalMemory());
			in.read(sequenceArr);
//			System.out.println("FreeMem3: " + Runtime.getRuntime().freeMemory());
//			System.out.println("TotalMem3: " + Runtime.getRuntime().totalMemory());

			sequence = ByteBuffer.wrap(sequenceArr).asReadOnlyBuffer();
//			System.out.println("FreeMem4: " + Runtime.getRuntime().freeMemory());
//			System.out.println("TotalMem4: " + Runtime.getRuntime().totalMemory());
			
//			this.original = fc.map(FileChannel.MapMode.READ_ONLY, 2*Integer.SIZE/Byte.SIZE + size, size);
//			byte[] originalArr = new byte[size];
//			System.out.println("FreeMem5: " + Runtime.getRuntime().freeMemory());
//			System.out.println("TotalMem5: " + Runtime.getRuntime().totalMemory());
//			in.read(originalArr);
//			System.out.println("FreeMem6: " + Runtime.getRuntime().freeMemory());
//			System.out.println("TotalMem6: " + Runtime.getRuntime().totalMemory());
//			original = ByteBuffer.wrap(originalArr).asReadOnlyBuffer();
//			System.out.println("FreeMem7: " + Runtime.getRuntime().freeMemory());
//			System.out.println("TotalMem7: " + Runtime.getRuntime().totalMemory());

			//this.original = new char[size];
			//ByteBuffer originalChars = fc.map(FileChannel.MapMode.READ_ONLY, 2*Integer.SIZE/Byte.SIZE + size, size);
			//for (int index=0; index<size; index++) this.original[index] = (char)originalChars.get(index);
			in.close();
			return id;
		}
		catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return 0;
	}



	/***** CONSTRUCTORS *****/
	/**
	 * Constructor. The alphabet will be created dynamically according from the 
	 * fasta file.
	 * @param filepath the path to the fasta file.
	 */
	public CompactFastaSequence(String filepath) {
		this(filepath, null);
	}

	/**
	 * Constructor using the specified alphabet set. If there is a letter not in
	 * the alphabet, it will be encoded as the TERMINATOR byte.
	 * @param filepath the path to the fasta file.
	 * @param alphabet the specifications alphabet string. This could take the 
	 *        predefined AminoAcid strings defined in this class or customized strings.
	 */
	public CompactFastaSequence(String filepath, String alphabet) {
		this(filepath, alphabet,  Constants.FILE_EXTENSION);
	}

	/**
	 * Constructor using the specified alphabet set. If there is a letter not in
	 * the alphabet, it will be encoded as the TERMINATOR byte.
	 * @param filepath the path to the fasta file.
	 * @param alphabet the specifications alphabet string. This could take the 
	 *        predefined AminoAcid strings defined in this class or customized strings.
	 * @param seqExtension the extension to use for the sequence file.
	 */
	public CompactFastaSequence(String filepath, String alphabet, String seqExtension) {

		this.seqExtension = seqExtension;

		String[] tokens = filepath.split("\\.");
		String extension = tokens[tokens.length-1];
		String basepath = filepath.substring(0, filepath.length()-extension.length()-1);

		this.baseFilepath = basepath;
		if(!extension.equalsIgnoreCase("fasta") && !extension.equalsIgnoreCase("fa")) {
			System.err.println("Input error: not a fasta file");
			System.exit(-1);
		}

		String metaFile = basepath + this.seqExtension + Constants.ANNO_FILE_SUFFIX;
		String sequenceFile = basepath + seqExtension;
		if(!new File(metaFile).exists() || !new File(sequenceFile).exists()) {      
			if(alphabet != null)  createObjectFromRawFile(filepath, alphabet);
			else                  createObjectFromRawFile(filepath);

		}

		int metaId = readMetaInfo();
		int seqId = readSequence();
		
		if(metaId == seqId) {
			initializeAlphabet(this.alphabetString);
			this.id = metaId;
		}
		else {
			System.err.println("The files " + metaFile + " and " + sequenceFile + " have different ids.");
			System.err.println("The problem can be solved by recreating the files");
			System.exit(-1);
		}
	}

	/***** CLASS METHODS *****/  
	@Override
	public Set<Byte> getAlphabetAsBytes() {
		return this.byte2alpha.keySet();
	}

	@Override
	public Collection<Character> getAlphabet() {
		ArrayList<Character> results = new ArrayList<Character>();
		for (char c : this.byte2alpha.values()) 
			if (c!='_') results.add(c);
		return results;
	}

	@Override
	public boolean isTerminator(long position) {
		return getByteAt(position)==Constants.TERMINATOR;
	}

	@Override
	public char toChar(byte b) {
		if (byte2alpha.containsKey(b)) return byte2alpha.get(b);
		return '?';
	}

	@Override
	public int getAlphabetSize() {
		return this.byte2alpha.size();
	}

	@Override
	public long getSize() {
		return this.size;
	}

	@Override
	public byte getByteAt(long position) {
		// forget boundary check for faster access
		if(position >= this.size) return Constants.TERMINATOR;
		return this.sequence.get((int)position);
	}

	@Override
	public String getSubsequence(long start, long end) {
		if(start >= end || end > this.size)     return null;
		char[] seq = new char[(int)(end-start)];
		for(long i=start; i<end; i++) {
			seq[(int)(i-start)] = toChar(this.sequence.get((int)i));
//			seq[(int)(i-start)] = (char)this.original.get((int)i);
		}
		return new String(seq);
	}

	@Override
	public char getCharAt(long position) {
		//return toChar(getByteAt(position));
		//return this.original[(int)position];
		return toChar(this.sequence.get((int)position));
	}

	@Override
	public String toString(byte[] sequence) {
		String retVal = "";
		for(byte item : sequence) {
			Character c = byte2alpha.get(item);
			if(c != null)     retVal += c;
			else              retVal += '?';
		}
		return retVal;
	}

	@Override
	public byte toByte(char c) {
		return alpha2byte.get(c);
	}

	@Override
	public byte[] getBytes(int start, int end) {
		byte[] result = new byte[end-start];
		for (int i = start; i <end; i++) {
			result[i-start] = getByteAt(i);
		}
		return result;
	}

	@Override
	public boolean isInAlphabet(char c) {
		return alpha2byte.containsKey(c);
	}

	@Override
	public boolean isValid(long position) {
		if (isTerminator(position)) return false;
		if (isInAlphabet(getCharAt(position))) return true;
		return false;
	}

	@Override
	public int getId() {
		return this.id;
	}

	@Override
	public String getAnnotation(long position) {
//		Entry<Integer, String> entry = annotations.higherEntry((int)position);
//		if(entry != null)
//			return entry.getValue();
//		else
//			return null;
		return null;
	}

	@Override
	public long getStartPosition(long position) {
		Long startPos = boundaries.floor(position);
		if(startPos==null)
			return 0;
		return startPos;
	}

	@Override
	public String getMatchingEntry(long position) {
//		Integer start = annotations.floorKey((int)position);	 // always "_" at start
//		Integer end = annotations.higherKey((int)position);	   // exclusive
//		if (start==null)    start = 0;
//		if (end==null)      end = (int)this.getSize();
//		while (!isValid(end-1)) end--;     // ensure that the last character is valid (exclusive)
//		return this.getSubsequence(start+1, end);
		Long start = boundaries.floor(position);
		Long end = boundaries.higher(position);
		if(start == null)	start = 0L;
		if(end == null)	end = this.getSize();
		while (!isValid(end-1)) end--;     // ensure that the last character is valid (exclusive)
		return this.getSubsequence(start+1, end);
	}

	@Override
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
}