package edu.ucsd.msjava.parser;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;

public class FullyBufferedLineReader implements LineReader {
	private int pointer;
	private byte[] buffer;

	private final byte CR = (byte)'\r';
	private final byte NL = (byte)'\n';
	int startIndex;
	
	public FullyBufferedLineReader(String fileName)
	{
		FileInputStream fin = null;
		try {
			fin = new FileInputStream(fileName);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		// load file into memory
		FileChannel in = fin.getChannel();
		ByteBuffer tempBuffer = null;
		try {
//			System.out.println(Integer.MAX_VALUE + "\t" + in.size() + "\t" + (int)in.size());
			tempBuffer = ByteBuffer.allocate((int)in.size());	// file size must be smaller than 2^32
			in.read(tempBuffer);
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		
		buffer = tempBuffer.array();
		pointer = 0;
		startIndex = 0;
	}
	
	public String readLine()	// line terminating char: \n or \r\n
	{
		if(pointer >= buffer.length)
			return null;
		while(pointer < buffer.length)
		{
			if(buffer[pointer] != NL)
				pointer++;
			else
			{
				String str;
				if(pointer > 0 && buffer[pointer-1] == CR)
					str = new String(buffer, startIndex, (pointer - startIndex - 1));
				else
					str = new String(buffer, startIndex, (pointer - startIndex));
				pointer++;
				startIndex = pointer;
				return str;
			}
		}
		String str = new String(buffer, startIndex, (pointer - startIndex));
		startIndex = pointer;
		return str;
	}
	
	public int getPosition()	{ return pointer; }
	public void seek(int position)	{ pointer = position; startIndex = pointer;}
	public void reset()		{ pointer = 0; startIndex = 0;}
	public int size()	{ return buffer.length; }
}
