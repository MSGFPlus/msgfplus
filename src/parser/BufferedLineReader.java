package parser;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

public class BufferedLineReader extends BufferedReader implements LineReader {

	public BufferedLineReader(String fileName) throws FileNotFoundException 
	{
            super(new FileReader(fileName));
	}

	@Override
	public String readLine() {
            try {
                    return super.readLine();
            } catch (IOException e) {
                    e.printStackTrace();
            }
            return null;
	}
}
