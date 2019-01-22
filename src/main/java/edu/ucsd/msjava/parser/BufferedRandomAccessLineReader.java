package edu.ucsd.msjava.parser;

import net.pempek.unicode.UnicodeBOMInputStream;
import javafx.util.Pair;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;

public class BufferedRandomAccessLineReader implements LineReader {
    private static final int DEFAULT_BUFFER_SIZE = 1 << 16;
    private long pointer;
    private byte[] buffer;
    int bufPointer;
    private int bufLength = -1;
    long bufStartingPos;

    private final byte CR = (byte) '\r';
    private final byte NL = (byte) '\n';
    private final FileChannel in;
    private long fileSize;
    int startIndex;
    int bufSize;
    int bomLength;

    public BufferedRandomAccessLineReader(String fileName) {
        this(fileName, DEFAULT_BUFFER_SIZE);
    }

    public BufferedRandomAccessLineReader(String fileName, int bufSize) {
        FileInputStream fin = null;
        try {
            fin = new FileInputStream(fileName);
        } catch (FileNotFoundException e1) {
            e1.printStackTrace();
        }

        in = fin.getChannel();
        try {
            fileSize = in.size();
        } catch (IOException e) {
            e.printStackTrace();
        }

        this.bufSize = bufSize;
        pointer = 0;
        fillBuffer();
    }

    /**
     * Compare the bytes in buf to the bytes associated with the given Byte Order Mark class
     * @param buf
     * @param bomType
     * @return True if the bytes match, otherwise false
     */
    private static boolean bytesMatchBOM(byte[] buf, UnicodeBOMInputStream.BOM bomType) {
        byte[] bomBytes = bomType.getBytes();
        int matchCount = 0;

        if (buf.length < bomBytes.length)
            return false;

        for (int i = 0; i < bomBytes.length; i++) {
            if (buf[i] == bomBytes[i])
                matchCount++;
            else
                break;
        }

        return (matchCount == bomBytes.length);
    }

    /**
     * Check for a byte order mark at the start of str
     * Returns the updated string with the byte order mark, if present
     * @param str
     * @return
     */
    public static String stripBOM(String str) {
        Pair<String, Integer> result = stripBOMAndGetLength(str);
        return result.getKey();
    }

    /**
     * Check for a byte order mark at the start of str
     * If found, remove it
     * @param str
     * @return Key/value pair where the key is the updated string and the value is the byte order mark length
     */
    public static Pair<String, Integer> stripBOMAndGetLength(String str) {
        // Check for byte order marks
        byte[] buf = str.getBytes();
        int copyOffset = 0;

        if (buf.length >= 4) {
            if (bytesMatchBOM(buf, net.pempek.unicode.UnicodeBOMInputStream.BOM.UTF_32_LE)) {
                copyOffset = 4;
            } else if (bytesMatchBOM(buf, net.pempek.unicode.UnicodeBOMInputStream.BOM.UTF_32_BE)) {
                copyOffset = 4;
            }
        }

        if (copyOffset == 0 && buf.length >= 3) {
            if (bytesMatchBOM(buf, net.pempek.unicode.UnicodeBOMInputStream.BOM.UTF_8)) {
                copyOffset = 3;
            }
        }

        if (copyOffset == 0 && buf.length >= 2) {
            if (bytesMatchBOM(buf, net.pempek.unicode.UnicodeBOMInputStream.BOM.UTF_16_LE)) {
                copyOffset = 2;
            } else if (bytesMatchBOM(buf, net.pempek.unicode.UnicodeBOMInputStream.BOM.UTF_16_BE)) {
                copyOffset = 2;
            }
        }

        if (copyOffset > 0) {
            str = new String(java.util.Arrays.copyOfRange(buf, copyOffset, buf.length));
        }

        return new Pair<>(str, copyOffset);
    }

    private int fillBuffer() {
        ByteBuffer tempBuffer = null;
        int bytesRead = -1;
        try {
            tempBuffer = ByteBuffer.allocate(bufSize);
            bytesRead = in.read(tempBuffer);
        } catch (IOException e1) {
            if (!Thread.currentThread().isInterrupted()) {
                e1.printStackTrace();
            }
        }

        buffer = tempBuffer.array();
        bufLength = bytesRead;
        startIndex = 0;
        bufPointer = 0;
        bufStartingPos = pointer;

        return bytesRead;
    }

    public String readLine() {
        if (pointer >= fileSize)
            return null;

        Boolean startOfFile = (pointer == 0);

        String str = readLineFromBuffer();

        if (startOfFile) {
            // Check for a byte order mark
            Pair<String, Integer> result = stripBOMAndGetLength(str);

            bomLength = result.getValue();
            if (bomLength > 0) {
                str = result.getKey();
            }
        }

        if (bufPointer == bufLength && bufLength == bufSize) {
            fillBuffer();
            str = str + readLine();
        } else if (pointer < fileSize) {
            bufPointer++;
            pointer++;
            startIndex = bufPointer;
        }
        return str;
    }

    private String readLineFromBuffer()    // line terminating char: \n or \r\n
    {
        if (pointer >= fileSize)
            return null;
        while (pointer < fileSize && bufPointer < bufLength) {
            if (buffer[bufPointer] != NL) {
                bufPointer++;
                pointer++;
            } else
                break;
        }

        String str;
        try {
            if (bufPointer > 0 && buffer[bufPointer - 1] == CR)
                str = new String(buffer, startIndex, (bufPointer - startIndex - 1));
            else
                str = new String(buffer, startIndex, (bufPointer - startIndex));

            return str;

        } catch (java.lang.ArrayIndexOutOfBoundsException e) {
            System.out.println("bufPointer " + bufPointer + " is larger than the buffer array, length " + buffer.length);
            throw e;
        }


    }

    /**
     * Byte order mark length: non-zero for Unicode files with a byte order mark
     * See https://en.wikipedia.org/wiki/Byte_order_mark
     * @return
     */
    public int getBOMLength() {
        return bomLength;
    }

    public long getPosition() {
        return pointer;
    }

    public void seek(long position) {
        pointer = position;
        if (position >= bufStartingPos && position < bufStartingPos + bufSize) {
            startIndex = bufPointer = (int) (position - bufStartingPos);
        } else {
            try {
                in.position(pointer);
            } catch (IOException e) {
                if (!Thread.currentThread().isInterrupted()) {
                    e.printStackTrace();
                }
            }
            fillBuffer();
        }
    }

    public void reset() {
        pointer = 0;
        startIndex = 0;
    }

    public int size() {
        return buffer.length;
    }

    public void close() throws IOException {
        in.close();
    }

    public static void main(String argv[]) throws Exception {
        long time = System.currentTimeMillis();
        String fileName = "/home/sangtaekim/Research/Data/ABRF/2011/UniProt.Yeast.NFISnr.contamsS48.fasta";
        BufferedRandomAccessLineReader in = new BufferedRandomAccessLineReader(fileName, 1 << 16);
//		BufferedReader in = new BufferedReader(new FileReader(fileName));
//		RandomAccessFile in = new RandomAccessFile(fileName, "r");
        String s;
        int lineNum = 0;
        long pos = 0;
        while ((s = in.readLine()) != null) {
            lineNum++;
            if (lineNum == 48232)
                System.out.println(lineNum + " " + s + " " + (pos = in.getPosition()));
        }
        in.seek(pos);
        System.out.println(in.readLine());
        System.out.println("Time: " + (System.currentTimeMillis() - time));
    }

}
