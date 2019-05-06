package org.systemsbiology.jrap.stax;
import java.io.*;
import java.nio.*;
import java.nio.channels.*;
import java.util.*;

public class ByteBufferIterator implements Iterator
{

    private int INITIAL_BUFFERSIZE = 10000;
    private int bufferSize = INITIAL_BUFFERSIZE;
    private FileInputStream fis = null;
    private FileChannel fc = null;
    private String fPath;
    private long fSize;
    private ByteBuffer bb = null;
    private long totBytesRead = 0;

    public void setBufferSize (int b) {
	      this.bufferSize = b;
          //bb = ByteBuffer.allocate(bufferSize);
    }

    public int getBufferSize() {
	      return this.bufferSize;
    }

    public long getFileSize() {
       return this.fSize;
    }

    public String getPath() {
       return this.fPath;
    }

    public long getFilePos() {
       return totBytesRead;
    }

    public ByteBufferIterator(String fN) throws IOException {
	   fPath = fN;
       fis = new FileInputStream(fN);
       fc = fis.getChannel();
       fSize = fc.size();
	}

    public ByteBufferIterator(String fN, int buflen) throws Exception {
	   fPath = fN;
       fis = new FileInputStream(fN);
       fc = fis.getChannel();
       fSize = fc.size();
       bufferSize = buflen;
    }

    public boolean hasNext() {
	   return totBytesRead < fSize;
    }
/*
    public ByteBuffer next() {
       try {
          if(bb == null) bb = ByteBuffer.allocate(bufferSize);
          bb.rewind();
		  int bytesRead = fc.read(bb);
          if(bytesRead > 0){
             totBytesRead += bytesRead;
             bb.limit(bytesRead);
          } else {
  		     fis.close();
          }
          bb.rewind(); 
          //System.out.println("read "+bytesRead+" bytes, current total is "+totBytesRead+"; and filesize is "+fSize);
       } catch (Exception e) {
	      System.err.println("Problem in ByteBufferIterator.next(): "+e);
          e.printStackTrace();
          return null;
       }
       return bb;
    }
*/
    public ByteBuffer next() {
        try {
            //dhmay 20100223, fixing issue with small files in which you can't try to read the full buffer size
            //on the last scan
            long numBytesToRead = Math.min(bufferSize, fSize-totBytesRead);
            bb = fc.map(FileChannel.MapMode.READ_ONLY, totBytesRead, numBytesToRead);
            int bytesRead = bb.capacity();
            if(bytesRead > 0){
                totBytesRead += bytesRead;
            } else {
                fis.close();
            }
            bb.rewind();
            //System.out.println("read "+bytesRead+" bytes, current total is "+totBytesRead+"; and filesize is "+fSize);
        } catch (Exception e) {
            System.err.println("Problem in ByteBufferIterator.next(): "+e);
            e.printStackTrace();
            return null;
        }
        return bb;
    }

    public void remove() {}

    protected void finalize() throws Throwable {
        try {fis.close();} catch (Throwable t) {};
        super.finalize();
    }
 }
