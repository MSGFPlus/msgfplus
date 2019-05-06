package org.systemsbiology.jrap.stax;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamReader;
import java.util.Iterator;
import java.util.regex.Pattern;
import java.util.regex.Matcher;
import java.io.StringReader;
import java.io.IOException;


/**
 * Created by IntelliJ IDEA.
 * User: tholzman
 * Date: Nov 16, 2009
 * Time: 3:30:12 PM
 * To change this template use File | Settings | File Templates.
 */

//Iterates through long strings within a file that begin with one
//pattern and end with another.  Useful for cutting out "<scan>...</scan>"
//etc.

public class EndPatternStringIterator implements Iterator {

    private static XMLInputFactory inputfactory = XMLInputFactory.newInstance();

    private LineIterator li;

    public LineIterator getLi() {
        return li;
    }

    public void setLi(LineIterator li) {
        this.li = li;
    }

    Pattern leftPat = null;
    Pattern rightPat = null;

    private String leftPatStr;
    public String getLeftPatStr() {
        return leftPatStr;
    }
    public void setLeftPatStr(String leftPatStr) {
        this.leftPatStr = leftPatStr;
        this.leftPat = Pattern.compile(leftPatStr);
    }

    private String rightPatStr;
    public String getRightPatStr() {
        return rightPatStr;
    }
    public void setRightPatStr(String rightPatStr) {
        this.rightPatStr = rightPatStr;
        this.rightPat = Pattern.compile(rightPatStr);
    }

    private long filePos = 0;
    public long getFilePos() {
       return filePos;
    }
    private int firstLineNo = 0;
    public int getFirstLineNo() {
       return this.firstLineNo;
    }

    public EndPatternStringIterator(String leftPat, String rightPat, LineIterator li) {
        setLi(li);
        setLeftPatStr(leftPat);
        setRightPatStr(rightPat);
    }

    public EndPatternStringIterator(String leftPat, String rightPat, String path) throws IOException {
        setLi(new LineIterator(new ByteBufferIterator(path)));
        setLeftPatStr(leftPat);
        setRightPatStr(rightPat);
    }


    StringBuilder curBuf = new StringBuilder();
    StringBuilder curLine = null;

    boolean noMore = false;
    public boolean hasNext() {
        return !noMore;
    }

    public StringBuilder next() {
        curBuf.setLength(0);
        //look for left pattern;
        //any part of a line left from last time?
        for(;;) {
           if(curLine == null || curLine.length() == 0) {
               if(li.hasNext()) {
                   curLine = li.next();
               } else {
                   noMore = true;
                   return curBuf;
               }
           }
           Matcher leftMatch = leftPat.matcher(curLine);
           if(!leftMatch.find()) {
               curLine.setLength(0);
               continue;
           }
           //set filepos and start concatenating
           int leftStartIndex = leftMatch.start();
           int leftEndIndex = leftMatch.end();
           filePos = li.getFilePos()+leftStartIndex;
           firstLineNo = li.getLineNum();
           curBuf.append(curLine.subSequence(leftStartIndex,leftEndIndex));
           //trim left moiety of curLine
           curLine.delete(0,leftEndIndex);
           break;
        }
        //find right pattern.  If there's an EOF, well, return the current working buffer.
        for(;;) {
           if(curLine.length() == 0) {
               curBuf.append(' ');
               if(li.hasNext()) {
                   curLine = li.next();
               } else {
                   noMore = true;
                   return curBuf;
               }
           }
           Matcher rightMatch = rightPat.matcher(curLine);
           if(!rightMatch.find())  {
               curBuf.append(curLine);
               curLine.setLength(0);
           } else {
               int rightEndIndex = rightMatch.end();
               curBuf.append(curLine.subSequence(0,rightEndIndex));
               curLine.delete(0,rightEndIndex);
               break;
           }
        }
        return curBuf;
    }

    public XMLStreamReader xmlsrNext() throws IOException {
        StringBuilder cursb = next();
        XMLStreamReader retVal = null;
        try {
            retVal = inputfactory.createXMLStreamReader(new StringBuilderReader(cursb));
        } catch (Exception e) {
            throw new IOException(e);
        }
        return retVal;
    }

    public XMLStreamReader xmlsrCur() throws Exception {
       //return inputfactory.createXMLStreamReader(new StringReader(new String(curBuf)));
       return inputfactory.createXMLStreamReader(new StringBuilderReader(curBuf));
    }

    public void remove() {}

    public static void main(String argv[]) {
       try {
          EndPatternStringIterator epsi =
            new EndPatternStringIterator(argv[0],argv[1],new LineIterator(new ByteBufferIterator(argv[2])));
          while(epsi.hasNext()) {
             StringBuilder sb = epsi.next();
             int lineno = epsi.getFirstLineNo();
             long filePos = epsi.getFilePos();
             System.out.println("Line: "+lineno+" filePos: "+ filePos+" "+sb);   
          }
       } catch (Exception e) {
           System.err.println(e);
           e.printStackTrace();
       }
    }
}
