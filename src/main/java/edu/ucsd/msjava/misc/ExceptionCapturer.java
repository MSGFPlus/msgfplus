/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.ucsd.msjava.misc;

/**
 * For use with Runnable implementations and ThreadPoolExecutorWithExceptions,
 * to allow throwing checked exceptions and then seeing them in the 
 * ThreadPoolExecutorWithExceptions to trigger thread pool shutdown.
 * @author Bryson
 */
public interface ExceptionCapturer {
    boolean hasException();
    Throwable getException();
}
