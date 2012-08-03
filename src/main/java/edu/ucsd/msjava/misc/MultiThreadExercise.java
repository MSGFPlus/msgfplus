package edu.ucsd.msjava.misc;

public class MultiThreadExercise implements Runnable {
	Thread t;
	int id;
	static int childNum = 0;
	MultiThreadExercise() {
		// Create a new, second thread
		t = new Thread(this, "Demo Thread");
		System.out.println("Child thread: " + t);
		t.start(); // Start the thread
		id = childNum++;
	}

	// This is the entry point for the second thread.
	public void run() {
		try {
			for(int i = 5; i > 0; i--) {
				System.out.println("Child " + id + " Thread: " + i);
				// Let the thread sleep for a while.
				Thread.sleep(500);
			}
		} catch (InterruptedException e) {
			System.out.println("Child " + id + " interrupted.");
		}
		System.out.println("Exiting child " + id + " thread.");
	}

	public static void main(String args[]) {
		System.out.println("NumProcessors: " + Runtime.getRuntime().availableProcessors());
		System.out.println("AvailableMem: " + Runtime.getRuntime().freeMemory());
		System.out.println("MaxMem: " + Runtime.getRuntime().maxMemory());
		
		new MultiThreadExercise(); // create a new thread
		new MultiThreadExercise(); // create a new thread
		try {
			for(int i = 5; i > 0; i--) {
				System.out.println("Main Thread: " + i);
				Thread.sleep(1000);
			}
		} catch (InterruptedException e) {
			System.out.println("Main thread interrupted.");
		}
		System.out.println("Main thread exiting.");
	}
}
