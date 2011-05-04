package misc;

public class MultiThreadExercise implements Runnable {
	Thread t;
	MultiThreadExercise() {
		// Create a new, second thread
		t = new Thread(this, "Demo Thread");
		System.out.println("Child thread: " + t);
		t.start(); // Start the thread
	}

	// This is the entry point for the second thread.
	public void run() {
		try {
			for(int i = 5; i > 0; i--) {
				System.out.println("Child Thread: " + i);
				// Let the thread sleep for a while.
				Thread.sleep(500);
			}
		} catch (InterruptedException e) {
			System.out.println("Child interrupted.");
		}
		System.out.println("Exiting child thread.");
	}

	public static void main(String args[]) {
		System.out.println("NumProcessors: " + Runtime.getRuntime().availableProcessors());
		System.out.println("AvailableMem: " + Runtime.getRuntime().freeMemory());
		System.out.println("MaxMem: " + Runtime.getRuntime().maxMemory());
		
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
