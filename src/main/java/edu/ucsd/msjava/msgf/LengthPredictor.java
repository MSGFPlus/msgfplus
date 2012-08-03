package edu.ucsd.msjava.msgf;

public class LengthPredictor {
	public static int	getLength(float parentMass)
	{
		int length = 0;
		
		if(parentMass < 882)
			length = 7;
		else if(parentMass  < 978)
			length = 8;
		else if(parentMass < 1080)
			length = 9;
		else if(parentMass < 1234)
			length = 10;
		else if(parentMass < 1447)
			length = 12;
		else if(parentMass < 1656)
			length = 14;
		else if(parentMass < 1870)
			length = 16;
		else if(parentMass < 2082)
			length = 18;
		else
			length = 20;
		return length;
	}
}
