package edu.ucsd.msjava.msutil;


/**
 * A mass object.
 * @author jung
 *
 */
public class Mass extends Matter {
  
  // holds the mass
  private float mass;
  
  // holds the nominal mass
  private int nominalMass;
  
  /**
   * Constructor.
   * @param mass the mass of this object.
   */
  public Mass(float mass) {
    this.mass = mass;
    this.nominalMass = Math.round(mass*Constants.INTEGER_MASS_SCALER);
  }
  
  public Mass(float mass, int nominalMass) {
	  this.mass = mass;
	  this.nominalMass = nominalMass;
  }

  /**
   * NominalMass setter
   * @param nominalMass
   */
  public void setNominalMass(int nominalMass)
  {
	  this.nominalMass = nominalMass;
  }
  
  /**
   * Gets the mass of this object. This is the mono isotopic mass.
   * @return
   */
  public float getMass()               { return mass; }
  
  /**
   * Gets the nominal mass of this object. 
   * @return nominal mass of this object.
   */
  public int getNominalMass()
  {
	  return nominalMass;
  }

	public boolean equals(Object obj) {
		if(!(obj instanceof Mass))
			return false;
		Mass m = (Mass)obj;
		return (this.compareTo(m)==0);
	}
}
