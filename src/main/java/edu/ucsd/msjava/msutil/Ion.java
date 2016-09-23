package edu.ucsd.msjava.msutil;

public class Ion {
    public Ion(float mass, int charge) {
        this.mass = mass;
        this.charge = charge;
    }

    public float getMz() {
        return (mass + charge * (float) Composition.ChargeCarrierMass()) / charge;
    }

    public float getMass() {
        return mass;
    }

    public int getCharge() {
        return charge;
    }

    private float mass;
    private int charge;
}
