package edu.ucsd.msjava.msgf;

import edu.ucsd.msjava.msutil.Peptide;
import edu.ucsd.msjava.msutil.Spectrum;

public class MSGFResult {
	public MSGFResult(Spectrum spec, Peptide annotation, GeneratingFunction<?> gf)
	{
		this.spec = spec;
		this.annotation = annotation;
		this.gf = gf;
	}
	
	public Spectrum getSpec() {
		return spec;
	}

	public Peptide getAnnotation() {
		return annotation;
	}

	public GeneratingFunction<?> getGf() {
		return gf;
	}

	public ProfileGF<?> getProfGF() {
		return profGF;
	}

	private Spectrum spec;
	private Peptide annotation;
	private GeneratingFunction<?> gf;
	private ProfileGF<?> profGF;
}
