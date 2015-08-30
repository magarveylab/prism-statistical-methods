package data;

import java.util.ArrayList;
import java.util.List;

public class SmilesPackage {

	private String real;
	private List<String> prism = new ArrayList<String>();
	private List<String> npsearcher = new ArrayList<String>();
	private List<String> antismash = new ArrayList<String>();
	
	public String real() {
		return real;
	}
	
	public void setReal(String real) {
		this.real = real;
	}
	
	public List<String> antismash() {
		return antismash;
	}
	
	public void setAntismash(List<String> antismash) {
		this.antismash = antismash;
	}

	public List<String> npsearcher() {
		return npsearcher;
	}

	public void setNpsearcher(List<String> npsearcher) {
		this.npsearcher = npsearcher;
	}

	public List<String> prism() {
		return prism;
	}

	public void setPrism(List<String> prism) {
		this.prism = prism;
	}
	
}
