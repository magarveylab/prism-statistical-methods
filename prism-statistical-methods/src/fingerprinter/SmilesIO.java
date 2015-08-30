package fingerprinter;

import java.io.IOException;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

public class SmilesIO {

	public static IAtomContainer molecule(String smiles) throws InvalidSmilesException, IOException {
		IAtomContainer mol = null;
		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		SmilesParser parser = new SmilesParser(builder);
		parser.setPreservingAromaticity(true);
		mol = parser.parseSmiles(smiles);
		return mol;
	}
	
}
