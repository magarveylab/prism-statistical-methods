package fingerprinter;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.BitSet;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.CircularFingerprinter;
import org.openscience.cdk.fingerprint.IFingerprinter;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.similarity.Tanimoto;

public class Ecfp6Fingerprinter {
	
	private static String file;
	private static String smilesFile;
	private static String outFile;
	private static String delimiter = "\t";
	private static String smiles;
	private static IFingerprinter fingerprinter = new CircularFingerprinter(4);
	
	public static void main(String[] args) throws IOException, CDKException {
		if (args.length > 1) {
			file = args[0];
			smilesFile = args[1];
			if (args.length > 2)
				outFile = args[2];
		} else {
			if (args.length == 0)
				throw new IllegalArgumentException("Error: input file must be specified!");
			if (args.length == 1)
				throw new IllegalArgumentException("Error: true SMILES must be specified!");
		}
		
		System.out.println("Reading input file " + file);
		System.out.println("Reading true SMILES file " + smilesFile);
		if (outFile != null)
			System.out.println("Printing output to file " + outFile);
		
		// create output writer
		BufferedWriter bw = null;
		if (outFile != null) 
			bw = new BufferedWriter(new FileWriter(outFile));
		
		// read smiles
		BufferedReader br = new BufferedReader(new FileReader(smilesFile));
		String line;
		while ((line = br.readLine()) != null) 
			smiles = line;
		br.close();
		
		// read library
		br = new BufferedReader(new FileReader(file));
		int lineIdx = 1;
		try {
			while ((line = br.readLine()) != null) {
				if (lineIdx == 1) {
					// detect the delimiter 
					if (line.split("\t").length == 2) {
						delimiter = "\t";
					} else if (line.split(",").length == 2) {
						delimiter = ",";
					}
				}
				
				// split the line
				String[] split = line.split(delimiter);
				if (split.length < 2)
					throw new IOException("Error: could not parse input on line " 
							+ lineIdx + ": not enough columns!");
				
				// get name and SMILES
				String name = split[0];
				String predictedSmiles = split[1];
				
				// calculate Tanimoto coefficient
				IAtomContainer queryMol = SmilesIO.molecule(smiles);
				BitSet queryFp = fingerprinter.getBitFingerprint(queryMol).asBitSet();
				IAtomContainer predictedMol = SmilesIO.molecule(predictedSmiles);
				BitSet predictedFp = fingerprinter.getBitFingerprint(predictedMol).asBitSet();
				Float tc = Tanimoto.calculate(queryFp, predictedFp);

				// output
				if (bw != null) {
					bw.append(name + "\t" + tc + "\t" + predictedSmiles + "\n");
				} else {
					System.out.println(name + "\t" + tc + "\t" + predictedSmiles);
				}
				
				// continue
				lineIdx++;
			}
		} finally {
			br.close();
		}
		
		System.out.println("Done.");
		System.exit(0);
	}

}
