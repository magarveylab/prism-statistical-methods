package accuracy;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVPrinter;
import org.apache.commons.math.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.inference.TTest;
import org.openscience.cdk.exception.CDKException;

import data.Ecfp6Fingerprint;
import data.SmilesPackage;
import data.TanimotoPackage;
import util.Files;
import util.Strings;
import exception.BadSmilesToFingerprinterException;
import exception.NoBitsetException;
import fingerprinter.SmilesIO;

public class AccuracyTest {

	private static String runDir = "/Users/michaelskinnider/Desktop/ecfp6/v4/";
	private static String parentDir = "/Users/michaelskinnider/Desktop/phylo_analysis_3/";
	private static String directory = "";
	
	public static void main(String[] args) throws IOException,
			BadSmilesToFingerprinterException, InterruptedException,
			CDKException, NoBitsetException {
		if (args.length > 0)
			directory = args[0];
		
		Test test = new Test();
		test.run();
	}

	private static class Test {

		public void run() throws IOException,
				BadSmilesToFingerprinterException, InterruptedException,
				CDKException, NoBitsetException {
			String[] subdirectories = Files.getAllDirectories(parentDir);
			if (subdirectories != null) {
				for (String subdirectory : subdirectories) {
					Map<String,SmilesPackage> results = new HashMap<String,SmilesPackage>();
					List<TanimotoPackage> tanimoto = new ArrayList<TanimotoPackage>();

					directory = parentDir + File.separator + subdirectory;
					loadKeys(results, tanimoto);

					parseAntismashOutput(results, tanimoto);
					parseNpsearcherOutput(results, tanimoto);
					parsePrismOutput(results, tanimoto);

					calculateTanimotoCoefficients(results, tanimoto);
					calculatePackageScores(results, tanimoto);
				}
			}
		}

		/**
		 * Load the keys for the hashmap of results, by instantiating a new
		 * hashmap entry for each cluster.
		 * 
		 * @param directory
		 *            the cluster directory
		 * @throws IOException
		 */
		public void loadKeys(Map<String,SmilesPackage> results, List<TanimotoPackage> tanimoto) throws IOException {
			BufferedReader br = new BufferedReader(new FileReader(directory + File.separator + "smiles.txt"));
			String line = null;
			while ((line = br.readLine()) != null) {
				String[] split = line.split("\t");
				String key = split[0];
				String smiles = split[1];
				SmilesPackage sp = new SmilesPackage();
				sp.setReal(smiles);
				results.put(key, sp);
			}
			System.out.println("Analyzing " + results.size() + " clusters in master SMILES file");
			br.close();
		}

		public void parseAntismashOutput(Map<String,SmilesPackage> results, List<TanimotoPackage> tanimoto) throws IOException {
			String extension = ".smi";
			File[] files = Files.getDirectoryFiles(directory, extension);
			for (File file : files) {
				String name = Strings.name(file.getAbsolutePath());
				String key = name.replace(extension, "");

				SmilesPackage sp = results.get(key);
				if (sp == null)
					throw new NullPointerException("Error: No entry for key " + key + "!");

				BufferedReader br = new BufferedReader(new FileReader(file));
				String line = null;
				int i = 0;
				while ((line = br.readLine()) != null && i < 1) {
					String smiles = line.replaceAll("\\[R\\d\\]", "[*]");
					List<String> antismash = new ArrayList<String>();
					antismash.add(smiles);
					sp.setAntismash(antismash);
					i++;
				}
				br.close();
			}
		}

		public void parseNpsearcherOutput(Map<String,SmilesPackage> results, List<TanimotoPackage> tanimoto) throws IOException {
			String extension = ".smiles";
			File[] files = Files.getDirectoryFiles(directory, extension);
			for (File file : files) {
				String name = Strings.name(file.getAbsolutePath());
				String key = name.replace(extension, "");

				SmilesPackage sp = results.get(key);
				if (sp == null)
					throw new NullPointerException("Error: No entry for key " + key + "!");

				List<String> npsearcher = sp.npsearcher();
				BufferedReader br = new BufferedReader(new FileReader(file));
				String line = null;
				while ((line = br.readLine()) != null) {
					if (line.indexOf(">") == -1 && npsearcher.size() <= 50) {
						String smiles = line.replace("X", "*");
						
						try {
							SmilesIO.molecule(smiles);
						} catch (CDKException e) {
							System.out.println("Could not parse NP.searcher SMILES " + smiles);
							continue;
						}
						
						npsearcher.add(smiles);
					}
				}

				br.close();
			}
		}

		public void parsePrismOutput(Map<String,SmilesPackage> results, List<TanimotoPackage> tanimoto) throws IOException {
			String extension = ".txt";
			File[] files = Files.getDirectoryFiles(directory, extension);
			for (File file : files) {
				List<String> smiles = new ArrayList<String>();
				BufferedReader br = new BufferedReader(new FileReader(file));
				String line = null;
				while ((line = br.readLine()) != null) {
					smiles.add(line.split("\t")[1]);
				}
				br.close();

				String key = Strings.name(file.getAbsolutePath()).replace(extension, "");
				if (!key.equals("smiles")) {
					if (results.containsKey(key)) {
						SmilesPackage sp = results.get(key);
						sp.setPrism(smiles);
					} else {
						System.out.println("Error: no key in results set for " + key);
					}
				}
			}
		}

		public void calculateTanimotoCoefficients(Map<String,SmilesPackage> results, List<TanimotoPackage> tanimoto) throws IOException,
				BadSmilesToFingerprinterException, InterruptedException,
				CDKException {
			for (Map.Entry<String, SmilesPackage> entry : results.entrySet()) {
				String key = entry.getKey();
				SmilesPackage sp = entry.getValue();
				TanimotoPackage tp = new TanimotoPackage();

				List<String> antismash = sp.antismash();
				List<String> npsearcher = sp.npsearcher();
				List<String> prism = sp.prism();
				String query = sp.real();
				
				List<Ecfp6Fingerprint> antismashFps = GetEcfp6Fingerprints.get(
						antismash, query, runDir + "antismash_tmp.txt", runDir);
				List<Ecfp6Fingerprint> prismFps = GetEcfp6Fingerprints.get(
						prism, query, runDir + "prism_tmp.txt", runDir);
				List<Ecfp6Fingerprint> npsearcherFps = GetEcfp6Fingerprints.get(
						npsearcher, query, runDir + "npsearcher_tmp.txt", runDir);
				
				float antismashAvg = median(antismashFps);
				float npsearcherAvg = median(npsearcherFps);
				float prismAvg = median(prismFps);

				if (Float.isNaN(antismashAvg))
					antismashAvg = 0.0f;
				if (Float.isNaN(npsearcherAvg))
					npsearcherAvg = 0.0f;
				if (Float.isNaN(prismAvg))
					prismAvg = 0.0f;				
					
				tp.antismashMed = antismashAvg;
				tp.npsearcherMed = npsearcherAvg;
				tp.prismMed = prismAvg;
				
				tp.antismashTop = top(antismashFps);
				tp.npsearcherTop = top(npsearcherFps);
				tp.prismTop = top(prismFps);
				
				tanimoto.add(tp);

				System.out.println("Calculated tanimoto scores for " + key + "\n"
						+ "\tPRISM:\t\t" + prismFps.size() + " Tanimoto scores\t" 
						+ "Avg: " + tp.prismMed + "\tTop: " + tp.prismTop
						+ "\n\tNP.searcher:\t" + npsearcherFps.size() + " Tanimoto scores\t" 
						+ "Avg: " + tp.npsearcherMed + "\tTop: " + tp.npsearcherTop
						+ "\n\tantiSMASH:\t" + antismashFps.size() + " Tanimoto scores\t" 
						+ "Avg: " + tp.antismashMed + "\tTop: " + tp.antismashTop);
			}
		}

		public void calculatePackageScores(Map<String,SmilesPackage> results, List<TanimotoPackage> tanimoto) throws IOException {
			int prismSize = 0;
			int npsearcherSize = 0;
			int antismashSize = 0;
			
			double[] prismAvg = new double[tanimoto.size()];
			double[] npsearcherAvg = new double[tanimoto.size()];
			double[] antismashAvg = new double[tanimoto.size()];

			double[] prismTop = new double[tanimoto.size()];
			double[] npsearcherTop = new double[tanimoto.size()];
			double[] antismashTop = new double[tanimoto.size()];

			for (int i = 0; i < tanimoto.size(); i++) {
				TanimotoPackage tp = tanimoto.get(i);
				if (tp.prismMed >= 0)
					prismSize++;
				prismAvg[i] = tp.prismMed;
				prismTop[i] = tp.prismTop;
				
				if (tp.npsearcherMed >= 0) 
					npsearcherSize++;
				npsearcherAvg[i] = tp.npsearcherMed;
				npsearcherTop[i] = tp.npsearcherTop;
				
				if (tp.antismashMed >= 0) 
					antismashSize++;
				antismashAvg[i] = tp.antismashMed;
				antismashTop[i] = tp.antismashTop;
			}
			
			// do statistical testing
			System.out.println("Executing t test on " + tanimoto.size() + " samples");
			TTest ttest = new TTest();
			double prism_antismash = ttest.tTest(prismAvg, antismashAvg);
			double prism_npsearcher = ttest.tTest(prismAvg, npsearcherAvg);
			double antismash_npsearcher = ttest.tTest(antismashAvg, npsearcherAvg);
			
			StandardDeviation sd = new StandardDeviation();
			double prismSd = sd.evaluate(prismAvg);
			double npsearcherSd = sd.evaluate(npsearcherAvg);
			double antismashSd = sd.evaluate(antismashAvg);

			double prismAvgAvg = average(prismAvg);
			double npsearcherAvgAvg = average(npsearcherAvg);
			double antismashAvgAvg = average(antismashAvg);
			
			System.out.println("----------------------------------");
			System.out.println("----------------------------------");
			System.out.println("Analysis completed.");
			System.out.println("Generated " + prismSize + " PRISM predictions\n\t"
					+ "Average median Tanimoto coefficient: " + prismAvgAvg);
			System.out.println("Generated " + antismashSize + " antiSMASH predictions\n\t"
					+ "Average median Tanimoto coefficient: " + antismashAvgAvg);
			System.out.println("Generated " + npsearcherSize + " NP.searcher predictions\n\t"
					+ "Average median Tanimoto coefficient: " + npsearcherAvgAvg);
			System.out.println("----------------------------------");
			System.out.println("----------------------------------");
			System.out.println("Standard distribution for PRISM predictions: " + prismSd);
			System.out.println("Standard distribution for antiSMASH predictions: " + antismashSd);
			System.out.println("Standard distribution for NP.searcher predictions: " + npsearcherSd);
			System.out.println("----------------------------------");
			System.out.println("----------------------------------");
			System.out.println("p-value for PRISM/antiSMASH comparison: " + prism_antismash);
			System.out.println("p-value for PRISM/npsearcher comparison: " + prism_npsearcher);
			System.out.println("p-value for npsearcher/antiSMASH comparison: " + antismash_npsearcher);
			
			// write 
			String out = runDir + new File(directory).getName() + ".csv";
			CSVFormat csvFileFormat = CSVFormat.EXCEL;
			CSVPrinter cw = new CSVPrinter(new FileWriter(out), csvFileFormat);

			List<String> header = new ArrayList<String>();
			header.add("Software");
			header.add("Tc");
			cw.printRecord(header);
			
			// add PRISM Tc's
			for (int i = 0; i < tanimoto.size(); i++) {
				List<String> row = new ArrayList<String>();
				TanimotoPackage tp = tanimoto.get(i);
				row.add("PRISM");
				row.add(tp.prismMed + "");
				cw.printRecord(row);
			}

			// add antiSMASH Tc's
			for (int i = 0; i < tanimoto.size(); i++) {
				List<String> row = new ArrayList<String>();
				TanimotoPackage tp = tanimoto.get(i);
				row.add("antiSMASH");
				row.add(tp.antismashMed + "");
				cw.printRecord(row);
			}

			// add PRISM Tc's
			for (int i = 0; i < tanimoto.size(); i++) {
				List<String> row = new ArrayList<String>();
				TanimotoPackage tp = tanimoto.get(i);
				row.add("NP.searcher");
				row.add(tp.npsearcherMed + "");
				cw.printRecord(row);
			}
			
			cw.flush();
			cw.close();
		}
		
		public Float median(List<Ecfp6Fingerprint> fingerprints) {
			Float median;
			if (fingerprints.size() == 0)
				return 0.0f;
			
			// sort by Tc 
			Collections.sort(fingerprints, new Comparator<Ecfp6Fingerprint>() { 
				public int compare(Ecfp6Fingerprint f1, Ecfp6Fingerprint f2) {
					return Float.compare(f1.tc, f2.tc);
				}
			});
			
			// pick median
			if (fingerprints.size() % 2 == 0) {
				median = (fingerprints.get(fingerprints.size()/2).tc 
						+ fingerprints.get(fingerprints.size()/2 - 1).tc) / 2; 
			} else {
				median = fingerprints.get(fingerprints.size() / 2).tc;
			}
			
			return median;
		}
		
		public Float top(List<Ecfp6Fingerprint> fingerprints) {
			Float top = 0.0f;
			for (Ecfp6Fingerprint f : fingerprints)
				if (f.tc > top)
					top = f.tc;
			return top;
		}
		
		public double average(double[] array) {
			double sum = 0.0d;
			for (double d : array)
				sum+=d;
			return sum / array.length;
		}
		
	}

}
