package accuracy;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import data.Ecfp6Fingerprint;

public class GetEcfp6Fingerprints {
		
	public static List<Ecfp6Fingerprint> get(List<String> querySmiles,
			String trueSmiles, String outPath, String dir) throws IOException,
			InterruptedException {
		if (querySmiles.size() == 0) {
			List<Ecfp6Fingerprint> singleton = new ArrayList<Ecfp6Fingerprint>();
			Ecfp6Fingerprint dummy = new Ecfp6Fingerprint();
			dummy.tc = 0.0f;
			dummy.name = "";
			dummy.smiles = "";
			singleton.add(dummy);
			return singleton;
		}
		
		String querySmilesPath = dir + "querysmiles_tmp.txt";
		BufferedWriter bw = new BufferedWriter(new FileWriter(querySmilesPath));
		int i = 1;
		for (String s : querySmiles) {
			String name = "Scaffold_" + i;
			bw.append(name + "\t" + s + "\n");
			i++;
		}
		bw.close();

		return get(querySmilesPath, trueSmiles, outPath, dir);
	}
	
	public static List<Ecfp6Fingerprint> get(String querySmilesPath,
			String trueSmiles, String outPath, String dir) throws IOException, InterruptedException {
		List<Ecfp6Fingerprint> fingerprints = new ArrayList<Ecfp6Fingerprint>();
		
		// write true smiles to file
		String trueSmilesPath = dir + "smiles_tmp.txt";
		BufferedWriter bw = new BufferedWriter(new FileWriter(trueSmilesPath));
		bw.append(trueSmiles);
		bw.close();
		
		// run ecfp6.jar
		String jarFile = dir + "ecfp6.jar";
		String[] cmd = { "java", "-jar", jarFile, querySmilesPath, trueSmilesPath };
		ProcessBuilder pb = new ProcessBuilder();
		pb.redirectErrorStream(true);
		pb.command(cmd);
		Process p = pb.start();
//		Process p = Runtime.getRuntime().exec(cmd);
		p.waitFor();

		// parse output
		BufferedReader in = new BufferedReader(new InputStreamReader(p.getInputStream()));
		BufferedWriter bw2 = new BufferedWriter(new FileWriter(outPath));
		String s;
		while ((s = in.readLine()) != null) {
			bw2.append(s);
			
			String[] split = s.split("\t");
			if (split.length < 3)
				continue;
			
			String name = split[0];
			String tcRaw = split[1];
			String smiles = split[2];
			Float tc = Float.parseFloat(tcRaw);
			Ecfp6Fingerprint fp = new Ecfp6Fingerprint();
			fp.name = name;
			fp.tc = tc;
			fp.smiles = smiles;
			fingerprints.add(fp);			
			System.out.println("Added Tanimoto coefficient " + tc + " for molecule " + name);
		}
		in.close();
		bw2.close();
		
		return fingerprints;
	}
	
}
