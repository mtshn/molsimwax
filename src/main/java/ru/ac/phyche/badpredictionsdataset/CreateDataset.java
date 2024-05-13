package ru.ac.phyche.badpredictionsdataset;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import javax.naming.spi.DirStateFactory.Result;

public class CreateDataset {

	private static String[] loadNonEmptyLines(String filename) throws IOException {
		ArrayList<String> result = new ArrayList<String>();
		BufferedReader br = new BufferedReader(new FileReader(filename));
		String s = br.readLine();
		while (s != null) {
			if (!s.trim().equals("")) {
				result.add(s.trim());
			}
			s = br.readLine();
		}
		br.close();
		return result.toArray(new String[result.size()]);
	}

	private static String[] groupDataSetLinesByCanonicalSmilesMedianValues(String[] a) {
		HashMap<String, ArrayList<Float>> ref = new HashMap<String, ArrayList<Float>>();
		HashMap<String, ArrayList<Float>> pred = new HashMap<String, ArrayList<Float>>();
		String[] can = new String[a.length];
		Arrays.stream(ArUtls.intsrnd(a.length)).parallel().forEach(i -> {
			String smiles = a[i].split("\\s+")[0];
			try {
				String smilesCan = ChemUtils.canonical(smiles, false);
				smilesCan = ChemUtils.canonical(smilesCan, false);
				can[i] = a[i].replace(smiles, smilesCan);
			} catch (Exception e) {
				e.printStackTrace();
				throw new RuntimeException(e.getLocalizedMessage());
			}
		});
		for (int i = 0; i < can.length; i++) {
			String[] sp = can[i].split("\\s+");
			String smiles = sp[0];
			float reff = Float.parseFloat(sp[1]);
			float predf = Float.parseFloat(sp[2]);
			if (ref.get(smiles) == null) {
				ArrayList<Float> al = new ArrayList<Float>();
				al.add(reff);
				ref.put(smiles, al);
			} else {
				ArrayList<Float> al = ref.get(smiles);
				al.add(reff);
				ref.put(smiles, al);
			}
			if (pred.get(smiles) == null) {
				ArrayList<Float> al = new ArrayList<Float>();
				al.add(predf);
				pred.put(smiles, al);
			} else {
				ArrayList<Float> al = pred.get(smiles);
				al.add(predf);
				pred.put(smiles, al);
			}
		}
		if (pred.size() != ref.size()) {
			new RuntimeException("Unknown parsing error");
		}
		String[] result = new String[pred.size()];
		int j = 0;
		String[] smilesArray = pred.keySet().toArray(new String[pred.keySet().size()]);
		Arrays.sort(smilesArray);
		for (String s : smilesArray) {
			float reff = ArUtls.median(ArUtls.toFloatArray(ref.get(s)));
			float predf = ArUtls.median(ArUtls.toFloatArray(pred.get(s)));
			result[j] = s + " " + reff + " " + predf;
			j++;
		}
		return result;
	}

	public static DifferenceEntry[] loadCVResultsWAX(String path) throws IOException {
		ArrayList<DifferenceEntry> result = new ArrayList<DifferenceEntry>();
		for (int i = 0; i < 5; i++) {
			String[] cnn = loadNonEmptyLines(path + "/cnn" + i + ".txt");
			String[] mlp = loadNonEmptyLines(path + "/mlp" + i + ".txt");
			String[] av = loadNonEmptyLines(path + "/average" + i + ".txt");
			for (int j = 0; j < cnn.length; j++) { // testing of files
				if (!cnn[j].split("\\s+")[0].equals(mlp[j].split("\\s+")[0])) {
					throw (new IOException("Files cnn and mlp do not match each other"));
				}
				float cnn2 = Float.parseFloat(cnn[j].split("\\s+")[2]);
				float mlp2 = Float.parseFloat(mlp[j].split("\\s+")[2]);
				float av2 = Float.parseFloat(av[j].split("\\s+")[2]);
				if (Math.abs(av2 * 2 - mlp2 - cnn2) > 2) {
					throw (new IOException("Files cnn and mlp do not match each other"));
				}
				float cnn1 = Float.parseFloat(cnn[j].split("\\s+")[1]);
				float mlp1 = Float.parseFloat(mlp[j].split("\\s+")[1]);
				if (Math.abs(mlp1 - cnn1) > 1) {
					throw (new IOException("Files cnn and mlp do not match each other"));
				}
			}
			cnn = groupDataSetLinesByCanonicalSmilesMedianValues(cnn);
			mlp = groupDataSetLinesByCanonicalSmilesMedianValues(mlp);
			for (int j = 0; j < cnn.length; j++) { // testing of arrays
				if (!cnn[j].split("\\s+")[0].equals(mlp[j].split("\\s+")[0])) {
					throw (new IOException("Arrays cnn and mlp do not match each other"));
				}
			}
			if (cnn.length != mlp.length) {
				throw (new IOException("Arrays cnn and mlp do not match each other"));
			}
			for (int j = 0; j < cnn.length; j++) {
				String[] cnnsp = cnn[j].split("\\s+");
				String[] mlpsp = mlp[j].split("\\s+");
				String smiles = cnnsp[0];
				float ref = Float.parseFloat(cnnsp[1]);
				float cnnf = Float.parseFloat(cnnsp[2]);
				float mlpf = Float.parseFloat(mlpsp[2]);
				DifferenceEntry de = new DifferenceEntry(ref, cnnf, mlpf, i, null, smiles);
				result.add(de);
			}
		}
		DifferenceEntry[] resultA = result.toArray(new DifferenceEntry[result.size()]);
		FingerprintsSimilarity.calculateAllFingerprints(resultA);
		return resultA;
	}

	public static void main(String[] args) throws IOException {
		DifferenceEntry[] de = loadCVResultsWAX("./data");
		int[][] similarMolecules = FingerprintsSimilarity.mostSimilarMoleculesFromOtherSubsets(de, 100);
		FileWriter fw = new FileWriter("tmp");
		for (int i=0;i<de.length;i++) {
			fw.write(de[i].getSmiles() + " " + de[i].getRiRef() + " " + de[i].getRiCNN() + " " + de[i].getRiMLP() + " " + de[i].getSubset());
			for(int j=0;j<100;j++) {
				fw.write(" "+similarMolecules[i][j]);
			}
			fw.write("\n");
		}
		fw.close();
	}

}
