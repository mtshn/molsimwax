package ru.ac.phyche.badpredictionsdataset;

import java.util.Arrays;

import org.openscience.cdk.exception.CDKException;

public class FingerprintsSimilarity {
	private static int[] toInts(float fp[]) {
		int[] r = new int[fp.length];
		for (int i = 0; i < fp.length; i++) {
			r[i] = (int) Math.round(fp[i]);
		}
		return r;
	}

	private static int[] fp(String smiles) throws CDKException {
		float[] f = ChemUtils.fingerprints(smiles, ChemUtils.FingerprintsType.CIRCULAR_6_1024);
		return toInts(f);
	}

	private static float tanimoto(int[] fp1, int[] fp2) {
		int all1 = 0;
		int all2 = 0;
		int both12 = 0;
		for (int i = 0; i < fp1.length; i++) {
			if (fp1[i] > 0) {
				all1++;
			}
			if (fp2[i] > 0) {
				all2++;
			}
			if ((fp1[i] > 0) && (fp2[i] > 0)) {
				both12++;
			}
		}
		return ((float) (both12)) / ((float) all1 + (float) all2 - ((float) (both12)));
	}

	public static void calculateAllFingerprints(DifferenceEntry[] a) {
		Arrays.stream(ArUtls.intsrnd(a.length)).parallel().forEach(i -> {
			try {
				a[i].setFingerprints(fp(a[i].getSmiles()));
			} catch (CDKException e) {
				e.printStackTrace();
				throw new RuntimeException(e.getMessage());
			}
		});
	}

	private static class ComparableResult implements Comparable {
		float similarity;
		int number;

		@Override
		public int compareTo(Object o) {
			ComparableResult o1 = (ComparableResult) o;
			return (-1) * ((Float) similarity).compareTo(o1.similarity);
		}
	}

	private static int[] mostSimilarMoleculesFromOtherSubsets(DifferenceEntry query, DifferenceEntry[] dataSet, int n) {
		ComparableResult[] x = new ComparableResult[dataSet.length];
		for (int i = 0; i < dataSet.length; i++) {
			x[i] = new ComparableResult();
			x[i].number = i;
			x[i].similarity = tanimoto(query.getFingerprints(), dataSet[i].getFingerprints());
			if (query.getSubset() == dataSet[i].getSubset()) {
				x[i].similarity = -1;
			}
		}
		Arrays.sort(x);
		int[] result = new int[n];
		for (int i = 0; i < n; i++) {
			result[i] = x[i].number;
		}
		return result;
	}

	public static int[][] mostSimilarMoleculesFromOtherSubsets(DifferenceEntry[] dataSet, int n) {
		int[][] results = new int[dataSet.length][];
		Arrays.stream(ArUtls.intsrnd(dataSet.length)).parallel().forEach(i -> {
			results[i] = mostSimilarMoleculesFromOtherSubsets(dataSet[i], dataSet, n);
		});
		return results;
	}

}
