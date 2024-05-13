package ru.ac.phyche.badpredictionsdataset;

public class DifferenceEntry {
	private float riRef;
	private float riCNN;
	private float riMLP;
	private int subset;
	private int[] fingerprints;
	private String smiles;

	public DifferenceEntry(float riRef, float riCNN, float riMLP, int subset, int[] fingerprints, String smiles) {
		this.riRef = riRef;
		this.riCNN = riCNN;
		this.riMLP = riMLP;
		this.subset = subset;
		this.fingerprints = fingerprints;
		this.smiles = smiles;
	}

	public float getRiRef() {
		return riRef;
	}

	public void setRiRef(float riRef) {
		this.riRef = riRef;
	}

	public float getRiCNN() {
		return riCNN;
	}

	public void setRiCNN(float riCNN) {
		this.riCNN = riCNN;
	}

	public float getRiMLP() {
		return riMLP;
	}

	public void setRiMLP(float riMLP) {
		this.riMLP = riMLP;
	}

	public int getSubset() {
		return subset;
	}

	public void setSubset(int subset) {
		this.subset = subset;
	}

	public int[] getFingerprints() {
		return fingerprints;
	}

	public void setFingerprints(int[] fingerprints) {
		this.fingerprints = fingerprints;
	}

	public String getSmiles() {
		return smiles;
	}

	public void setSmiles(String smiles) {
		this.smiles = smiles;
	}

}
