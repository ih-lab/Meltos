package lineage;

import java.util.ArrayList;
import cc.mallet.types.Dirichlet;

public class NodeStats {
	
	private double[] means;
	private double[] vars;
	private int numSamples;
	private double[] alphas;
	private Dirichlet dir;
	
	
	public NodeStats(ArrayList<SVEntry> entries, int numS){
		numSamples=numS;
		setMeans(entries);
		setVariances(entries);
		//setDirich(entries);
	}
	
	public NodeStats(double[] artiMeans, double[] artiVars, int numS){
		numSamples = numS;
		means = artiMeans;
		vars = artiVars;
	}
	
	
	public NodeStats(double[] mean, double[] stdev, int numS, int[] indices){
		means = new double[indices.length];
		vars = new double[indices.length];
		int ind = 0;
		for(int x = 0; x<means.length; x++){
			if(indices[ind]>x||ind==mean.length){
				means[x]=0;
				vars[x]=Parameters.VAF_ERROR_MARGIN/4;//Arbitrary.
			}
			else{
				means[x]=mean[ind];
				vars[x]=Math.pow(stdev[ind],2);
				ind++;
			}
		}
	}
	
	public void setDirich(ArrayList<SVEntry> entries){
		alphas = new double[numSamples];
		for(SVEntry sv: entries)
			for(int x = 0; x<sv.approxVAF.length; x++)
				alphas[x]+=sv.approxVAF[x];
		dir=new Dirichlet(alphas);
	}
	
	public void setMeans(ArrayList<SVEntry> entries){
		means = new double[numSamples];
		if(entries.isEmpty())
			return;
		for(SVEntry entry: entries)
			for(int x = 0; x<numSamples; x++)
				means[x]+=entry.getVAF(x)/entries.size();
	}
	
	
	public void setVariances(ArrayList<SVEntry> entries){
		vars = new double[numSamples];
		if(entries.isEmpty())
			return;
		for(SVEntry entry: entries)
			for(int x = 0; x<numSamples; x++)
				vars[x]+=(entry.getVAF(x)-means[x])*(entry.getVAF(x)-means[x])/(entries.size()-1);
            

	}
	
	public double probGenerate(NodeStats node2){
		double prob = 1.0;
		
		for(int x = 0; x<numSamples; x++){
			double v1 = pdf(node2.means[x], means[x], Math.sqrt(vars[x]));
			double v2 = pdf(means[x], node2.means[x], Math.sqrt(node2.vars[x]));
			if(!Double.isNaN(v1))
				prob = prob*v1;
			if(!Double.isNaN(v2))
				prob = prob*v2;
			
		}
		
		return prob;
	}
	
	public double probGenerateSV(SVEntry sv){
		double prob = 1.0;
		for(int x = 0; x<numSamples; x++){
			double v1 = pdf(sv.approxVAF[x], means[x], Math.sqrt(vars[x]));
			if(!Double.isNaN(v1))
				prob = prob*pdf(sv.approxVAF[x], means[x], Math.sqrt(vars[x]));
		}
		return prob;
	}
	
	public static double pdf(double x) {
        return Math.exp(-x*x / 2) / Math.sqrt(2 * Math.PI);
    }

    // return pdf(x, mu, signma) = Gaussian pdf with mean mu and stddev sigma
    public static double pdf(double x, double mu, double sigma) {
        return pdf((x - mu) / sigma) / sigma;
    }

}
