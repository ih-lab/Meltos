package lineage;

import java.util.ArrayList;


public class SampleStats {

	
	protected double mean;
	protected double std;
	public SampleStats(ArrayList<Double> VAFs){
		setMean(VAFs);
		setVariance(VAFs);
	}
	
	public SampleStats(ArrayList<SVEntry> VAFs, int index){
		int count=0;
		ArrayList<Double> sample = new ArrayList<Double>();
		for(SVEntry s: VAFs)
			sample.add(s.getVAF(index));
		setMean(sample);
		setVariance(sample);
	}
	
	public void setMean(ArrayList<Double> VAFs){
		mean = 0.0;
		int count = 0;
		if(VAFs.isEmpty())
			return;
		for(Double vaf: VAFs)
			if(vaf!=0.0){
				mean+=vaf;
				count++;
			}
		if(count!=0)
			mean = mean/count;
	}
	
	public void setVariance(ArrayList<Double> VAFs){
		std = 0.0;
		int count = 0;
		if(VAFs.isEmpty())
			return;
        for(Double a :VAFs)
        	if(a!=0){
        		std += (a-mean)*(a-mean);
        		count++;
        	}
        if(count>1)
			std = Math.sqrt(std/(count-1));
        else
        	std = 0.0;
	}
	
	public double convertVAF(double VAF, SampleStats template){
		if(template.std==0)
			return VAF;
		else
			return (VAF-template.mean)/template.std*std+mean;
		
	}
}
