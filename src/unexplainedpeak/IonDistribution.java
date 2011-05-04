package unexplainedpeak;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;

import msutil.IonType;
import msutil.Peptide;

public class IonDistribution {
	int numBins;
	HashMap<IonType, HashMap<Integer, float[]>> map;
		
	public IonDistribution(int numBins){
		this.numBins = numBins;
		map = new HashMap<IonType, HashMap<Integer, float[]>>();
	}
	
	public void add(IonType ion, Peptide pep, float mz){
		//if(ion instanceof IonType.PrecursorIon) return;
		int pepLen = pep.size();
		float ratio = mz/pep.getParentMass();
		
		HashMap<Integer, float[]> dists;
		if(map.containsKey(ion))
			dists = map.get(ion);
		else
			dists = new HashMap<Integer, float[]>();
			
		float[] dist;
		if(dists.containsKey(pepLen))
			dist = dists.get(pepLen);
		else 
			dist = new float[numBins+1];
		
		dist[(int)(numBins * ratio)]++;
		
		dists.put(pepLen, dist);
		map.put(ion, dists);			
	}
	
	
	public float[] getDist(IonType ion, int pepLen){
		return map.get(ion).get(pepLen);
	}
	
	public HashMap<IonType, Float> getDist(float parentMass, float peakmz){	
		return getDist((int)(parentMass/121.16), (int)(numBins*(peakmz/parentMass)));
	}
	
	public HashMap<IonType, Float> getDist(int pepLen, int position){
		if(position > numBins) return null;
		HashMap<IonType, Float> dist = new HashMap<IonType, Float>();
		float sum = 0;
		for(IonType ion : map.keySet()){
			if(!map.get(ion).containsKey(pepLen)) continue;
			float toAdd = map.get(ion).get(pepLen)[position];
			sum += toAdd;
			if(dist.containsKey(ion))
				toAdd += dist.get(ion);
			dist.put(ion, toAdd);
		}
		if (sum == 0) return null;
		for(IonType ion : map.keySet()){
			if(dist.containsKey(ion))
			dist.put(ion, map.get(ion).get(pepLen)[position]/sum);
		}
		
		return dist;
	}
	
	
	public void write(String origfilename) throws IOException{
		String filename =  origfilename.substring(0,origfilename.lastIndexOf('.')) + "_dist.m";
		BufferedWriter out = new BufferedWriter(new FileWriter(filename));
		
		out.write(this.toString());
		
		out.close();
	}
	
	public String toString(){ // very slow!! do not call often 
		String s = "";
		int minlen = 7, maxlen = 30;
		int i=1;
		for(IonType ion : map.keySet()){
			s+="IonOrder(" + i + ",:)=";
			s+="'"+Util.getIonString(ion)+"';\n";
			i++;
		}
		
		int inOrder = 1;
		for(IonType ion : map.keySet()){
			
			int sum = 0;
			for(int len = minlen ;len<maxlen;len++){
				if(!map.get(ion).containsKey(len)) continue;
				s+="length(:," + inOrder + ", " + len + ")=[\t";
				for(float f: getDist(ion, len)){
					s += f+"\t";
					sum += f;
				}
				s+="];\n";
				
			}
			//s+="%total ="  + sum + "\n";
			inOrder++;
			
		}
		

		/*
		 rankLimit
%rank1
s= size(dist);
i=1:s(2);

t=zeros(9,s(2));
for l=1:9
t(l,:) = t(l,:) + l*ones(1,s(2));
end

%ions = [20,14,11,2,18,5,8,25,26];
ions = [26,25,8,5,18,2,11,14,20];
for j=7:4:s(3)
k=1:9;
createfigure1(i/s(2), t', dist(ions,:,j));hold all;
title(['Probability that each ion type appears in relative position - length : ' , num2str(j)],'FontSize',20)
grid on
end


%rank3
s= size(dist);
i=1:s(2);

t=zeros(9,s(2));
for l=1:9
t(l,:) = t(l,:) + l*ones(1,s(2));
end

ions = [33,9,5,21,15,11,32,2,23];
for j=7:4:s(3)
k=1:9;
createfigure2(i/s(2), t', dist(ions,:,j));hold all;
title(['Probability that each ion type appears in relative position - length : ' , num2str(j)],'FontSize',20)
grid on
end

%rank5
s= size(dist);
i=1:s(2);

t=zeros(9,s(2));
for l=1:9
t(l,:) = t(l,:) + l*ones(1,s(2));
end

ions = [11,10,5,24,16,12,37,1,26];
for j=7:4:s(3)
k=1:9;
createfigure3(i/s(2), t', dist(ions,:,j));hold all;
title(['Probability that each ion type appears in relative position - length : ' , num2str(j)],'FontSize',20)
grid on
end


		 * */
		
		for(int len = minlen ;len<maxlen;len++){
			for(int p=0;p<=numBins;p++){
				if(getDist(len, p) == null) continue;
				s+="dist(:, " + p + ", " + len + " )=[\t";
				for(IonType ion : map.keySet()){
					if(getDist(len, p).get(ion) != null)
						s += getDist(len, p).get(ion) + "\t";
					else s += "0\t";
				}
				s += "];\n";
			}
		}
		
		
		
		return s;
	}
}
