/*
This program accepts a folder of .dta files, percentage data 
to be retained and the output folder as arguments. The output 
folder contains .ms2 file which can be directly used with 
Crux-Tide search software.

java msreduce ./inputFolder xx ./outputFolder
xx = percentage data to be retained e.g. 10, 20, 30.

Copyright (C) Muaaz Gul Awan and Fahad Saeed 

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/


import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import static java.lang.Math.floor;
import java.util.ArrayList;
import java.util.Random;
import java.text.DecimalFormat;

public class msreduce {

    public static String firstLine;
    static float sample_amount;
    static double dataSizeA = 0,dataSizeN=0,timeStart,timeEnd,totalTime=0,timeStartWithSort,timeEndWithSort,totalTimeWithSort,timeStartWidth,timeEndWidth,totalTimeWidth;
    static  String targetFolder;
   
    
    public static void main(String[] args) throws FileNotFoundException, IOException {
        targetFolder =  args[0];
        sample_amount = Float.parseFloat(args[1]);
        final File folder = new File(targetFolder);
        String temp;
		String outPutFolder = args[2];
        int fSize=0;
        float bias_val = (float) 10;
        float avgWidth=0;
        avgWidth = getAvgWidth(folder, targetFolder);
        
        PrintWriter writer = new PrintWriter(outPutFolder +sample_amount+ ".ms2", "UTF-8");
        int myScans = 0;
        for (final File fileEntry : folder.listFiles()) {
            if (fileEntry.isDirectory()) {
                listFilesForFolder(fileEntry);
            } else {
                temp = fileEntry.getName();
                
                ArrayList<String> vals = sampleReader(targetFolder + temp);
		//System.out.println("file: "+temp);
		myScans++;
		if(vals.size() >15){
                sampleWriter(temp, writer, sample_amount, vals, bias_val, fSize,avgWidth, myScans);
             //  System.out.println("size:"+vals.size());
               }
            }
        }
        System.out.println("Percentage of Data Retained: "+(((dataSizeN/dataSizeA)*100)));
        //System.out.println("time taken: "+totalTime);
		//System.out.println("peaks retained:"+dataSizeN);
		//System.out.println("peaks total:"+dataSizeA);
        writer.close();
    }

    public static ArrayList<String> sampleReader(String File_Name) throws FileNotFoundException, IOException {

        BufferedReader reader = new BufferedReader(new FileReader(File_Name));
        ArrayList<String> vals = new ArrayList<String>();
        String str;
        int i = 0;
        while ((str = reader.readLine()) != null) {
            if (i == 0) {
                firstLine = str;
                i++;
            } else {
                vals.add(str);
            }
        }

        reader.close();

        return vals;
    }

    public static void sampleWriter(String source, PrintWriter writer, float sample_val, ArrayList<String> data, float bias_val, int fSize, float avgWidth, int myScan) throws FileNotFoundException, UnsupportedEncodingException {
        
        String[] headsName = source.split("\\.");
        String[] headsFile = firstLine.split(" ");
        float prec_mass_peptide;
        ArrayList<Peak> dataSpectrum = new ArrayList<Peak>();
        ArrayList<ArrayList<Peak>> quants;
        float spectrumWidth;
        dataSizeA += data.size();
        dataSpectrum = stringToSpectrum(data);
        ArrayList<Peak> recursList = new ArrayList<Peak>();
        spectrumWidth = getWidth(dataSpectrum);
        float intensitySpread = (spectrumWidth/avgWidth)*100;
        float peaksReq = (sample_val/100)*data.size();

        timeStart = System.currentTimeMillis();        
        //below is the classification code
        //5 levels
        if (intensitySpread < 25){
            quants = quantization(dataSpectrum, 5);
            dataSpectrum = dynamicTune(quants, recursList, peaksReq, 5, 4);
        }
            
        
        //7 levels
        else if ((intensitySpread >=25) && (intensitySpread  < 50)){
            quants = quantization(dataSpectrum, 7);
            dataSpectrum = dynamicTune(quants, recursList, peaksReq, 5, 6);
        }
           
        //9 levels
        else if ((intensitySpread >=50) && (intensitySpread  < 75)){
            quants = quantization(dataSpectrum, 9);
             dataSpectrum = dynamicTune(quants, recursList, peaksReq, 5, 8);
        }
        //11 levels
        else if ((intensitySpread)>=75){
            quants = quantization(dataSpectrum, 11);
            dataSpectrum = dynamicTune(quants, recursList, peaksReq, 5, 10);
        }
        
      
        
        dataSpectrum = sortSpectraInsertion(dataSpectrum);

        timeEnd = System.currentTimeMillis();

        totalTime += timeEnd - timeStart;
        dataSizeN +=dataSpectrum.size();
        //converting protonated mass value of precursor ion from .dta file to m/z of the precursor ion
        prec_mass_peptide = (float) ((Float.parseFloat(headsFile[0]) + (Float.parseFloat(headsFile[1]) - 1))) / (Float.parseFloat(headsFile[1]));
        if(headsName.length > 2)
           writer.println("S\t" + headsName[2] + "\t" + headsName[2] + "\t" + String.format("%.6f",prec_mass_peptide));
        else
           writer.println("S\t" + " " + "\t" + " " + "\t" + String.format("%.6f",prec_mass_peptide));  
        writer.println("Z\t" + headsFile[1] + "\t" + headsFile[0]);
//writer.flush();
        //writing everything to files now
        for (int i = 0; i < dataSpectrum.size(); i++) {
            writer.println(String.format("%.6f",dataSpectrum.get(i).m_z) + "\t" + String.format("%.6f",dataSpectrum.get(i).intensity));
        writer.flush();
		}
		

    }

	//lists all the files in a given folder
    public static void listFilesForFolder(final File folder) throws FileNotFoundException, UnsupportedEncodingException {
        PrintWriter writer = new PrintWriter("files.txt", "UTF-8");
        for (final File fileEntry : folder.listFiles()) {
            if (fileEntry.isDirectory()) {
                listFilesForFolder(fileEntry);
            } else {
                writer.println(fileEntry.getName());
            }
        }
    }

	//sorts spectra 
    public static ArrayList<Integer> sort(ArrayList<Integer> lst) {

        int temp;
        int largest = 0;
        for (int j = 0; j < lst.size(); j++) {
            for (int i = 0; i < lst.size() - 1; i++) {
                largest = lst.get(i);
                temp = lst.get(i + 1);
                if (temp < largest) {
                    lst.set(i, temp);
                    lst.set(i + 1, largest);

                }
            }
        }
        return lst;
    }

	//sorts spectra wrt m_z values
    public static ArrayList<Peak> sortSpectra(ArrayList<Peak> input) {
        ArrayList<Peak> result = new ArrayList<Peak>();
        Peak temp = new Peak();
        Peak max = new Peak();
        if (input.size() > 0) {
            max.m_z = input.get(0).m_z;
        }
        for (int j = 0; j < input.size(); j++) {
            for (int i = 0; i < input.size() - 1; i++) {
                max = input.get(i);
                temp = input.get(i + 1);
                if (temp.m_z < max.m_z) {
                    input.set(i, temp);
                    input.set(i + 1, max);
                }
            }
        }
        return input;
    }
    
    //sorts spectra wrt m_z values using insertion sort
      public static ArrayList<Peak> sortSpectraInsertion(ArrayList<Peak> input) {
        int j;
        for(int i = 1; i < input.size(); i++){
            j=i;
            while(j>0 && input.get(j-1).m_z > input.get(j).m_z){
                Peak temp = input.get(j-1);
                input.set(j-1, input.get(j));
                input.set(j, temp);
                j--;
            }
            
            
        }
        
        
        return input;
    }
    
      
      //sorts spectra by intensity
    public static ArrayList<Peak> sortSpectraByIntensity(ArrayList<Peak> input) {
        Peak temp = new Peak();
        Peak max = new Peak();
        if (input.size() > 0) {
            max.intensity = input.get(0).intensity;
        }
        for (int j = 0; j < input.size(); j++) {
            for (int i = 0; i < input.size() - 1; i++) {
                max = input.get(i);
                temp = input.get(i + 1);
                if (temp.intensity < max.intensity) {
                    input.set(i, temp);
                    input.set(i + 1, max);
                }
            }
        }
        return input;
    }
    
    //sorts spectra by intensity using isnertion sort
    public static ArrayList<Peak> sortSpectraByIntensityInsertion(ArrayList<Peak> input) {
        int j;
        for(int i = 1; i < input.size(); i++){
            j=i;
            while(j>0 && input.get(j-1).intensity > input.get(j).intensity){
                Peak temp = input.get(j-1);
                input.set(j-1, input.get(j));
                input.set(j, temp);
                j--;
            }
            
            
        }
        return input;
    }
    
   //performs skewed sampling on spectra 
    public static ArrayList<Peak> skewSampling(ArrayList<Peak> input, float threshold, Peak maxIntensity) {
        ArrayList<Peak> output = new ArrayList<Peak>();
        float thresholdVal = threshold;
        float thresholdIntensity;
        thresholdIntensity = thresholdVal * maxIntensity.intensity;
        for (int i = 0; i < input.size(); i++) {
            if (input.get(i).intensity > thresholdIntensity) {
                output.add(input.get(i));
            }
        }

        return output;
    }

	//finds maximum wrt intensity
    public static Peak maxPeak(ArrayList<Peak> input) {
        Peak max = new Peak();
        max = input.get(0);
        for (int i = 0; i < input.size(); i++) {
            if (max.intensity < input.get(i).intensity) {
                max = input.get(i);
            }

        }
        return max;
    }
    //finds smallest peak wrt intensity
    public static Peak minPeak(ArrayList<Peak> input) {
        Peak min = new Peak();
        min = input.get(0);
        for (int i = 0; i < input.size(); i++) {
            if (min.intensity > input.get(i).intensity) {
                min = input.get(i);
            }

        }
        return min;
    }

	//randomly samples a spectrum based on the sampling rate provided
    public static ArrayList<String> randomSampler(ArrayList<String> data, float sample_val) {
        Random random = new Random();
        ArrayList<Integer> indexList = new ArrayList<Integer>();
        ArrayList<String> outList = new ArrayList<String>();
        float percentage = (float) sample_val / 100;
        int ind = 0;

        while (indexList.size() < floor(data.size() * percentage)) {
            ind = random.nextInt(data.size());
            if (!indexList.contains(ind)) {
                indexList.add(ind);
            }
        }

        for (int i = 0; i < floor(data.size() * percentage); i++) {
            outList.add(data.get(indexList.get(i)));
        }

        return outList;
    }

	/*
	converts string-spectrum to spectrum-class
	*/
    public static ArrayList<Peak> stringToSpectrum(ArrayList<String> input) {
        ArrayList<Peak> output = new ArrayList<Peak>();
        Peak tempSpectrum;
        String temp;
        String[] temp2;

        for (int i = 0; i < input.size(); i++) {
            tempSpectrum = new Peak();
            temp = input.get(i);
            temp2 = temp.split("\\ ");
            tempSpectrum.intensity = Float.parseFloat(temp2[1]);
            tempSpectrum.m_z = Float.parseFloat(temp2[0]);
            output.add(i, tempSpectrum);
            
        }

        return output;
    }
	
    //creates fsets
    public static ArrayList<Fset> createFsets(ArrayList<Peak> data, int fSize) {
        ArrayList<Fset> Fsets = new ArrayList<Fset>();

        for (int i = 0; i < (data.size() - (fSize - 1)); i++) {
            Fset temp = new Fset();
            for (int j = i; j < i + fSize; j++) {
                temp.sets.add(data.get(j));
            }
            Fsets.add(temp);
        }

        return Fsets;
    }

	//randomly samples Fsets
	
    public static ArrayList<Fset> randomSamplerFset(ArrayList<Fset> data, float sample_val) {
        Random random = new Random();
        ArrayList<Integer> indexList = new ArrayList<Integer>();
        ArrayList<Fset> outList = new ArrayList<Fset>();
        float percentage = (float) sample_val / 100;
        int ind = 0;

        while (indexList.size() < floor(data.size() * percentage)) {
            ind = random.nextInt(data.size());
            if (!indexList.contains(ind)) {
                indexList.add(ind);
            }
        }

        for (int i = 0; i < floor(data.size() * percentage); i++) {
            outList.add(data.get(indexList.get(i)));
        }

        return outList;
    }

	// merges Fsets
    public static ArrayList<Peak> FsetsMergeToSpectrum(ArrayList<Fset> Fsets) {
        ArrayList<Peak> spectrum = new ArrayList<Peak>();

        for (int i = 0; i < Fsets.size(); i++) {
            for (int j = 0; j < Fsets.get(i).sets.size(); j++) {
                if (!spectrum.contains(Fsets.get(i).sets.get(j))) {
                    spectrum.add(Fsets.get(i).sets.get(j));
                }
            }
        }

        return spectrum;
    }

	/*
	   filters fsets based on a threshold 
	  */
    public static Fset filterFset(Fset setIn, int opt, float threshold) {
        Fset setOut = new Fset();
        Peak maxPeak,minPeak;
        float cutOffVal;
        if (opt == 1) {
            maxPeak = maxPeak(setIn.sets);
            cutOffVal = maxPeak.intensity * (threshold / 100);
            for (int i = 0; i < setIn.sets.size(); i++) {
                if (setIn.sets.get(i).intensity > cutOffVal) {
                    setOut.sets.add(setIn.sets.get(i));
                }
            }
        } else if (opt == 0) {
            minPeak = minPeak(setIn.sets);
            cutOffVal = minPeak.intensity * (threshold / 100);
            for (int i = 0; i < setIn.sets.size(); i++) {
                if (setIn.sets.get(i).intensity > cutOffVal) {
                    setOut.sets.add(setIn.sets.get(i));
                }
            }
        } else if (opt == 2) {
            float sum = 0;
            for (int i = 0; i < setIn.sets.size(); i++) {
                sum += setIn.sets.get(i).intensity;
            }

            cutOffVal = (sum / setIn.sets.size()) * (threshold / 100);

            for (int i = 0; i < setIn.sets.size(); i++) {
                if (setIn.sets.get(i).intensity > cutOffVal) {
                    setOut.sets.add(setIn.sets.get(i));
                }
            }
        }
        
        else
            System.out.println("ERRORR!!!! Wrong Selection for Filter");

        return setOut;
    }
    
	 /*
	   Fset related calculations
	  */
    public static ArrayList<Fset> filterFsetList(ArrayList<Fset> listIn, float threshold, int opt){
        ArrayList<Fset> listOut = new ArrayList<Fset>();
        for( int i = 0; i < listIn.size(); i++)
            listOut.add(filterFset(listIn.get(i),opt , threshold));
        
        return listOut;
    }
    
	
	  /*
	   calculates the average intensity of largest 10 peaks
	  */
	 
      public static float avgMax10Peaks(ArrayList<Peak> unsortedIn){
        float avg=0,sum = 0;
        ArrayList<Peak> sortedIn = sortSpectraByIntensityInsertion(unsortedIn);
        ArrayList<Peak> max10 = new ArrayList<Peak>();
        if(sortedIn.size() > 14){
        for( int i = sortedIn.size()-1; i > sortedIn.size()-11 ; i--)
            max10.add(sortedIn.get(i));
        
        for (int i = 0; i < max10.size(); i++)
            sum += max10.get(i).intensity;
            
            }
        avg = sum/10;
        return avg;
    }
   
       /*
	  calculates the average intensity of largest 3 peaks
	 */
      public static float avgMax3Peaks(ArrayList<Peak> unsortedIn){
        float avg=0,sum = 0;
        ArrayList<Peak> sortedIn = sortSpectraByIntensityInsertion(unsortedIn);
        ArrayList<Peak> max3 = new ArrayList<Peak>();
        
        for( int i = sortedIn.size()-1; i > sortedIn.size()-4 ; i--)
            max3.add(sortedIn.get(i));
        
        for (int i = 0; i < max3.size(); i++)
            sum += max3.get(i).intensity;
        avg = sum/3;
        return avg;
    }
   
     /*
	  calculates the average intensity of smallest 10 peaks
	 */
       public static float avgMin10Peaks(ArrayList<Peak> unsortedIn){
           float avg=0,sum = 0;
           ArrayList<Peak> sortedIn = sortSpectraByIntensityInsertion(unsortedIn);
           ArrayList<Peak> min10 = new ArrayList<Peak>();
           if(unsortedIn.size() > 14){
           for( int i = 0; i <10 ; i++)
               min10.add(sortedIn.get(i));
        
           for (int i = 0; i < min10.size(); i++)
               sum += min10.get(i).intensity;
               
               }
           avg = sum/10;
           return avg;
    }
	
      //the output list is unsorted.
	        
	  	/**
	     quantizes given spectra 
	
	    @param ArrayList<Peak> listIn
	     input spectra
	  
	    @param int numSamples
	     number of quanta
		 
	    **/
	  
	  
      public static ArrayList<ArrayList<Peak>> quantization(ArrayList<Peak> listIn,int numSamples){
          ArrayList<ArrayList<Peak>> outList = new ArrayList<ArrayList<Peak>>();
          ArrayList<Peak> tempList = new ArrayList<Peak>();
          float maxAvg, minAvg,testVar=0;
          float jumpInc = (float)1/numSamples,jump=0,jumpLag=0;
          float refVal = avgMax10Peaks(listIn);
          
          for (int i = 0; i < numSamples; i++){
              if(i == 0){
              jumpLag = 0;
              jump=jumpInc;
              }
              for (int j = 0; j < listIn.size(); j++){
                  if(i == numSamples-1){
                      if((listIn.get(j).intensity >= (refVal*jumpLag))){
                         tempList.add(listIn.get(j));
                      }
                      else{
                      }
                  }
                  else{
                      if((listIn.get(j).intensity >= (refVal*jumpLag)) && (listIn.get(j).intensity < (refVal*jump))){
                      tempList.add(listIn.get(j));
                      }        
                  }
              }
              outList.add( new ArrayList<Peak>(tempList));
              tempList.clear();
              
                      jump=jump+jumpInc;
                      jumpLag +=jumpInc;
          }
          
       
      
          return outList;
      }
      
	  	/**
	     randomly samples give spectra
	
	    @param ArrayList<Peak> data
	     input spectra
	  
	    @param float sample_val 
	     sampling rate
		 
	 **/
	  
      public static ArrayList<Peak> randomSamplerPeaks(ArrayList<Peak> data, float sample_val) {
        Random random = new Random();
        ArrayList<Integer> indexList = new ArrayList<Integer>();
        ArrayList<Peak> outList = new ArrayList<Peak>();
        float percentage = (float) sample_val / 100;
        int ind = 0;

        while (indexList.size() < floor(data.size() * percentage)) {
            ind = random.nextInt(data.size());
            if (!indexList.contains(ind)) {
                indexList.add(ind);
            }
        }

        for (int i = 0; i < floor(data.size() * percentage); i++) {
            outList.add(data.get(indexList.get(i)));
        }

        return outList;
    }
      
	  
      public static float getWidth(ArrayList<Peak> listIn){
          float width,maxAvg,minAvg;
          
          maxAvg = avgMax10Peaks(listIn);
          minAvg = avgMin10Peaks(listIn);
          width = maxAvg-minAvg;
          
      return width;
      }
      
	    /**
	     recursively evaluates sampling weights for each quanta.
	
	    @param ArrayList<ArrayList<Peak>> quantsIn 
	    quantized spectra 
	  
	    @param String ArrayList<Peak> sampledList
	     list containing the sampled spectra
		 
		@param float numOfPeaksReq
		 percentage of peaks required to be sampled
		 
		 
		@param int tolPeaks
		 tolerance in number of peaks.
		 
		@param int index 
		 starting index
	 **/
      public static ArrayList<Peak> dynamicTune(ArrayList<ArrayList<Peak>> quantsIn, ArrayList<Peak> sampledList,float numOfPeaksReq, int tolPeaks,int index){
          
          if ((quantsIn.get(index).size() >= numOfPeaksReq-tolPeaks) && (quantsIn.get(index).size() <= numOfPeaksReq+tolPeaks)){
              for(int i=0; i < quantsIn.get(index).size(); i++)
                  sampledList.add(quantsIn.get(index).get(i));
              return sampledList;
          }
          else if(quantsIn.get(index).size()> numOfPeaksReq+tolPeaks){
              ArrayList<Peak> tempoHold = new ArrayList<Peak>();
              tempoHold = randomSamplerPeaks(quantsIn.get(index), (numOfPeaksReq/quantsIn.get(index).size())*100);
              for (int i = 0; i < tempoHold.size(); i++)
                  sampledList.add(tempoHold.get(i));
              return sampledList;
          }
          
          else{
              for(int i=0; i < quantsIn.get(index).size(); i++)
                  sampledList.add(quantsIn.get(index).get(i));
              numOfPeaksReq = numOfPeaksReq - quantsIn.get(index).size();
              index = index -1;
              return dynamicTune(quantsIn, sampledList, numOfPeaksReq, tolPeaks, index);
          }
      }
      
    /**
	gets average width for spectra  
	
	  @param File folder 
	  folder object containing the spectra files for which width is to be calculated.
	  
	  @param String targetFolder
	  address to the folder which contains spectra files
	 **/
      public static float getAvgWidth(File folder,String targetFolder) throws FileNotFoundException, UnsupportedEncodingException, IOException{
          float avgWidth=0,width=0,sum=0,total=0;
          String temp;
          ArrayList<Peak> testing;
                 
        for (final File fileEntry : folder.listFiles()) {
            if (fileEntry.isDirectory()) {
                listFilesForFolder(fileEntry);
            } else {
                temp = fileEntry.getName();
                ArrayList<String> vals = sampleReader(targetFolder + temp);
		//System.out.println("file: "+temp);
               testing = stringToSpectrum(vals);
               
               timeStartWidth = System.currentTimeMillis();
               width = getWidth(testing);
               sum+=width;
               total++;
              timeEndWidth = System.currentTimeMillis();
               
               totalTimeWidth += timeEndWidth-timeStartWidth  ;
            }
        }
          avgWidth = sum/total;
          return avgWidth;
      }
}

