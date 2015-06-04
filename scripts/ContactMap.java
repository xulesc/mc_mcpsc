

import java.util.StringTokenizer;
import java.io.*;
import java.awt.*;

public class ContactMap{


    /* ###################################################### */
    /* ###################################################### */
    /* ###################################################### */
    /*  HERE IT STARTS ContactMap CLASS */
    /* ###################################################### */
    /* ###################################################### */
    /* ###################################################### */

    static int MAXCAATOMS = 1000;

    /* In myContacts[i][] we keep a list of residues to which i is connected*/
    private int myContacts[][];
    /* In numberContactsPerResidue[i] we keep the number of residues that are in contact with i*/
    private int numberContactsPerResidue[];
    private int residuesNumber;
    private String proteinId;
    private int totalContactsNumber;
    private double threshold;
    private PrintWriter fileToWrite;
    
    
    public ContactMap(String fileName){
     	readContactMap(fileName);
	threshold = -1;
	
    }

    public ContactMap()
    {

    }

    public void setName(String name)
    {
	proteinId=name;
    }

    /* The matrix is built in such a way that it is a triangular superior */
    public void addContact(Integer hI,Integer tI)
    {
	int h=hI.intValue();
	int t=tI.intValue();
	int tmpInteger;
	/* the pairs (key,value) in the hashtable are stored in such a way that key<=value*/
	if(h>t)
	    {
		tmpInteger=t;
		t= h;
		h=tmpInteger;
	    }
	myContacts[h][numberContactsPerResidue[h]]=t;
	numberContactsPerResidue[h]=numberContactsPerResidue[h]+1;
	totalContactsNumber ++;
    }




    
    public void print()
    {
	int i;
	int j;

	for(i=0;i<residuesNumber;i++)
	    {
		if(numberContactsPerResidue[i]>0)
		    {
			System.out.println("Residue "+i+" contacts:");
			for(j=0;j<numberContactsPerResidue[i];j++)
			    {
				System.out.println("                       -->"+myContacts[i][j]);
			    }
		    }

	    }

    }


    public boolean isInContact(int r, int contRes)
    {
	int i;
	boolean res = false;

	for(i=0;(i<this.getContactsNumber(r))&&(!res);i++)
	    {
		res = getContactResidueId(r,i) == contRes;
	    }

	return res;
    }


    public int getResiduesNumber()
    {
	return residuesNumber;
    }

    public int getContactResidueId(int r,int pos)
    {
	if (r>residuesNumber) return -1;
	if (pos>numberContactsPerResidue[r]) return -1;

	return myContacts[r][pos];
    }

    /** Gives the number of residues to which i is in contact with*/
    public int getContactsNumber(int i)
    {
	if(i>residuesNumber) return -1;

	return numberContactsPerResidue[i];
    }


/* We assume that the format for a contact map files is:
 nro_residues  # Number of Residues
 nro_contacts  # Number of Contacts
 resA1 resB1
 resA2 resB2
 ...........
 resAn resBn

where resAi, resBi are numbers.
*/
    private  void readContactMap(String fileName)
    {
	
	
	
	/* Variable definitions to read the instance file.*/ 
	File              instanceFile;
	FileInputStream   fileInputStream;
	InputStreamReader inputStreamReader;
	BufferedReader    bufferedReader;
	String            fileLine;
	StringTokenizer   aLineTokenizer;
	String            aLineToken;
	int               i;
	int               j;
	Integer           headResidue;
        Integer           tailResidue;
	int               numberOfContacts;
	
	try
	    {
		instanceFile      = new File(fileName);
		fileInputStream   = new FileInputStream(instanceFile);
		inputStreamReader = new InputStreamReader(fileInputStream);
		bufferedReader    = new BufferedReader(inputStreamReader); /* This is _actually_ from
									      where we are gona read */
		

		
		// read protein size in residues number
		fileLine = bufferedReader.readLine();
		aLineTokenizer = new StringTokenizer(fileLine," ");
		aLineToken     = (aLineTokenizer.nextToken()).trim();
		residuesNumber = (new Integer(aLineToken)).intValue();
	        myContacts     =  new int[residuesNumber][residuesNumber];
		numberContactsPerResidue = new int[residuesNumber];
		for(i=0;i<residuesNumber;i++)
		    {
			for(j=0;j<residuesNumber;j++)
			    {
				myContacts[i][j]=0;
			    }
			numberContactsPerResidue[i]=0;
		    }

		// read in the number of contacts
		fileLine = bufferedReader.readLine();
		aLineTokenizer = new StringTokenizer(fileLine," ");
		aLineToken     = (aLineTokenizer.nextToken()).trim();
		numberOfContacts = (new Integer(aLineToken)).intValue();
		// set protein id
		this.setName(fileName+" - total contacts:"+numberOfContacts);
		
		totalContactsNumber = 0;
		do
		    {
			// reads a pair of residues
			fileLine = bufferedReader.readLine();
			aLineTokenizer = new StringTokenizer(fileLine," ");
			aLineToken = (aLineTokenizer.nextToken()).trim();
			headResidue = (new Integer(aLineToken));
			aLineToken = (aLineTokenizer.nextToken()).trim();
			tailResidue = (new Integer(aLineToken));
			this.addContact(headResidue,tailResidue);
			
		    } while (bufferedReader.ready() == true);
		bufferedReader.close();
		if(totalContactsNumber!=numberOfContacts)
		    {
			System.out.println("WARNING - number of contacts declared in file "+fileName+" differs from contacts specified within the file.");
		    }
	    }
	catch (Throwable e)
	    {
		System.out.println("ERROR - exception "+e+" was generated in readInstance() in ContactMap.java");
	    }
	

     

    }



    public void createContactMapFromModelPDB(String fileName, double aThreshold)
    {
	/* Variable definitions to read the instance file.*/ 
	File              instanceFile;
	FileInputStream   fileInputStream;
	InputStreamReader inputStreamReader;
	BufferedReader    bufferedReader;
	String            fileLine;
	StringTokenizer   aLineTokenizer;
	String            aLineToken1;
	String            aLineToken2;
	String            aLineToken3;
	String            aLineToken4;
	String            aLineToken5;
	String            aLineToken6;
	double            aLineToken7;
	double            aLineToken8;
	double            aLineToken9;
	int               i;
	int               j;
	Integer           headResidue;
        Integer           tailResidue;
	int               numberOfContacts;
	boolean           found;
	double CAAtomsCoordinates[][];
	int numberCAs;
	double            x1,y1,z1,x2,y2,z2;


	CAAtomsCoordinates = new double[MAXCAATOMS][3];
	numberCAs = 0;
	threshold = aThreshold;
	try
	    {
		instanceFile      = new File(fileName);
		fileInputStream   = new FileInputStream(instanceFile);
		inputStreamReader = new InputStreamReader(fileInputStream);
		bufferedReader    = new BufferedReader(inputStreamReader); /* This is _actually_ from where we are gona read */
		found = false;
		// read protein size in residues number
		System.out.println("Reading pdb Model in"+fileName+" :");
		while(!found)
		    {

			fileLine = bufferedReader.readLine();
			System.out.println(fileLine);
			found = fileLine.indexOf("MOLECULE")!=-1;
			
		    }
		
		/* here we know we reach the juici part of the file */
		found = false;
		while((!found)&&(bufferedReader.ready()))
		    {
			aLineToken1 ="";
			aLineToken2 ="";
			aLineToken3 ="";
			aLineToken4 ="";
			aLineToken5 ="";
			aLineToken6 ="";
			aLineToken7 =0.0;
			aLineToken8 =0.0;
			aLineToken9 =0.0;

			fileLine = bufferedReader.readLine();
			if(fileLine.length()>0)
			    {
				aLineTokenizer = new StringTokenizer(fileLine," ");
				if(aLineTokenizer.hasMoreTokens())
				    aLineToken1     = (aLineTokenizer.nextToken()).trim();
				if(aLineTokenizer.hasMoreTokens())
				    aLineToken2     = (aLineTokenizer.nextToken()).trim();
				if(aLineTokenizer.hasMoreTokens())
				    aLineToken3     = (aLineTokenizer.nextToken()).trim();
				if(aLineTokenizer.hasMoreTokens())
				    aLineToken4     = (aLineTokenizer.nextToken()).trim();
				if(aLineTokenizer.hasMoreTokens())
				    aLineToken5     = (aLineTokenizer.nextToken()).trim();
				if(aLineTokenizer.hasMoreTokens())
				    aLineToken6     = (aLineTokenizer.nextToken()).trim();
				if(aLineTokenizer.hasMoreTokens())
				    aLineToken7     = new Double((aLineTokenizer.nextToken()).trim()).doubleValue();
				if(aLineTokenizer.hasMoreTokens())
				    aLineToken8     = new Double((aLineTokenizer.nextToken()).trim()).doubleValue();
				if(aLineTokenizer.hasMoreTokens())
				    aLineToken9     = new Double((aLineTokenizer.nextToken()).trim()).doubleValue();
			    
			CAAtomsCoordinates[numberCAs][0]=aLineToken7;
			CAAtomsCoordinates[numberCAs][1]=aLineToken8;
			CAAtomsCoordinates[numberCAs][2]=aLineToken9;
			numberCAs++;
			    }
			found = fileLine.indexOf("END")!=-1;
		    }
		bufferedReader.close();
			
	    }
	catch (Throwable e)
	    {
		System.out.println("ERROR - exception "+e+" was generated in createContactMapFromModelPDB() in ContactMap.java");
	    }
	
	/* we initialize this contact map */
	residuesNumber = numberCAs;
	myContacts     =  new int[residuesNumber][residuesNumber];
	numberContactsPerResidue = new int[residuesNumber];
	for(i=0;i<residuesNumber;i++)
	    {
		for(j=0;j<residuesNumber;j++)
		    {
			myContacts[i][j]=0;
		    }
		numberContactsPerResidue[i]=0;
	    }
	numberOfContacts = 0;
	for(i=0;i<numberCAs;i++)
	    {
		for(j=i+2;j<numberCAs;j++)
		    {
			x1 = CAAtomsCoordinates[i][0];
			y1 = CAAtomsCoordinates[i][1];
			z1 = CAAtomsCoordinates[i][2];
			x2 = CAAtomsCoordinates[j][0];
			y2 = CAAtomsCoordinates[j][1];
			z2 = CAAtomsCoordinates[j][2];
			if(distance3D(x1,y1,z1,x2,y2,z2)<=threshold)
			    {
				addContact(new Integer(i),new Integer(j));
				numberOfContacts++;
			    }
		    }
	    }
	// set protein id
	this.setName(fileName+" - total contacts:"+numberOfContacts);
		

    }



    public void saveContactMap(String fileName)
    {

	
	int i,j;

      /* we open the file to write */
      
      try
	{
	  fileToWrite   = new PrintWriter(new BufferedWriter(new FileWriter(fileName)));
	  writeToFile(residuesNumber     +"       # Number of Residues\n");
	  writeToFile(totalContactsNumber+"       # Number of Contacts at "+threshold+" Angstroms \n");
	  for(i=0;i<residuesNumber;i++)
	      {
		  for (j=0;j<numberContactsPerResidue[i];j++)
		      {
			  writeToFile(i+"    "+myContacts[i][j]+"\n");
		      }
	      }
	} 
      catch (IOException e)
	{
	  System.out.println("ERROR - couldn't open file for writting "+fileName+ " in ContactMap.java");
	  //	  System.exit(1);
	}

      closeFile();
    }

    protected void writeToFile(String aString)
    {
	fileToWrite.print(aString);
	fileToWrite.flush();
    }
    protected void closeFile()
    {
	fileToWrite.close();
    }
    
    protected double distance3D(double x1, double y1, double z1, double x2, double y2, double z2)
    {
	double d=0.0;

	d = Math.sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2) );
	//	if (d<20) System.out.println(d);

	return (d); 

    }

}
