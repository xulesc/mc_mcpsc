



import java.util.*;
import java.awt.*;



public class BuildContactMapFromPDB {

  public static void  main(String []argv) throws CloneNotSupportedException
   {
     

     /* Problem and Individual Variables */
       ContactMap aCM;
       ContactMapDisplay aCMD;
       Frame aF;
       int i;
       boolean display= (Boolean.valueOf(argv[2])).booleanValue();;

       aCM = new ContactMap();
       aCM.createContactMapFromModelPDB(argv[0],new Double(argv[1]).doubleValue());
       if(display)
	   {
	       aF = new Frame("Displaying Contact Map");
	       aF.setSize(670,520);
	       aF.show();
	       
	       aCMD = new ContactMapDisplay(aF,0,0,670,520,Color.black,aCM,aCM,10,240,10,280);
	       aCMD.display(20,20);
	   }
       aCM.saveContactMap(argv[0]+".cm");


     for(i=0;display;i+=0);



   }
}
