

import java.io.*;
import java.awt.*;
import java.util.*;

public class ContactMapDisplay extends Canvas
    {
	static Color alignmentColor  = Color.yellow;
	static Color residueColor    = Color.red;
	static Color lineColor       = Color.blue;
	static Color backgroundColor = Color.black;
	protected int radius;
	protected Frame myFrame;
	protected int xPos; /* where to print */
	protected int yPos; /* where to print */
	protected int iniX; /* where to place the canvas within myFrame */
	protected int iniY; /* where to place the canvas within myFrame */
	protected ContactMap cm1;
	protected ContactMap cm2;
	protected int xcm1,ycm1;
	protected int xcm2,ycm2;
	protected double xScale;


	public ContactMapDisplay (Frame aFrame, int xIni, int yIni, int w, int h, Color aC, ContactMap acm1, ContactMap acm2, int axcm1,int aycm1, int axcm2, int aycm2)
	{
	    
	    myFrame= aFrame;
	    this.setSize(w,h);
	    this.setLocation(xIni,yIni);
	    backgroundColor = aC;
	    this.setBackground(backgroundColor);
	    myFrame.add(this);
	    iniX = xIni;
	    iniY = yIni;
	    cm1  = acm1;
	    cm2  = acm2;
	    xcm1 = axcm1;
	    xcm2 = axcm2;
	    ycm1 = aycm1;
	    ycm2 = aycm2;
	    if(cm1.getResiduesNumber()>cm2.getResiduesNumber())
		{
		    xScale = (w-2*xPos)/(cm1.getResiduesNumber()+1);
		}
	    else
		{
		    xScale = (w-2*xPos)/(cm2.getResiduesNumber()+1);
		}
	    radius=(int)xScale/4;
	}

	
	
	public void display( int posX, int posY)
	{

	    xPos = posX;
	    yPos = posY;
	    Graphics g= this.getGraphics();
	    update(g);
	    
	}

	public void paint(Graphics g)
	{

	    paintBackground(g);
	    paintContactMaps(g);

	}

	public void update(Graphics g)
	{
	    paintBackground(g);
	    //	    System.out.println("In update()");
	    paintContactMaps(g);
	}

	public void close()
	{
	    myFrame.setVisible(false);
	    myFrame.dispose();

	}

	public void deleteMyCanvasBackground()
	{
	    paintBackground(this.getGraphics());
	}


	public void plotAlignment(int alignment[], double fitness)
	{/* the alignment vector length is always the same as that of the number of residues in the first contact map as we keep the shortest contact map in
	    aCM1. */

	    int i;
	    int l=cm1.getResiduesNumber();
	    int x1,y1,x2,y2;

	    x1=x2=y1=y2=-1000;
            update(this.getGraphics());
            if (alignment==null) return;
	  
	    for(i=0;i<l;i++)
		{
		    if(alignment[i]>-1)
			{
			    x1 =(int)(i*xScale+xPos+xcm1);
			    y1 = yPos+ycm1;
			    x2 =(int)(alignment[i]*xScale+xPos+xcm2);
			    y2 = yPos+ycm2;
			    myDrawLine(this.getGraphics(),x1,y1,x2,y2,alignmentColor);

			}
		}
	    drawText(this.getGraphics(),"Fitness: "+fitness,(int)( cm2.getResiduesNumber()*xScale+xPos+xcm2+25 ), yPos+ycm1,Color.white);

	}

	protected void paintBackground(Graphics g)
	{
	    g.setColor(backgroundColor);
	    //  System.out.println("X="+this.getLocation().x+", Y="+this.getLocation().y);
	    g.fillRect(this.getLocation().x,this.getLocation().y,this.getLocation().x+this.getSize().width,this.getLocation().y+this.getSize().height);
	}

	protected void paintContactMaps(Graphics g)
	{
	    int i,x1,y1;

	    paintVertices(cm1,xcm1,ycm1,g);
	    paintEdges(cm1,xcm1,ycm1,g,0);
	    paintVertices(cm2,xcm2,ycm2,g);
	    paintEdges(cm2,xcm2,ycm2,g,1);

	 
	}

	private void paintVertices(ContactMap aCM, int xcm,int ycm, Graphics g)
	{
	    int numRes;
	    int i;
	    int x;
	    int y;

	    numRes = aCM.getResiduesNumber();
	    for(i=0;i<numRes;i++)
		{
		    x =(int)(i*xScale+xPos+xcm);
		    y = yPos+ycm;
		    plotResidue(g,x,y,residueColor);
		}
	    
	}


	private void paintEdges(ContactMap aCM, int xcm,int ycm, Graphics g, int topDown)
	{
	    int resNum;
	    int contactsPerRes;
	    int i,j,f;
	    double x,y,arcRadius;
	    double h,w;

	   
	    resNum=aCM.getResiduesNumber();
	    for(i=0;i<resNum;i++)
		{
		    contactsPerRes=aCM.getContactsNumber(i);
		    for(f=0;f<contactsPerRes;f++)
			{
			   

			    w= (aCM.getContactResidueId(i,f)*xScale+xPos+xcm)-(i*xScale+xPos+xcm);
			    h=w;
			    y=yPos+ycm-(h/2); /* the same height where the residue lies */
			    x=i*xScale+xPos+xcm;
			    if(topDown==0)
				{
				    myDrawArc(g,(int)Math.round(x),(int)Math.round(y),(int)Math.round(w),(int)Math.round(h),0,180,lineColor);   
				}
			    else
				{
				    myDrawArc(g,(int)Math.round(x),(int)Math.round(y),(int)Math.round(w),(int)Math.round(h),0,-180,lineColor);   
				}

			}



		}

	}

	private void drawText(Graphics aG, String s, int x, int y, Color aC)
	{
	    Color tmp;
	    
	    
	    tmp = aG.getColor();
	    aG.setColor(aC);
	    aG.drawString(s,x,y);
	    aG.setColor(tmp);
	}

	private void plotResidue(Graphics aG,double x,double y, Color aC)
	{
	    Color tmp;
	    
	    
	    tmp = aG.getColor();
	    aG.setColor(aC);
	    aG.fillOval((int)x,(int)y,radius,radius);
	    aG.setColor(tmp);
	}
	
	private void myDrawLine(Graphics aG,double x1, double y1, double x2, double y2, Color aC)
	{
	    Color tmp;
	    
	    
	    tmp = aG.getColor();
	    aG.setColor(aC);
	    aG.drawLine((int)x1,(int)y1,(int)x2,(int)y2);
	    aG.setColor(tmp);
	}
	
	

	private void myDrawArc(Graphics aG,int x, int y, int w, int h,int sa, int ea, Color aC)
	{
	    Color tmp;
	    
	    
	    tmp = aG.getColor();
	    aG.setColor(aC);
	    aG.drawArc(x,y,w,h,sa,ea);
	    aG.setColor(tmp);
	}
	
	





}
