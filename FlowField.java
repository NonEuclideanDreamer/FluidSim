import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Random;

import javax.imageio.ImageIO;

public class FlowField 
{
	static Random rand=new Random();

	public static boolean fourier=false,//Couldn't get Fourier to work, might change in the future
						  viscosity=false, //viscosity only implemented for vist<=1/delt
						  obs=true; //Are there obstacles with boundary conditions present?
	public static int[]scale= {1080,1080}; // Dimensions of the screen= Toroidal world
	public static String name="testfour", //name of files
			type="png";					//format type of saved picture
	public static double delt=0.1,delx=1, //discrete steps in time&space
			visc=2,						//visvosity of medium. Must be <1/delt
			sqrt2=Math.sqrt(2); 
	public static int qterm=0,			//needed for the relaxation scheme when solvin Poisson
			x=400, r=5;				//For the Gust force
	public static int maxit=10000;	//breakoff of the simulation
	public static double[][]normal;
	public static DecimalFormat df=new DecimalFormat("0000");
	public static PotentialField[] q= {new PotentialField(), new PotentialField()};

	//FlowFiel consisting of velocity field, color field and current time
	public ForceField[] velocity,
				color;	
	public int time;	
	
	//Constructor, initialize color field, there could be an initial velocity field too
	public FlowField() 
	{
		double[] r= {256,0,0},y= {256,256,0},g= {0,256,0},c= {0,256,256},b= {0,0,256},m= {256,0,256},w= {256,256,256},s= {0,0,0},p= {128,0,256},o={256,128,0},
				l= {128,256,0}, mint=  {0,256,128},a= {0,85,170},wine= {170,0,85},gold= {255,215,0}, silver= {192,192,192};
		time=0;
		velocity=new ForceField[] {new ForceField(2),new ForceField(2)};
		
		//Choose an initial condition for the color field
		ForceField initcol=//ForceField.stripes(r,y,o,g,c,mint,b,m,p});
							//ForceField.rgbw();
							//ForceField.ryb();
							ForceField.multistripes(new double[][] {s,s,s,o,s,s,o,o,s,o,o,o,p,o,o,o,p,p,o,o,p,p,p,o,p,p,p,s,p,p,s,s,p,s,s,s});
		color=new ForceField[] {initcol,new ForceField()};////colors(wine,offset )),ForceField.stripes(colors(wine, offset)) };//
	}
	
	
	public static void main(String[] args)
	{
		FlowField euler=new FlowField();
		BufferedImage canvas=new BufferedImage(scale[0],scale[1],BufferedImage.TYPE_3BYTE_BGR);
		
		//initialize Obstacle(array of pixel coordinates, must be ordered&accompanied by array of normal vectors
		
		int[][]obstacle=
			{};//only if obs==false
		//ObstacleCircle(x,270,40,30);
		//ObstacleDiag();
		//ObstacleHor();
		//ObstacleBows(460,480,2);//
		
		//needed for gust-force
		x=rand.nextInt(scale[0]);
		r=(int) Math.abs(rand.nextGaussian(54,18))+1;
		//	System.out.println("x="+x+", r="+r+", strength="+(1600/Math.sqrt(r)));
	
		while(euler.time<maxit)
		{
			euler.colordraw(canvas);
			
			for(int i=0;i<5;i++)
			{

				euler.update(obstacle);//
			
			}
			euler.time++;
			System.out.println(euler.time);
		}
	}

	//********************************************************************************************
	// Making a time step: adding the force, applying pressure, advecting, and viscosity if chosen
	//********************************************************************************************
	public void update(int[][]obstacle)
	{
		
		/*	if(time%200==101)//for a gust force
			{
				x=rand.nextInt(scale[0]);
				r=(int) Math.abs(rand.nextGaussian(54,18))+1;
				System.out.println("x="+x+", r="+r+", strength="+(800/Math.sqrt(r)));
			}*/
		int next=(time+1)%2, now=time%2;
	
		//choose a force pattern, vary parameters for random gusts uncomment above if-statement
		velocity[now].
					stir(time*0.01, 20, 1, 2);
					//diamondstir(time*0.01,20,1);
					//diamond stir
					//comb(time*0.01,100,4,20,1);
					//gust(time/100.0, 30, 10,1);
					//levelcomb(time*0.02,scale[1]/8,6);

		int k=0;
		
		if(fourier) //not working properly
		{
		
			System.out.println("Fourier");
			velocity[now]=velocity[next].Fourier(1024);
		
			//viscuous term 
			if(viscosity)
			{
				System.out.println("visc");
				velocity[now].divide(velocity[now],visc*delt);	
				System.out.println("pressure");
			}
			
			//project//pressure
			velocity[next].project(velocity[now]);
			
			System.out.println("back");
			velocity[now]=velocity[next].invFourier(1024);
		}
		
		else
		{
			//Pressureterm; Making the velocity gradientfree
			if(obs)		qterm=relax(velocity[now],obstacle);
			else 	qterm=relax(velocity[now]);
			velocity[now].subtract(q[qterm].gradient());
			
			if(viscosity)
			{
				for(int i=0;i<scale[0];i++)
					for(int j=0;j<scale[1];j++)
					{  
						double[]del2=velocity[now].del2(i,j);
						for(int k1=0;k1<2;k1++)
							velocity[next].vectorfield[i][j][k1]=velocity[now].vectorfield[i][j][k1]-visc*del2[k1];
					}
				velocity[now].vectorfield=velocity[next].vectorfield.clone();
			}
		}
		//advection of u-field
		for(int i=0;i<scale[0];i++)
			for(int j=0;j<scale[1];j++)
				{  
					int ip=(i+1)%scale[0],im=(i+scale[0]-1)%scale[0],jp=(j+1)%scale[1],jm=(j+scale[1]-1)%scale[1],id=200,jd=200;//torus	
					double[]loc= {((i-velocity[now].vectorfield[i][j][0]*delt/delx)%scale[0]+scale[0])%scale[0],((j-velocity[now].vectorfield[i][j][1]*delt/delx)%scale[1]+scale[1])%scale[1]};
					velocity[next].vectorfield[i][j]=average(velocity[now].vectorfield,loc);
					color[next].vectorfield[i][j]=average(color[now].vectorfield,loc);
				}
	}
	
	//****************
	// Drawing Methods
	//****************
	
	//drawing the color field
	private void colordraw(BufferedImage canvas )
	{
		int now=time%2;
		for(int i=0;i<scale[0];i++)
			for(int j=0;j<scale[1];j++)
			{
				
				double	red=color[now].vectorfield[i][j][0],
					green=color[now].vectorfield[i][j][1],
					blue=color[now].vectorfield[i][j][2];
				int	x=sRGB(red,green,blue);
					//new Color(red,green,blue).getRGB();
				canvas.setRGB(i, j, x);
			}

		File outputfile = new File(name+df.format(time)+"."+type);
		try 
		{  
			ImageIO.write(canvas, type, outputfile);
		} 
		catch (IOException e) 		
		{
			System.out.println("IOException");
			e.printStackTrace();
		}
		
	}
	
	//drawing the color field and the obstacle in black
	private void colordraw(BufferedImage canvas , int[][]obstacle)
	{
		int now=time%2;
		for(int i=0;i<scale[0];i++)
			for(int j=0;j<scale[1];j++)
			{
				
				double	red=color[now].vectorfield[i][j][0],
					green=color[now].vectorfield[i][j][1],
					blue=color[now].vectorfield[i][j][2];
				int	x=sRGB(red,green,blue);//new Color(red,green,blue).getRGB();
				canvas.setRGB(i, j, x);
			}
		
		for(int i=0;i<obstacle.length;i++) canvas.setRGB(obstacle[i][0],obstacle[i][1],Color.black.getRGB());
		
		File outputfile = new File(name+df.format(time/10)+"."+type);
		try 
		{  
			ImageIO.write(canvas, type, outputfile);
		} 
		catch (IOException e) 		
		{
			System.out.println("IOException");
			e.printStackTrace();
		}
		
	}
	
	//color by velocity field(I only used it before adding the color field...)
	private void draw(BufferedImage canvas) {
		int now=time%2;
		for(int i=0;i<scale[0];i++)
			for(int j=0;j<scale[1];j++)
			{
				
				int	red=0,
					blue=(int)Math.max(0, Math.min(255, (velocity[now].vectorfield[i][j][0]*64+128))),
					green=(int)Math.max(0, Math.min(255, velocity[now].vectorfield[i][j][1]*64+128)),
					x=new Color(red,green,blue).getRGB();
				canvas.setRGB(i, j, x);
			}
		
		File outputfile = new File(name+df.format(time)+"."+type);
		try 
		{  
			ImageIO.write(canvas, type, outputfile);
		} 
		catch (IOException e) 		
		{
			System.out.println("IOException");
			e.printStackTrace();
		}
		
	}

	//converting the color to srgb, uncomment stuff to actually do that...
	private static int sRGB(double r, double g, double b) 
	{
		double[] y=new double[3];
		 y[0]=r;//(3.2406*r-1.5372*g-.4986*b);
		y[1]=g;//(-.9689*r+1.8758*g+.0415*b);
		y[2]=b;//(.0557*r-.204*g+1.057*b);
		
		int[]x=new int[3];
		/*for(int i=0;i<3;i++)
		{if(y[i]>.8015)x[i]=(int)(Math.pow(y[i],1/2.4)*26.795-14.08);
		else x[i]=(int) (12.92*y[i]);}
	*/
		for(int i=0;i<3;i++)
		x[i]=(int) Math.min(255, Math.max(0, y[i]));
		
		return new Color(x[0],x[1],x[2]).getRGB();
	}
	
	//color list with fixed color fallback and rainboish in between
	public static double[][]colors(double[]fallback,double d)
	{
		int base=2;
		double[][]out=new double[12][3];
		for(int i=1;i<12;i+=2)
			out[i]=fallback.clone();
	
		out[0]=new double[] {(1-2*d)*fallback[0]+d*256,256-fallback[1],fallback[2]};
		out[2]=new double[] {256-fallback[0],(2*d-1)*fallback[1]+(1-d)*256,fallback[2]};
		out[8]=new double[] {fallback[0],(1-2*d)*fallback[1]+d*256,256-fallback[2]};
		out[10]=new double[] {fallback[0],256-fallback[1],(2*d-1)*fallback[2]+(1-d)*256};
		out[4]=new double[] {256-fallback[0],fallback[1],(1-2*d)*fallback[2]+d*256};
		out[6]=new double[] {(2*d-1)*fallback[0]+(1-d)*256,fallback[1],256-fallback[2]};
		return out;
	}
	
	//***************************************************
	// Relaxation scheme for solving the Poisson Equation
	//***************************************************
	public static int relax (ForceField f)
	{
		q[0].set(0);
		PotentialField div=f.divergence(delx);
		int l=0, now=qterm,next=1-now;
		double max=10,bound=0.08;
		while(l<100&&max>bound)//600
		{
			max=q[next].poisson(q[now],div,delx);
			l++;
			now=1-now;
			next=1-next;
		}
		return now;
	}
	public static int relax (ForceField f, int[][]obstacle)
	{
		q[0].set(0);
		PotentialField div=f.divergence(delx);
		int l=0, now=0,next=1-now;
		double max=10,bound=0.05;
		//System.out.println("i="+f.vectorfield[0][0].length);
		while(l<1000)//500l<400&&max>bound&&
		{
			max=q[next].poisson(q[now],div,delx,obstacle, normal	,f);
			l++;
			now=1-now;
			next=1-next;
		}
	//	System.out.println(q[now].potential[109][50]);
		return now;
	}
	
	//********************************************************************************
	// Taking the weighted average of the field value of the for cells surrounding loc
	//********************************************************************************
	public static double[] average(double[][][] vectorfield, double[] loc) 
	{
		int x=(int)loc[0], y=(int)loc[1], xp=(x+1)%scale[0], yp=(y+1)%scale[1];
		double dx=loc[0]%1, dy=loc[1]%1;
		double[]out=new double[vectorfield[x][y].length];
		for(int i=0;i<vectorfield[x][y].length;i++)
		{
			out[i]=dx*(dy*vectorfield[xp][yp][i]+(1-dy)*vectorfield[xp][y][i])+(1-dx)*(dy*vectorfield[x][yp][i]+(1-dy)*vectorfield[x][y][i]);
		}
		return out;
	}
	public static double average(double[][] field, double[] loc) 
	{
		int x=(int)loc[0], y=(int)loc[1], xp=(x+1)%scale[0], yp=(y+1)%scale[1];
		double dx=loc[0]%1, dy=loc[1]%1;
		double out=dx*(dy*field[xp][yp]+(1-dy)*field[xp][y])+(1-dx)*(dy*field[x][yp]+(1-dy)*field[x][y]);
		return out;
	}
	
	
	//****************************************************
	// printing field values to the terminal for debugging
	//****************************************************
	static void print(double[][][]tensor)
	{
		System.out.print("{");
		for(int i=0;i<tensor.length;i++)
		{
			System.out.print("{");
			for(int j=0;j<tensor[i].length;j++)
			{
				System.out.print(" {");
				for(int k=0;k<tensor[i][j].length;k++)
					System.out.print(Math.round(tensor[i][j][k]*100)*.01+", ");
				System.out.print("},");
			}
			System.out.println("},");
		}
		System.out.println("}");
	}
	static void print(double[][] m) 
	{
		System.out.print("{");
		for(int i=0;i<m.length;i++)
		{
			print(m[i]);
		}
		System.out.println("}");
	}
	static void print(double[] m) 
	{
		System.out.print("{");
		for(int i=0;i<m.length;i++)
		{
			System.out.print(Math.round(m[i]*100)*0.01+", ");	
		}
		System.out.println("},");
	}
	
	//*********************************************************
	// Obstacles and their normal arrays for boundary behavior
	//*********************************************************

	//vertical boun on left & right side
	public static int[][] ObstacleVert()
	{
		int[][]out= new int[2*scale[1]][2];
		normal=new double[out.length][2];
		for(int i=0;i<scale[1];i++)
		{
				out[i][0]=0;
				out[i][1]=i;
				normal[i][0]=2;
				normal[i][1]=0;
				
				out[i+scale[1]][0]=scale[0]-1;
				out[i+scale[1]][1]=i;
				normal[i+scale[1]][0]=-2;
				normal[i+scale[1]][1]=0;

		}
		return out;
	}
	
	//diagonal bound
	public static int[][] ObstacleDiag()
	{
		int[][]out= new int[3*scale[1]][2];
		normal=new double[out.length][2];
		for(int i=0;i<scale[1]-1;i++)
		{
				out[i][0]=i;
				out[i][1]=i+1;
				normal[i][0]=-(2);
				normal[i][1]=(2);
				
				out[i+scale[1]][0]=i;
				out[i+scale[1]][1]=i;
				normal[i+scale[1]][0]=(2);
				normal[i+scale[1]][1]=-(2);		
		}
		return out;
	}
	
	//horizontal bound on top & bottom
	public static int[][] ObstacleHor()
	{
		int[][]out= new int[2*scale[0]][2];
		normal=new double[out.length][2];
		for(int i=0;i<scale[0];i++)
		{
				out[i][1]=0;
				out[i][0]=i;
				normal[i][1]=2;
				normal[i][0]=0;
				
				out[i+scale[0]][1]=scale[1]-1;
				out[i+scale[0]][0]=i;
				normal[i+scale[0]][1]=-2;
				normal[i+scale[0]][0]=0;

		}
		return out;
	}
	
	// circular obstacle
	public static int[][] ObstacleCircle(int x, int y, int r)
	{
		int[][]out= new int[(int)(Math.PI*r*r)][2];
		normal=new double[out.length][2];
		int k=0;
		for(int i=-r;i<r;i++)
		{
			double s=Math.sqrt(r*r-i*i);
			for(int j=(int) -s;j<s;j++)
			{
				out[k][0]=i+x;
				out[k][1]=j+y;
				normal[k][0]=i;
				normal[k][1]=j;
				double norm=Math.sqrt(normal[k][0]*normal[k][0]+normal[k][1]*normal[k][1]);
				if(norm!=0) {
				normal[k][0]*=2/norm;
				normal[k][1]*=2/norm;}
				k++;
			}
		}
		return out;
	}
	
	// ringshaped obstacle, r2<r
	public static int[][] ObstacleCircle(int x, int y, int r,double r2)
	{
		int[][]out= new int[(int)(Math.PI*r*r)][2];
		normal=new double[out.length][2];
		int k=0;
		double rmid=(r+r2)/2;
		for(int i=-r;i<r;i++)
		{
			double s=Math.sqrt(r*r-i*i);
			for(int j=(int) -s;j<s;j++)
			{
				out[k][0]=i+x;
				out[k][1]=j+y;
				normal[k][0]=i;
				normal[k][1]=j;
				double norm=Math.sqrt(normal[k][0]*normal[k][0]+normal[k][1]*normal[k][1]);
				if(norm>r2) 
				{
					if(norm<rmid)norm*=-1;
				normal[k][0]*=2/norm;
				normal[k][1]*=2/norm;
				k++;
				}
			}
		}
		return out;
	}
	
	
	public static int[][] nestedSemiCircles(int n)
	{
		int[][] out =new int[10000][2];
		normal=new double[out.length][2];

		int k=0;
		for(int i=0;i<scale[0];i++)
			for(int j=0;j<scale[1];j++)
			{
				int x=scale[0]-i, y=scale[1]-j;
				double r=Math.sqrt(i*i+j*j), rlog=Math.log(r)%2;
				if(rlog<.2)
				{
					out[k][0]=i;
					out[k][1]=j;
					normal[k][0]=-x/r*2;
					normal[k][1]=-y/r*2;
					k++;
				}
				else if(rlog>1.8)
				{
					out[k][0]=i;
					out[k][1]=j;
					normal[k][0]=x/r*2;
					normal[k][1]=y/r*2;
					k++;
				}
			}return out;
	}
	
	//r2>r1, 2*r2*n<=scale[0]
	public static int[][] ObstacleBows(double r1,double r2, int n )
	{
		int[][]out= new int[(int)(n*Math.PI*(r1+r2)*(r2-r1))*2][2];
		normal=new double[out.length][2];
		int k=0;
		double step=1.0*scale[0]/n, gap=step-2*r2, midr=Math.sqrt((r1+r2)*(r1+r2)/4);
		for(int i=0;i<scale[0];i++)
		{
			double x=i%step;
			if (x>gap/2&&x<step-gap/2)
			{
				double xcent=(int)(i/step)*step+step/2;
				double bound=scale[1]/2,s2=r1*r1-Math.pow(xcent-i, 2);
				if(s2>0)bound-=Math.sqrt(s2);
				for( int j=scale[1]/2-(int)(Math.sqrt(r2*r2-Math.pow(xcent-i, 2)));j<bound;j++)
				{
					out[k][0]=i;
					out[k][1]=j;
					normal[k][0]=i-xcent;
					normal[k][1]=j-scale[1]/2;
					double norm=Math.sqrt(normal[k][0]*normal[k][0]+normal[k][1]*normal[k][1]);
					if(norm!=0)
					if(norm<midr)norm*=-1;
					{
						normal[k][0]*=2/norm;
						normal[k][1]*=2/norm;
					}
					System.out.println("("+i+","+j+")");
					k++;
				}
			}
			if(x<r2||x>step-r2)
			{
				double xcent=(int)((i+step/2)/step)*step;
				double bound=scale[1]/2,s2=r1*r1-Math.pow(xcent-i, 2);
				if(s2>0)bound+=Math.sqrt(s2);
				for( int j=(int) bound+1;j<scale[1]/2+(int)(Math.sqrt(r2*r2-Math.pow(xcent-i, 2)));j++)
				{
					out[k][0]=i;
					out[k][1]=j;
					normal[k][0]=i-xcent;
					normal[k][1]=j-scale[1]/2;
					double norm=Math.sqrt(normal[k][0]*normal[k][0]+normal[k][1]*normal[k][1]);
					if(norm!=0)
					if(norm<midr)norm*=-1;
					{
						normal[k][0]*=2/norm;
						normal[k][1]*=2/norm;
					}
					System.out.println("("+i+","+j+")");

					k++;
				}
			}
		}
		return out;
	}
}
