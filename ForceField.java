//******************************************************************************
// FluidSim/ForceField.java
// author: Non-Euclidean Dreamer
// gives methods to work with vectorfields, I also use them for the color field
//******************************************************************************

import java.util.Random;

public class ForceField 
{
	static int dim=2;//I only implemented 2d-stuff, but that might change in the future
	public static int[] scale=FlowField.scale;
	
	public double[][][]vectorfield;
	
	//***********
	//Contructors
	//***********
	
	public ForceField()
	{
		vectorfield=new double[scale[0]][scale[1]][dim];
	}
	public ForceField(int d)
	{
		vectorfield=new double[scale[0]][scale[1]][d];
	}
	
	public ForceField(double x, double y) //uniform field
	{
		vectorfield=new double[scale[0]][scale[1]][dim];	
		
		for(int i=0;i<scale[0];i++)
			for(int j=0;j<scale[1];j++)
				vectorfield[i][j]=new double[] {x,y};
	}

	public ForceField(double[][][] vf) 
	{
		vectorfield=vf;
	}

	//chaotic an thus highly unstable initial state
	public static ForceField random(double bound) 
	{
		
		
		Random rand=new Random();
		double[][][]vectorfield=new double[scale[0]][scale[1]][dim];	
		for(int i=0;i<scale[0];i++)
			for(int j=0;j<scale[1];j++)
				for(int k=0;k<dim;k++) vectorfield[i][j][k]=rand.nextGaussian(0, bound);
		
		return new ForceField(vectorfield);
	}
	
	//***************************************
	// Initial Conditions for the Color Field
	//***************************************
	public static ForceField rgbw() 
	{
		ForceField out=new ForceField(3);
		for(int i=0;i<scale[0]/2;i++)
		{	
			for(int j=0;j<scale[1]/2;j++) {out.vectorfield[i][j][0]=255.999;out.vectorfield[i][j][1]=255.999;}
			for(int j=scale[1]/2;j<scale[1];j++) {out.vectorfield[i][j][2]=255.999;}
		}
		for(int i=scale[0]/2;i<scale[0];i++)
		{	
		
			for(int j=scale[1]/2;j<scale[1];j++)
				{
					out.vectorfield[i][j][0]=255.999;
					
				}
		}		
		
		return out;
	}

	
	public static ForceField ryb()
	{
		ForceField out=new ForceField(3);
		for(int i=0;i<scale[0]/3;i++)
		{	
			for(int j=0;j<scale[1];j++) 
			{
				out.vectorfield[i][j][1]=255.999;//	out.vectorfield[i][j][0]=255.999;
				out.vectorfield[i+scale[0]/3][j][0]=255.999;
				out.vectorfield[i+scale[0]*2/3][j][2]=255.999;

			}

		}
		return out;
	}
	public static ForceField stripes(double[][]cl)
	{
		int n=cl.length,width=scale[0]/n;
		ForceField out=new ForceField(3);
		for(int k=0;k<n;k++)
			for(int l=0;l<3;l++)
		for(int i=0;i<scale[0]/n;i++)
		{	
			for(int j=0;j<scale[1];j++) 
			{
				out.vectorfield[i+k*width][j][l]=cl[k][l];
			}
		}
		return out;
	}
	public static ForceField multistripes(double[][]cl)
	{
		int n=cl.length,width=scale[0]/n;double offset;
		ForceField out=new ForceField(3);
		Random rand=new Random();
		for(int k=0;k<n;k++)
			for(int l=0;l<3;l++)
			{ offset=0;
			for(int j=0;j<scale[1];j++) 
			{
				offset+=3.6*rand.nextDouble()-1.8;
		for(int i=0;i<scale[0]/n;i++)
		{	
			
				out.vectorfield[i+k*width][j][l]=cl[k][l]+offset;
			}
		}
			}
		return out;
	}
	

	//***********************************************
	// Methods working with the vectorfield structure
	//***********************************************
	
	//Divergence of the field, where delx is the cellsize
	public PotentialField divergence(double delx) 
	{
		PotentialField out=new PotentialField();
		for(int i=0;i<scale[0];i++)
			for(int j=0;j<scale[1];j++)
			{
				int im=(i+scale[0]-1)%scale[0],jm=(j+scale[1]-1)%scale[1],ip=(i+1)%scale[0],jp=(j+1)%scale[1];//neigbours in toroidal universe
				out.potential[i][j]=0.5*(vectorfield[ip][j][0]-vectorfield[im][j][0]+vectorfield[i][jp][1]-vectorfield[i][jm][1])/delx;		
			}
		//	System.out.println("div=");FlowField.print(out.potential);
		return out;
	}
	
	//Gives the second derivative at (x|y)
	public double[] del2(int x, int y)
	{
		int xp=(x+1)%scale[0], xm=(x+scale[0]-1)%scale[0],yp=(y+1)%scale[1], ym=(y+scale[1]-1)%scale[1];
		double[] out=new double[2];
		out[0]=(vectorfield[xp][y][0]-4*vectorfield[x][y][0]+vectorfield[xm][y][0]+vectorfield[x][yp][0]+vectorfield[x][ym][0])/4;
		out[1]=(vectorfield[xp][y][1]-4*vectorfield[x][y][1]+vectorfield[xm][y][1]+vectorfield[x][yp][1]+vectorfield[x][ym][1])/4;
		return out;
	}
	
	//*******************************
	// Modifying the Vector Field
	//*******************************
	
	//action of a force stirring in a sinusoidal curve, mult decides the shape of the curve, for integer it's closed
	public void stir(double t,int r,double strength, double mult)
	{
		
		double f0=-strength*Math.sin(t),
				f1=strength*mult*Math.cos(mult*t)*scale[1]/scale[0];
		for(int i=-r;i<r;i++) 
		{
			double s=Math.sqrt(r*r-i*i);
			for(int j=(int) -s;j<s;j++)
			{
				double factor=0.5*(1-Math.cos(Math.PI*Math.sqrt(s*s-j*j)/r));
				vectorfield[(int) (i+scale[0]/2+(scale[0]/2-r)*Math.cos(t))][(int) (j+scale[1]/2+(scale[1]/2-r)*Math.sin(mult*t))][0]+=f0*factor;
				vectorfield[(int) (i+scale[0]/2+(scale[0]/2-r)*Math.cos(t))][(int) (j+scale[1]/2+(scale[1]/2-r)*Math.sin(mult*t))][1]+=f1*factor;
			}
		}
	}
	
	//action of force stirrin in a diamond shape
	public void diamondstir(double t,int r,double strength)
	{
		double f0=strength*Math.cos(t),
				f1=-strength*Math.sin(t),x,y;
		if((t/Math.PI*2)%4<1) {x=scale[0]-r;y=scale[1]-r;f0*=-1;}
		else if((t/Math.PI*2)%4<2) {x=r;y=scale[1]-r;f1*=-1;}
		else if((t/Math.PI*2)%4<3) {x=r;y=r;f0*=-1;}
		else  {x=scale[0]-r;y=r;f1*=-1;}
		for(int i=-r;i<r;i++) 
		{
			double s=Math.sqrt(r*r-i*i);
			for(int j=(int) -s;j<s;j++)
			{
				double factor=0.5*(1-Math.cos(Math.PI*Math.sqrt(s*s-j*j)/r));
				vectorfield[(int) (i+x-(scale[0]/2-r)*Math.cos(t))][(int) (j+y-(scale[1]/2-r)*Math.sin(t))][0]+=f0*factor;
				vectorfield[(int) (i+x-(scale[0]/2-r)*Math.cos(t))][(int) (j+y-(scale[1]/2-r)*Math.sin(t))][1]+=f1*factor;
			}
		}
	}

	//combing force acting on vectorfield
	public void comb(double t, int w,int n, double r, double str) 
	{
		int x;
		//System.out.println("wh");
		int c=0, start=(scale[1]-1-(n-1)*w)/2;
		for(int k=start;k<scale[1]-start;k+=w)
		{
			double u=t+k*2.0/scale[1];
			if(u%2<1)
			x= (int) (r+(u%1)*(scale[0]-2*r));
		else	x=(int) ((scale[0]-2*r)*(1-(u%1))+r);
			for(int i=-(int)r;i<r;i++) 
			{
				double s=Math.sqrt(r*r-i*i);
				for(int j=(int) -s;j<s;j++)
			{
				double factor=0.5*(1-Math.cos(Math.PI*Math.sqrt(s*s-j*j)/r));
				//if(c%2==0)
				vectorfield[(x+i)][j+k][0]+=str*factor*Math.signum(1-(u%2));
						//	else vectorfield[scale[0]-1-(x+i)][j+k][0]-=str*factor*Math.signum(1-(u%2));	
			}
		}c++;
		}
	}
	
	//combing force getting ever finer
	public void levelcomb(double t, double r, double str)
	{
		int level=(int)t+1;
		//double scal=Math.pow(level, 1.0/3);
		comb(t,scale[1]/level,level,(r/level),str);
	}
	public void gust(double t,double x, int r, double str)
	{
		if(t%2<1)
		{
			int y=1;
			double sign;
			if(t%4<2)sign=1;
			else {sign=-1;y=scale[1]-2;}
			System.out.println("swoosh");
			int xi;
			for(int i=-r;i<r;i++)
			{
				xi=(int) ((x+i+scale[0])%scale[0]);
				double factor=0.5*(1+Math.cos(Math.PI*Math.abs(i)/r));
				vectorfield[xi][y][1]+=factor*str*sign;
			}
		}
	}
	
	//subtracts vectorfield s from this
	public void subtract(ForceField s) 
	{
		for(int i=0;i<scale[0];i++)
			for(int j=0;j<scale[1];j++)
				for(int k=0;k<2;k++)
				{
					vectorfield[i][j][k]-=s.vectorfield[i][j][k];		
				}
	}

	//adds vectorfield s to this
	public void add(ForceField s)
	{
		for(int i=0;i<scale[0];i++)
			for(int j=0;j<scale[1];j++)
				for(int k=0;k<2;k++)
				{
					vectorfield[i][j][k]+=s.vectorfield[i][j][k];		
				}
	}
	
	public void add(double[] ds, double delt) 
	{
		for(int i=0;i<scale[0];i++)
			for(int j=0;j<scale[1];j++)
				for(int m=0;m<2;m++)
			{
				vectorfield[i][j][m]+=delt*ds[m];
			}
	}
	
	//divides each comp of f by 1+d*|i| overwriting this (Made sense in my head), needed for Fourier
	public void divide(ForceField f, double d) 
	{
		double norm;
		for(int k=0;k<scale[0];k++)
			for(int l=0;l<scale[1];l++)
			{
				norm=k*k+l*l;
				for (int m=0;m<4;m++)
				vectorfield[k][l][m]=f.vectorfield[k][l][m]/(1+d*norm);
			}
	}
	
	//***********************************
	// Trying to figure ouf Fourier, but I couldn't get it to work fast...
	//*************************************
	public ForceField Fourier(int n)
	{
		
		ForceField out=new ForceField(fft(n,1,0,0));//((0r,0i.1r,1i))
	/*	
	// Non-Fast Fourier Transform. Technically works, but slooooow
			for(int j=0;j<scale[1];j++)	
				for(int k=0;k<scale[0];k++)
					for(int l=0;l<scale[1];l++)
					{
						argument=-Math.PI*(2.0*k*i/scale[0]+2.0*l*j/scale[1]);
						for(int m=0;m<2;m++)
						{
							out.vectorfield[k][l][2*m]+=Math.cos(argument)*vectorfield[i][j][m];
							out.vectorfield[k][l][2*m+1]+=Math.sin(argument)*vectorfield[i][j][m];
						}
						/*out.vectorfield[k][l][0]+=Math.cos(argument)*vectorfield[i][j][0]-Math.sin(argument)*vectorfield[i][j][1];
						out.vectorfield[k][l][1]+=Math.sin(argument)*vectorfield[i][j][0]+Math.cos(argument)*vectorfield[i][j][1];
					}
				*/			
		return out;	
		
	}
	public ForceField invFourier(int n)
	{
	//	double argument,sc=scale[0]*scale[1];
		ForceField out=new ForceField(invfft(n,1,0,0));
	/*	for(int i=0;i<scale[0];i++)
			for(int j=0;j<scale[1];j++)	
				for(int k=0;k<scale[0];k++)
					for(int l=0;l<scale[1];l++)
					{
						argument=Math.PI*(2.0*k*i/scale[0]+2.0*l*j/scale[1]);
						for(int m=0;m<2;m++)
						{
							out.vectorfield[k][l][m]+=(Math.cos(argument)*vectorfield[i][j][2*m]-Math.sin(argument)*vectorfield[i][j][2*m+1])/sc;
						//	out.vectorfield[k][l][2*m+1]+=Math.sin(argument)*vectorfield[i][j][m];
						}
					}
				*/	
		return out;	
	}
	

	//projecting to gradient free part in fourier
	public void project(ForceField f) 
	{
		double real,im,norm;
		for(int k=0;k<scale[0];k++)
			for(int l=0;l<scale[1];l++)
			{
				norm=(k*k+l*l);
				real=k*f.vectorfield[k][l][0]+l*f.vectorfield[k][l][2]/norm;
				im=k*f.vectorfield[k][l][1]+l*f.vectorfield[k][l][3]/norm;
				vectorfield[k][l][0]=f.vectorfield[k][l][0]-k*real;
				vectorfield[k][l][1]=f.vectorfield[k][l][1]-k*im;
				vectorfield[k][l][2]=f.vectorfield[k][l][2]-l*real;
				vectorfield[k][l][3]=f.vectorfield[k][l][3]-l*im;
			}
	}

	//fast fourier transform. not fast or not working
	public double[][][] fft(int n, int step, int istart, int jstart)
	{
		double[][][]out=new double[n][n][4];
		
		if(n==1)
		{
			out[0][0][0]=vectorfield[istart][jstart][0];
			out[0][0][2]=vectorfield[istart][jstart][1];
		}
		else
		for(int k=0;k<n/2;k++)
			for(int l=0;l<n/2;l++)
			{
				double[][][]eo=fft(n/2,2*step,istart,jstart+step),ee=fft(n/2,2*step,istart,jstart),oe=fft(n/2,2*step,istart+step,jstart),oo=fft(n/2,2*step,istart+step,jstart+step);
				double cosk=Math.cos(-2*Math.PI*k/n), sink=Math.sin(-2*Math.PI*k/n), cosl= Math.cos(-2*Math.PI*l/n),sinl=Math.sin(-2*Math.PI*l/n),coskl=Math.cos(-2*Math.PI*(k+l)/n), sinkl=Math.sin(-2*Math.PI*(k+l)/n);
				for(int m=0;m<2;m++)
				{
				out[k][l][2*m]=ee[k][l][2*m]+cosl*eo[k][l][2*m]+cosk*oe[k][l][2*m]+coskl*oo[k][l][2*m]-(sinl*eo[k][l][2*m+1]+sink*oe[k][l][2*m+1]+sinkl*oo[k][l][2*m+1]);
				out[k][l][2*m+1]=ee[k][l][2*m+1]+cosl*eo[k][l][2*m+1]+cosk*oe[k][l][2*m+1]+coskl*oo[k][l][2*m+1]+sinl*eo[k][l][2*m]+sink*oe[k][l][2*m]+sinkl*oo[k][l][2*m];
				
				out[k+n/2][l][2*m]=ee[k][l][2*m]+cosl*eo[k][l][2*m]-cosk*oe[k][l][2*m]-coskl*oo[k][l][2*m]-(sinl*eo[k][l][2*m+1]-sink*oe[k][l][2*m+1]-sinkl*oo[k][l][2*m+1]);
				out[k+n/2][l][2*m+1]=ee[k][l][2*m+1]+cosl*eo[k][l][2*m+1]-cosk*oe[k][l][2*m+1]-coskl*oo[k][l][2*m+1]+sinl*eo[k][l][2*m]-sink*oe[k][l][2*m]-sinkl*oo[k][l][2*m];
				
				out[k][l+n/2][2*m]=ee[k][l][2*m]-cosl*eo[k][l][2*m]+cosk*oe[k][l][2*m]-coskl*oo[k][l][2*m]-(-sinl*eo[k][l][2*m+1]+sink*oe[k][l][2*m+1]-sinkl*oo[k][l][2*m+1]);
				out[k][l+n/2][2*m+1]=ee[k][l][2*m+1]-cosl*eo[k][l][2*m+1]+cosk*oe[k][l][2*m+1]-coskl*oo[k][l][2*m+1]-sinl*eo[k][l][2*m]+sink*oe[k][l][2*m]-sinkl*oo[k][l][2*m];
				
				out[k+n/2][l+n/2][2*m]=ee[k][l][2*m]-cosl*eo[k][l][2*m]-cosk*oe[k][l][2*m]+coskl*oo[k][l][2*m]-(-sinl*eo[k][l][2*m+1]-sink*oe[k][l][2*m+1]+sinkl*oo[k][l][2*m+1]);
				out[k+n/2][l+n/2][2*m+1]=ee[k][l][2*m+1]-cosl*eo[k][l][2*m+1]-cosk*oe[k][l][2*m+1]+coskl*oo[k][l][2*m+1]-sinl*eo[k][l][2*m]-sink*oe[k][l][2*m]+sinkl*oo[k][l][2*m];
				}
				
			}
		return out;
	}
	//inverse fast fourier	transform, not functional
	public double[][][] invfft(int n, int step, int istart, int jstart)
	{
		double[][][]out=new double[n][n][4];
		
		if(n==1)
		{
			out[0][0][0]=vectorfield[istart][jstart][0];
			out[0][0][2]=vectorfield[istart][jstart][1];
		}
		else
		for(int k=0;k<n/2;k++)
			for(int l=0;l<n/2;l++)
			{
				double[][][]eo=fft(n/2,2*step,istart,jstart+step),ee=fft(n/2,2*step,istart,jstart),oe=fft(n/2,2*step,istart+step,jstart),oo=fft(n/2,2*step,istart+step,jstart+step);
				double cosk=Math.cos(2*Math.PI*k/n), sink=Math.sin(2*Math.PI*k/n), cosl= Math.cos(2*Math.PI*l/n),sinl=Math.sin(2*Math.PI*l/n),coskl=Math.cos(2*Math.PI*(k+l)/n), sinkl=Math.sin(2*Math.PI*(k+l)/n);
				for(int m=0;m<2;m++)
				{
				out[k][l][2*m]=(ee[k][l][2*m]+cosl*eo[k][l][2*m]+cosk*oe[k][l][2*m]+coskl*oo[k][l][2*m]-(sinl*eo[k][l][2*m+1]+sink*oe[k][l][2*m+1]+sinkl*oo[k][l][2*m+1]))/4;
				out[k][l][2*m+1]=(ee[k][l][2*m+1]+cosl*eo[k][l][2*m+1]+cosk*oe[k][l][2*m+1]+coskl*oo[k][l][2*m+1]+sinl*eo[k][l][2*m]+sink*oe[k][l][2*m]+sinkl*oo[k][l][2*m])/4;
				
				out[k+n/2][l][2*m]=(ee[k][l][2*m]+cosl*eo[k][l][2*m]-cosk*oe[k][l][2*m]-coskl*oo[k][l][2*m]-(sinl*eo[k][l][2*m+1]-sink*oe[k][l][2*m+1]-sinkl*oo[k][l][2*m+1]))/4;
				out[k+n/2][l][2*m+1]=(ee[k][l][2*m+1]+cosl*eo[k][l][2*m+1]-cosk*oe[k][l][2*m+1]-coskl*oo[k][l][2*m+1]+sinl*eo[k][l][2*m]-sink*oe[k][l][2*m]-sinkl*oo[k][l][2*m])/4;
				
				out[k][l+n/2][2*m]=(ee[k][l][2*m]-cosl*eo[k][l][2*m]+cosk*oe[k][l][2*m]-coskl*oo[k][l][2*m]-(-sinl*eo[k][l][2*m+1]+sink*oe[k][l][2*m+1]-sinkl*oo[k][l][2*m+1]))/4;
				out[k][l+n/2][2*m+1]=(ee[k][l][2*m+1]-cosl*eo[k][l][2*m+1]+cosk*oe[k][l][2*m+1]-coskl*oo[k][l][2*m+1]-sinl*eo[k][l][2*m]+sink*oe[k][l][2*m]-sinkl*oo[k][l][2*m])/4;
				
				out[k+n/2][l+n/2][2*m]=(ee[k][l][2*m]-cosl*eo[k][l][2*m]-cosk*oe[k][l][2*m]+coskl*oo[k][l][2*m]-(-sinl*eo[k][l][2*m+1]-sink*oe[k][l][2*m+1]+sinkl*oo[k][l][2*m+1]))/4;
				out[k+n/2][l+n/2][2*m+1]=(ee[k][l][2*m+1]-cosl*eo[k][l][2*m+1]-cosk*oe[k][l][2*m+1]+coskl*oo[k][l][2*m+1]-sinl*eo[k][l][2*m]-sink*oe[k][l][2*m]+sinkl*oo[k][l][2*m])/4;
				}
				
			}
		return out;
	}

}
