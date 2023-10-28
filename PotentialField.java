//******************************************************************************
// FluidSim/PotentialField.java
// author: Non-Euclidean Dreamer
// Class of Scalar Fields, mostly used as a PotentialField for solving Poisson
//******************************************************************************

public class PotentialField 
{
	static int[] scale=FlowField.scale;
	public double[][]potential;
	
	
	//************
	//Constructors
	//************
	public PotentialField()
	{
		potential=new double[scale[0]][scale[1]];
	}
	
	//uniformly equal to p
	public PotentialField(double p) 
	{
		potential=new double[scale[0]][scale[1]];
		for(int i=0;i<scale[0];i++)
			for(int j=0;j<scale[1];j++)
				potential[i][j]=p;
	}	

	

	//Overwrite uniformly
	public void set(double p)
	{
		for(int i=0;i<scale[0];i++)
			for(int j=0;j<scale[1];j++)
				potential[i][j]=p;
	}
	
	//********************************************
	// Methods working with scalar field structure
	//********************************************
	
	// One step in a relaxation scheme to solve the Poisson equation del2 out=div, with p being the preceding solution
	public double poisson(PotentialField p, PotentialField div, double delx) 
	{

		double out=0;
		for(int i=0;i<scale[0];i++)
			for(int j=0;j<scale[1];j++)
			{
				int im=(i+scale[0]-1)%scale[0],jm=(j+scale[1]-1)%scale[1],ip=(i+1)%scale[0],jp=(j+1)%scale[1];
				potential[i][j]=(p.potential[im][j]+p.potential[ip][j]+p.potential[i][jm]+p.potential[i][jp]-delx*delx*div.potential[i][j])/4;
				out=Math.max(out, Math.abs(potential[i][j]-p.potential[i][j]));
			}
				return out;
	}
	
	//Dito, but with boundary condition on obstacle enforcing the gradient to equal f normal to obstacle
	public double poisson(PotentialField p, PotentialField div, double strength, int[][]obstacle,double[][]norm,ForceField f) 
	{
		double out=0; int k=0;
		for(int i=0;i<scale[0];i++)
		{	for(int j=0;j<scale[1];j++)
			{
				if(obstacle[k][0]==i&&obstacle[k][1]==j)
				{
					if(norm[k][0]==0&&norm[k][1]==0)potential[i][j]=0;else {
					double[] loc= {(i+norm[k][0]+scale[0])%scale[0],(j+norm[k][1]+scale[1])%scale[1]};
					double[]mid= {(i+norm[k][0]/2+scale[0])%scale[0],(j+norm[k][1]/2+scale[1])%scale[1]};
					double qloc=FlowField.average(p.potential, loc);
					double[]umid=FlowField.average(f.vectorfield, mid);
			
					potential[i][j]=qloc-(umid[0]*norm[k][0]+umid[1]*norm[k][1])*strength;
					}
					k++;
				}
				else 
				{
					int im=(i+scale[0]-1)%scale[0],jm=(j+scale[1]-1)%scale[1],ip=(i+1)%scale[0],jp=(j+1)%scale[1];
					potential[i][j]=(p.potential[im][j]+p.potential[ip][j]+p.potential[i][jm]+p.potential[i][jp]-strength*strength*div.potential[i][j])/4;
				}
				out=Math.max(out, Math.abs(potential[i][j]-p.potential[i][j]));
				}
		}
				return out;
	}
	
	public ForceField gradient() 
	{
		ForceField out=new ForceField();
		for(int i=0;i<scale[0];i++)
			for(int j=0;j<scale[1];j++)
			{
				int im=(i+scale[0]-1)%scale[0],jm=(j+scale[1]-1)%scale[1],ip=(i+1)%scale[0],jp=(j+1)%scale[1];
				out.vectorfield[i][j][0]=(potential[ip][j]-potential[im][j])/2;
				out.vectorfield[i][j][1]=(potential[i][jp]-potential[i][jm])/2;
			}
		return out;
	}
	
}
