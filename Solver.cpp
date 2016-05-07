#include "Solver.h"
#include "Grid.h"
#include "Grid_int.h"
#include "Stencil.h"
#include <fstream>

Solver::Solver(int l_level, int n_Vcycle)
{
    level= l_level;
    Vcycle = n_Vcycle;

    ngp_= pow(2,level)+1;
}

Solver::Solver(int l_level, int n_Vcycle, std::vector<double> &u)
{

level= l_level;
Vcycle = n_Vcycle;

ngp_= pow(2,level)+1; 
u_initial = u;

for(int i= level;i>0;--i)
	{
	Grid *g =  new Grid(i);
	lev_Vec.push_back(g);
	}
lev_Vec[0]->u_app = u_initial;
}

//------------------------------getu_app----------------------------------------------------------//

std::vector<double> Solver::get_u_app(int l_level)
{
 int x=level- l_level;
 std::vector<double> z= lev_Vec[x]->u_app;
 return z;
}

///*********************************MAPPING FUNCTION**********************************************//

int map(int i,int j,int ngl)
{
return  (i*ngl)+j;
}

///*********************************PRINT FUNCTIONS**********************************************//

void Solver::display_u()
{
std::cout<< "level "<< level << std::endl; 
//ngp_= pow(2,level)+1; 
for(size_t i=0; i<ngp_; ++i)
	{
			{	
			for(size_t j=0; j<ngp_;++j)
			std::cout<< u_initial[i*ngp_+j] << "\t";
			}
	std::cout<<std::endl;
	}
}


void Solver::display_u_app(int l_level)
{
std::cout<< '\n' <<"U approximation at level " << l_level << std::endl;
int x = level-l_level;

double ngl_=lev_Vec[x]->get_ngpValue();

//double ngl_= pow(2,l_level)+1;

for(size_t i=0; i<ngl_; ++i)
	{
		{	
		for(size_t j=0; j<ngl_;++j)
		std::cout<< lev_Vec[x]->u_app[i*ngl_+j] << "\t";
		}
	std::cout<<std::endl;
	}
}



void Solver::display_frc(int l_level)
{
std::cout<< '\n'<< "Force vector at level " << l_level << std::endl;

int x = level-l_level;
double ngl_= lev_Vec[x]->Grid::get_ngpValue();


for(size_t i=0; i<ngl_; ++i)
	{
		{	
		for(size_t j=0; j<ngl_;++j)
		std::cout<< lev_Vec[x]->frc[i*ngl_+j] << "\t";
		}
	std::cout<<std::endl;
	}
}


void Solver::display_res(int l_level)
{
std::cout<< '\n' << "Residual vector at level " << l_level << std::endl;

int x = level-l_level;
double ngl_= lev_Vec[x]->Grid::get_ngpValue();


for(size_t i=0; i<ngl_; ++i)
	{
		{	
		for(size_t j=0; j<ngl_;++j)
		std::cout<< lev_Vec[x]->res[i*ngl_+j] << "\t";
		}
	std::cout<<std::endl;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
///********************************* RED-BLACK G.S. *********************************************///
////////////////////////////////////////////////////////////////////////////////////////////////////

void Solver::RBGS(int l_level)
{
int x=level-l_level;
double ngl_= lev_Vec[x]->Grid::get_ngpValue();
double h2_= 1.0/((ngl_-1)*(ngl_-1));

///---------------------------------- RED UPDATE -------------------------------------------------//

	 for(int i=1;i<ngl_-1;++i){
	
		if(i & 1)
		{
            for(int j=1;j<ngl_-1;j+=2){

			lev_Vec[x]->u_app[map(i,j,ngl_)] = 0.25*(lev_Vec[x]->u_app[map(i,j-1,ngl_)]	+
								lev_Vec[x]->u_app[map(i,j+1,ngl_)]	+
								lev_Vec[x]->u_app[map(i-1,j,ngl_)]	+
								lev_Vec[x]->u_app[map(i+1,j,ngl_)]	+
                                (h2_ * lev_Vec[x]->frc[map(i,j,ngl_)]));
		}
	    }
		else{

			for(int j=2;j<ngl_-1;j+=2){
		
			lev_Vec[x]->u_app[map(i,j,ngl_)] = 0.25*(lev_Vec[x]->u_app[map(i,j-1,ngl_)]	+
								lev_Vec[x]->u_app[map(i,j+1,ngl_)]	+
								lev_Vec[x]->u_app[map(i-1,j,ngl_)]	+
								lev_Vec[x]->u_app[map(i+1,j,ngl_)]	+
                                (h2_ * lev_Vec[x]->frc[map(i,j,ngl_)]));
		}
		}

	 }

///---------------------------------- BLACK UPDATE -----------------------------------------------//

	for(int i=1;i<ngl_-1;++i){
	
		if(i & 1)
        {
	  	for(int j=2;j<ngl_-1;j+=2)
			{
			lev_Vec[x]->u_app[map(i,j,ngl_)] = 0.25*(lev_Vec[x]->u_app[map(i,j-1,ngl_)]	+
								lev_Vec[x]->u_app[map(i,j+1,ngl_)]	+
								lev_Vec[x]->u_app[map(i-1,j,ngl_)]	+
								lev_Vec[x]->u_app[map(i+1,j,ngl_)]	+
								(h2_*lev_Vec[x]->frc[map(i,j,ngl_)]));
			}
	    	}
		else
		{

			for(int j=1;j<ngl_-1;j+=2)
			{
		
			lev_Vec[x]->u_app[map(i,j,ngl_)] = 0.25*(lev_Vec[x]->u_app[map(i,j-1,ngl_)]	+
								lev_Vec[x]->u_app[map(i,j+1,ngl_)]	+
								lev_Vec[x]->u_app[map(i-1,j,ngl_)]	+
								lev_Vec[x]->u_app[map(i+1,j,ngl_)]	+
								(h2_*lev_Vec[x]->frc[map(i,j,ngl_)]));
			}
		}

	 }
}

///***************************** SMOOTHING FUNCTIONS *********************************************//

void Solver::pre_smoothing(int l_level)
{
    for(int i=0;i<2;++i)
        this->RBGS(l_level);
}

void Solver::post_smoothing(int l_level)
{
	this->RBGS(l_level);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
///*********************************** RESIDUAL **************************************************//
////////////////////////////////////////////////////////////////////////////////////////////////////

void Solver::residual(int l_level)
{

int x=level-l_level;
double ngl_= lev_Vec[x]->Grid::get_ngpValue();
double h2_ = ((ngl_-1)*(ngl_-1)); // h2_ =(h^2)
//std::cout<< "h sqr "<< h2_ <<std::endl ;
double norm=0;
double temp=0;
for(int i=1;i<ngl_-1;++i)
	{
	for(int j=1;j<ngl_-1;++j)
		{
		//std::cout<< "i "<< i <<" j " << j << " to map " << map(i,j,ngl_) << std::endl ;
		lev_Vec[x]->res[map(i,j,ngl_)]=	lev_Vec[x]->frc[map(i,j,ngl_)]+	(h2_*(lev_Vec[x]->u_app[map(i,j-1,ngl_)]		+
										lev_Vec[x]->u_app[map(i,j+1,ngl_)]			+
										lev_Vec[x]->u_app[map(i-1,j,ngl_)]			+
										lev_Vec[x]->u_app[map(i+1,j,ngl_)]			-
										(4*lev_Vec[x]->u_app[map(i,j,ngl_)])));
		/*if(x!=0)
		{
		std::cout<< "\n neighbour points " << lev_Vec[x]->u_app[map(i,j-1,ngl_)] << "\n"
						<< lev_Vec[x]->u_app[map(i,j+1,ngl_)] << "\n"
						<<  lev_Vec[x]->u_app[map(i-1,j,ngl_)] <<"\n"
						<<lev_Vec[x]->u_app[map(i+1,j,ngl_)] << std::endl;
		std::cout<< "\n center point " << lev_Vec[x]->u_app[map(i,j,ngl_)] <<std::endl;
		}
		std::cout<< lev_Vec[x]->u_app[map(i,j,ngl_)] << '\t' ;*/
		}
	//std::cout<<std::endl;
	}
	
///////Residual Norm
for(int i=1;i<ngl_-1;++i)
	{
	for(int j=1;j<ngl_-1;++j)
		{
			temp=lev_Vec[x]->res[map(i,j,ngl_)];
			norm+=(temp*temp);
		}
	}
norm= norm/((ngl_-2)*(ngl_-2));
norm = sqrt(norm);
std::cout<< "the residual norm at level " << l_level << " is " << norm << std::endl;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
///*********************************** RESTRICTION ***********************************************//
////////////////////////////////////////////////////////////////////////////////////////////////////


void Solver::restriction(int l_level)
{
int x=level-l_level;
double ngl_= lev_Vec[x]->get_ngpValue();
double ngl_1dwn = lev_Vec[x+1]->get_ngpValue();
Stencil_rest r;


for(int i=0;i<ngl_;++i)
{
	if(i%2==0)
	{
	for(int j=0;j<ngl_;++j)
	{
		if(j%2==0)
		{	
			/*if((j==0 || j==ngl_-1))
			{
			//std::cout<< "i " << i << " j " <<j<< " map " << map(i,j,ngl_) << " to " << map(i/2,j/2,ngl_1dwn) << std::endl;
			lev_Vec[x+1]->u_app[map(i/2,j/2,ngl_1dwn)] = lev_Vec[x]->u_app[map(i,j,ngl_)];
			}	
			else if((i==0 || i==ngl_-1))
			{
			//std::cout<< "i " << i << " j " <<j<< " map " << map(i,j,ngl_) << " to " << map(i/2,j/2,ngl_1dwn) << std::endl;
			lev_Vec[x+1]->u_app[map(i/2,j/2,ngl_1dwn)] = lev_Vec[x]->u_app[map(i,j,ngl_)];
			}
			else */
			lev_Vec[x+1]->u_app[map(i/2,j/2,ngl_1dwn)]=0;	
		}
	}
}
}


/*for(int i=2;i<ngl_-2;i+=2)
{
	for(int j=2;j<ngl_-2;j+=2)
	{
	lev_Vec[x+1]->u_app[map(i/2,j/2,ngl_1dwn)]=   r.W*lev_Vec[x]->u_app[map(i-1,j,ngl_)] 		+ 
                		               		r.C*lev_Vec[x]->u_app[map(i,j,ngl_)]   		+ 
                                			r.E*lev_Vec[x]->u_app[map(i+1,j,ngl_)]  	+ 
				                	r.SW*lev_Vec[x]->u_app[map(i-1,j-1,ngl_)] 	+ 
							r.S*lev_Vec[x]->u_app[map(i,j-1,ngl_)] 		+ 
							r.SE*lev_Vec[x]->u_app[map(i-1,j+1,ngl_)] 	+ 
							r.N*lev_Vec[x]->u_app[map(i,j+1,ngl_)] 		+ 
							r.NW*lev_Vec[x]->u_app[map(i+1,j-1,ngl_)] 	+ 
							r.NE*lev_Vec[x]->u_app[map(i+1,j+1,ngl_)];
	}
}*/

for(int i=2;i<ngl_-2;i+=2)
{
	for(int j=2;j<ngl_-2;j+=2)
	{
	lev_Vec[x+1]->frc[map(i/2,j/2,ngl_1dwn)]=    	r.W*lev_Vec[x]->res[map(i-1,j,ngl_)] 		+ 
							r.C*lev_Vec[x]->res[map(i,j,ngl_)]   		+ 
							r.E*lev_Vec[x]->res[map(i+1,j,ngl_)]  		+ 
							r.SW*lev_Vec[x]->res[map(i-1,j-1,ngl_)] 	+ 
							r.S*lev_Vec[x]->res[map(i,j-1,ngl_)] 		+ 
							r.SE*lev_Vec[x]->res[map(i-1,j+1,ngl_)] 	+ 
							r.N*lev_Vec[x]->res[map(i,j+1,ngl_)] 		+ 
							r.NW*lev_Vec[x]->res[map(i+1,j-1,ngl_)] 	+ 
							r.NE*lev_Vec[x]->res[map(i+1,j+1,ngl_)];
	}
}

}


////////////////////////////////////////////////////////////////////////////////////////////////////
///*********************************** PROLONGATION **********************************************//
////////////////////////////////////////////////////////////////////////////////////////////////////

void Solver::prolongation(int l_level)
{
int x=level-l_level; 
double ngl_= lev_Vec[x]->get_ngpValue();   //coarse grid
double ngl_up = lev_Vec[x-1]->get_ngpValue();  //fine grid


Stencil_prol p;

for(size_t i=1;i<ngl_-1;i++)
	{
	for(size_t j=1;j<ngl_-1;j++)
		{

		lev_Vec[x-1]->u_app[map(2*i-1,2*j-1,ngl_up)]	+= 	p.SW*lev_Vec[x]->u_app[map(i,j,ngl_)];
		lev_Vec[x-1]->u_app[map(2*i-1,2*j,ngl_up)]	+= 	p.S*lev_Vec[x]->u_app[map(i,j,ngl_)];
		lev_Vec[x-1]->u_app[map(2*i-1,2*j+1,ngl_up)]  	+= 	p.SE*lev_Vec[x]->u_app[map(i,j,ngl_)];
		lev_Vec[x-1]->u_app[map(2*i,2*j-1,ngl_up)] 	+= 	p.W*lev_Vec[x]->u_app[map(i,j,ngl_)];
		lev_Vec[x-1]->u_app[map(2*i,2*j,ngl_up)]         	+= 	p.C*lev_Vec[x]->u_app[map(i,j,ngl_)];
		lev_Vec[x-1]->u_app[map(2*i,2*j+1,ngl_up)] 	+= 	p.E*lev_Vec[x]->u_app[map(i,j,ngl_)];
		lev_Vec[x-1]->u_app[map(2*i+1,2*j-1,ngl_up)] 	+=    	p.NW*lev_Vec[x]->u_app[map(i,j,ngl_)];
		lev_Vec[x-1]->u_app[map(2*i+1,2*j,ngl_up)] 	+= 	p.N*lev_Vec[x]->u_app[map(i,j,ngl_)];
		lev_Vec[x-1]->u_app[map(2*i+1,2*j+1,ngl_up)] 	+= 	p.NE*lev_Vec[x]->u_app[map(i,j,ngl_)];
		}
	}

}

void store(double ngp_,std::vector<double>u)
{
double hgl_ = 1.0 / (ngp_-1);

std::ofstream mg;
mg.open("multigrid.txt");
for(size_t i=0; i<ngp_; ++i)
{
    {
    for(size_t j=0; j<ngp_;++j)
    mg<< i*hgl_<<"\t"<< j*hgl_<<"\t" <<u[i*ngp_+j] << "\n";
    }
}
mg.close();
}


//////////////////////////////Simulation///////////////////////////////////////////////////

void Solver::Simulation()
{
int l_Level = this-> level;
int n_Vcycle = this-> Vcycle;

/// Initialise grid.
Grid_int v(l_Level);
//std::cout<< "\nInitial level = " << l_Level << " number of grid points = " << v.get_ngpValue()*v.get_ngpValue() << std::endl;
//v.display_grid_int();

/// Apply Boundary conditions.
v.boundary_con();
//std::cout<< "\nAfter applying BC "<< std::endl;
//v.display_grid_int();

std::vector<double> u = v.get_Xvalue();

for(int i=1; i<=n_Vcycle; ++i)
    {
  std::cout<< "Current V-cycle " << i <<std::endl;

    Solver S(l_Level,i,u);

    //std::cout<< "\n Restriction starts from here " << std::endl;
        for (int j =l_Level; j>0; --j) // Pre Smoothning -> U Print -> Residual -> Residual print -> Restriction -> Force Print
        {
        S.pre_smoothing(j); //2*S.RBGS(j);

       // std::cout<< "\nAfter 2 RBGS "<< std::endl;
       //S.display_u_app(j);


        S.residual(j);
        //S.display_frc(j);
        //std::cout<< "\n Residual "<< std::endl;
        //S.display_res(j);

        if(j!=1)
        {
        S.restriction(j);
        //S.display_u_app(j-1);
        //std::cout<< "\n Force "<< std::endl;
        //S.display_frc(j);
        //S.display_frc(j-1);
        }
    }

        S.post_smoothing(1);
        //std::cout<< "After 1 RBGS u at coarsest level"<<std::endl;
        //S.display_u_app(1);

    //std::cout<< "\n Prolongation starts from here " << std::endl;

    for (int k =1; k<l_Level; ++k) // Prolongation -> U Print +1 Level -> Residual -> Residual print -> Restriction -> Force Print
    {

        S.prolongation(k);
        //if(k==1)
      //  std::cout<< "\nProlongated u "<<std::endl;
        //S.display_u_app(k+1);

        if(k!=1)
        {
            S.post_smoothing(k+1);
        //    std::cout<< "After 1 RBGS u"<<std::endl;
            //S.display_u_app(k+1);

            /*if(k==l_Level-1)
            S.display_u_app(k+1);*/
        }
    }

    //std::cout<<" \n u after one v cycle " << std::endl;
    //S.display_u_app(l_Level);
    //Reassign u for next V-cycle.
    u=S.get_u_app(l_Level);


}

//std::cout<<" final u after " << n_Vcycle << " v cycle " << std::endl;
/*for(size_t i=0; i<ngp_; ++i)
{
    {
    for(size_t j=0; j<ngp_;++j)
    std::cout<< u[i*ngp_+j] << "\t";
    }
std::cout<<std::endl;
}*/

store(ngp_,u);

}



Solver::~Solver()
{
//for(int i= 0;i<level;i++)
//delete lev_Vec[i];
}

