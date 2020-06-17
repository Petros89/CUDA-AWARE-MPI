/// this program solves the 3D heat equation on a 3D structured, cartesian grid using MPI.


/**
 * The equation we want to solve can be expressed in the following way:
 
 *
 * T_t = Dx * T_xx + Dy * T_yy + Dz * T_zz,
 *
 *
 * We use a second order accurate central scheme for the space derivatives, i.e. we have (in 1D):


 * d^2 T(x) / dx^2 = T_xx ~= (T[i+1] - 2*T[i] + T[i-1]) / (dx^2)
 *
 * which we can apply in each coordinate direction equivalently. dx is the spacing between to adjacent cells, i.e. the
 * distance from one cell to its neighbors. It can be different for the y and z direction, however, within the same
 * direction it is always constant. For the time derivative, we use a first order Euler time integration scheme like so:
 *
 * dT(x) / dt = T_t ~= (T[n+1] - T[n]) / dt
 *
 * Here, n is the timestep from the previous solution and n+1 is the timestep for the next solution. In this way we can
 * integrate our solution in time. Combining the two above approximations, we could write (dor a 1D equation)
 * T_t = Dx * T_xx =>
 * (T[n+1] - T[n]) / dt = Dx * (T[i+1] - 2*T[i] + T[i-1]) / (dx^2)
 *
 * We can solve this for T[n+1] to yield:
 * T[n+1] = T[n] + (dt * Dx / (dx^2)) * (T[i+1] - 2*T[i] + T[i-1])
 *
 * We have the information of the right hand side available, thus we can calculate T[n+1] for each i.
 * For i=0 or i=iend we need to specify boundary conditions and for all T[n] we need to specify initial conditions.
 * With those information available, we can loop over time and calculate an updated solution until the solution between
 * two consequtive time steps does not change more than a user-defined convergence threshold.
 *
 * For more information on the heat equation, you may check the following link:
 * https://www.uni-muenster.de/imperia/md/content/physik_tp/lectures/ws2016-2017/num_methods_i/heat.pdf
 */

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <array>
#include <fstream>
#include <limits>
#include <cmath>
#include <chrono>
#include <cassert>
#include "mpi.h"
#include<string>
#include <cuda.h>
#include <malloc.h>

#define DIM_THREAD_BLOCK_X 2
#define DIM_THREAD_BLOCK_Y 1
using namespace std;


//method to build 3D pointer
double*** CreateGrid(int m,int n,int t)
{
    int i = 0;
    int k = 0;
    double*** result = NULL; 
    if((m > 0) && (n > 0) && (t > 0))
    {
        double** pp = NULL;
        double* p = NULL;
        result = (double***)malloc(m * sizeof(double**));     
        pp = (double**)malloc(m * n * sizeof(double*));      
        p = (double*)malloc(m * n * t * sizeof(double));     
        if((result != NULL) && (pp != NULL) && (p != NULL))
        {
            for(i = 0;i < m;i++)
            {
                result[i] = pp + i * n; 
                for (k = 0;k < n;k++)
                {
                    result[i][k] = p + k * t; 
                }
                p = p + n*t;
            }
        }
        else
        {
            free(result);
            free(pp);
            free(p);
            result = NULL;
            pp = NULL;
            p = NULL;
        }
    }
    return result;
}
void FreeGrid(double*** p)
{
    if(*p != NULL)
    {
        if(**p != NULL)
        {
            free(**p);
            **p = NULL;
        }
        free(*p);
        *p = NULL;
    }
    free(p);
    p = NULL;
}

/*********************************************************************************************
                                                   GPU kernel methods
**********************************************************************************************/


__global__  void computeT(double*** TBegin, double ***TEnd, double ***Tres_gpu,int numX, int numY, int numZ, double Dx, double Dy, double Dz) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	
	if(i>1&&i<(numX-2)){

	//Distribute data by X direction
	if(i<(numX)){
   
        	for (unsigned j = 1; j < numY - 1; ++j){
            		for (unsigned k = 1; i < numZ - 1; ++k) {
                		TEnd[i][j][k] = TBegin[i][j][k] +
                   		 Dx * (TBegin[i+1][j][k] - 2.0 * TBegin[i][j][k] + TBegin[i-1][j][k] )+
                   		 Dy * (TBegin[i][j+1][k] - 2.0 * TBegin[i][j][k] + TBegin[i][j-1][k]) +
                   		 Dz * (TBegin[i][j][k+1] - 2.0 * TBegin[i][j][k] + TBegin[i][j][k-1]);
				Tres_gpu[i][j][k]=TEnd[i][j][k];


            		}

		}

	}
	}

}






// based on compiler flag, use either floats or doubles for floating point operations


using floatT = double;

#define MPI_FLOAT_T MPI_DOUBLE


/// enum used to index over the respective coordinate direction
enum COORDINATE { X = 0, Y, Z };

/// enum used to access the respective direction on each local processor

/**
 *  0: LEFT
 *  1: RIGHT
 *  2: BOTTOM
 *  3: TOP
 *  4: BACK
 *  5: FRONT
 */
enum DIRECTION { LEFT = 0, RIGHT, BOTTOM, TOP, BACK, FRONT };

/// the number of physical dimensions, here 3 as we have a 3D domain

#define NUMBER_OF_DIMENSIONS 3

int main(int argc, char** argv)
{
    /// if USE_MPI is defined (see makefile), execute the following code


  /// default ranks and size (number of processors), will be rearranged by cartesian topology
 
    int rankDefaultMPICOMM, sizeDefaultMPICOMM;

    /// status and requests for non-blocking communications, i.e. MPI_IAllreduce(...) and MPI_IRecv(...)
   

    MPI_Status  status[NUMBER_OF_DIMENSIONS * 2];
    MPI_Status  postStatus[NUMBER_OF_DIMENSIONS];
    MPI_Request request[NUMBER_OF_DIMENSIONS * 2];
    MPI_Request reduceRequest;

    /// buffers into which we write data that we want to send and receive using MPI
  
    /**
     * sendbuffer will be received into receivebuffer\
    
     */
    std::array<std::vector<floatT>, NUMBER_OF_DIMENSIONS * 2> sendBuffer;
    std::array<std::vector<floatT>, NUMBER_OF_DIMENSIONS * 2> receiveBuffer;

    /// initialise MPI and get default ranks and size
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankDefaultMPICOMM);
    MPI_Comm_size(MPI_COMM_WORLD, &sizeDefaultMPICOMM);

    /// new MPI communicator for cartesian topologies
   
    MPI_Comm MPI_COMM_CART;

    /// new rank and size for cartesian topology
    int       rank, size;

    /// tag used later during MPI_Send(...)
   
    int       tagSend[NUMBER_OF_DIMENSIONS * 2];

    /// tag used later during MPI_IRecv(...)
    int       tagReceive[NUMBER_OF_DIMENSIONS * 2];

    /// the dimensions are equivalent to how we want our domain to be partitioned.


    int       dimension3D[NUMBER_OF_DIMENSIONS] = { 0, 0, 0 };

    /// the coordinate in the current cartesian topology for the sub processor
    

    int       coordinates3D[NUMBER_OF_DIMENSIONS];

    /// flags to indicate if we have period boundary conditions
    
    const int periods3D[NUMBER_OF_DIMENSIONS] = { false, false, false };

    /// neighbors hold the rank of the neighboring processors and are accessed with the DIRECTION enum
   

    int       neighbors[NUMBER_OF_DIMENSIONS * 2];

    /// MPI tries to find the best possible partition of our domain and stores that in dimension3D
 

    MPI_Dims_create(sizeDefaultMPICOMM, NUMBER_OF_DIMENSIONS, dimension3D);

    /// based on the partition, we create a new cartesian topology which simplifies communication
    
    MPI_Cart_create(MPI_COMM_WORLD, NUMBER_OF_DIMENSIONS, dimension3D, periods3D, true, &MPI_COMM_CART);

    /// These calls will find the direct neighbors for each processors and return MPI_PROC_NULL if no neighbor is found.
    
    MPI_Cart_shift(MPI_COMM_CART, COORDINATE::X, 1, &neighbors[DIRECTION::LEFT], &neighbors[DIRECTION::RIGHT]);
    MPI_Cart_shift(MPI_COMM_CART, COORDINATE::Y, 1, &neighbors[DIRECTION::BOTTOM], &neighbors[DIRECTION::TOP]);
    MPI_Cart_shift(MPI_COMM_CART, COORDINATE::Z, 1, &neighbors[DIRECTION::BACK], &neighbors[DIRECTION::FRONT]);

    /// get the new rank and size for the cartesian topology
   

    MPI_Comm_rank(MPI_COMM_CART, &rank);
    MPI_Comm_size(MPI_COMM_CART, &size);

    /// get the coordinates inside our cartesian topology
 
    MPI_Cart_coords(MPI_COMM_CART, rank, NUMBER_OF_DIMENSIONS, coordinates3D);

    /// if USE_SEQUENTIAL is defined (see makefile), execute the following code
     



    /// check that we have the right number of input arguments
    
    /**
     * this is the order in which we need to pass in the command line argument:


     *
     * argv[0]: name of compiled program
     * argv[1]: number of cells in the x direction
     * argv[2]: number of cells in the y direction
     * argv[3]: number of cells in the z direction
     * argv[4]: maximum number of iterations to be used by time loop    
     * argv[5]: convergence criterion to be used to check if a solution has converged 
     */
    if (rank == 0) {
        if (argc != 6) {
            std::cout << "Incorrect number of command line arguments specified, use the following syntax:\n" << std::endl;
            std::cout << "bin/HeatEquation3D NUM_CELLS_X NUM_CELLS_Y NUM_CELLS_Z ITER_MAX EPS" << std::endl;
            std::cout << "\nor, using MPI, use the following syntax:\n" << std::endl;
            std::cout << "mpirun -n NUM_PROCS bin/HeatEquation3D NUM_CELLS_X NUM_CELLS_Y NUM_CELLS_Z ITER_MAX EPS" << std::endl;
            std::cout << "\nSee source code for additional informations!" << std::endl;
            std::abort();
        }
        else {
            std::cout << "Runnung HeatEquation3D with the following arguments: " << std::endl;
            std::cout << "executable:               " << argv[0] << std::endl;
            std::cout << "number of cells in x:     " << std::stoi(argv[1]) << std::endl;
            std::cout << "number of cells in y:     " << std::stoi(argv[2]) << std::endl;
            std::cout << "number of cells in z:     " << std::stoi(argv[3]) << std::endl;
            std::cout << "max number of iterations: " << std::stoi(argv[4]) << std::endl;


            std::cout << "convergence threshold:    " << std::stod(argv[5]) << "\n" << std::endl;

        }
    }

    /// maximum number of iterations to perform in time loop
   
    const unsigned iterMax = std::stoi(argv[4]);

    /// convergence criterion, which, once met, will terminate the calculation
    


    const floatT eps = std::stod(argv[5]);


    /// both variables are used to calculate the convergence and normalise the result.
         /**
     * We have two normalisation factors as we have to perform a reduction first (if we use MPI) to have a globally
     * available normalisation factor

          */
    floatT globalNorm = 1.0;
    floatT norm = 1.0;

    /// the break conditions used for checking of convergence has been achieved and the simulation should be stopped.
       int breakCondition = false;
    int globalBreakCondition = false;

    /// number of points (in total, not per processor) in x, y and z.
       unsigned numCells[NUMBER_OF_DIMENSIONS];
    numCells[COORDINATE::X] = std::stoi(argv[1]);
    numCells[COORDINATE::Y] = std::stoi(argv[2]);
    numCells[COORDINATE::Z] = std::stoi(argv[3]);

    /// length of the domain in x, y and z.
        floatT domainLength[NUMBER_OF_DIMENSIONS];
    domainLength[COORDINATE::X] = 1.0;
    domainLength[COORDINATE::Y] = 1.0;
    domainLength[COORDINATE::Z] = 1.0;

    /// thermal conductivity parameter. 。

    const floatT alpha = 1.0;

    /// The courant fridrich levy number                    
    const floatT CFL = 0.4;

    /// the distance between cells in the x, y and z direction.
    
    floatT spacing[NUMBER_OF_DIMENSIONS];
    spacing[COORDINATE::X] = domainLength[COORDINATE::X] / static_cast<floatT>(numCells[COORDINATE::X] - 1.0);
    spacing[COORDINATE::Y] = domainLength[COORDINATE::Y] / static_cast<floatT>(numCells[COORDINATE::Y] - 1.0);
    spacing[COORDINATE::Z] = domainLength[COORDINATE::Z] / static_cast<floatT>(numCells[COORDINATE::Z] - 1.0);

    /// the timestep to be used in the time integration.
   
    const floatT dt = CFL * 1.0 / (NUMBER_OF_DIMENSIONS * 2) *
        std::pow(std::min({ spacing[COORDINATE::X], spacing[COORDINATE::Y], spacing[COORDINATE::Z] }), 2.0) / alpha;

    /// thermal diffusivity strength in the x, y and z direction.
   

    const floatT Dx = dt * alpha / (std::pow(spacing[COORDINATE::X], 2.0));
    const floatT Dy = dt * alpha / (std::pow(spacing[COORDINATE::Y], 2.0));
    const floatT Dz = dt * alpha / (std::pow(spacing[COORDINATE::Z], 2.0));

    /// numer of iterations taken to converge solution. will be set once simulation has converged.
       unsigned finalNumIterations = 0;


    /// assure that the partition given to use by MPI can be used to partition our domain in each direction
  
    assert((numCells[COORDINATE::X] - 1) % dimension3D[COORDINATE::X] == 0 &&
        "Can not partition data for given number of processors in x!");
    assert((numCells[COORDINATE::Y] - 1) % dimension3D[COORDINATE::Y] == 0 &&
        "Can not partition data for given number of processors in y!");
    assert((numCells[COORDINATE::Z] - 1) % dimension3D[COORDINATE::Z] == 0 &&
        "Can not partition data for given number of processors in z!");

    /// chunck contains the number of cells in the x, y and z direction for each sub domain.
   

    const unsigned chunck[NUMBER_OF_DIMENSIONS] = {
      ((numCells[COORDINATE::X] - 1) / dimension3D[COORDINATE::X]) + 1,
      ((numCells[COORDINATE::Y] - 1) / dimension3D[COORDINATE::Y]) + 1,
      ((numCells[COORDINATE::Z] - 1) / dimension3D[COORDINATE::Z]) + 1
    };


    /// Create a solution vector

    std::vector<std::vector<std::vector<floatT>>> T, T0;

    /// resize both T and T0 for each sub-domain
      T.resize(chunck[COORDINATE::X]);
    T0.resize(chunck[COORDINATE::X]);
    for (unsigned i = 0; i < chunck[COORDINATE::X]; ++i) {
        T[i].resize(chunck[COORDINATE::Y]);
        T0[i].resize(chunck[COORDINATE::Y]);
        for (unsigned j = 0; j < chunck[COORDINATE::Y]; ++j) {
            T[i][j].resize(chunck[COORDINATE::Z]);
            T0[i][j].resize(chunck[COORDINATE::Z]);
        }
    }

    /// initialise each solution vector on each sub-domain with zero everywhere
       for (unsigned i = 0; i < chunck[COORDINATE::X]; ++i)
        for (unsigned j = 0; j < chunck[COORDINATE::Y]; ++j)
            for (unsigned k = 0; k < chunck[COORDINATE::Z]; ++k)
                T[i][j][k] = 0.0;

    /// apply boundary conditions on the top of the domain
    
    if (neighbors[DIRECTION::TOP] == MPI_PROC_NULL)

        for (unsigned i = 0; i < chunck[COORDINATE::X]; ++i)
            for (unsigned k = 0; k < chunck[COORDINATE::Z]; ++k)
                T[i][chunck[COORDINATE::Y] - 1][k] = 1.0;

    /// apply boundary conditions on the left-side of the domain
     

    if (neighbors[DIRECTION::LEFT] == MPI_PROC_NULL)

        for (unsigned j = 0; j < chunck[COORDINATE::Y]; ++j)
            for (unsigned k = 0; k < chunck[COORDINATE::Z]; ++k)
                T[0][j][k] = (coordinates3D[COORDINATE::Y] * (chunck[COORDINATE::Y] - 1) + j) * spacing[COORDINATE::Y];

    /// apply boundary conditions on the right-side of the domain
      
    if (neighbors[DIRECTION::RIGHT] == MPI_PROC_NULL)

        for (unsigned j = 0; j < chunck[COORDINATE::Y]; ++j)
            for (unsigned k = 0; k < chunck[COORDINATE::Z]; ++k)
                T[chunck[COORDINATE::X] - 1][j][k] = (coordinates3D[COORDINATE::Y] * (chunck[COORDINATE::Y] - 1) + j) * spacing[COORDINATE::Y];

    /// apply boundary conditions on the back-side of the domain
    
    if (neighbors[DIRECTION::BACK] == MPI_PROC_NULL)

        for (unsigned i = 0; i < chunck[COORDINATE::X]; ++i)
            for (unsigned j = 0; j < chunck[COORDINATE::Y]; ++j)
                T[i][j][0] = (coordinates3D[COORDINATE::Y] * (chunck[COORDINATE::Y] - 1) + j) * spacing[COORDINATE::Y];

    /// apply boundary conditions on the front-side of the domain
     
    if (neighbors[DIRECTION::FRONT] == MPI_PROC_NULL)

        for (unsigned i = 0; i < chunck[COORDINATE::X]; ++i)
            for (unsigned j = 0; j < chunck[COORDINATE::Y]; ++j)
                T[i][j][chunck[COORDINATE::Z] - 1] = (coordinates3D[COORDINATE::Y] * (chunck[COORDINATE::Y] - 1) + j) * spacing[COORDINATE::Y];

    /// if we use MPI, make sure that our send and recieve buffers are correctly allocated
    

  /// allocate storage for left-side send- and recievebuffer
  
    if (neighbors[DIRECTION::LEFT] != MPI_PROC_NULL) {
        sendBuffer[DIRECTION::LEFT].resize((chunck[COORDINATE::Y] - 1) * (chunck[COORDINATE::Z] - 1));
        receiveBuffer[DIRECTION::LEFT].resize((chunck[COORDINATE::Y] - 1) * (chunck[COORDINATE::Z] - 1));
    }
    else {
        sendBuffer[DIRECTION::LEFT].resize(1);
        receiveBuffer[DIRECTION::LEFT].resize(1);
    }

    /// allocate storage for right-side send- and recievebuffer
   
    if (neighbors[DIRECTION::RIGHT] != MPI_PROC_NULL) {
        sendBuffer[DIRECTION::RIGHT].resize((chunck[COORDINATE::Y] - 1) * (chunck[COORDINATE::Z] - 1));
        receiveBuffer[DIRECTION::RIGHT].resize((chunck[COORDINATE::Y] - 1) * (chunck[COORDINATE::Z] - 1));
    }
    else {
        sendBuffer[DIRECTION::RIGHT].resize(1);
        receiveBuffer[DIRECTION::RIGHT].resize(1);
    }

    /// allocate storage for bottom-side send- and recievebuffer
   
    if (neighbors[DIRECTION::BOTTOM] != MPI_PROC_NULL) {
        sendBuffer[DIRECTION::BOTTOM].resize((chunck[COORDINATE::X] - 1) * (chunck[COORDINATE::Z] - 1));
        receiveBuffer[DIRECTION::BOTTOM].resize((chunck[COORDINATE::X] - 1) * (chunck[COORDINATE::Z] - 1));
    }
    else {
        sendBuffer[DIRECTION::BOTTOM].resize(1);
        receiveBuffer[DIRECTION::BOTTOM].resize(1);
    }

    /// allocate storage for top-side send- and recievebuffer
    
    if (neighbors[DIRECTION::TOP] != MPI_PROC_NULL) {
        sendBuffer[DIRECTION::TOP].resize((chunck[COORDINATE::X] - 1) * (chunck[COORDINATE::Z] - 1));
        receiveBuffer[DIRECTION::TOP].resize((chunck[COORDINATE::X] - 1) * (chunck[COORDINATE::Z] - 1));
    }
    else {
        sendBuffer[DIRECTION::TOP].resize(1);
        receiveBuffer[DIRECTION::TOP].resize(1);

    }

    /// allocate storage for back-side send- and recievebuffer
    
    if (neighbors[DIRECTION::BACK] != MPI_PROC_NULL) {
        sendBuffer[DIRECTION::BACK].resize((chunck[COORDINATE::X] - 1) * (chunck[COORDINATE::Y] - 1));
        receiveBuffer[DIRECTION::BACK].resize((chunck[COORDINATE::X] - 1) * (chunck[COORDINATE::Y] - 1));
    }
    else {
        sendBuffer[DIRECTION::BACK].resize(1);
        receiveBuffer[DIRECTION::BACK].resize(1);
    }

    /// allocate storage for front-side send- and recievebuffer
   
    if (neighbors[DIRECTION::FRONT] != MPI_PROC_NULL) {
        sendBuffer[DIRECTION::FRONT].resize((chunck[COORDINATE::X] - 1) * (chunck[COORDINATE::Y] - 1));
        receiveBuffer[DIRECTION::FRONT].resize((chunck[COORDINATE::X] - 1) * (chunck[COORDINATE::Y] - 1));
    }
    else {
        sendBuffer[DIRECTION::FRONT].resize(1);
        receiveBuffer[DIRECTION::FRONT].resize(1);
    }

    /// start timing (we don't want any setup time to be included, thus we start it just before the time loop)
   
    auto start = MPI_Wtime();


    /// main time loop
    /**
     * this is where we solve the actual partial differential equation and do the communication among processors.
          */



   



    for (unsigned time = 0; time < iterMax; ++time)
    {
        /// copy the solution from the previous timestep into T, which holds the solution of the last iteration
         
        for (unsigned i = 0; i < chunck[COORDINATE::X]; ++i)
            for (unsigned j = 0; j < chunck[COORDINATE::Y]; ++j)
                for (unsigned k = 0; k < chunck[COORDINATE::Z]; ++k)
                    T0[i][j][k] = T[i][j][k];

        // HALO communication step



  /// preparing the send buffer (the data we want to send to the left neighbor), if a neighbor exists


  /**
   * for simplicity, we write the 2D array (the face on the boundary) into a 1D array which we can easily send.
   * It is important that once we receive the it we are aware that the array containing the data is 1D now.

   */
        unsigned counter = 0;
        if (neighbors[DIRECTION::LEFT] != MPI_PROC_NULL)
            for (unsigned j = 1; j < chunck[COORDINATE::Y] - 1; ++j)
                for (unsigned k = 1; k < chunck[COORDINATE::Z] - 1; ++k)
                    sendBuffer[DIRECTION::LEFT][counter++] = T0[1][j][k];

        /// preparing the send buffer (the data we want to send to the right neighbor), if a neighbor exists
       
        counter = 0;
        if (neighbors[DIRECTION::RIGHT] != MPI_PROC_NULL)
            for (unsigned j = 1; j < chunck[COORDINATE::Y] - 1; ++j)
                for (unsigned k = 1; k < chunck[COORDINATE::Z] - 1; ++k)
                    sendBuffer[DIRECTION::RIGHT][counter++] = T0[chunck[COORDINATE::X] - 2][j][k];

        /// preparing the send buffer (the data we want to send to the bottom neighbor), if a neighbor exists
                counter = 0;
        if (neighbors[DIRECTION::BOTTOM] != MPI_PROC_NULL)
            for (unsigned i = 1; i < chunck[COORDINATE::X] - 1; ++i)
                for (unsigned k = 1; k < chunck[COORDINATE::Z] - 1; ++k)
                    sendBuffer[DIRECTION::BOTTOM][counter++] = T0[i][1][k];

        /// preparing the send buffer (the data we want to send to the top neighbor), if a neighbor exists
      

        counter = 0;
        if (neighbors[DIRECTION::TOP] != MPI_PROC_NULL)
            for (unsigned i = 1; i < chunck[COORDINATE::X] - 1; ++i)
                for (unsigned k = 1; k < chunck[COORDINATE::Z] - 1; ++k)
                    sendBuffer[DIRECTION::TOP][counter++] = T0[i][chunck[COORDINATE::Y] - 2][k];

        /// preparing the send buffer (the data we want to send to the back neighbor), if a neighbor exists
                counter = 0;
        if (neighbors[DIRECTION::BACK] != MPI_PROC_NULL)
            for (unsigned i = 1; i < chunck[COORDINATE::X] - 1; ++i)
                for (unsigned j = 1; j < chunck[COORDINATE::Y] - 1; ++j)
                    sendBuffer[DIRECTION::BACK][counter++] = T0[i][j][1];

        /// preparing the send buffer (the data we want to send to the front neighbor), if a neighbor exists
       
        counter = 0;
        if (neighbors[DIRECTION::FRONT] != MPI_PROC_NULL)
            for (unsigned i = 1; i < chunck[COORDINATE::X] - 1; ++i)
                for (unsigned j = 1; j < chunck[COORDINATE::Y] - 1; ++j)
                    sendBuffer[DIRECTION::FRONT][counter++] = T0[i][j][chunck[COORDINATE::Z] - 2];


       

        /// prepare the tags we need to append to the send message for each send (in each direction) and receive
        

        for (unsigned index = 0; index < NUMBER_OF_DIMENSIONS * 2; ++index) {
            tagSend[index] = 100 + neighbors[index];
            tagReceive[index] = 100 + rank;
        }

        /// send the prepared send buffer to the neighbors using non-blocking MPI_Isend(...)
                MPI_Isend(&sendBuffer[DIRECTION::LEFT][0], (chunck[COORDINATE::Y] - 1) * (chunck[COORDINATE::Z] - 1),
            MPI_FLOAT_T, neighbors[DIRECTION::LEFT], tagSend[DIRECTION::LEFT], MPI_COMM_CART,
            &request[DIRECTION::LEFT]);

        MPI_Isend(&sendBuffer[DIRECTION::RIGHT][0], (chunck[COORDINATE::Y] - 1) * (chunck[COORDINATE::Z] - 1),
            MPI_FLOAT_T, neighbors[DIRECTION::RIGHT], tagSend[DIRECTION::RIGHT], MPI_COMM_CART,
            &request[DIRECTION::RIGHT]);

        MPI_Isend(&sendBuffer[DIRECTION::BOTTOM][0], (chunck[COORDINATE::X] - 1) * (chunck[COORDINATE::Z] - 1),
            MPI_FLOAT_T, neighbors[DIRECTION::BOTTOM], tagSend[DIRECTION::BOTTOM], MPI_COMM_CART,
            &request[DIRECTION::BOTTOM]);

        MPI_Isend(&sendBuffer[DIRECTION::TOP][0], (chunck[COORDINATE::X] - 1) * (chunck[COORDINATE::Z] - 1),
            MPI_FLOAT_T, neighbors[DIRECTION::TOP], tagSend[DIRECTION::TOP], MPI_COMM_CART,
            &request[DIRECTION::TOP]);

        MPI_Isend(&sendBuffer[DIRECTION::BACK][0], (chunck[COORDINATE::X] - 1) * (chunck[COORDINATE::Y] - 1),
            MPI_FLOAT_T, neighbors[DIRECTION::BACK], tagSend[DIRECTION::BACK], MPI_COMM_CART,
            &request[DIRECTION::BACK]);

        MPI_Isend(&sendBuffer[DIRECTION::FRONT][0], (chunck[COORDINATE::X] - 1) * (chunck[COORDINATE::Y] - 1),
            MPI_FLOAT_T, neighbors[DIRECTION::FRONT], tagSend[DIRECTION::FRONT], MPI_COMM_CART,
            &request[DIRECTION::FRONT]);


        /*****************************************************************************************************************
                                                          GPU BEGIN
   ****************************************************************************************************************/
        // compute internal domain (no halos required)
        
	//Assign GPU to each process
	int deviceCount;

  	cudaGetDeviceCount(&deviceCount);
  	int device_id = rank%deviceCount;
  	cudaSetDevice(device_id); 
	
        int numX = chunck[COORDINATE::X];
        int numY = chunck[COORDINATE::Y];
        int numZ = chunck[COORDINATE::X];
	int num = numX*numY*numZ;
        double ***TBegin = CreateGrid(numX, numY, numZ);
	double ***TEnd = CreateGrid(numX, numY, numZ);
	double ***Tres = CreateGrid(numX, numY, numZ);
	double ***Tres_gpu = CreateGrid(numX, numY, numZ);
	double ***Thost = CreateGrid(numX, numY, numZ);
	double ***Thost0 = CreateGrid(numX, numY, numZ);
       
	for (unsigned i = 0; i < chunck[COORDINATE::X]; ++i){
            for (unsigned j = 0; j < chunck[COORDINATE::Y]; ++j){
                for (unsigned k = 0; k < chunck[COORDINATE::Z]; ++k){
                    *(*(*(Thost0 + i) + j) + k)=T0[i][j][k];
			*(*(*(Thost + i) + j) + k)=T[i][j][k];
	//printf("Thost0: %f", Thost0[i][j][k]);
	//printf("\n");


		}
	    }
	}
	
	//malloc memory 
        cudaMalloc((void**)&TBegin, sizeof(double) * num);
        cudaMalloc((void**)&TEnd, sizeof(double) * num);
	cudaMalloc((void**)&Tres_gpu, sizeof(double) * num);


        cudaMemcpy(TBegin, Thost0, sizeof(double) * num, cudaMemcpyHostToDevice);
        cudaMemcpy(TEnd, Thost, sizeof(double) * num, cudaMemcpyHostToDevice);
	
	//set grid and block size
        dim3 block(DIM_THREAD_BLOCK_X, DIM_THREAD_BLOCK_Y);
        dim3 grid((size_t)ceil((double)(numX-2)/ ((double)DIM_THREAD_BLOCK_X)), 1);
	
	//gpu compute kernel
        computeT <<<grid,block>>> (TBegin, TEnd, Tres_gpu,numX, numY, numZ, Dx, Dy,Dz);
        
	//pass data from gpu to cpu
        cudaMemcpy(Tres, Tres_gpu, sizeof(double) *num, cudaMemcpyDeviceToHost);
	cudaThreadSynchronize();

        for (unsigned i = 0; i < chunck[COORDINATE::X]; ++i){
            for (unsigned j = 0; j < chunck[COORDINATE::Y]; ++j){
                for (unsigned k = 0; k < chunck[COORDINATE::Z]; ++k){
                    T[i][j][k]=Tres[i][j][k];
			
		}
	    }
	}
	FreeGrid(TBegin);
	FreeGrid(TEnd);
	FreeGrid(Tres);
	FreeGrid(Tres_gpu);
	FreeGrid(Thost);
	FreeGrid(Thost0);
	cudaFree(TBegin);
	cudaFree(TEnd);
	cudaFree(Tres_gpu);

        /// now work on the halo cells
       

  /// receive the halo information from each neighbor, if exists.


        MPI_Recv(&receiveBuffer[DIRECTION::LEFT][0], (chunck[COORDINATE::Y] - 1) * (chunck[COORDINATE::Z] - 1),
            MPI_FLOAT_T, neighbors[DIRECTION::LEFT], tagReceive[DIRECTION::LEFT], MPI_COMM_CART,
            &status[DIRECTION::LEFT]);

        MPI_Recv(&receiveBuffer[DIRECTION::RIGHT][0], (chunck[COORDINATE::Y] - 1) * (chunck[COORDINATE::Z] - 1),
            MPI_FLOAT_T, neighbors[DIRECTION::RIGHT], tagReceive[DIRECTION::RIGHT], MPI_COMM_CART,
            &status[DIRECTION::RIGHT]);

        MPI_Recv(&receiveBuffer[DIRECTION::BOTTOM][0], (chunck[COORDINATE::X] - 1) * (chunck[COORDINATE::Z] - 1),
            MPI_FLOAT_T, neighbors[DIRECTION::BOTTOM], tagReceive[DIRECTION::BOTTOM], MPI_COMM_CART,
            &status[DIRECTION::BOTTOM]);

        MPI_Recv(&receiveBuffer[DIRECTION::TOP][0], (chunck[COORDINATE::X] - 1) * (chunck[COORDINATE::Z] - 1),
            MPI_FLOAT_T, neighbors[DIRECTION::TOP], tagReceive[DIRECTION::TOP], MPI_COMM_CART,
            &status[DIRECTION::TOP]);

        MPI_Recv(&receiveBuffer[DIRECTION::BACK][0], (chunck[COORDINATE::X] - 1) * (chunck[COORDINATE::Y] - 1),
            MPI_FLOAT_T, neighbors[DIRECTION::BACK], tagReceive[DIRECTION::BACK], MPI_COMM_CART,
            &status[DIRECTION::BACK]);

        MPI_Recv(&receiveBuffer[DIRECTION::FRONT][0], (chunck[COORDINATE::X] - 1) * (chunck[COORDINATE::Y] - 1),
            MPI_FLOAT_T, neighbors[DIRECTION::FRONT], tagReceive[DIRECTION::FRONT], MPI_COMM_CART,
            &status[DIRECTION::FRONT]);

        /// make sure that all communications have been executed
        
        /**
         * even though we use a blocking receive here, since we used a non-blocking send, we have to wait for all
         * communications to have finished before continuing.
        
         */
        MPI_Waitall(NUMBER_OF_DIMENSIONS * 2, request, status);

        /// now that we have the halo cells, we update the boundaries using information from other processors
      

        if (neighbors[DIRECTION::LEFT] != MPI_PROC_NULL) {
            const auto& THalo = receiveBuffer[DIRECTION::LEFT];
            unsigned i = 0;
            unsigned counter = 0;

            for (unsigned j = 1; j < chunck[COORDINATE::Y] - 1; ++j)
                for (unsigned k = 1; k < chunck[COORDINATE::Z] - 1; ++k) {
                    T[i][j][k] = T0[i][j][k] +
                        Dx * (T0[i + 1][j][k] - 2.0 * T0[i][j][k] + THalo[counter++]) +
                        Dy * (T0[i][j + 1][k] - 2.0 * T0[i][j][k] + T0[i][j - 1][k]) +
                        Dz * (T0[i][j][k + 1] - 2.0 * T0[i][j][k] + T0[i][j][k - 1]);
                }
        }

        /// do the same as above, this time for the right neighbor halo data
      

        if (neighbors[DIRECTION::RIGHT] != MPI_PROC_NULL) {
            const auto& THalo = receiveBuffer[DIRECTION::RIGHT];
            unsigned i = chunck[COORDINATE::X] - 1;
            unsigned counter = 0;

            for (unsigned j = 1; j < chunck[COORDINATE::Y] - 1; ++j)
                for (unsigned k = 1; k < chunck[COORDINATE::Z] - 1; ++k) {
                    T[i][j][k] = T0[i][j][k] +
                        Dx * (THalo[counter++] - 2.0 * T0[i][j][k] + T0[i - 1][j][k]) +
                        Dy * (T0[i][j + 1][k] - 2.0 * T0[i][j][k] + T0[i][j - 1][k]) +
                        Dz * (T0[i][j][k + 1] - 2.0 * T0[i][j][k] + T0[i][j][k - 1]);
                }
        }

        /// do the same as above, this time for the bottom neighbor halo data
       
        if (neighbors[DIRECTION::BOTTOM] != MPI_PROC_NULL) {
            const auto& THalo = receiveBuffer[DIRECTION::BOTTOM];
            unsigned j = 0;
            unsigned counter = 0;

            for (unsigned i = 1; i < chunck[COORDINATE::X] - 1; ++i)
                for (unsigned k = 1; k < chunck[COORDINATE::Z] - 1; ++k) {
                    T[i][j][k] = T0[i][j][k] +
                        Dx * (T0[i + 1][j][k] - 2.0 * T0[i][j][k] + T0[i - 1][j][k]) +
                        Dy * (T0[i][j + 1][k] - 2.0 * T0[i][j][k] + THalo[counter++]) +
                        Dz * (T0[i][j][k + 1] - 2.0 * T0[i][j][k] + T0[i][j][k - 1]);
                }
        }

        /// do the same as above, this time for the top neighbor halo data
      
        if (neighbors[DIRECTION::TOP] != MPI_PROC_NULL) {
            const auto& THalo = receiveBuffer[DIRECTION::TOP];
            unsigned j = chunck[COORDINATE::Y] - 1;
            unsigned counter = 0;

            for (unsigned i = 1; i < chunck[COORDINATE::X] - 1; ++i)
                for (unsigned k = 1; k < chunck[COORDINATE::Z] - 1; ++k) {
                    T[i][j][k] = T0[i][j][k] +
                        Dx * (T0[i + 1][j][k] - 2.0 * T0[i][j][k] + T0[i - 1][j][k]) +
                        Dy * (THalo[counter++] - 2.0 * T0[i][j][k] + T0[i][j - 1][k]) +
                        Dz * (T0[i][j][k + 1] - 2.0 * T0[i][j][k] + T0[i][j][k - 1]);
                }
        }

        /// do the same as above, this time for the back neighbor halo data
        
        if (neighbors[DIRECTION::BACK] != MPI_PROC_NULL) {
            const auto& THalo = receiveBuffer[DIRECTION::BACK];
            unsigned k = 0;
            unsigned counter = 0;

            for (unsigned i = 1; i < chunck[COORDINATE::X] - 1; ++i)
                for (unsigned j = 1; j < chunck[COORDINATE::Y] - 1; ++j) {
                    T[i][j][k] = T0[i][j][k] +
                        Dx * (T0[i + 1][j][k] - 2.0 * T0[i][j][k] + T0[i - 1][j][k]) +
                        Dy * (T0[i][j + 1][k] - 2.0 * T0[i][j][k] + T0[i][j - 1][k]) +
                        Dz * (T0[i][j][k + 1] - 2.0 * T0[i][j][k] + THalo[counter++]);
                }
        }

        /// do the same as above, this time for the front neighbor halo data
       
        if (neighbors[DIRECTION::FRONT] != MPI_PROC_NULL) {
            const auto& THalo = receiveBuffer[DIRECTION::FRONT];
            unsigned k = chunck[COORDINATE::Z] - 1;
            unsigned counter = 0;

            for (unsigned i = 1; i < chunck[COORDINATE::X] - 1; ++i)
                for (unsigned j = 1; j < chunck[COORDINATE::Y] - 1; ++j) {
                    T[i][j][k] = T0[i][j][k] +
                        Dx * (T0[i + 1][j][k] - 2.0 * T0[i][j][k] + T0[i - 1][j][k]) +
                        Dy * (T0[i][j + 1][k] - 2.0 * T0[i][j][k] + T0[i][j - 1][k]) +
                        Dz * (THalo[counter++] - 2.0 * T0[i][j][k] + T0[i][j][k - 1]);
                }
        }
        /************************************************************************************************************

                                                                GPU      END

       ***********************************************************************************************************/
        /// update edges of halo elements
       
        if (neighbors[DIRECTION::LEFT] != MPI_PROC_NULL) {
            if (neighbors[DIRECTION::BOTTOM] != MPI_PROC_NULL) {
                unsigned i = 0;
                unsigned j = 0;
                for (unsigned k = 1; k < chunck[COORDINATE::Z] - 1; ++k)
                    T[i][j][k] = 2.0 * T[i + 1][j][k] - T[i + 2][j][k];
            }
            if (neighbors[DIRECTION::TOP] != MPI_PROC_NULL) {
                unsigned i = 0;
                unsigned j = chunck[COORDINATE::Y] - 1;
                for (unsigned k = 1; k < chunck[COORDINATE::Z] - 1; ++k)
                    T[i][j][k] = 2.0 * T[i + 1][j][k] - T[i + 2][j][k];
            }
            if (neighbors[DIRECTION::BACK] != MPI_PROC_NULL) {
                unsigned i = 0;
                unsigned k = 0;
                for (unsigned j = 1; j < chunck[COORDINATE::Y] - 1; ++j)
                    T[i][j][k] = 2.0 * T[i + 1][j][k] - T[i + 2][j][k];
            }
            if (neighbors[DIRECTION::FRONT] != MPI_PROC_NULL) {
                unsigned i = 0;
                unsigned k = chunck[COORDINATE::Z] - 1;
                for (unsigned j = 1; j < chunck[COORDINATE::Y] - 1; ++j)
                    T[i][j][k] = 2.0 * T[i + 1][j][k] - T[i + 2][j][k];
            }
        }

        if (neighbors[DIRECTION::RIGHT] != MPI_PROC_NULL) {
            if (neighbors[DIRECTION::BOTTOM] != MPI_PROC_NULL) {
                unsigned i = chunck[COORDINATE::X] - 1;
                unsigned j = 0;
                for (unsigned k = 1; k < chunck[COORDINATE::Z] - 1; ++k)
                    T[i][j][k] = 2.0 * T[i - 1][j][k] - T[i - 2][j][k];
            }
            if (neighbors[DIRECTION::TOP] != MPI_PROC_NULL) {
                unsigned i = chunck[COORDINATE::X] - 1;
                unsigned j = chunck[COORDINATE::Y] - 1;
                for (unsigned k = 1; k < chunck[COORDINATE::Z] - 1; ++k)
                    T[i][j][k] = 2.0 * T[i - 1][j][k] - T[i - 2][j][k];
            }
            if (neighbors[DIRECTION::BACK] != MPI_PROC_NULL) {
                unsigned i = chunck[COORDINATE::X] - 1;
                unsigned k = 0;
                for (unsigned j = 1; j < chunck[COORDINATE::Y] - 1; ++j)
                    T[i][j][k] = 2.0 * T[i - 1][j][k] - T[i - 2][j][k];
            }
            if (neighbors[DIRECTION::FRONT] != MPI_PROC_NULL) {
                unsigned i = chunck[COORDINATE::X] - 1;
                unsigned k = chunck[COORDINATE::Z] - 1;
                for (unsigned j = 1; j < chunck[COORDINATE::Y] - 1; ++j)
                    T[i][j][k] = 2.0 * T[i - 1][j][k] - T[i - 2][j][k];
            }
        }

        if (neighbors[DIRECTION::BACK] != MPI_PROC_NULL) {
            if (neighbors[DIRECTION::BOTTOM] != MPI_PROC_NULL) {
                unsigned j = 0;
                unsigned k = 0;
                for (unsigned i = 1; i < chunck[COORDINATE::X] - 1; ++i)
                    T[i][j][k] = 2.0 * T[i][j][k + 1] - T[i][j][k + 2];
            }
            if (neighbors[DIRECTION::TOP] != MPI_PROC_NULL) {
                unsigned j = chunck[COORDINATE::Y] - 1;
                unsigned k = 0;
                for (unsigned i = 1; i < chunck[COORDINATE::X] - 1; ++i)
                    T[i][j][k] = 2.0 * T[i][j][k + 1] - T[i][j][k + 2];
            }
        }

        if (neighbors[DIRECTION::FRONT] != MPI_PROC_NULL) {
            if (neighbors[DIRECTION::BOTTOM] != MPI_PROC_NULL) {
                unsigned j = 0;
                unsigned k = chunck[COORDINATE::Z] - 1;
                for (unsigned i = 1; i < chunck[COORDINATE::X] - 1; ++i)
                    T[i][j][k] = 2.0 * T[i][j][k - 1] - T[i][j][k - 2];
            }
            if (neighbors[DIRECTION::TOP] != MPI_PROC_NULL) {
                unsigned j = chunck[COORDINATE::Y] - 1;
                unsigned k = chunck[COORDINATE::Z] - 1;
                for (unsigned i = 1; i < chunck[COORDINATE::X] - 1; ++i)
                    T[i][j][k] = 2.0 * T[i][j][k - 1] - T[i][j][k - 2];
            }
        }
        /// finished with halo edges extrapolation
       

        /// at last, we update the boundary points through weighted averages
      
        if ((neighbors[DIRECTION::LEFT] != MPI_PROC_NULL) && (neighbors[DIRECTION::BOTTOM] != MPI_PROC_NULL) &&
            (neighbors[DIRECTION::BACK] != MPI_PROC_NULL)) {
            unsigned i = 0;
            unsigned j = 0;
            unsigned k = 0;
            T[i][j][k] = 1.0 / 3.0 * (T[i + 1][j][k] + T[i][j + 1][k] + T[i][j][k + 1]);
        }

        if ((neighbors[DIRECTION::LEFT] != MPI_PROC_NULL) && (neighbors[DIRECTION::BOTTOM] != MPI_PROC_NULL) &&
            (neighbors[DIRECTION::FRONT] != MPI_PROC_NULL)) {
            unsigned i = 0;
            unsigned j = 0;
            unsigned k = chunck[COORDINATE::Z] - 1;
            T[i][j][k] = 1.0 / 3.0 * (T[i + 1][j][k] + T[i][j + 1][k] + T[i][j][k - 1]);
        }

        if ((neighbors[DIRECTION::LEFT] != MPI_PROC_NULL) && (neighbors[DIRECTION::TOP] != MPI_PROC_NULL) &&
            (neighbors[DIRECTION::BACK] != MPI_PROC_NULL)) {
            unsigned i = 0;
            unsigned j = chunck[COORDINATE::Y] - 1;
            unsigned k = 0;
            T[i][j][k] = 1.0 / 3.0 * (T[i + 1][j][k] + T[i][j - 1][k] + T[i][j][k + 1]);
        }

        if ((neighbors[DIRECTION::LEFT] != MPI_PROC_NULL) && (neighbors[DIRECTION::TOP] != MPI_PROC_NULL) &&
            (neighbors[DIRECTION::FRONT] != MPI_PROC_NULL)) {
            unsigned i = 0;
            unsigned j = chunck[COORDINATE::Y] - 1;
            unsigned k = chunck[COORDINATE::Z] - 1;
            T[i][j][k] = 1.0 / 3.0 * (T[i + 1][j][k] + T[i][j - 1][k] + T[i][j][k - 1]);
        }

        if ((neighbors[DIRECTION::RIGHT] != MPI_PROC_NULL) && (neighbors[DIRECTION::BOTTOM] != MPI_PROC_NULL) &&
            (neighbors[DIRECTION::BACK] != MPI_PROC_NULL)) {
            unsigned i = chunck[COORDINATE::X] - 1;
            unsigned j = 0;
            unsigned k = 0;
            T[i][j][k] = 1.0 / 3.0 * (T[i - 1][j][k] + T[i][j + 1][k] + T[i][j][k + 1]);
        }

        if ((neighbors[DIRECTION::RIGHT] != MPI_PROC_NULL) && (neighbors[DIRECTION::BOTTOM] != MPI_PROC_NULL) &&
            (neighbors[DIRECTION::FRONT] != MPI_PROC_NULL)) {
            unsigned i = chunck[COORDINATE::X] - 1;
            unsigned j = 0;
            unsigned k = chunck[COORDINATE::Z] - 1;
            T[i][j][k] = 1.0 / 3.0 * (T[i - 1][j][k] + T[i][j + 1][k] + T[i][j][k - 1]);
        }

        if ((neighbors[DIRECTION::RIGHT] != MPI_PROC_NULL) && (neighbors[DIRECTION::TOP] != MPI_PROC_NULL) &&
            (neighbors[DIRECTION::BACK] != MPI_PROC_NULL)) {
            unsigned i = chunck[COORDINATE::X] - 1;
            unsigned j = chunck[COORDINATE::Y] - 1;
            unsigned k = 0;
            T[i][j][k] = 1.0 / 3.0 * (T[i - 1][j][k] + T[i][j - 1][k] + T[i][j][k + 1]);
        }

        if ((neighbors[DIRECTION::RIGHT] != MPI_PROC_NULL) && (neighbors[DIRECTION::TOP] != MPI_PROC_NULL) &&
            (neighbors[DIRECTION::FRONT] != MPI_PROC_NULL)) {
            unsigned i = chunck[COORDINATE::X] - 1;
            unsigned j = chunck[COORDINATE::Y] - 1;
            unsigned k = chunck[COORDINATE::Z] - 1;
            T[i][j][k] = 1.0 / 3.0 * (T[i - 1][j][k] + T[i][j - 1][k] + T[i][j][k - 1]);
        }
        /// finished with halo corner points
     


/// calculate the difference between the current and previous (last time step) solution.
  

        floatT res = std::numeric_limits<floatT>::min();
        for (unsigned i = 1; i < chunck[COORDINATE::X] - 1; ++i)
            for (unsigned j = 1; j < chunck[COORDINATE::Y] - 1; ++j)
                for (unsigned k = 1; k < chunck[COORDINATE::Z] - 1; ++k)
                    if (std::fabs(T[i][j][k] - T0[i][j][k]) > res)
                        res = std::fabs(T[i][j][k] - T0[i][j][k]);

        /// if it is the first time step, store the residual as the normalisation factor
      


        if (time == 0)
            if (res != 0.0)
                norm = res;

        /// For MPI, we have to communicate the norm by selecting the lowest among all processors
       

        if (time == 0) {
            MPI_Iallreduce(&norm, &globalNorm, 1, MPI_FLOAT_T, MPI_MIN, MPI_COMM_CART, &reduceRequest);
            MPI_Wait(&reduceRequest, MPI_STATUS_IGNORE);
        }


        /// if we want to debug, it may be useful to see the residuals. Turned of for release builds for performance.
          
//#if defined(USE_DEBUG)
//        if (rank == 0) {
//            std::cout << "time: " << std::setw(10) << time;
//            std::cout << std::scientific << std::setw(15) << std::setprecision(5) << ", residual: ";
//            std::cout << res / norm << std::endl;
//        }
//#endif

        /// check if the current residual has dropped below our defined convergence threshold "eps"
         
        if (res / norm < eps)
            breakCondition = true;

        /// Again, for MPI we need to among all processors if we can break from the loop
       


        MPI_Iallreduce(&breakCondition, &globalBreakCondition, 1, MPI_INT, MPI_MAX, MPI_COMM_CART, &reduceRequest);
        MPI_Wait(&reduceRequest, MPI_STATUS_IGNORE);


        /// final check if we can break, the above was just preparation for this check.
         


        if (globalBreakCondition) {
            finalNumIterations = time;
            break;
        }
    }
    /// done with the time loop
   

    /// output the timing information to screen.
    

    auto end = MPI_Wtime();
    if (rank == 0) {
        std::cout << "Computational time (parallel): " << std::fixed << (end - start) << "\n" << std::endl;
        if (globalBreakCondition) {
            std::cout << "Simulation has converged in " << finalNumIterations << " iterations";
            std::cout << " with a convergence threshold of " << std::scientific << eps << std::endl;
        }
        else
            std::cout << "Simulation did not converge within " << iterMax << " iterations." << std::endl;
    }


    /// calculate the error we have made against the analytic solution
   

    double globalError = 0.0;
    double error = 0.0;
    for (unsigned k = 1; k < chunck[COORDINATE::Z] - 1; ++k)
        for (unsigned j = 1; j < chunck[COORDINATE::Y] - 1; ++j)
            for (unsigned i = 1; i < chunck[COORDINATE::X] - 1; ++i)
                error += std::sqrt(std::pow(T[i][j][k] - (coordinates3D[COORDINATE::Y] * (chunck[COORDINATE::Y] - 1) + j) * spacing[COORDINATE::Y], 2.0));
    error /= ((chunck[COORDINATE::X] - 2) * (chunck[COORDINATE::Y] - 2) * (chunck[COORDINATE::Z] - 2));
    MPI_Iallreduce(&error, &globalError, 1, MPI_FLOAT_T, MPI_SUM, MPI_COMM_CART, &reduceRequest);
    MPI_Wait(&reduceRequest, MPI_STATUS_IGNORE);
    if (rank == 0)
        std::cout << "L2-norm error: " << std::fixed << std::setprecision(4) << 100 * error << " %" << std::endl;


    /// output the solution in a format readable by a post processor, such as paraview.
     

    std::vector<floatT> receiveBufferPostProcess;
    receiveBufferPostProcess.resize(chunck[COORDINATE::X] * chunck[COORDINATE::Y] * chunck[COORDINATE::Z]);
    if (rank > 0 && size != 1)
    {
        int counter = 0;
        for (unsigned k = 0; k < chunck[COORDINATE::Z]; ++k)
            for (unsigned j = 0; j < chunck[COORDINATE::Y]; ++j)
                for (unsigned i = 0; i < chunck[COORDINATE::X]; ++i)
                    receiveBufferPostProcess[counter++] = T[i][j][k];

        MPI_Send(&receiveBufferPostProcess[0], chunck[COORDINATE::X] * chunck[COORDINATE::Y] * chunck[COORDINATE::Z], MPI_FLOAT_T, 0, 200 + rank, MPI_COMM_CART);
        MPI_Send(&coordinates3D[0], NUMBER_OF_DIMENSIONS, MPI_INT, 0, 300 + rank, MPI_COMM_CART);
    }
    if (rank == 0 && size != 1)
    {
        std::ofstream out("output/out.dat");
        out << "TITLE=\"out\"" << std::endl;
        out << "VARIABLES = \"X\", \"Y\", \"Z\", \"T\", \"rank\"" << std::endl;
        out << "ZONE T = \"" << rank << "\", I=" << chunck[COORDINATE::X] << ", J=" << chunck[COORDINATE::Y] << ", K=" << chunck[COORDINATE::Z] << ", F=POINT" << std::endl;
        for (unsigned k = 0; k < chunck[COORDINATE::Z]; ++k)
            for (unsigned j = 0; j < chunck[COORDINATE::Y]; ++j)
                for (unsigned i = 0; i < chunck[COORDINATE::X]; ++i)
                {
                    out << std::scientific << std::setprecision(5) << std::setw(15) << (coordinates3D[COORDINATE::X] * (chunck[COORDINATE::X] - 1) + i) * spacing[COORDINATE::X];
                    out << std::scientific << std::setprecision(5) << std::setw(15) << (coordinates3D[COORDINATE::Y] * (chunck[COORDINATE::Y] - 1) + j) * spacing[COORDINATE::Y];
                    out << std::scientific << std::setprecision(5) << std::setw(15) << (coordinates3D[COORDINATE::Z] * (chunck[COORDINATE::Z] - 1) + k) * spacing[COORDINATE::Z];
                    out << std::scientific << std::setprecision(5) << std::setw(15) << T[i][j][k];
                    out << std::fixed << std::setw(5) << rank << std::endl;
                }

        for (int recvRank = 1; recvRank < size; ++recvRank)
        {
            int coordinates3DFromReceivedRank[NUMBER_OF_DIMENSIONS];
            MPI_Recv(&receiveBufferPostProcess[0], chunck[COORDINATE::X] * chunck[COORDINATE::Y] * chunck[COORDINATE::Z], MPI_FLOAT_T, recvRank, 200 + recvRank, MPI_COMM_CART, &postStatus[0]);
            MPI_Recv(&coordinates3DFromReceivedRank[0], NUMBER_OF_DIMENSIONS, MPI_INT, recvRank, 300 + recvRank, MPI_COMM_CART, &postStatus[1]);

            out << "ZONE T = \"" << rank << "\", I=" << chunck[COORDINATE::X] << ", J=" << chunck[COORDINATE::Y] << ", K=" << chunck[COORDINATE::Z] << ", F=POINT" << std::endl;
            int counter = 0;
            for (unsigned k = 0; k < chunck[COORDINATE::Z]; ++k)
                for (unsigned j = 0; j < chunck[COORDINATE::Y]; ++j)
                    for (unsigned i = 0; i < chunck[COORDINATE::X]; ++i)
                    {
                        out << std::scientific << std::setprecision(5) << std::setw(15) << (coordinates3DFromReceivedRank[COORDINATE::X] * (chunck[COORDINATE::X] - 1) + i) * spacing[COORDINATE::X];
                        out << std::scientific << std::setprecision(5) << std::setw(15) << (coordinates3DFromReceivedRank[COORDINATE::Y] * (chunck[COORDINATE::Y] - 1) + j) * spacing[COORDINATE::Y];
                        out << std::scientific << std::setprecision(5) << std::setw(15) << (coordinates3DFromReceivedRank[COORDINATE::Z] * (chunck[COORDINATE::Z] - 1) + k) * spacing[COORDINATE::Z];
                        out << std::scientific << std::setprecision(5) << std::setw(15) << receiveBufferPostProcess[counter++];
                        out << std::fixed << std::setw(5) << recvRank << std::endl;
                    }
        }
        out.close();
    }
    if (size == 1)
    {
        std::ofstream out("output/out.dat");
        out << "TITLE=\"out\"" << std::endl;
        out << "VARIABLES = \"X\", \"Y\", \"Z\", \"T\"" << std::endl;
        out << "ZONE T = \"" << rank << "\", I=" << chunck[COORDINATE::X] << ", J=" << chunck[COORDINATE::Y] << ", K=" << chunck[COORDINATE::Z] << ", F=POINT" << std::endl;
        for (unsigned k = 0; k < chunck[COORDINATE::Z]; ++k)
            for (unsigned j = 0; j < chunck[COORDINATE::Y]; ++j)
                for (unsigned i = 0; i < chunck[COORDINATE::X]; ++i)
                {
                    out << std::scientific << std::setprecision(5) << std::setw(15) << (coordinates3D[COORDINATE::X] * (chunck[COORDINATE::X] - 1) + i) * spacing[COORDINATE::X];
                    out << std::scientific << std::setprecision(5) << std::setw(15) << (coordinates3D[COORDINATE::Y] * (chunck[COORDINATE::Y] - 1) + j) * spacing[COORDINATE::Y];
                    out << std::scientific << std::setprecision(5) << std::setw(15) << (coordinates3D[COORDINATE::Z] * (chunck[COORDINATE::Z] - 1) + k) * spacing[COORDINATE::Z];
                    out << std::scientific << std::setprecision(5) << std::setw(15) << T[i][j][k] << std::endl;
                }
        out.close();
    }



    MPI_Finalize();

    return 0;
}