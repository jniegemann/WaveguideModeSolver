#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"

#include "WaveguideMode.h"
#include "WaveguideDefinitionsFields2.h"
#include "ClockTickCounter.h"

//#include "Eigensolver/ArpackEigensolver.h"

#include <Eigen/SparseCore>
#include <unsupported/Eigen/SparseExtra>
#include <Eigen/Eigenvalues>

#include <boost/program_options.hpp>

#include <stdio.h>
#include <memory>
#include <chrono>
#include <iostream> 
#include <iomanip> 

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace boost::program_options;


// Define the parameters and material constants of our system
const double a = 1e-6;      // Length scale (here 1um)
const double c = 299792458; // Speed of light

/**
 * Routine which takes a mesh and a list of materials and computes eigenvalues for a given fixed order, closest to targetKz.
 */
template<typename T>
WaveguideMode<T> findMode(MeshSharedPtr mesh, std::vector<Material<T>> materialList, unsigned int order, double k0, T targetKz, bool extractEigenvectors)
{
  // Add dispersive (k0-dependent materials to the list)
  const T i(0.0,1.0);
  const T epsInf = 5.0;
  const T fp = 2175e12;
  const T gamma = 4.35e12;
  const T f = k0/(2*M_PI)*c/a;
  const T epsAg = epsInf - fp*fp/( f*(f+i*gamma));		// TODO: Remove hardcoded material properties

  materialList.push_back(Material<T>("Ag",epsAg, epsAg, epsAg, 1.0, 1.0, 1.0));

  const double sigma = 2.0;

  // Now, for each material, we also add the corresponding PML materials
  const std::vector<Material<T>> tempMaterialList(materialList);
  for( const auto& curMat : tempMaterialList )
  {
    const T sx(1,sigma/k0 );
    const T sy(1,sigma/k0 );

    materialList.push_back(Material<T>(curMat.name + "PmlX" , curMat.epsXX/sx,    curMat.epsYY*sx,    curMat.epsZZ*sx,    curMat.invMuXX*sx,    curMat.invMuYY/sx,    curMat.invMuZZ/sx));
    materialList.push_back(Material<T>(curMat.name + "PmlY" , curMat.epsXX*sy,    curMat.epsYY/sy,    curMat.epsZZ*sy,    curMat.invMuXX/sy,    curMat.invMuYY*sy,    curMat.invMuZZ/sy));
    materialList.push_back(Material<T>(curMat.name + "PmlXY", curMat.epsXX*sy/sx, curMat.epsYY*sx/sy, curMat.epsZZ*sx*sy, curMat.invMuXX*sx/sy, curMat.invMuYY*sy/sx, curMat.invMuZZ/(sx*sy)));
  }

  // Now, we have to make sure that we only have materials that actually occur in the mesh. This seems a new requirement in hermes2d 3.0
  std::vector<Material<T>> reducedMaterialList;
  for( const auto& curMaterial : materialList )
  {
	  if( mesh->get_element_markers_conversion().get_internal_marker(curMaterial.name).valid )
	  {
		  reducedMaterialList.push_back( curMaterial );
	  	  std::cout << "Adding material: " << curMaterial.name << std::endl;
	  }
  }

  // Initialize boundary conditions (PEC, i.e. homogeneous Dirichlet conditions)
  DefaultEssentialBCConst<T> bcEt("OuterBdy", 0.0);
  EssentialBCs<T> bcsEt(&bcEt);
  DefaultEssentialBCConst<T> bcEz("OuterBdy", 0.0);
  EssentialBCs<T> bcsEz(&bcEz);

  // Create an Hcurl space and an H1 space with default shapeset.
  SpaceSharedPtr<T>	spaceEt{ new HcurlSpace<T>(mesh, &bcsEt, order-1)};
  SpaceSharedPtr<T>  spaceEz{ new H1Space<T>(mesh, &bcsEz, order)};
  std::vector<SpaceSharedPtr<T>> spaceWaveguide({spaceEt, spaceEz});

  // Initialize mesh solutions for both field components
  MeshFunctionSharedPtr<T> slnEt{ new Solution<T>()};
  MeshFunctionSharedPtr<T> slnEz{ new Solution<T>()};

  // Build a combined space for the two component
  const unsigned int ndof = Space<T>::get_num_dofs(spaceWaveguide);
  std::cout << "Number of dofs: " << ndof << std::endl;

  // Assemble the left-hand matrix
  std::cout << "Starting Matrix Assembly..."; std::cout.flush();
  const auto assemblyStartTime = monotonic_seconds();
  WeakFormWaveguideLeft<T> wfWaveguideLeft(k0, reducedMaterialList);
  DiscreteProblem<T> dpWaveguideLeft(&wfWaveguideLeft, spaceWaveguide);
  auto matrixLeft = std::make_shared<CSCMatrix<T>>();
  dpWaveguideLeft.assemble(matrixLeft.get());

  // Assemble the right-hand matrix
  WeakFormWaveguideRight<T> wfWaveguideRight(k0, reducedMaterialList);
  DiscreteProblem<T> dpWaveguideRight(&wfWaveguideRight, spaceWaveguide);
  auto matrixRight = std::make_shared<CSCMatrix<T>>();
  dpWaveguideRight.assemble(matrixRight.get());
  const auto assemblyComputeTime = monotonic_seconds();

  std::cout << "done (took " << assemblyComputeTime - assemblyStartTime << "s)!" << std::endl;

  // Map the CSC-Matrix to an Eigen3-Matrix so we can use the convenient methods from Eigen
  Eigen::MappedSparseMatrix<T> L(matrixLeft->get_size(), matrixLeft->get_size(), matrixLeft->get_nnz(), matrixLeft->get_Ap(), matrixLeft->get_Ai(), matrixLeft->get_Ax());
  Eigen::MappedSparseMatrix<T> R(matrixRight->get_size(), matrixRight->get_size(), matrixRight->get_nnz(), matrixRight->get_Ap(), matrixRight->get_Ai(), matrixRight->get_Ax());

  // For debugging: Dump the matrices in matrix market format!
  // THOSE MATRICES ARE NOT CORRECT!!!
  Eigen::saveMarket(L,"Left.smat");
  Eigen::saveMarket(R,"Right.smat");
//
//  // Call iterative eigensolver
//  std::cout << "Starting Eigensolver..."; std::cout.flush();
//  const auto eigensolverStartTime = monotonic_seconds();
//  Math::ArpackSupport<Eigen::MappedSparseMatrix<T>> eigensolver(L, R, extractEigenvectors);
//  const auto eigensolverSetupTime = monotonic_seconds();
//  eigensolver.compute(-targetKz*targetKz, 4);
//  const auto eigensolverComputeTime = monotonic_seconds();
//  std::cout << "done (took " << eigensolverComputeTime - eigensolverStartTime << "s)!" << std::endl;
//
//  std::cout << "Starting eigenvalue extraction..."; std::cout.flush();
//  Eigen::VectorXcd eigenvaluesRaw = -eigensolver.eigenvalues();
//  std::cout << "done!" << std::endl;
//
//  // Iterate over all eigenvalues and keep only the interesting ones
//  std::cout << "Starting eigenvalue filtering..."; std::cout.flush();
//  std::vector<std::pair<std::complex<double>,unsigned int>> evList;
//  for(unsigned int curEvIdx=0; curEvIdx<eigenvaluesRaw.size(); ++curEvIdx)
//  {
//    const std::complex<double> curEv = eigenvaluesRaw[curEvIdx];
//    const T curKz = std::sqrt(curEv);
//    const T curNeff = curKz/k0;
//
//      // Remove eigenvalues close to zero or with either real or imaginary part negative!
//      if( (std::norm(curEv) > 0.001) && (std::real(curEv) > 0.0) && (std::imag(curNeff) > 0.0)  )
//      {
//        evList.push_back( std::make_pair(curEv, curEvIdx) );
//      }
//    }
//  std::cout << "done!" << std::endl;
//
//  std::cout << "Starting eigenvalue sorting:" << std::endl;
//    // Sort with closest to target!
//    std::sort(evList.begin(),evList.end(),[&](const std::pair<std::complex<double>,unsigned int>& a, const std::pair<std::complex<double>,unsigned int>& b) -> bool { return std::norm(a.first-targetKz*targetKz) < std::norm(b.first-targetKz*targetKz);} );
//
//    // Print all filtered modes!
//    for(auto& curEv : evList)
//    {
//      const std::complex<double> kz = std::sqrt(curEv.first);
//
//      std::cout << "n_eff = " << std::setprecision(16) << std::real(kz)/k0 << " " << std::imag(kz)/k0 << " Error: " << std::abs(std::real(kz) - std::real(targetKz))/std::abs(std::real(targetKz)) << ", " << std::abs(std::imag(kz) - std::imag(targetKz))/std::abs(std::imag(targetKz)) << " -- Target: " << targetKz << " Result: " << std::real(kz) << " " << std::imag(kz) << std::endl;
//    }
//
//    const unsigned int modeIdx = 0;
//
//    // Pick the first mode and convert to a hermes solution (i.e. split the components)
//    const std::complex<double> kz = std::sqrt(evList[modeIdx].first);
////    targetKz = kz;
//
//  if( extractEigenvectors )
//  {
//    Eigen::Matrix<T, Eigen::Dynamic, 1> curEigenvector = eigensolver.eigenvectors().col(evList[modeIdx].second);
//    Solution<T>::vector_to_solutions((T*)curEigenvector.data(), spaceWaveguide, Hermes::vector<MeshFunctionSharedPtr<T>>(slnEt, slnEz));
//  }

  // Just some fake eigenvalue so the code works without too many changes
  const std::complex<double> kz(0,0);

  return WaveguideMode<T>(k0, kz, slnEt, slnEz);
}


int main(int argc, char* argv[])
{
  typedef std::complex<double> T;	// We have complex matrices due to the PMLs and/or the metals!

  std::string meshFile;				///< Contains the filename of the mesh file
  unsigned int order = 4;			///< Default order
  unsigned int initRefinement = 1;	///< Default refinement
  bool plotResults = false;

  // Define the allowed command line parameters and parse them
  try
  {
    options_description description("Allowed options:");

    description.add_options()	("help,h", "Shows a list of available options")
				("mesh,m", value<std::string>(&meshFile)->required(), "Specifies the mesh file")
                ("order,p", value<unsigned int>(&order)->default_value(4), "Specifies the spatial order of the basis functions")
                ("show,s", value<bool>(&plotResults)->default_value(false), "Plots the waveguide mode with the mesh")
				("refinement,r", value<unsigned int>(&initRefinement)->default_value(1), "Specifies the initial uniform refinement!");

    // Evaluate the command line arguments using the boost::program_options library
    variables_map vm;
    store( command_line_parser(argc, argv).options(description).run(), vm);

    // If the user asked for help, show the list of options and leave
    if( vm.count("help") )
    {
      std::cout << description << std::endl;
      return EXIT_SUCCESS;
    }

    notify(vm);
  }
  catch(std::exception& e)
  {
    std::cerr << "Error: " << e.what() << "\n";
    abort();
  }

  // Load the initial mesh
  MeshSharedPtr meshEt(new Mesh);
  MeshReaderH2D mloader;
  mloader.load(meshFile.c_str(), meshEt);

  // Perform initial mesh refinements (optional).
  for (int i = 0; i < initRefinement; i++)
  {
    meshEt->refine_all_elements();
  } 

  // Prepare the materials (including PMLs)
  const double k0 = 4.053667940115862; // 1.55 um;

  // Here, we push a number of (hard-coded) material properties
  std::vector<Material<T>> materialList;
  const T epsSubstrate = 1.444*1.444;
  const T metalEps = T(-129.53,3.465);
  const T SiEps = T(12,0);
  const T HfO2Eps = T(3.5,0);
  const T tcoEps = T(3.84,0.0089);
  const T tcoOffEps = T(-0.29,0.62);
  materialList.push_back(Material<T>("Metal", metalEps, metalEps, metalEps, 1.0, 1.0, 1.0));
  materialList.push_back(Material<T>("Dielectric", SiEps, SiEps, SiEps, 1.0, 1.0, 1.0));
  materialList.push_back(Material<T>("HighK",HfO2Eps, HfO2Eps, HfO2Eps, 1.0, 1.0, 1.0));
  materialList.push_back(Material<T>("TcoOn", tcoEps, tcoEps, tcoEps, 1.0, 1.0, 1.0));
  materialList.push_back(Material<T>("TcoOff", tcoOffEps, tcoOffEps, tcoOffEps, 1.0, 1.0, 1.0));
  materialList.push_back(Material<T>("Substrate", epsSubstrate, epsSubstrate, epsSubstrate, 1.0, 1.0, 1.0));

  // Initial target (this is already a very accurate solution which also serves as a reference)!
  T kzTarget1 = T(13.03494294296519, 0.01380173768411812);

  auto wgMode1 = findMode(meshEt, materialList, order, k0, kzTarget1, plotResults);

  // View the solution if requested!
  if( plotResults )
  {
	  wgMode1.plot();
  }

  return 0;
};
