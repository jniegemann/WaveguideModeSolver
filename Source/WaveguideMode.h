#ifndef _WAVEGUIDEMODE_H
#define _WAVEGUIDEMODE_H

#include "function/solution.h"
#include "function/filter.h"

#include "views/vector_view.h"
#include "views/scalar_view.h"

/**
 *
 */
template<typename Scalar>
class WaveguideMode
{
  public:

	/**
	 * Constructor
	 */
    WaveguideMode(double in_k0, Scalar in_kz, Hermes::Hermes2D::MeshFunctionSharedPtr<Scalar> in_slnEt, Hermes::Hermes2D::MeshFunctionSharedPtr<Scalar> in_slnEz)
      : k0(in_k0),kz(in_kz), slnEt(in_slnEt), slnEz(in_slnEz)
    {

    }

    /**
     * Calculates the energy flux of this mode in z-direction
     */
    Scalar calculateEnergyFlux();

    /**
     * Display the real part of the tangential and the normal component of the mode
     */
    void plot();

    Scalar getEffectiveIndex() const {return kz/k0;}
    Scalar getKz() const {return kz;}
    double getK0() const {return k0;}

  protected:
    double k0;
    Scalar kz;

    Hermes::Hermes2D::MeshFunctionSharedPtr<Scalar> slnEt;
    Hermes::Hermes2D::MeshFunctionSharedPtr<Scalar> slnEz;
};


/**
 * Normalizes the mode such that the flux of energy in z-direction is unity!
 */
template<typename Scalar>
Scalar WaveguideMode<Scalar>::calculateEnergyFlux()
{
  using namespace Hermes::Hermes2D;

  Scalar overlapIntegral = 0.0;

//  // Generate a union mesh for the integration
//  const unsigned int num = 2;
//
//  const Mesh* meshes[num];
//  meshes[0] = slnEt->get_mesh();
//  meshes[1] = slnEz->get_mesh();
//
//  Mesh* unionMesh = new Mesh;
//  Traverse trav;
//  trav.begin(num, meshes);
//  UniData** unidata = trav.construct_union_mesh(unionMesh);
//  trav.finish();
//  const unsigned int order = 4;
//
//  // element for the element loop
//  Element *e;
//
//  // refmap for computing Jacobian
//  RefMap rm;
//  rm.set_quad_2d(&g_quad_2d_std);
//
//  // Iterate over each element of the union mesh
//  for_all_active_elements(e, unionMesh)
//  {
//	// Get the respective element
//	slnEt->set_quad_2d(&g_quad_2d_std);
//	slnEt->set_active_element(unidata[0][e->id].e);
//	slnEt->set_transform(unidata[0][e->id].idx);
//	slnEt->set_quad_order(order);
//
//	slnEz->set_quad_2d(&g_quad_2d_std);
//	slnEz->set_active_element(unidata[1][e->id].e);
//	slnEz->set_transform(unidata[1][e->id].idx);
//	slnEz->set_quad_order(order);
//
//	rm.set_active_element(e);
//
//	// get the quadrature points
//	int np = g_quad_2d_std.get_num_points(order, e->get_mode());
//	double3 *pt = g_quad_2d_std.get_points(order, e->get_mode());
//
//	Scalar* Ex = slnEt->get_fn_values(0);
//	Scalar* Ey = slnEt->get_fn_values(1);
//	Scalar* dEzDx = slnEz->get_dx_values(0);
//	Scalar* dEzDy = slnEz->get_dy_values(0);
//	double* jac = rm.get_jacobian(order);
//
//	// loop over points and integrate the energy
//	for( int j = 0; j < np; ++j )
//	{
//	  const double w = pt[j][2];
////	  std::cout << "Pts: " << jac[j] << " " << w << ", " << Ex[j] << ", " << Ey[j] << ", " << dEzDx[j] << " " << dEzDy[j] << std::endl;
//
//	  overlapIntegral += w*jac[j]*( kz*(Ex[j]*Ex[j] + Ey[j]*Ey[j]) - (Ex[j]*dEzDx[j] + Ey[j]*dEzDy[j]));
//	}
//  }

  return overlapIntegral;
}

template<typename Scalar>
void WaveguideMode<Scalar>::plot()
{
    // Initialize views.
	Hermes::Hermes2D::Views::VectorView s_view_0("E_t", new Hermes::Hermes2D::Views::WinGeom(0, 0, 440, 350));
	Hermes::Hermes2D::Views::ScalarView s_view_1("E_z", new Hermes::Hermes2D::Views::WinGeom(880, 0, 440, 350));

	Hermes::Hermes2D::RealFilter slnEtReal(slnEt);
    s_view_0.show(&slnEtReal);
//    o_view_0.show(spaceEt);

    Hermes::Hermes2D::RealFilter slnEzReal(slnEz);
    s_view_1.show(&slnEzReal);
//    o_view_1.show(spaceEz);

    Hermes::Hermes2D::Views::View::wait(Hermes::Hermes2D::Views::HERMES_WAIT_KEYPRESS);
}


#endif
